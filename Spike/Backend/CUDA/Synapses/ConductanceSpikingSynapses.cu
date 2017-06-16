// -*- mode: c++ -*-
#include "Spike/Backend/CUDA/Synapses/ConductanceSpikingSynapses.hpp"

SPIKE_EXPORT_BACKEND_TYPE(CUDA, ConductanceSpikingSynapses);

namespace Backend {
  namespace CUDA {
    // ConductanceSpikingSynapses Destructor
    ConductanceSpikingSynapses::~ConductanceSpikingSynapses() {
      CudaSafeCall(cudaFree(biological_conductance_scaling_constants_lambda));
      CudaSafeCall(cudaFree(synapse_decay_id));
      CudaSafeCall(cudaFree(neuron_wise_conductance_trace));
      CudaSafeCall(cudaFree(decay_term_values));
      CudaSafeCall(cudaFree(circular_spikenum_buffer));
      CudaSafeCall(cudaFree(spikeid_buffer));
      free(h_decay_term_values);
      free(h_synapse_decay_id);
      free(h_neuron_wise_conductance_trace);
    }

    void ConductanceSpikingSynapses::prepare() {
      SpikingSynapses::prepare();
      buffersize = (frontend()->maximum_axonal_delay_in_timesteps + 2);

      // Set up tau and reversal potential values and ids (Host-Side)
      h_synapse_decay_id = (int*)realloc(h_synapse_decay_id, frontend()->total_number_of_synapses*sizeof(int));
      // Prepare the tau synaptic conductance host-side vars
      for (int syn_id = 0; syn_id < frontend()->total_number_of_synapses; syn_id++){
        float tau_g = frontend()->decay_terms_tau_g[syn_id];
        float reversal_pot = frontend()->reversal_potentials_Vhat[syn_id];
        int id = -1;
        // Get Tau ID
        for (int i = 0; i < num_decay_terms; i++){
          // If this combination exists already, find it and assign the correct id
          if ((tau_g == h_decay_term_values[i]) && (reversal_pot == h_reversal_values[i]))
            id = i;       
        }
        if (id < 0){
	  // If this combination of tau/reversal potential doesn't exist, add it:
          num_decay_terms += 1;
          h_decay_term_values = (float*)realloc(h_decay_term_values, (num_decay_terms)*sizeof(float));
          h_reversal_values = (float*)realloc(h_reversal_values, (num_decay_terms)*sizeof(float));
          id = num_decay_terms - 1;
          h_decay_term_values[id] = tau_g;
          h_reversal_values[id] = reversal_pot;
        } 
        // Set this tau id
        h_synapse_decay_id[syn_id] = id;
      }
      // Set up per neuron conductances
      conductance_trace_length = frontend()->neuron_pop_size*num_decay_terms;
      h_neuron_wise_conductance_trace = (float*)realloc(h_neuron_wise_conductance_trace, conductance_trace_length*sizeof(float));
      for (int id = 0; id < conductance_trace_length; id++)
        h_neuron_wise_conductance_trace[id] = 0.0f;

      // Carry out remaining device actions
      allocate_device_pointers();
      copy_constants_and_initial_efficacies_to_device();
    }

    void ConductanceSpikingSynapses::reset_state() {
      SpikingSynapses::reset_state();
      CudaSafeCall(cudaMemcpy(
        neuron_wise_conductance_trace,
        h_neuron_wise_conductance_trace,
        sizeof(float)*conductance_trace_length, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemset(
	circular_spikenum_buffer,
	0, sizeof(int)*(buffersize)));
    }


    void ConductanceSpikingSynapses::allocate_device_pointers() {
      CudaSafeCall(cudaMalloc((void **)&biological_conductance_scaling_constants_lambda, sizeof(float)*frontend()->total_number_of_synapses));
      CudaSafeCall(cudaMalloc((void **)&active_synapse_indices, sizeof(int)*frontend()->total_number_of_synapses));
      CudaSafeCall(cudaMalloc((void **)&synapse_decay_id, sizeof(int)*frontend()->total_number_of_synapses));
      CudaSafeCall(cudaMalloc((void **)&neuron_wise_conductance_trace, sizeof(float)*conductance_trace_length));
      CudaSafeCall(cudaMalloc((void **)&decay_term_values, sizeof(float)*num_decay_terms));
      CudaSafeCall(cudaMalloc((void **)&reversal_values, sizeof(float)*num_decay_terms));
      CudaSafeCall(cudaMalloc((void **)&circular_spikenum_buffer, sizeof(int)*(buffersize)));
      CudaSafeCall(cudaMalloc((void **)&spikeid_buffer, sizeof(int)*(frontend()->maximum_axonal_delay_in_timesteps * frontend()->maximum_number_of_afferent_synapses)));
    }

    void ConductanceSpikingSynapses::copy_constants_and_initial_efficacies_to_device() {
      CudaSafeCall(cudaMemcpy(biological_conductance_scaling_constants_lambda,
        frontend()->biological_conductance_scaling_constants_lambda,
        sizeof(float)*frontend()->total_number_of_synapses, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemcpy(
        synapse_decay_id,
        h_synapse_decay_id,
        sizeof(int)*frontend()->total_number_of_synapses, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemcpy(
        neuron_wise_conductance_trace,
        h_neuron_wise_conductance_trace,
        sizeof(float)*conductance_trace_length, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemcpy(
        decay_term_values,
        h_decay_term_values,
        sizeof(float)*num_decay_terms, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemcpy(
        reversal_values,
        h_reversal_values,
        sizeof(float)*num_decay_terms, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemset(
	circular_spikenum_buffer,
	0, sizeof(int)*buffersize));

    }



    /* STATE UPDATE */
    void ConductanceSpikingSynapses::state_update
    (::SpikingNeurons* neurons,
     ::SpikingNeurons* input_neurons,
     float current_time_in_seconds, float timestep) {

      // Setting up access to neuron backends
      ::Backend::CUDA::SpikingNeurons* neurons_backend =
        dynamic_cast<::Backend::CUDA::SpikingNeurons*>(neurons->backend());
      assert(neurons_backend);
      ::Backend::CUDA::SpikingNeurons* input_neurons_backend =
        dynamic_cast<::Backend::CUDA::SpikingNeurons*>(input_neurons->backend());
      assert(input_neurons_backend);

      // First, get the current location in our circular buffer
      int bufferloc = (int)(round(current_time_in_seconds / timestep)) % buffersize;
      
        update_conductance_synapses_kernel<<<neurons_backend->number_of_neuron_blocks_per_grid, threads_per_block>>>(
	        postsynaptic_neuron_indices,
		biological_conductance_scaling_constants_lambda,
		synaptic_efficacies_or_weights,
                neurons_backend->per_neuron_efferent_synapse_count,
                neurons_backend->per_neuron_efferent_synapse_total,
                neurons_backend->per_neuron_efferent_synapse_indices,
                input_neurons_backend->per_neuron_efferent_synapse_count,
                input_neurons_backend->per_neuron_efferent_synapse_total,
                input_neurons_backend->per_neuron_efferent_synapse_indices,
                delays,
                spikes_travelling_to_synapse,
		time_of_last_spike_to_reach_synapse,
                neurons_backend->last_spike_time_of_each_neuron,
                input_neurons_backend->last_spike_time_of_each_neuron,
                current_time_in_seconds,
		circular_spikenum_buffer,
		spikeid_buffer,
		bufferloc,
		buffersize,
                timestep,
		input_neurons_backend->frontend()->total_number_of_neurons,
                (neurons_backend->frontend()->total_number_of_neurons + input_neurons_backend->frontend()->total_number_of_neurons),
     decay_term_values,
     reversal_values,
     num_decay_terms,
     synapse_decay_id,
     neuron_wise_conductance_trace,
     neurons_backend->current_injections,
     neurons_backend->membrane_potentials_v);
        CudaCheckError();
    }


    /* KERNELS BELOW */
    __global__ void update_conductance_synapses_kernel(
		int* postsynaptic_neuron_indices,
		float* d_biological_conductance_scaling_constants_lambda,
		float* d_synaptic_efficacies_or_weights,
		int* d_per_neuron_efferent_synapse_count,
        	int* d_per_neuron_efferent_synapse_total,
                int* d_per_neuron_efferent_synapse_indices,
		int* d_per_input_neuron_efferent_synapse_count,
        	int* d_per_input_neuron_efferent_synapse_total,
                int* d_per_input_neuron_efferent_synapse_indices,
                int* d_delays,
                int* d_spikes_travelling_to_synapse,
		float* d_time_of_last_spike_to_reach_synapse,
                float* d_last_spike_time_of_each_neuron,
                float* d_last_spike_time_of_each_input_neuron,
                float current_time_in_seconds,
		int* circular_spikenum_buffer,
		int* spikeid_buffer,
		int bufferloc,
		int buffersize,
                float timestep,
		int num_input_neurons,
                size_t total_number_of_neurons,
                  float* decay_term_values,
                  float* reversal_values,
                  int num_decay_terms,
                  int* synapse_decay_values,
                  float* neuron_wise_conductance_traces,
                  float* d_neurons_current_injections,
                  float * d_membrane_potentials_v){

      // Kernel to find activated neurons and to include synapse to list
      int indx = threadIdx.x + blockIdx.x * blockDim.x;
      while (indx < total_number_of_neurons){
   	int idx = indx - (num_input_neurons); 

	// For all neurons, check if spikes have occurred
        bool presynaptic_is_input = PRESYNAPTIC_IS_INPUT(idx);
	bool corr_idx = CORRECTED_PRESYNAPTIC_ID(idx, presynaptic_is_input);
        float effecttime = presynaptic_is_input ? d_last_spike_time_of_each_input_neuron[corr_idx] : d_last_spike_time_of_each_neuron[corr_idx];

        // Check if spike occurred within the last timestep    
        if (fabs(effecttime - current_time_in_seconds) < 0.5*timestep){
	  int synapse_count = presynaptic_is_input ? d_per_input_neuron_efferent_synapse_count[corr_idx] : d_per_neuron_efferent_synapse_count[corr_idx];  
          // For each of this neuron's efferent synapses
          for (int i = 0; i < synapse_count; i++){
            int synapse_id = presynaptic_is_input ? d_per_input_neuron_efferent_synapse_indices[d_per_input_neuron_efferent_synapse_total[corr_idx] - i - 1] : d_per_neuron_efferent_synapse_indices[d_per_neuron_efferent_synapse_total[corr_idx] - i - 1];
	    int targetloc = (bufferloc + d_delays[synapse_id]) % buffersize;
	    int pos = atomicAdd(&circular_spikenum_buffer[targetloc], 1);
	    spikeid_buffer[buffersize*targetloc + pos] = synapse_id;
          }
        }
        __syncthreads();
        indx += blockDim.x * gridDim.x;
      }

      // Kernel to decay conductance values and use new values to update current injections
      indx = threadIdx.x + blockIdx.x * blockDim.x;
      while (indx < (total_number_of_neurons - num_input_neurons)){
	int idx = indx;
	// For all non-input neurons, update the current injection
	float membrane_potential_v = d_membrane_potentials_v[idx];

	for (int decay_id = 0; decay_id < num_decay_terms; decay_id++){
	  if (decay_id == 0)
		  d_neurons_current_injections[idx] = 0.0f;
	  float synaptic_conductance_g = neuron_wise_conductance_traces[idx + decay_id*total_number_of_neurons];
	  // First decay the conductance values as required
	  synaptic_conductance_g *= expf(- timestep / decay_term_values[decay_id]);
	  neuron_wise_conductance_traces[idx + decay_id*total_number_of_neurons] = synaptic_conductance_g;
	  d_neurons_current_injections[idx] += synaptic_conductance_g * (reversal_values[decay_id] - membrane_potential_v);
	}
	indx += blockDim.x * gridDim.x;
      }

      
      // Kernel to create impact from awaiting synapses
      indx = threadIdx.x + blockIdx.x * blockDim.x;
      while (indx < circular_spikenum_buffer[bufferloc]){
	if (indx == 0){
		int last_address = bufferloc - 1;
		if (last_address < 0)
			last_address = buffersize - 1;
		circular_spikenum_buffer[last_address] = 0;
	}
        int synapse_id = spikeid_buffer[buffersize*bufferloc + indx];
	d_time_of_last_spike_to_reach_synapse[synapse_id] = current_time_in_seconds;
        int postsynaptic_neuron_id = postsynaptic_neuron_indices[synapse_id];
	int trace_id = synapse_decay_values[synapse_id];
	float synaptic_efficacy = d_biological_conductance_scaling_constants_lambda[synapse_id] * d_synaptic_efficacies_or_weights[synapse_id];
	atomicAdd(&neuron_wise_conductance_traces[total_number_of_neurons*trace_id + postsynaptic_neuron_id], synaptic_efficacy);
	indx += blockDim.x * gridDim.x;
      }
    }

  }
}
