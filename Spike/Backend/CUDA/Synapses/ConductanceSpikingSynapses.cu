// -*- mode: c++ -*-
#include "Spike/Backend/CUDA/Synapses/ConductanceSpikingSynapses.hpp"

SPIKE_EXPORT_BACKEND_TYPE(CUDA, ConductanceSpikingSynapses);

namespace Backend {
  namespace CUDA {
    __device__ injection_kernel conductance_device_kernel = conductance_spiking_current_injection_kernel;

    // ConductanceSpikingSynapses Destructor
    ConductanceSpikingSynapses::~ConductanceSpikingSynapses() {
      CudaSafeCall(cudaFree(neuron_wise_conductance_trace));
      CudaSafeCall(cudaFree(d_synaptic_data));
      CudaSafeCall(cudaFree(d_decay_factors_g));
      CudaSafeCall(cudaFree(d_reversal_potentials_Vhat));
      free(h_neuron_wise_conductance_trace);
    }

    void ConductanceSpikingSynapses::prepare() {
      SpikingSynapses::prepare();

      // Set up per neuron conductances
      conductance_trace_length = frontend()->neuron_pop_size*frontend()->num_syn_labels;
      h_neuron_wise_conductance_trace = (float*)realloc(h_neuron_wise_conductance_trace, conductance_trace_length*sizeof(float));
      for (int id = 0; id < conductance_trace_length; id++)
        h_neuron_wise_conductance_trace[id] = 0.0f;

      // Carry out remaining device actions
      allocate_device_pointers();
      copy_constants_and_initial_efficacies_to_device();

      conductance_spiking_synapses_data_struct temp_synaptic_data;
      memcpy(&temp_synaptic_data, synaptic_data, sizeof(spiking_synapses_data_struct));
      free(synaptic_data);
      synaptic_data = new conductance_spiking_synapses_data_struct();
      memcpy(synaptic_data, &temp_synaptic_data, sizeof(spiking_synapses_data_struct));
      conductance_spiking_synapses_data_struct* this_synaptic_data = static_cast<conductance_spiking_synapses_data_struct*>(synaptic_data); 
      this_synaptic_data->decay_factors_g = d_decay_factors_g;
      this_synaptic_data->reversal_potentials_Vhat = d_reversal_potentials_Vhat;
      this_synaptic_data->neuron_wise_conductance_trace = neuron_wise_conductance_trace;
      this_synaptic_data->synapse_type = CONDUCTANCE;
      CudaSafeCall(cudaMemcpy(
        d_synaptic_data,
        synaptic_data,
        sizeof(conductance_spiking_synapses_data_struct), cudaMemcpyHostToDevice));

    }

    void ConductanceSpikingSynapses::reset_state() {
      SpikingSynapses::reset_state();
      CudaSafeCall(cudaMemcpy(
        neuron_wise_conductance_trace,
        h_neuron_wise_conductance_trace,
        sizeof(float)*conductance_trace_length, cudaMemcpyHostToDevice));

    }


    void ConductanceSpikingSynapses::allocate_device_pointers() {
      CudaSafeCall(cudaMalloc((void **)&neuron_wise_conductance_trace, sizeof(float)*conductance_trace_length));
      CudaSafeCall(cudaMalloc((void **)&d_decay_factors_g, sizeof(float)*frontend()->num_syn_labels));
      CudaSafeCall(cudaMalloc((void **)&d_reversal_potentials_Vhat, sizeof(float)*frontend()->num_syn_labels));
      CudaSafeCall(cudaFree(d_synaptic_data));
      CudaSafeCall(cudaMalloc((void **)&d_synaptic_data, sizeof(conductance_spiking_synapses_data_struct)));
      CudaSafeCall(cudaMemcpyFromSymbol(
            &host_injection_kernel,
            conductance_device_kernel,
            sizeof(injection_kernel)));
    }

    void ConductanceSpikingSynapses::copy_constants_and_initial_efficacies_to_device() {
      CudaSafeCall(cudaMemcpy(
        neuron_wise_conductance_trace,
        h_neuron_wise_conductance_trace,
        sizeof(float)*conductance_trace_length, cudaMemcpyHostToDevice));
      vector<float> decay_vals_g;
      for (int syn_label_indx = 0; syn_label_indx < frontend()->num_syn_labels; syn_label_indx++)
        decay_vals_g.push_back((expf(-frontend()->model->timestep / frontend()->decay_terms_tau_g[syn_label_indx])));
      CudaSafeCall(cudaMemcpy(
        d_decay_factors_g,
        decay_vals_g.data(),
        sizeof(float)*frontend()->num_syn_labels, cudaMemcpyHostToDevice));
      CudaSafeCall(cudaMemcpy(
        d_reversal_potentials_Vhat,
        &(frontend()->reversal_potentials_Vhat[0]),
        sizeof(float)*frontend()->num_syn_labels, cudaMemcpyHostToDevice));
    }



    /* STATE UPDATE */
    void ConductanceSpikingSynapses::state_update
    (::SpikingNeurons* neurons,
     ::SpikingNeurons* input_neurons,
     float current_time_in_seconds, float timestep) {
      SpikingSynapses::state_update(neurons, input_neurons, current_time_in_seconds, timestep);
    }


    /* KERNELS BELOW */
    __device__ float conductance_spiking_current_injection_kernel(
        spiking_synapses_data_struct* in_synaptic_data,
        spiking_neurons_data_struct* neuron_data,
        float multiplication_to_volts,
        float current_membrane_voltage,
        float current_time_in_seconds,
        float timestep,
        int idx,
        int g){
      
      conductance_spiking_synapses_data_struct* synaptic_data = (conductance_spiking_synapses_data_struct*) in_synaptic_data;
        
      int bufferloc = (((int)roundf(current_time_in_seconds / timestep) + g) % synaptic_data->neuron_inputs.temporal_buffersize)*synaptic_data->neuron_inputs.input_buffersize;
        
      float total_current = 0.0f;
      for (int syn_label = 0; syn_label < synaptic_data->num_syn_labels; syn_label++){
        float decay_factor = synaptic_data->decay_factors_g[syn_label];
        float reversal_value = synaptic_data->reversal_potentials_Vhat[syn_label];
        float synaptic_conductance_g = synaptic_data->neuron_wise_conductance_trace[syn_label + idx*synaptic_data->num_syn_labels];
        // Update the synaptic conductance
        synaptic_conductance_g = decay_factor*synaptic_conductance_g + synaptic_data->neuron_inputs.circular_input_buffer[bufferloc + syn_label + idx*synaptic_data->num_syn_labels];
        // Reset the conductance update
        synaptic_data->neuron_inputs.circular_input_buffer[bufferloc + syn_label + idx*synaptic_data->num_syn_labels] = 0.0f;
        total_current += synaptic_conductance_g*(reversal_value - current_membrane_voltage);
    
        synaptic_data->neuron_wise_conductance_trace[syn_label + idx*synaptic_data->num_syn_labels] = synaptic_conductance_g;
    
      }
      return total_current*multiplication_to_volts;
    };

  }
}
