#include "Spike/Backend/CUDA/Neurons/PoissonInputSpikingNeurons.hpp"

namespace Backend {
  namespace CUDA {
    PoissonInputSpikingNeurons::~PoissonInputSpikingNeurons() {
      CudaSafeCall(cudaFree(rates));
    }

    void PoissonInputSpikingNeurons::allocate_device_pointers(int maximum_axonal_delay_in_timesteps, bool high_fidelity_spike_storage) {
      InputSpikingNeurons::allocate_device_pointers(maximum_axonal_delay_in_timesteps, high_fidelity_spike_storage);
      CudaSafeCall(cudaMalloc((void **)&d_rates, sizeof(float)*total_number_of_neurons));
    }

    void PoissonInputSpikingNeurons::copy_constants_to_device() {
      InputSpikingNeurons::copy_constants_to_device();

      if (rates != NULL) {
        CudaSafeCall(cudaMemcpy(d_rates, rates, sizeof(float)*total_number_of_neurons, cudaMemcpyHostToDevice));
      }
    }

    void PoissonInputSpikingNeurons::update_membrane_potentials(float timestep, float current_time_in_seconds) {

      poisson_update_membrane_potentials_kernel<<<random_state_manager->block_dimensions, random_state_manager->threads_per_block>>>
        (random_state_manager->d_states,
         d_rates,
         d_membrane_potentials_v,
         timestep,
         d_thresholds_for_action_potential_spikes,
         total_number_of_neurons,
         current_stimulus_index);

	CudaCheckError();
    }

    __global__ void poisson_update_membrane_potentials_kernel(curandState_t* d_states,
                                                              float *d_rates,
                                                              float *d_membrane_potentials_v,
                                                              float timestep,
                                                              float * d_thresholds_for_action_potential_spikes,
                                                              size_t total_number_of_input_neurons,
                                                              int current_stimulus_index) {

	 
      int t_idx = threadIdx.x + blockIdx.x * blockDim.x;
      int idx = t_idx;
      while (idx < total_number_of_input_neurons){

        int rate_index = (total_number_of_input_neurons * current_stimulus_index) + idx;

        float rate = d_rates[rate_index];

        if (rate > 0.1) {
          // Creates random float between 0 and 1 from uniform distribution
          // d_states effectively provides a different seed for each thread
          // curand_uniform produces different float every time you call it
          float random_float = curand_uniform(&d_states[t_idx]);
			
          // if the randomnumber is less than the rate
          if (random_float < (rate * timestep)) {
            // Puts membrane potential above default spiking threshold
            d_membrane_potentials_v[idx] = d_thresholds_for_action_potential_spikes[idx] + 0.02;
          } 
        }

        idx += blockDim.x * gridDim.x;

      }
      __syncthreads();
    }
  }
}
