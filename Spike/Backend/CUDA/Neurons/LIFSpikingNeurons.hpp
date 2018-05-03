#pragma once

#include "Spike/Neurons/LIFSpikingNeurons.hpp"
#include "Spike/Neurons/SpikingNeurons.hpp"
#include "Spike/Backend/CUDA/CUDABackend.hpp"
#include "SpikingNeurons.hpp"

#include <cuda.h>
#include <vector_types.h>

namespace Backend {
  namespace CUDA {
    struct lif_spiking_neurons_data_struct: spiking_neurons_data_struct {
	float* membrane_time_constants_tau_m;
	float* membrane_resistances_R;

    };

    class LIFSpikingNeurons : public virtual ::Backend::CUDA::SpikingNeurons,
                              public virtual ::Backend::LIFSpikingNeurons {
    public:
      float * membrane_time_constants_tau_m = nullptr;
      float * membrane_resistances_R = nullptr;

      ~LIFSpikingNeurons() override;
      SPIKE_MAKE_BACKEND_CONSTRUCTOR(LIFSpikingNeurons);
      using ::Backend::LIFSpikingNeurons::frontend;

      lif_spiking_neurons_data_struct* neuron_data;

      void prepare() override;
      void reset_state() override;

      void copy_constants_to_device(); // Not virtual
      void allocate_device_pointers(); // Not virtual

      void state_update(float current_time_in_seconds, float timestep) override;
    };

    __global__ void lif_update_membrane_potentials(float *d_membrane_potentials_v,
                                                   float * d_last_spike_time_of_each_neuron,
                                                   float * d_membrane_resistances_R,
                                                   float * d_membrane_time_constants_tau_m,
                                                   float * d_resting_potentials,
                                                   float* d_current_injections,
						   float* d_total_current_conductance,
						   float* d_thresholds_for_action_potential_spikes,
                                                   float background_current,
                                                   float timestep,
						   int timestep_grouping,
                                                   float current_time_in_seconds,
                                                   float refactory_period_in_seconds,
                                                   size_t total_number_of_neurons);
  }
}
