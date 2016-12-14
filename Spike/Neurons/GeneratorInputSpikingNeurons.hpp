#ifndef GeneratorInputSpikingNeurons_H
#define GeneratorInputSpikingNeurons_H

#include "InputSpikingNeurons.hpp"

struct generator_input_spiking_neuron_parameters_struct : input_spiking_neuron_parameters_struct {
	generator_input_spiking_neuron_parameters_struct() { input_spiking_neuron_parameters_struct(); }
};

class GeneratorInputSpikingNeurons; // forward definition

namespace Backend {
  class GeneratorInputSpikingNeurons : public virtual InputSpikingNeurons {
  public:
    ADD_FRONTEND_GETTER(GeneratorInputSpikingNeurons);
  };
} 

#include "Spike/Backend/Dummy/Neurons/GeneratorInputSpikingNeurons.hpp"
#ifdef SPIKE_WITH_CUDA
#include "Spike/Backend/CUDA/Neurons/GeneratorInputSpikingNeurons.hpp"
#endif

class GeneratorInputSpikingNeurons : public InputSpikingNeurons {
public:
  // Constructor/Destructor
  GeneratorInputSpikingNeurons();
  ~GeneratorInputSpikingNeurons();

  ADD_BACKEND_GETSET(GeneratorInputSpikingNeurons, InputSpikingNeurons);
  
  // Variables
  int length_of_longest_stimulus;

  // Host Pointers
  int* number_of_spikes_in_stimuli = nullptr;
  int** neuron_id_matrix_for_stimuli = nullptr;
  float** spike_times_matrix_for_stimuli = nullptr;

  void update_membrane_potentials(float timestep, float current_time_in_seconds) override;

  void AddStimulus(int spikenumber, int* ids, float* spiketimes);

private:
  ::Backend::GeneratorInputSpikingNeurons* _backend = nullptr;
};

#endif
