#include "PoissonInputSpikingNeurons.hpp"
#include "Spike/Helpers/TerminalHelpers.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <algorithm> // For random shuffle
#include <cassert>

using namespace std;

PoissonInputSpikingNeurons::~PoissonInputSpikingNeurons() {
  free(random_state_manager);
  free(rates);
}


int PoissonInputSpikingNeurons::AddGroup(neuron_parameters_struct * group_params){
  int new_group_id = InputSpikingNeurons::AddGroup(group_params);
  poisson_input_spiking_neuron_parameters_struct * poisson_input_spiking_group_params = (poisson_input_spiking_neuron_parameters_struct*)group_params;
  rate = poisson_input_spiking_group_params->rate;
  set_up_rates();
  return new_group_id;
}


void PoissonInputSpikingNeurons::set_up_rates() {
  rates = (float*)realloc(rates, sizeof(float)*total_number_of_neurons);
  for (int i = 0; i < total_number_of_neurons; i++) {
    rates[i] = rate;
  }
}


void PoissonInputSpikingNeurons::init_random_state(bool force) {
  assert(backend() && "Backend needs to have been prepared before calling this!");
  if (force || !random_state_manager) {
    random_state_manager = new RandomStateManager();
    random_state_manager->init_backend(backend()->context);
  }
}

void PoissonInputSpikingNeurons::prepare_backend_early() {
  set_up_rates();
  InputSpikingNeurons::prepare_backend_early();
  init_random_state();
}
  
void PoissonInputSpikingNeurons::select_stimulus(int stimulus_index){
  InputSpikingNeurons::select_stimulus(stimulus_index);
  if (_backend) reset_state();
}

void PoissonInputSpikingNeurons::state_update(float current_time_in_seconds, float timestep) {
  backend()->state_update(current_time_in_seconds, timestep);
}

SPIKE_MAKE_INIT_BACKEND(PoissonInputSpikingNeurons);
