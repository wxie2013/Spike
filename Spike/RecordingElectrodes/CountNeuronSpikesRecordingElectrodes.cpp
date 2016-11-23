#include "CountNeuronSpikesRecordingElectrodes.hpp"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "../Helpers/TerminalHelpers.hpp"
#include <string>
#include <time.h>
using namespace std;

// CountNeuronSpikesRecordingElectrodes Constructor
CountNeuronSpikesRecordingElectrodes::CountNeuronSpikesRecordingElectrodes(SpikingNeurons * neurons_parameter, SpikingSynapses * synapses_parameter, string full_directory_name_for_simulation_data_files_param, const char * prefix_string_param) 
		: RecordingElectrodes(neurons_parameter, synapses_parameter, full_directory_name_for_simulation_data_files_param, prefix_string_param) { }

CountNeuronSpikesRecordingElectrodes::~CountNeuronSpikesRecordingElectrodes() {
  // TODO
}

void CountNeuronSpikesRecordingElectrodes::initialise_count_neuron_spikes_recording_electrodes() {
  // allocate_pointers_for_spike_count();
  reset_state();
  // prepare_backend();
  backend()->prepare();
}

void CountNeuronSpikesRecordingElectrodes::reset_state() {
  backend()->reset_state();
  // NB: RecordingElectrodes::reset_state is pure virtual at the moment ::
  // RecordingElectrodes::reset_state();
}

void CountNeuronSpikesRecordingElectrodes::add_spikes_to_per_neuron_spike_count(float current_time_in_seconds) {
  backend()->add_spikes_to_per_neuron_spike_count(this, current_time_in_seconds);
}

MAKE_PREPARE_BACKEND(CountNeuronSpikesRecordingElectrodes);