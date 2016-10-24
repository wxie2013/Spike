#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm> // For random shuffle
#include <time.h>
#include <string>
#include <sys/stat.h>

#include "Simulator.h"

#include "../Helpers/CUDAErrorCheckHelpers.h"
#include "../Helpers/TerminalHelpers.h"
#include "../Helpers/TimerWithMessages.h"
#include "../Helpers/RandomStateManager.h"

using namespace std;

//string full_directory_name_for_simulation_data_files ("output/");


// Constructor
Simulator::Simulator(){

	full_directory_name_for_simulation_data_files = "output/";

	count_neuron_spikes_recording_electrodes = NULL;
	input_count_neuron_spikes_recording_electrodes = NULL;
	collect_neuron_spikes_recording_electrodes = NULL;
	input_collect_neuron_spikes_recording_electrodes = NULL;
	network_state_archive_recording_electrodes = NULL;

	high_fidelity_spike_storage = false;


	number_of_simulations_run = 0;

}


// Destructor
Simulator::~Simulator(){

}

void Simulator::CreateDirectoryForSimulationDataFiles(string directory_name_for_simulation_data_files) {
	if (mkdir(("output/"+directory_name_for_simulation_data_files).c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH)==0)
		printf("\nDirectory created\n");
	else
		print_message_and_exit("\nERROR: You must set a different experiment name to avoid overwriting the results\n");
	full_directory_name_for_simulation_data_files = "output/"+directory_name_for_simulation_data_files+"/";
}


void Simulator::SetSpikingModel(SpikingModel * spiking_model_parameter) {
	spiking_model = spiking_model_parameter;
}


void Simulator::prepare_recording_electrodes(Simulator_Recording_Electrodes_To_Use_Struct * recording_electrodes_to_use_struct_param, Collect_Neuron_Spikes_Optional_Parameters * collect_neuron_spikes_optional_parameters) {


	TimerWithMessages * timer = new TimerWithMessages("Setting up recording electrodes...");

	recording_electrodes_to_use_struct = recording_electrodes_to_use_struct_param;


	if (recording_electrodes_to_use_struct->count_neuron_spikes_recording_electrodes_bool) {
		count_neuron_spikes_recording_electrodes = new CountNeuronSpikesRecordingElectrodes(spiking_model->spiking_neurons, spiking_model->spiking_synapses, full_directory_name_for_simulation_data_files, "Neurons");
		count_neuron_spikes_recording_electrodes->initialise_count_neuron_spikes_recording_electrodes();
	}

	if (recording_electrodes_to_use_struct->input_count_neuron_spikes_recording_electrodes_bool) {
		input_count_neuron_spikes_recording_electrodes = new CountNeuronSpikesRecordingElectrodes(spiking_model->input_spiking_neurons, spiking_model->spiking_synapses, full_directory_name_for_simulation_data_files, "Input_Neurons");
		input_count_neuron_spikes_recording_electrodes->initialise_count_neuron_spikes_recording_electrodes();
	}

	if (recording_electrodes_to_use_struct->collect_neuron_spikes_recording_electrodes_bool) {
		collect_neuron_spikes_recording_electrodes = new CollectNeuronSpikesRecordingElectrodes(spiking_model->spiking_neurons, spiking_model->spiking_synapses, full_directory_name_for_simulation_data_files, "Neurons");
		collect_neuron_spikes_recording_electrodes->initialise_collect_neuron_spikes_recording_electrodes(collect_neuron_spikes_optional_parameters);
	}

	if (recording_electrodes_to_use_struct->input_collect_neuron_spikes_recording_electrodes_bool) {
		input_collect_neuron_spikes_recording_electrodes = new CollectNeuronSpikesRecordingElectrodes(spiking_model->input_spiking_neurons, spiking_model->spiking_synapses, full_directory_name_for_simulation_data_files, "Input_Neurons");
		input_collect_neuron_spikes_recording_electrodes->initialise_collect_neuron_spikes_recording_electrodes(collect_neuron_spikes_optional_parameters);
	}

	if (recording_electrodes_to_use_struct->network_state_archive_recording_electrodes_bool) {
		network_state_archive_recording_electrodes = new NetworkStateArchiveRecordingElectrodes(spiking_model->spiking_neurons, spiking_model->spiking_synapses, full_directory_name_for_simulation_data_files, "Synapses");
	}


	timer->stop_timer_and_log_time_and_message("Recording electrodes setup.", true);


}

void Simulator::reset_all_recording_electrodes() {

	if (recording_electrodes_to_use_struct->count_neuron_spikes_recording_electrodes_bool) {
		count_neuron_spikes_recording_electrodes->reset_pointers_for_spike_count();
	}

	if (recording_electrodes_to_use_struct->input_count_neuron_spikes_recording_electrodes_bool) {
		input_count_neuron_spikes_recording_electrodes->reset_pointers_for_spike_count();
	}

	if (recording_electrodes_to_use_struct->collect_neuron_spikes_recording_electrodes_bool) {
		collect_neuron_spikes_recording_electrodes->delete_and_reset_collected_spikes();
	}

	if (recording_electrodes_to_use_struct->input_collect_neuron_spikes_recording_electrodes_bool) {
		input_collect_neuron_spikes_recording_electrodes->delete_and_reset_collected_spikes();
	}

}



// void Simulator::RunSimulationToCountNeuronSpikes(float presentation_time_per_stimulus_per_epoch, bool collect_spikes, bool save_collected_spikes_and_states_to_file, SpikeAnalyser *spike_analyser, bool human_readable_storage, bool isTrained) {
// 	bool number_of_epochs = 1;
// 	bool apply_stdp_to_relevant_synapses = false;
// 	bool count_spikes_per_neuron = true;
// 	int stimulus_presentation_order_seed = 0; // Shouldn't be needed if stimuli presentation not random
// 	Stimuli_Presentation_Struct * stimuli_presentation_params = new Stimuli_Presentation_Struct();
// 	// stimuli_presentation_params->presentation_format = PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_OBJECTS;
// 	stimuli_presentation_params->presentation_format = PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_STIMULI;
// 	stimuli_presentation_params->object_order = OBJECT_ORDER_ORIGINAL;
// 	stimuli_presentation_params->transform_order = TRANSFORM_ORDER_ORIGINAL;
	
// 	if (!isTrained)
// 		recording_electrodes->write_initial_synaptic_weights_to_file(spiking_model->spiking_synapses, human_readable_storage);
	
// 	RunSimulation(presentation_time_per_stimulus_per_epoch, number_of_epochs, collect_spikes, save_collected_spikes_and_states_to_file, apply_stdp_to_relevant_synapses, count_spikes_per_neuron, stimuli_presentation_params, stimulus_presentation_order_seed, spike_analyser,human_readable_storage,isTrained);
	
// 	if (isTrained)
// 		recording_electrodes->write_network_state_to_file(spiking_model->spiking_synapses, human_readable_storage);

// }


// void Simulator::RunSimulationToCollectEvents(float presentation_time_per_stimulus_per_epoch, bool isTrained) {
// 	bool number_of_epochs = 1;
// 	bool apply_stdp_to_relevant_synapses = false;
// 	bool count_spikes_per_neuron = true;
// 	int stimulus_presentation_order_seed = 0; // Shouldn't be needed if stimuli presentation not random
// 	Stimuli_Presentation_Struct * stimuli_presentation_params = new Stimuli_Presentation_Struct();
// 	// stimuli_presentation_params->presentation_format = PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_OBJECTS;
// 	stimuli_presentation_params->presentation_format = PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_STIMULI;
// 	stimuli_presentation_params->object_order = OBJECT_ORDER_ORIGINAL;
// 	stimuli_presentation_params->transform_order = TRANSFORM_ORDER_ORIGINAL;

// 	bool collect_spikes = false;
// 	bool save_collected_spikes_and_states_to_file = false;

// 	SpikeAnalyser * spike_analyser = NULL;
// 	bool human_readable_storage = false;
	
// 	RunSimulation(presentation_time_per_stimulus_per_epoch, number_of_epochs, collect_spikes, save_collected_spikes_and_states_to_file, apply_stdp_to_relevant_synapses, count_spikes_per_neuron, stimuli_presentation_params, stimulus_presentation_order_seed, spike_analyser, human_readable_storage,isTrained);
	
// }

// void Simulator::RunSimulationToTrainNetwork(float presentation_time_per_stimulus_per_epoch, int number_of_epochs, Stimuli_Presentation_Struct * stimuli_presentation_params, int stimulus_presentation_order_seed) {

// 	bool apply_stdp_to_relevant_synapses = true;
// 	bool count_spikes_per_neuron = false;
// 	bool collect_spikes = false;
// 	bool save_collected_spikes_and_states_to_file = false;

// 	RunSimulation(presentation_time_per_stimulus_per_epoch, number_of_epochs, collect_spikes, save_collected_spikes_and_states_to_file, apply_stdp_to_relevant_synapses, count_spikes_per_neuron, stimuli_presentation_params, stimulus_presentation_order_seed, NULL, false, false);
// }




void Simulator::RunSimulation(Simulator_File_Storage_Options_Struct simulator_file_storage_options_struct, float presentation_time_per_stimulus_per_epoch, int number_of_epochs, bool save_collected_spikes_and_states_to_file, bool apply_stdp_to_relevant_synapses, Stimuli_Presentation_Struct * stimuli_presentation_params, int stimulus_presentation_order_seed, SpikeAnalyser *spike_analyser, bool human_readable_storage, bool isTrained) {

	check_for_epochs_and_begin_simulation_message(spiking_model->timestep, spiking_model->input_spiking_neurons->total_number_of_input_stimuli, number_of_epochs, collect_spikes, save_collected_spikes_and_states_to_file, spiking_model->spiking_neurons->total_number_of_neurons, spiking_model->input_spiking_neurons->total_number_of_neurons, spiking_model->spiking_synapses->total_number_of_synapses);
	// Should print something about stimuli_presentation_params as old stuff removed from check_for_epochs...
	TimerWithMessages * simulation_timer = new TimerWithMessages();

	// Set seed for stimulus presentation order
	srand(stimulus_presentation_order_seed);


	if (number_of_simulations_run > 0) reset_all_recording_electrodes();

	if (simulator_file_storage_options_struct->write_initial_synaptic_weights_to_file_bool) {
	
		network_state_archive_recording_electrodes->write_initial_synaptic_weights_to_file(spiking_model->spiking_synapses, human_readable_storage);
	
	}

	for (int epoch_number = 0; epoch_number < number_of_epochs; epoch_number++) {

		TimerWithMessages * epoch_timer = new TimerWithMessages();
		printf("Starting Epoch: %d\n", epoch_number);

		spiking_model->reset_model_activities();

		float current_time_in_seconds = 0.0f;

		int* stimuli_presentation_order = spiking_model->input_spiking_neurons->setup_stimuli_presentation_order(stimuli_presentation_params);
		for (int stimulus_index = 0; stimulus_index < spiking_model->input_spiking_neurons->total_number_of_input_stimuli; stimulus_index++) {

			perform_pre_stimulus_presentation_instructions();

			int number_of_timesteps_per_stimulus_per_epoch = presentation_time_per_stimulus_per_epoch / spiking_model->timestep;
		

			for (int timestep_index = 0; timestep_index < number_of_timesteps_per_stimulus_per_epoch; timestep_index++){
				
				spiking_model->perform_per_timestep_model_instructions(current_time_in_seconds, apply_stdp_to_relevant_synapses);

				perform_per_timestep_recording_electrode_instructions();

				current_time_in_seconds += float(spiking_model->timestep);

			}

			perform_post_stimulus_presentation_instructions();
			
		}

		perform_post_epoch_instructions();
		
	}
	
	perform_end_of_simulation_instructions();
	
}



void Simulator::perform_per_timestep_recording_electrode_instructions(float current_time_in_seconds, int timestep_index, int number_of_timesteps_per_stimulus_per_epoch) {

	// Do various recording electrode operations


	// JI PSEUDO CODE FOR COLLECTING EVENTS START

	// if synapse_spike_arrived

	// 	int neuron_spike_count_for_postsynaptic_neuron = recording_electrodes->d_per_neuron_spike_counts[postsynaptic_neuron_id];

	// 	if (current_time_in_seconds > (d_ordered_spike_times_data[d_beginning_spike_time_int_indices_per_neuron[neuron_index] + neuron_spike_count_for_postsynaptic_neuron] - window) {
	// 		d_events_as_bools_per_neuron_and_spike_data[d_beginning_event_bool_indices_per_neuron + neuron_spike_count_for_postsynaptic_neuron * neurons->per_neuron_afferent_synapse_count[postsynaptic_neuron_index] + synapses->d_synapse_postsynaptic_neuron_count_index[synapse_index]] = true;
	// }


	// JI PSEUDO CODE FOR COLLECTING EVENTS END



	if (count_neuron_spikes_recording_electrodes_to_use_struct->count_neuron_spikes_recording_electrodes_bool) {
	
		count_neuron_spikes_recording_electrodes->add_spikes_to_per_neuron_spike_count(current_time_in_seconds);
	
	}

	if (count_input_neuron_spikes_recording_electrodes_to_use_struct->count_neuron_spikes_recording_electrodes_bool) {
	
		count_input_neuron_spikes_recording_electrodes->add_spikes_to_per_neuron_spike_count(current_time_in_seconds);
	
	}

	if (recording_electrodes_to_use_struct->collect_neuron_spikes_recording_electrodes_bool){

		collect_neuron_spikes_recording_electrodes->collect_spikes_for_timestep(current_time_in_seconds);
		collect_neuron_spikes_recording_electrodes->copy_spikes_from_device_to_host_and_reset_device_spikes_if_device_spike_count_above_threshold(current_time_in_seconds, timestep_index, number_of_timesteps_per_stimulus_per_epoch );
	
	}

	if (recording_electrodes_to_use_struct->collect_input_neuron_spikes_recording_electrodes_bool) {
		
		collect_input_neuron_spikes_recording_electrodes->collect_spikes_for_timestep(current_time_in_seconds);
		collect_input_neuron_spikes_recording_electrodes->copy_spikes_from_device_to_host_and_reset_device_spikes_if_device_spike_count_above_threshold(current_time_in_seconds, timestep_index, number_of_timesteps_per_stimulus_per_epoch );

	}

}


void Simulator::perform_pre_stimulus_presentation_instructions() {

	printf("Stimulus: %d\n", stimuli_presentation_order[stimulus_index]);
	printf("stimuli_presentation_params->presentation_format: %d\n", stimuli_presentation_params->presentation_format);

	switch (stimuli_presentation_params->presentation_format) {
		case PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_STIMULI: case PRESENTATION_FORMAT_RANDOM_RESET_BETWEEN_EACH_STIMULUS:
		{
			spiking_model->reset_model_activities();

			// Reset time (Useful for Generator Neurons Specifically)
			current_time_in_seconds = 0.0f;

			break;
		}
		case PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_OBJECTS:
		{
			bool stimulus_is_new_object = spiking_model->input_spiking_neurons->stimulus_is_new_object_for_object_by_object_presentation(stimulus_index);
			(stimulus_is_new_object) ? printf("stimulus_is_new_object\n") : printf("stimulus_is_NOT_new_object\n");

			if (stimulus_is_new_object) {
				spiking_model->reset_model_activities();
			}

			break;
		}
		default:
			break;

	}


	spiking_model->input_spiking_neurons->current_stimulus_index = stimuli_presentation_order[stimulus_index];
	spiking_model->input_spiking_neurons->reset_neuron_activities();


}



void Simulator::perform_post_stimulus_presentation_instructions() {

	if (count_neuron_spikes_recording_electrodes_to_use_struct->count_neuron_spikes_recording_electrodes_bool && spike_analyser) {

		spike_analyser->store_spike_counts_for_stimulus_index(spiking_model->input_spiking_neurons->current_stimulus_index, recording_electrodes->d_per_neuron_spike_counts);
		count_neuron_spikes_recording_electrodes->reset_pointers_for_spike_count();

	}

}


void Simulator::perform_post_epoch_instructions() {

	printf("Epoch %d, Complete.\n", epoch_number);
	epoch_timer->stop_timer_and_log_time_and_message(" ", true);
	
	if (recording_electrodes_to_use_struct->collect_neuron_spikes_recording_electrodes_bool) printf(" Number of Spikes: %d\n", collect_neuron_spikes_recording_electrodes->h_total_number_of_spikes_stored_on_host);
	if (recording_electrodes_to_use_struct->collect_input_neuron_spikes_recording_electrodes_bool) printf(" Number of Input Spikes: %d\n", collect_input_neuron_spikes_recording_electrodes->h_total_number_of_spikes_stored_on_host);

	if (simulator_file_storage_options_struct->save_recorded_neuron_spikes_to_file) collect_neuron_spikes_recording_electrodes->write_spikes_to_file(epoch_number, human_readable_storage, isTrained);
	if (simulator_file_storage_options_struct->save_recorded_input_neuron_spikes_to_file) collect_input_neuron_spikes_recording_electrodes->write_spikes_to_file(epoch_number, human_readable_storage,isTrained);

}


void Simulator::perform_end_of_simulation_instructions() {

	simulation_timer->stop_timer_and_log_time_and_message("Simulation Complete!", true);

	if (recording_electrodes_to_use_struct->network_state_archive_recording_electrodes_bool) {

		network_state_archive_recording_electrodes->write_network_state_to_file(spiking_model->spiking_synapses, human_readable_storage);

	}

	number_of_simulations_run++;

}