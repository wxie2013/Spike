#include "aki_model.h"

/*
 *  Main function in which the network is created and run
 */

aki_model::aki_model()
{
    load_run_config_parameters();

    set_model_parameters();

    define_spiking_model();
}

aki_model::~aki_model()
{
    delete model;
    delete input_neurons;
    delete lif_spiking_neurons;
    delete conductance_spiking_synapses;
    delete spike_monitor;
    delete input_spike_monitor;
    delete adding_synapses_timer;
    delete G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
    delete adding_input_neurons_timer;
    delete adding_neurons_timer;
    delete EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS;
    delete INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS;
    delete STDP_PARAMS;
    delete evans_stdp;
}

//__ run the model ...
void aki_model::run_spiking_model(bool binary_output_only=true)
{
    // if it's already trained, load any weights before finalising the model
    if (starting_time != 0){
        load_weights(model, current_weight);
    }
    model->run(simtime, plasticity_on);

    //use binary mode. The text mode is too slow ..
    if(!binary_output_only) {
        spike_monitor->save_spikes_as_txt(neuron_dir);
        input_spike_monitor->save_spikes_as_txt(input_dir);
        model->spiking_synapses->save_connectivity_as_txt(synapse_dir);
    }

    spike_monitor->save_spikes_as_binary(neuron_dir);
    input_spike_monitor->save_spikes_as_binary(input_dir);
    model->spiking_synapses->save_connectivity_as_binary(synapse_dir);
}

//
// Defining the Spiking Model
void aki_model::define_spiking_model()
{
    // Create the SpikingModel
    model = new SpikingModel();

    // Set up the simulator with a timestep at which the neuron, synapse and STDP properties will be calculated
    model->SetTimestep(timestep);

    // Choose an input neuron type
    input_neurons = new ImagePoissonInputSpikingNeurons();

    // Choose your neuron type
    lif_spiking_neurons = new LIFSpikingNeurons();

    // Choose your synapse type
    conductance_spiking_synapses = new ConductanceSpikingSynapses();

    // Allocate your chosen components to the simulator
    model->input_spiking_neurons = input_neurons;
    model->spiking_neurons = lif_spiking_neurons;
    model->spiking_synapses = conductance_spiking_synapses;

    // setup STDP. According to the https://sites.google.com/view/spike-simulator/home/simple-example, add plastiscity rule need to be here
    setup_STDP();

    // ADD ANY ACTIVITY MONITORS OR PLASTICITY RULES YOU WISH FOR
    spike_monitor = new SpikingActivityMonitor(lif_spiking_neurons);
    input_spike_monitor = new SpikingActivityMonitor(input_neurons);
    model->AddActivityMonitor(spike_monitor);
    model->AddActivityMonitor(input_spike_monitor);

    // setup neuron groups
    setup_neuron_groups();

    // SETTING UP SYNAPSES
    setup_synapses();
}

//__
void aki_model::setup_synapses()
{
    // Creating a synapses parameter structure for connections from the input neurons to the excitatory neurons
    adding_synapses_timer = new TimerWithMessages("Adding Synapses...\n");

    define_synapses_parameters();
    make_synapses_connections();

    adding_synapses_timer->stop_timer_and_log_time_and_message("Synapses Added.", true);
}

//__
void aki_model::make_synapses_connections()
{
    for (int l=0; l<number_of_layers; l++){
        if(l==0)
            model->AddSynapseGroupsForNeuronGroupAndEachInputGroup(EXCITATORY_NEURONS[l], G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
        else{
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FF[l-1];
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[l-1];
            //.. E2E_FF has 2 connection per pre->post. One can directly set max_number_of_connections_per_pair = 2 in parameter setting, through there are some differences.  See Nashir's answer
           
            for (int connection_number = 0; connection_number < max_number_of_connections_per_pair; connection_number++){
                model->AddSynapseGroup(EXCITATORY_NEURONS[l-1], EXCITATORY_NEURONS[l], E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
            }

            if(E2E_FB_ON)  // FF has 2 connect pre-post. why FB does not have ???
                model->AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l-1], E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
        }

        E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_E2I_L[l];
        model->AddSynapseGroup(EXCITATORY_NEURONS[l], INHIBITORY_NEURONS[l], E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);

        if (inh_layer_on[l]){
            I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_I2E_L[l];
            model->AddSynapseGroup(INHIBITORY_NEURONS[l], EXCITATORY_NEURONS[l], I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
        }

        if(E2E_L_ON)
            model->AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l], E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
    }
}
//__ read synaptic information from existing binary output ___
void aki_model::read_synaptic_data()
{
    float weight;
    int delay, pre_ID, post_ID;

    string SynapticWeights_file = existing_synapse_dir+"SynapticWeights.bin";
    string SynapticDelays_file = existing_synapse_dir+"SynapticDelays.bin";
    string PresynapticIDs_file = existing_synapse_dir+"PresynapticIDs.bin";
    string PostsynapticIDs_file = existing_synapse_dir+"PostsynapticIDs.bin";

    ifstream SynapticWeights;
    ifstream SynapticDelays;
    ifstream PresynapticIDs;
    ifstream PostsynapticIDs;

    SynapticWeights.open(SynapticWeights_file, ios::in | ios::binary);
    SynapticDelays.open(SynapticDelays_file, ios::in | ios::binary);
    PresynapticIDs.open(PresynapticIDs_file, ios::in | ios::binary);
    PostsynapticIDs.open(PostsynapticIDs_file, ios::in | ios::binary);

    //..
    if(!SynapticWeights.good() || !SynapticDelays.good() || !PresynapticIDs.good() || !PostsynapticIDs.good()) {
        cout<<" !!!! following input binary data file does not exist,  exit !!!!"<<endl;
        cout<<SynapticWeights_file<<endl;
        cout<<SynapticDelays_file<<endl;
        cout<<PresynapticIDs_file<<endl;
        cout<<PostsynapticIDs_file<<endl;

        exit(0);
    }
    //.. map pre_ID with other information of a synapses
    while (SynapticWeights.good()) {
        SynapticWeights.read((char*)&weight, sizeof(float));
        SynapticDelays.read((char*)&delay, sizeof(int));
        PresynapticIDs.read((char*)&pre_ID, sizeof(int));
        PostsynapticIDs.read((char*)&post_ID, sizeof(int));

        if(SynapticWeights.eof()) {
            SynapticWeights.close();
            SynapticDelays.close();
            PresynapticIDs.close();
            PostsynapticIDs.close();
            break;
        }

        synapse_pre_ID_vec.push_back(pre_ID);
        synapse_post_ID_vec.push_back(post_ID);
        synapse_delay_vec.push_back(delay);
        synapse_weight_vec.push_back(weight);
    }
    cout<<" --input synapses size: "<<synapse_pre_ID_vec.size()<<endl;

    //.. now read the starting ID and ending ID of synapses in a group 
    string synapse_start_end_ID_in_group_file = existing_synapse_dir+"synapse_start_end_ID_in_group_file.txt";
    int start, end;
    ifstream fin;
    fin.open(synapse_start_end_ID_in_group_file, ifstream::in);
    while(!fin.good()) {
        fin >> start >> end;
        synapse_start_end_ID_in_group.push_back(make_pair(start, end));
    }
}
//_add synapses from existing out put_
void aki_model::AddSynapseGroup(int id1, int id2, conductance_spiking_synapse_parameters_struct* SYN_PARAMS, SpikingModel* Model, int start, int end) 
{
    vector<int> prevec, postvec;
    vector<float> weightvec;
    vector<float> delayvec; //.. in seconds

    for(int i = start; i<end; i++) {
        prevec.push_back(synapse_pre_ID_vec[i]);
        postvec.push_back(synapse_post_ID_vec[i]);
        weightvec.push_back(synapse_weight_vec[i]);
        delayvec.push_back(synapse_delay_vec[i]*timestep);
    }

    SYN_PARAMS->pairwise_connect_presynaptic = prevec;
    SYN_PARAMS->pairwise_connect_postsynaptic = postvec;
    SYN_PARAMS->pairwise_connect_weight = weightvec;
    SYN_PARAMS->pairwise_connect_delay = delayvec;
    SYN_PARAMS->connectivity_type = CONNECTIVITY_TYPE_PAIRWISE;

    Model->AddSynapseGroup(id1, id2, SYN_PARAMS);
}

//__ loading existing synapses connections. For 4-layer model, making synapstic connect is too slow. 
void aki_model::load_synapses_connections()
{
    //.. read exiting synaptic data ...
    read_synaptic_data();

    int groupID = 0; //.. synapse group ID
    int start = 0;
    int end = 0;

    for (int l=0; l<number_of_layers; l++) {

        if(l==0) {
            for(int i = 1; i<=8; i++) {//..8 Garbor type: -1, -2, .., -8, otained from SpikingModel::AddSynapseGroupsForNeuronGroupAndEachInputGroup printout
                start = synapse_start_end_ID_in_group[groupID].first;
                end = synapse_start_end_ID_in_group[groupID].second;
                AddSynapseGroup(-i, EXCITATORY_NEURONS[l], G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);
                groupID++;
            }
        } else{
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FF[l-1];
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[l-1];

            //.. E2E_FF has 2 connection per pre->post. One can directly set max_number_of_connections_per_pair = 2 in parameter setting, through there are some differences.  See Nashir's answer
            for (int connection_number = 0; connection_number < max_number_of_connections_per_pair; connection_number++){

                start = synapse_start_end_ID_in_group[groupID].first;
                end = synapse_start_end_ID_in_group[groupID].second;

                AddSynapseGroup(EXCITATORY_NEURONS[l-1], EXCITATORY_NEURONS[l], E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);

                groupID++;
            }

            if(E2E_FB_ON) {  // FF has 2 connect pre-post. why FB does not have ???

                start = synapse_start_end_ID_in_group[groupID].first;
                end = synapse_start_end_ID_in_group[groupID].second;

                AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l-1], E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);

                groupID++;
            }
        }

        E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_E2I_L[l];

        start = synapse_start_end_ID_in_group[groupID].first;
        end = synapse_start_end_ID_in_group[groupID].second;

        AddSynapseGroup(EXCITATORY_NEURONS[l], INHIBITORY_NEURONS[l], E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);

        groupID++;

        if (inh_layer_on[l]){

            start = synapse_start_end_ID_in_group[groupID].first;
            end = synapse_start_end_ID_in_group[groupID].second;

            I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_I2E_L[l];
            AddSynapseGroup(INHIBITORY_NEURONS[l], EXCITATORY_NEURONS[l], I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);

            groupID++;
        }

        if(E2E_L_ON) {
            start = synapse_start_end_ID_in_group[groupID].first;
            end = synapse_start_end_ID_in_group[groupID].second;

            AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l], E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS, model, start, end);

            groupID++;
        }
    }

    cout<<" ... number of synapse group: "<<groupID<<endl;
}

//___
void aki_model::define_synapses_parameters()
{

    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = timestep; 
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = timestep;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_G2E_FF;
    //from Brunel10K.cpp,  Biological Scaling factors (ensures that voltage is in mV)
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    // In aki's model, learning on this set of synapses was off because of no learning in input neurons. Remove the line below to match that. WX: the paper mentioned G2E is also modified. 
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->plasticity_vec.push_back(evans_stdp);
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = weight_range_bottom;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = weight_range_top;


    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_FF_minDelay;//5.0*timestep;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_FF_maxDelay;//3.0f*pow(10, -3);
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_FF;
    if (fanInCount_E2E_FF%max_number_of_connections_per_pair!=0)
    {
        printf("number of total syn connections has to be a multiple of max_number_of_connections_per_pair");
    }

    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_E2E_FF;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->plasticity_vec.push_back(evans_stdp);
    // this line is assigned later when assemble layers since each later the value is different
    //E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FF[0];
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_E2E_FF;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = weight_range_bottom;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = weight_range_top;


    E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    if(E2E_FB_ON){
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_FB_minDelay;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_FB_maxDelay;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_FB;

        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_E2E_FB;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->plasticity_vec.push_back(evans_stdp);
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FB;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_E2E_FB;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = weight_range_bottom;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = weight_range_top;
    }


    // no stdp for E2I ???
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2I_L_minDelay; //5.0*timestep;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2I_L_maxDelay; //3.0f*pow(10, -3);
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;  
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = 0.5;//weight_range_bottom;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = 0.5;//weight_range_top;

    // no stdp for I2E ???
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = I2E_L_minDelay;//5.0*timestep;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = I2E_L_maxDelay;//3.0f*pow(10, -3);
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = -70.0*pow(10, -3);
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = 0.5;//weight_range_bottom;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = 0.5;//weight_range_top;

    E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    if(E2E_L_ON){
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_L_minDelay;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_L_maxDelay;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_L;

        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_E2E_L;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
        if (E2E_L_STDP_ON)
            E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->plasticity_vec.push_back(evans_stdp);
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_L;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_E2E_L;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = weight_range_bottom;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = weight_range_top;
    }
}

//___
void aki_model::setup_neuron_groups()
{
    // ADD INPUT NEURONS 
    adding_input_neurons_timer = new TimerWithMessages("Adding Input Neurons...\n");

    // GaborFilter result: Need to include this into a this code
    input_neurons->set_up_rates("FileList.txt", "FilterParameters.txt", (source + filepath + inputs_for_test_name).c_str(), max_FR_of_input_Gabor);
    equalize_rates(input_neurons, 0.1f);

    image_poisson_input_spiking_group_params = new image_poisson_input_spiking_neuron_parameters_struct();
    image_poisson_input_spiking_group_params->rate = 30.0f;
    input_neurons->AddGroupForEachGaborType(image_poisson_input_spiking_group_params);

    adding_input_neurons_timer->stop_timer_and_log_time_and_message("Input Neurons Added.", true);

    // SETTING UP NEURON GROUPS
    // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
    adding_neurons_timer = new TimerWithMessages("Adding Neurons...\n");

    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS = new lif_spiking_neuron_parameters_struct();
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[0] = dim_excit_layer;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[1] = dim_excit_layer;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->resting_potential_v0 = -0.074f;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->threshold_for_action_potential_spike = -0.053f;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_capacitance_Cm = 500.0*pow(10, -12); // 500pF
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_leakage_conductance_g0 = 25.0*pow(10, -9);
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->absolute_refractory_period = absolute_refractory_period;

    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS = new lif_spiking_neuron_parameters_struct();
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[0] = dim_inhib_layer;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[1] = dim_inhib_layer;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->resting_potential_v0 = -0.082f;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->threshold_for_action_potential_spike = -0.053f;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_capacitance_Cm = 214.0*pow(10, -12);
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_leakage_conductance_g0 = 18.0*pow(10, -9);
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->absolute_refractory_period = absolute_refractory_period;

    // Create populations of excitatory and inhibitory neurons for the defined number of layers
    for (int l=0;l<number_of_layers;l++){
        EXCITATORY_NEURONS.push_back(model->AddNeuronGroup(EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
        INHIBITORY_NEURONS.push_back(model->AddNeuronGroup(INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
        cout<<"Neuron Group "<<EXCITATORY_NEURONS[l]<<": Excitatory layer "<<l<<endl;
        cout<<"Neuron Group "<<INHIBITORY_NEURONS[l]<<": Inhibitory layer "<<l<<endl;
    }

    adding_neurons_timer->stop_timer_and_log_time_and_message("Neurons Added.", true);
}

//__
void aki_model::setup_STDP()
{
    // setup STDP. According to the https://sites.google.com/view/spike-simulator/home/simple-example, add plastiscity rule need to be here
    STDP_PARAMS = new evans_stdp_plasticity_parameters_struct();
    STDP_PARAMS->decay_term_tau_C = decay_term_tau_C;
    STDP_PARAMS->decay_term_tau_D = decay_term_tau_D;
    STDP_PARAMS->model_parameter_alpha_D = 0.5;   
    STDP_PARAMS->synaptic_neurotransmitter_concentration_alpha_C = 0.5;  
    STDP_PARAMS->learning_rate_rho = learning_rate_rho;

    evans_stdp = new EvansSTDPPlasticity(conductance_spiking_synapses, lif_spiking_neurons, input_neurons, STDP_PARAMS);

    model->AddPlasticityRule(evans_stdp);
}

//__
void aki_model::set_model_parameters()
{
    //..
    original_timestep = 0.00002;      // This value is the timestep used in Aki's spiking investigations

    // Since the model can be run under different connectivity styles, these booleans turn them on/off
    E2E_FB_ON = true;
    E2E_L_ON = true;
    E2E_L_STDP_ON = true;

    // In order to set up a sensible set of FF exc and inh values, a set of booleans have been set up to turn on/off the values
    for(int i = 0; i<number_of_layers; i++)
        inh_layer_on[i] = true;


    /*
     *
     *  Visual Model General Settings
     *
     */
    // Network Parameters
    max_number_of_connections_per_pair = 2;  // The maximum number of connections refers to multiple synaptic contacts pre->post
    dim_excit_layer = 64;  // The dimension of the excitatory layers (grid with this width)
    dim_inhib_layer = 32;  // The dimension of the inhibitory layers (as above)

    // Measure of the radius of the Fan-in
    // G2E = Gabor to excitatory, E2E = excitatory to excitatory, E2I = excitatory to inhibitory, I2E = inhibitory to excitatory
    // FF = feed forward, L = Lateral, FB = Feedback
    gaussian_synapses_standard_deviation_G2E_FF = 1.0;//12.0;//12.0;

    float tmp[number_of_layers-1] = {8.0,12.0,16.0};// {12.0,18.0,24.0};//{6.0,9.0,12.0};//{8.0,12.0,16.0};//{12.0,18.0,18.0};
    for(int i = 0; i<number_of_layers-1; i++)
        gaussian_synapses_standard_deviation_E2E_FF[i] = tmp[i];

    gaussian_synapses_standard_deviation_E2E_FB = 8.0;//12.0;
    gaussian_synapses_standard_deviation_E2E_L = 4.0;
    gaussian_synapses_standard_deviation_E2I_L = 1.0;
    gaussian_synapses_standard_deviation_I2E_L = 8.0;

    // Fan-in Number
    fanInCount_G2E_FF = 30;
    fanInCount_E2E_FF = 100;
    fanInCount_E2E_FB = 10;//{0, 10} //means two scenarios, 0 or 10 
    fanInCount_E2E_L = 10; //{0, 10} //means two scenarios, 0 or 10 
    fanInCount_E2I_L = 30;
    fanInCount_I2E_L = 30;

    if (fanInCount_E2E_FF%max_number_of_connections_per_pair!=0){
        printf("total_number_of_new_synapses has to be a multiple of max_number_of_connections_per_pair");
    }

    // Synaptic Parameters
    // Range of axonal transmission delay (0.1 ms - 10 ms)
    // timestep is defined above
    min_delay = 0.0001; //0.1 ms
    max_delay = 0.01; // 10ms
    max_FR_of_input_Gabor = 100.0f;
    absolute_refractory_period = 0.002;

    //Synaptic Parameters
    weight_range_bottom = 0.0;
    weight_range_top = 1.0;
    learning_rate_rho = 0.1;

    // calculating different Connections
    E2E_FF_minDelay = min_delay;
    E2E_FF_maxDelay = max_delay;//3.0f*pow(10, -3);
    E2I_L_minDelay = min_delay;
    E2I_L_maxDelay = max_delay;//3.0f*pow(10, -3);
    I2E_L_minDelay = min_delay;
    I2E_L_maxDelay = max_delay;//3.0f*pow(10, -3);
    E2E_FB_minDelay = min_delay;
    E2E_FB_maxDelay = max_delay;
    E2E_L_minDelay = min_delay;
    E2E_L_maxDelay = max_delay;

    // Below are the decay rates of the variables for learning: Pre/Post synaptic activities C and D (See Ben Evans)
    // Aki's model tried 5, 25, 125 ms for both Tau_C and Tau_D. the shorter, the more PGs
    decay_term_tau_C = 0.003;//aki's model:0.005(In Ben's model, tau_C/tau_D = 0.003/0.005 v 0.015/0.025 v 0.075/0.125, and the first one produces the best result)
    decay_term_tau_D = 0.005;

    // Biological Scaling Constant = How much you multiply the weights up or down for realism/stability
    // If this value is roughly on the order of the Leakage Conductance, it will be close to one input spike -> one output spike (n.b. depends on syn tau)
    //float biological_conductance_scaling_constant_lambda_G2E_FF = 0.1 * 0.0001 * original_timestep;  //..0.2ns
    //float biological_conductance_scaling_constant_lambda_E2E_FF = 0.00005 * original_timestep; // 1ns
    //float biological_conductance_scaling_constant_lambda_E2E_FB = 0.1 * 0.0001 * original_timestep; //0.2ns
    //float biological_conductance_scaling_constant_lambda_E2E_L  = 0.000001 * original_timestep; //0.02ns
    //float biological_conductance_scaling_constant_lambda_E2I_L  = 0.001 * original_timestep; // 20ns
    //float biological_conductance_scaling_constant_lambda_I2E_L  = 0.005 * original_timestep; // 100ns

    biological_conductance_scaling_constant_lambda_G2E_FF = 0.2 * 0.0001 * original_timestep; //. below 0.3ns, no spiking. set to 0.4ns 
    biological_conductance_scaling_constant_lambda_E2E_FF = 0.00005 * original_timestep; 
    biological_conductance_scaling_constant_lambda_E2E_FB = 0.0001 * original_timestep; 
    biological_conductance_scaling_constant_lambda_E2E_L  = 0.000001 * original_timestep; 
    biological_conductance_scaling_constant_lambda_E2I_L  = 0.001 * original_timestep; 
    biological_conductance_scaling_constant_lambda_I2E_L  = 0.005 * original_timestep; 


    // is the re-adjust the scaling factor come from optimization?
    float tmp_1[number_of_layers-1] = {
        0.625f * biological_conductance_scaling_constant_lambda_E2E_FF,
        0.5f * biological_conductance_scaling_constant_lambda_E2E_FF,
        0.75f * biological_conductance_scaling_constant_lambda_E2E_FF};
    for(int i = 0; i<number_of_layers-1; i++)
        layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[i] = tmp_1[i];

    // Aki's model fixed at 40ns for all layers. Use these values for now and change to 40ns if needed
    float tmp_2[number_of_layers] = {
        1.1f * biological_conductance_scaling_constant_lambda_E2I_L,
        1.625f * biological_conductance_scaling_constant_lambda_E2I_L,
        0.875f * biological_conductance_scaling_constant_lambda_E2I_L,
        1.6f * biological_conductance_scaling_constant_lambda_E2I_L};
    for(int i = 0; i<number_of_layers; i++)
        layerwise_biological_conductance_scaling_constant_lambda_E2I_L[i] = tmp_2[i];

    // Aki's model fixed at 80 ns. Use these values for now and change to 80ns if needed
    float tmp_3[number_of_layers] = {
        0.04f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.375f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.2f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.325f * biological_conductance_scaling_constant_lambda_I2E_L};
    for(int i = 0; i<number_of_layers; i++)
        layerwise_biological_conductance_scaling_constant_lambda_I2E_L[i] = tmp_3[i];



    // Tau G = Synaptic Conductance Decay TIME CONSTANT for each synapse type (#1) (Seconds)
    // Most of these values are set to 150ms for trace-like learning. Other than Exc->Inh and Inh->Exc
    decay_term_tau_g_G2E_FF = 0.15;
    decay_term_tau_g_E2E_FF = 0.15;//0.15;//0.15;//0.002 v. 0.15 and 0.15 is better?
    decay_term_tau_g_E2E_FB = 0.15;
    decay_term_tau_g_E2E_L = 0.15;//0.002 v. 0.15 and 0.15 is better?
    decay_term_tau_g_E2I_L = 0.002;
    decay_term_tau_g_I2E_L = 0.025;//Aki's model 0.005;//In Ben's model, 0.005 v 0.025 and latter produced better result
}

//__
void aki_model::load_run_config_parameters()
{
    //.. loading run configuation parameters ...
    //
    starting_time = 0;  // if !=0, it's trained and just load the parameter and continue training. When train again, starting_time is the end of the last training time. 
    simtime = 2.0f;  //simulation time in seconds starting from the starting_time; Can be set higher than usual for any period for generally test the network behaviour.
    plasticity_on = false;  // turn on the plasticity or not
    timestep = 0.00002;  // can set to any value for testing. set it to original_timestep when run 

    // Files/Paths relevent to the input set
    source = "";
    filepath = "";
    inputs_for_test_name = "";
    current_weight= "";  //since default it no training yet, there's no current weight location

    configFile.open ("run_config.txt", ifstream::in);
    if(!configFile.good()) {
        cout<<" run configuation file does not exist. Using the default value"<<endl;
    } else {
        configFile >> starting_time;
        configFile >> simtime;
        configFile >> plasticity_on;
        configFile >> timestep;
        configFile >> source;
        configFile >> filepath;
        configFile >> inputs_for_test_name;
        configFile >> current_weight;
        configFile >> existing_synapse_dir;
    }

    // output file location ...
    current_weight = source + current_weight;
    output_location = source + "output/Start"+to_string(starting_time)+"End"+to_string(starting_time+simtime)+"/";
    neuron_dir  = output_location + "neuron_dir/";
    input_dir  = output_location + "input_dir/";
    synapse_dir  = output_location + "synapse_dir/";

    cout<<starting_time<<" "<<to_string(starting_time)<<endl;

    if( mkdir(output_location.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) || 
        mkdir(neuron_dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) || 
        mkdir(input_dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) || 
        mkdir(synapse_dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR)) {

        cout<<" !!! One or all of  following output directory can not be created, exit !!!"<<endl;
        cout<<output_location<<endl;
        cout<<neuron_dir<<endl;
        cout<<input_dir<<endl;
        cout<<synapse_dir<<endl;
        exit(0);
    }

    cout<<" -----------------------------------------------"<<endl;
    cout<<" starting_time: "<<starting_time<<endl;
    cout<<" new_sim_time: "<<simtime<<endl;
    cout<<" plasticity_on: "<<plasticity_on<<endl;
    cout<<" timestep: "<<timestep<<endl;
    cout<<" input file: "<<source + filepath + inputs_for_test_name <<endl;
    cout<<" current weight : "<<current_weight<<endl;
    cout<<" output location: "<<output_location<<endl;
    cout<<" exiting synaptic data location: "<<existing_synapse_dir<<endl;
    cout<<" -----------------------------------------------"<<endl;
}

/** Function to load weights from a file into a SpikingModel.
  save time. No need to run the model again
  note: model->spiking_synapses->load_weights_from_binary(current_weight);  does not work.
 */
void aki_model::load_weights( SpikingModel* Model,  std::string weightloc )  
{
    std::ifstream weightfile;
    std::vector<float> WeightsToLoad; // This vector should ultimately hold the list of replacement weights

    weightfile.open (weightloc, ios::in | ios::binary);
    if(!weightfile.good()) {
        cout<<" !!! Current weight file: "<< weightloc << "does not exist,  exit !!!"<<endl;
        exit(0);
    }
    while( weightfile.good() )
    {
        float currentweight;
        weightfile.read((char*)&currentweight, sizeof(float));
        if (weightfile.eof()) {
            weightfile.close();
            break;
        }
        WeightsToLoad.push_back(currentweight);
    }
    // Check if you have the correct number of weights
    if (WeightsToLoad.size() != Model->spiking_synapses->total_number_of_synapses){
        printf("%d, %d\n", (int)WeightsToLoad.size(), Model->spiking_synapses->total_number_of_synapses);
        printf("The number of weights being loaded is not equivalent to the model.");
        exit(2);
    }
    // If the previous test is passed, apply the weights
    for (int i=0; i < WeightsToLoad.size(); i++){
        Model->spiking_synapses->synaptic_efficacies_or_weights[i] = WeightsToLoad[i];
    }

    cout<<WeightsToLoad.size()<<" Weights Loaded from "<<weightloc<<endl;
}

/** Function to equalize the mean rate of the stimuli being presented to the network.
 *  Not strictly necessary if the stimuli are set up well.
 */
void aki_model::equalize_rates( ImagePoissonInputSpikingNeurons* input_neurons, float target)      
{
    // Rates can be altered here without much issue
    int num_rates_per_image = input_neurons->total_number_of_rates_per_image;
    int num_images = input_neurons->total_number_of_input_stimuli;

    for (int image_index = 0; image_index < num_images; image_index++){
        float meanval = 0.0f;
        for (int rate_index = 0; rate_index < num_rates_per_image; rate_index++){
            meanval += input_neurons->gabor_input_rates[image_index*num_rates_per_image + rate_index];
        }
        meanval /= float(num_rates_per_image);
        // printf("%f\n", meanval);

        float multiplication_factor = target / meanval;
        for (int rate_index = 0; rate_index < num_rates_per_image; rate_index++){
            input_neurons->gabor_input_rates[image_index*num_rates_per_image + rate_index] *= multiplication_factor;
        }
    }
}

