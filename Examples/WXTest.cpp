/*  try to reproduce the model for the paper Psychological Review, 125(4), 545-571.
 *  with the latest version of SPIKE
 *  reference: 
 1. branch aki_test/Experiment/ConductanceExperiment1.cpp.
 2. Gisi_Model.cpp in https://github.com/nasiryahm/LabVisionIntro
 */

#include "Spike/Spike.hpp"

#include <sys/stat.h>

using namespace std;

/** Function to load weights from a file into a SpikingModel.
  save time. No need to run the model again
 */
void load_weights(
        SpikingModel* Model,      /**< SpikingModel Pointer to the model which should load weights */
        std::string weightloc,    /**< String path to the file from which weights should be loaded */
        bool binaryfile)      /**< Boolean flag indicating if the file is a binary file */
{
    std::ifstream weightfile;
    std::vector<float> WeightsToLoad; // This vector should ultimately hold the list of replacement weights

    if (binaryfile){
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
    } else {
        weightfile.open(weightloc);
        if(!weightfile.good()) {
            cout<<" !!! Current weight file: "<< weightloc << " does not exist,  exit !!!"<<endl;
            exit(0);
        }
        while( weightfile.good() )
        {
            string fileline;
            getline( weightfile, fileline);
            if (weightfile.eof()) {
                weightfile.close();
                break;
            }
            // Put each line into the float vector
            WeightsToLoad.push_back( std::stof(fileline) );
        }
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
void equalize_rates(
        ImagePoissonInputSpikingNeurons* input_neurons, /**< ImagePoissonInputSpikingNeuron pointer to initialized input population */
        float target)                   /**< float value indicating desired mean FR */
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

/*
 *  Main function in which the network is created and run
 */

int main (int argc, char *argv[])
{
    /*
     *
     *  General Simulation Settings and Parameters
     *
     */

    //.. loading run configuation parameters ...
    //
    int starting_time = 0;  // if !=0, it's trained and just load the parameter and continue training. When train again, starting_time is the end of the last training time. 
    float simtime = 2.0f;  //simulation time in seconds starting from the starting_time; Can be set higher than usual for any period for generally test the network behaviour.
    bool plasticity_on = false;  // turn on the plasticity or not
    float timestep = 0.00002;  // can set to any value for testing. set it to original_timestep when run 

    // Files/Paths relevent to the input set
    string source = "/home/wxie/AI/Spike/master/Spike/";
    string filepath = source + "Data/MatlabGaborFilter/";
    string inputs_for_test_name = "simpleShapes";
    string current_weight= "";  //since default it no training yet, there's no current weight location
    
    ifstream configFile;
    configFile.open ("run_config.txt", ifstream::in);
    if(!configFile.good()) {
        cout<<" run configuation file does not exist. Using the default value"<<endl;
    } else {
        configFile >> starting_time;
        configFile >> simtime;
        configFile >> plasticity_on;
        configFile >> timestep;
        configFile >> filepath;
        configFile >> inputs_for_test_name;
        configFile >> current_weight;
    }

    // output file location ...
    filepath = source + filepath;
    current_weight = source + current_weight;
    string output_location = source + "output/Start"+to_string(starting_time)+"End"+to_string(starting_time+simtime)+"/";
    string neuron_dir  = output_location + "neuron_dir/";
    string input_dir  = output_location + "input_dir/";
    string synapse_dir  = output_location + "synapse_dir/";

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
    cout<<" input file: "<<filepath + inputs_for_test_name <<endl;
    cout<<" current weight : "<<current_weight<<endl;
    cout<<" output location: "<<output_location<<endl;
    cout<<" -----------------------------------------------"<<endl;

    //..
    float original_timestep = 0.00002;      // This value is the timestep used in Aki's spiking investigations

    // Since the model can be run under different connectivity styles, these booleans turn them on/off
    bool E2E_FB_ON = true;
    bool E2E_L_ON = true;
    bool E2E_L_STDP_ON = true;

    // In order to set up a sensible set of FF exc and inh values, a set of booleans have been set up to turn on/off the values
    bool inh_layer_on[] = {true, true, true, true};  //.. 4-layer network

    /*
     *
     *  Visual Model General Settings
     *
     */
    // Network Parameters
    const int number_of_layers = 4;  // This value is explicitly assumed in this model. Not recommended to change unless you understand what else may need changing in this file.
    int max_number_of_connections_per_pair = 2;  // The maximum number of connections refers to multiple synaptic contacts pre->post
    int dim_excit_layer = 64;  // The dimension of the excitatory layers (grid with this width)
    int dim_inhib_layer = 32;  // The dimension of the inhibitory layers (as above)

    // Measure of the radius of the Fan-in
    // G2E = Gabor to excitatory, E2E = excitatory to excitatory, E2I = excitatory to inhibitory, I2E = inhibitory to excitatory
    // FF = feed forward, L = Lateral, FB = Feedback
    float gaussian_synapses_standard_deviation_G2E_FF = 1.0;//12.0;//12.0;
    float gaussian_synapses_standard_deviation_E2E_FF[number_of_layers-1] = {8.0,12.0,16.0};// {12.0,18.0,24.0};//{6.0,9.0,12.0};//{8.0,12.0,16.0};//{12.0,18.0,18.0};
    float gaussian_synapses_standard_deviation_E2E_FB = 8.0;//12.0;
    float gaussian_synapses_standard_deviation_E2E_L = 4.0;
    float gaussian_synapses_standard_deviation_E2I_L = 1.0;
    float gaussian_synapses_standard_deviation_I2E_L = 8.0;

    // Fan-in Number
    int fanInCount_G2E_FF = 30;
    int fanInCount_E2E_FF = 100;
    int fanInCount_E2E_FB = 10;//{0, 10] //means two scenarios, 0 or 10 
    int fanInCount_E2E_L = 10; //{0, 10} //means two scenarios, 0 or 10 
    int fanInCount_E2I_L = 30;
    int fanInCount_I2E_L = 30;

    if (fanInCount_E2E_FF%max_number_of_connections_per_pair!=0){
        printf("total_number_of_new_synapses has to be a multiple of max_number_of_connections_per_pair");
        return 0;
    }

    // Synaptic Parameters
    // Range of axonal transmission delay
    // timestep is defined above
    float min_delay = 0.0001//5.0*timestep; // In timesteps
    float max_delay = 0.01; // In seconds (10ms)
    float max_FR_of_input_Gabor = 100.0f;
    float absolute_refractory_period = 0.002;

    //Synaptic Parameters
    float weight_range_bottom = 0.0;
    float weight_range_top = 1.0;
    float learning_rate_rho = 0.1;

    // calculating different Connections
    float E2E_FF_minDelay = min_delay;
    float E2E_FF_maxDelay = max_delay;//3.0f*pow(10, -3);
    float E2I_L_minDelay = min_delay;
    float E2I_L_maxDelay = max_delay;//3.0f*pow(10, -3);
    float I2E_L_minDelay = min_delay;
    float I2E_L_maxDelay = max_delay;//3.0f*pow(10, -3);
    float E2E_FB_minDelay = min_delay;
    float E2E_FB_maxDelay = max_delay;
    float E2E_L_minDelay = min_delay;
    float E2E_L_maxDelay = max_delay;

    // Below are the decay rates of the variables for learning: Pre/Post synaptic activities C and D (See Ben Evans)
    // Aki's model tried 5, 25, 125 ms for both Tau_C and Tau_D. the shorter, the more PGs
    float decay_term_tau_C = 0.003;//aki's model:0.005(In Ben's model, tau_C/tau_D = 0.003/0.005 v 0.015/0.025 v 0.075/0.125, and the first one produces the best result)
    float decay_term_tau_D = 0.005;

    // Biological Scaling Constant = How much you multiply the weights up or down for realism/stability
    // If this value is roughly on the order of the Leakage Conductance, it will be close to one input spike -> one output spike (n.b. depends on syn tau)
    float biological_conductance_scaling_constant_lambda_G2E_FF = 0.1 * 0.0001 * original_timestep;  //..0.2ns
    float biological_conductance_scaling_constant_lambda_E2E_FF = 0.00005 * original_timestep; // 1ns
    float biological_conductance_scaling_constant_lambda_E2E_FB = 0.1 * 0.0001 * original_timestep; //0.2ns
    float biological_conductance_scaling_constant_lambda_E2E_L  = 0.000001 * original_timestep; //0.02ns
    float biological_conductance_scaling_constant_lambda_E2I_L  = 0.001 * original_timestep; // 20ns
    float biological_conductance_scaling_constant_lambda_I2E_L  = 0.005 * original_timestep; // 100ns

    // is the re-adjust the scaling factor come from optimization? 
    float layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[number_of_layers-1] = {
        0.625f * biological_conductance_scaling_constant_lambda_E2E_FF,
        0.5f * biological_conductance_scaling_constant_lambda_E2E_FF,
        0.75f * biological_conductance_scaling_constant_lambda_E2E_FF};

    // Aki's model fixed at 40ns for all layers. Use these values for now and change to 40ns if needed
    float layerwise_biological_conductance_scaling_constant_lambda_E2I_L[number_of_layers] = {
        1.1f * biological_conductance_scaling_constant_lambda_E2I_L,
        1.625f * biological_conductance_scaling_constant_lambda_E2I_L,
        0.875f * biological_conductance_scaling_constant_lambda_E2I_L,
        1.6f * biological_conductance_scaling_constant_lambda_E2I_L};

    // Aki's model fixed at 80 ns. Use these values for now and change to 80ns if needed 
    float layerwise_biological_conductance_scaling_constant_lambda_I2E_L[number_of_layers] = {
        0.04f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.375f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.2f * biological_conductance_scaling_constant_lambda_I2E_L,
        0.325f * biological_conductance_scaling_constant_lambda_I2E_L};


    // Tau G = Synaptic Conductance Decay TIME CONSTANT for each synapse type (#1) (Seconds)
    // Most of these values are set to 150ms for trace-like learning. Other than Exc->Inh and Inh->Exc
    float decay_term_tau_g_G2E_FF = 0.15;
    float decay_term_tau_g_E2E_FF = 0.15;//0.15;//0.15;//0.002 v. 0.15 and 0.15 is better?
    float decay_term_tau_g_E2E_FB = 0.15;
    float decay_term_tau_g_E2E_L = 0.15;//0.002 v. 0.15 and 0.15 is better?
    float decay_term_tau_g_E2I_L = 0.002;
    float decay_term_tau_g_I2E_L = 0.025;//Aki's model 0.005;//In Ben's model, 0.005 v 0.025 and latter produced better result


    /*
     *
     *  Defining the Spiking Model
     *
     */

    // Create the SpikingModel
    SpikingModel* model = new SpikingModel();

    // Set up the simulator with a timestep at which the neuron, synapse and STDP properties will be calculated
    model->SetTimestep(timestep);

    // Choose an input neuron type
    ImagePoissonInputSpikingNeurons* input_neurons = new ImagePoissonInputSpikingNeurons();

    // Choose your neuron type
    LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();

    // Choose your synapse type
    ConductanceSpikingSynapses * conductance_spiking_synapses = new ConductanceSpikingSynapses();

    // Allocate your chosen components to the simulator
    model->input_spiking_neurons = input_neurons;
    model->spiking_neurons = lif_spiking_neurons;
    model->spiking_synapses = conductance_spiking_synapses;

    // setup STDP. According to the https://sites.google.com/view/spike-simulator/home/simple-example, add plastiscity rule need to be here
    evans_stdp_plasticity_parameters_struct * STDP_PARAMS = new evans_stdp_plasticity_parameters_struct();
    STDP_PARAMS->decay_term_tau_C = decay_term_tau_C;
    STDP_PARAMS->decay_term_tau_D = decay_term_tau_D;
    STDP_PARAMS->model_parameter_alpha_D = 0.5;   
    STDP_PARAMS->synaptic_neurotransmitter_concentration_alpha_C = 0.5;  
    STDP_PARAMS->learning_rate_rho = learning_rate_rho;

    EvansSTDPPlasticity* evans_stdp = new EvansSTDPPlasticity(conductance_spiking_synapses, lif_spiking_neurons, input_neurons, STDP_PARAMS);
    // shall I use this one? 
    //WeightNormSTDPPlasticity * evans_stdp = new WeightNormSTDPPlasticity((SpikingSynapses *) conductance_spiking_synapses, (SpikingNeurons *) lif_spiking_neurons, (SpikingNeurons *) input_neurons, (stdp_plasticity_parameters_struct *) STDP_PARAMS); // the paper used weight normalized rule (???)

    model->AddPlasticityRule(evans_stdp);

    // ADD ANY ACTIVITY MONITORS OR PLASTICITY RULES YOU WISH FOR
    SpikingActivityMonitor* spike_monitor = new SpikingActivityMonitor(lif_spiking_neurons);
    SpikingActivityMonitor* input_spike_monitor = new SpikingActivityMonitor(input_neurons);
    model->AddActivityMonitor(spike_monitor);
    model->AddActivityMonitor(input_spike_monitor);

    /*
       SETUP PROPERTIES AND CREATE NETWORK:
     * Note:
     * All Neuron, Synapse and STDP types have associated parameters structures.
     * These structures are defined in the header file for that class and allow us to set properties.
     */
    // ADD INPUT NEURONS 
    TimerWithMessages * adding_input_neurons_timer = new TimerWithMessages("Adding Input Neurons...\n");

    // GaborFilter result: Need to include this into a this code
    input_neurons->set_up_rates("FileList.txt", "FilterParameters.txt", ("../../Data/MatlabGaborFilter/"+inputs_for_test_name+"/").c_str(), max_FR_of_input_Gabor);
    equalize_rates(input_neurons, 0.1f);

    image_poisson_input_spiking_neuron_parameters_struct * image_poisson_input_spiking_group_params = new image_poisson_input_spiking_neuron_parameters_struct();
    image_poisson_input_spiking_group_params->rate = 30.0f;
    input_neurons->AddGroupForEachGaborType(image_poisson_input_spiking_group_params);

    adding_input_neurons_timer->stop_timer_and_log_time_and_message("Input Neurons Added.", true);

    // SETTING UP NEURON GROUPS
    // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
    TimerWithMessages * adding_neurons_timer = new TimerWithMessages("Adding Neurons...\n");

    lif_spiking_neuron_parameters_struct * EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS = new lif_spiking_neuron_parameters_struct();
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[0] = dim_excit_layer;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[1] = dim_excit_layer;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->resting_potential_v0 = -0.074f;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->threshold_for_action_potential_spike = -0.053f;
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_capacitance_Cm = 500.0*pow(10, -12); // 500pF
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_leakage_conductance_g0 = 25.0*pow(10, -9);
    EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS->absolute_refractory_period = absolute_refractory_period;

    lif_spiking_neuron_parameters_struct * INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS = new lif_spiking_neuron_parameters_struct();
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[0] = dim_inhib_layer;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->group_shape[1] = dim_inhib_layer;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->resting_potential_v0 = -0.082f;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->threshold_for_action_potential_spike = -0.053f;
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_capacitance_Cm = 214.0*pow(10, -12);
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->somatic_leakage_conductance_g0 = 18.0*pow(10, -9);
    INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS->absolute_refractory_period = absolute_refractory_period;

    // Create populations of excitatory and inhibitory neurons for the defined number of layers
    vector<int> EXCITATORY_NEURONS;
    vector<int> INHIBITORY_NEURONS;
    for (int l=0;l<number_of_layers;l++){
        EXCITATORY_NEURONS.push_back(model->AddNeuronGroup(EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
        INHIBITORY_NEURONS.push_back(model->AddNeuronGroup(INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
        cout<<"Neuron Group "<<EXCITATORY_NEURONS[l]<<": Excitatory layer "<<l<<endl;
        cout<<"Neuron Group "<<INHIBITORY_NEURONS[l]<<": Inhibitory layer "<<l<<endl;
    }

    adding_neurons_timer->stop_timer_and_log_time_and_message("Neurons Added.", true);

    // SETTING UP SYNAPSES
    // Creating a synapses parameter structure for connections from the input neurons to the excitatory neurons
    TimerWithMessages * adding_synapses_timer = new TimerWithMessages("Adding Synapses...\n");

    conductance_spiking_synapse_parameters_struct * G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = timestep; //???
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = timestep;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1; //??? set to 1.  later when forming layers, it AddSynapseGroup twice. Why not just set it to 2 here? 
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_G2E_FF;
    //from Brunel10K.cpp,  Biological Scaling factors (ensures that voltage is in mV)
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    // In aki's model, learning on this set of synapses was off because of no learning in input neurons. Remove the line below to match that.
    //G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.plasticity_vec.push_back(evans_stdp);
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_G2E_FF;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = weight_range_bottom;
    G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = weight_range_top;


    conductance_spiking_synapse_parameters_struct * E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_FF_minDelay;//5.0*timestep;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_FF_maxDelay;//3.0f*pow(10, -3);
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1;
    E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_FF;
    if (fanInCount_E2E_FF%max_number_of_connections_per_pair!=0)
    {
        printf("number of total syn connections has to be a multiple of max_number_of_connections_per_pair");
        return 0;
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


    conductance_spiking_synapse_parameters_struct * E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    if(E2E_FB_ON){
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_FB_minDelay;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_FB_maxDelay;
        E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1; 
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
    conductance_spiking_synapse_parameters_struct * E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2I_L_minDelay; //5.0*timestep;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2I_L_maxDelay; //3.0f*pow(10, -3);
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1; 
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = 0.0;  //??? should this be -70mv
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_E2I_L;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = 0.5;//weight_range_bottom;
    E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = 0.5;//weight_range_top;

    // no stdp for I2E ???
    conductance_spiking_synapse_parameters_struct * I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = I2E_L_minDelay;//5.0*timestep;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = I2E_L_maxDelay;//3.0f*pow(10, -3);
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1; 
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_per_postsynaptic_neuron = fanInCount_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = biological_conductance_scaling_constant_lambda_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->reversal_potential_Vhat = -70.0*pow(10, -3);
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->decay_term_tau_g = decay_term_tau_g_I2E_L;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[0] = 0.5;//weight_range_bottom;
    I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_range[1] = 0.5;//weight_range_top;

    conductance_spiking_synapse_parameters_struct * E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS = new conductance_spiking_synapse_parameters_struct();
    if(E2E_L_ON){
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[0] = E2E_L_minDelay;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->delay_range[1] = E2E_L_maxDelay;
        E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->max_number_of_connections_per_pair = 1; 
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

    for (int l=0; l<number_of_layers; l++){
        if(l==0)
            model->AddSynapseGroupsForNeuronGroupAndEachInputGroup(EXCITATORY_NEURONS[l], G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
        else{
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FF[l-1];
            E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS->weight_scaling_constant = layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[l-1];
            //.. E2E_FF has 2 connection per pre->post
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

    adding_synapses_timer->stop_timer_and_log_time_and_message("Synapses Added.", true);

    /*
     *
     *  Finalize and Run the model
     *
     */

    // if it's already trained, load any weights before finalising the model
    if (starting_time != 0){
        load_weights(model, current_weight, true);
    }
    model->run(simtime, plasticity_on);

    //use binary mode. The text mode is too slow ..
    spike_monitor->save_spikes_as_binary(neuron_dir);
    input_spike_monitor->save_spikes_as_binary(input_dir);
    model->spiking_synapses->save_connectivity_as_binary(synapse_dir);
}
