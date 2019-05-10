/*  try to reproduce the model for the paper Psychological Review, 125(4), 545-571.
 *  with the latest version of SPIKE
 *  reference: 
 1. branch aki_test/Experiment/ConductanceExperiment1.cpp.
 2. Gisi_Model.cpp in https://github.com/nasiryahm/LabVisionIntro
 */

#include "Spike/Spike.hpp"

#include <sys/stat.h>

using namespace std;

const int number_of_layers = 4;  // This value is explicitly assumed in this model. Not recommended to change unless you understand what else may need changing in this file.

class aki_model
{
    private:
        float starting_time;
        float simtime;
        bool plasticity_on;
        float timestep;

        // Files/Paths relevent to the input set
        string source;
        string filepath;
        string inputs_for_test_name;
        string current_weight;
        string output_location;
        string neuron_dir;
        string input_dir;
        string synapse_dir;

    private:
        ifstream configFile;

        //..
        float original_timestep;

        // Since the model can be run under different connectivity styles, these booleans turn them on/off
        bool E2E_FB_ON;
        bool E2E_L_ON;
        bool E2E_L_STDP_ON;

        // In order to set up a sensible set of FF exc and inh values, a set of booleans have been set up to turn on/off the values
        bool inh_layer_on[number_of_layers];

        /*
         *
         *  Visual Model General Settings
         *
         */
        // Network Parameters
        int max_number_of_connections_per_pair;
        int dim_excit_layer;
        int dim_inhib_layer;

        // Measure of the radius of the Fan-in
        // G2E = Gabor to excitatory, E2E = excitatory to excitatory, E2I = excitatory to inhibitory, I2E = inhibitory to excitatory
        // FF = feed forward, L = Lateral, FB = Feedback
        float gaussian_synapses_standard_deviation_G2E_FF;
        float gaussian_synapses_standard_deviation_E2E_FF[number_of_layers-1];
        float gaussian_synapses_standard_deviation_E2E_FB;
        float gaussian_synapses_standard_deviation_E2E_L;
        float gaussian_synapses_standard_deviation_E2I_L;
        float gaussian_synapses_standard_deviation_I2E_L;

        // Fan-in Number
        int fanInCount_G2E_FF;
        int fanInCount_E2E_FF;
        int fanInCount_E2E_FB;
        int fanInCount_E2E_L;
        int fanInCount_E2I_L;
        int fanInCount_I2E_L;


        // Synaptic Parameters
        // Range of axonal transmission delay (0.1 ms - 10 ms)
        // timestep is defined above
        float min_delay;
        float max_delay;
        float max_FR_of_input_Gabor;
        float absolute_refractory_period;

        //Synaptic Parameters
        float weight_range_bottom;
        float weight_range_top;
        float learning_rate_rho;

        // calculating different Connections
        float E2E_FF_minDelay;
        float E2E_FF_maxDelay;
        float E2I_L_minDelay;
        float E2I_L_maxDelay;
        float I2E_L_minDelay;
        float I2E_L_maxDelay;
        float E2E_FB_minDelay;
        float E2E_FB_maxDelay;
        float E2E_L_minDelay;
        float E2E_L_maxDelay;

        // Below are the decay rates of the variables for learning: Pre/Post synaptic activities C and D (See Ben Evans)
        // Aki's model tried 5, 25, 125 ms for both Tau_C and Tau_D. the shorter, the more PGs
        float decay_term_tau_C;
        float decay_term_tau_D;

        // Biological Scaling Constant = How much you multiply the weights up or down for realism/stability
        // If this value is roughly on the order of the Leakage Conductance, it will be close to one input spike -> one output spike (n.b. depends on syn tau)
        //float biological_conductance_scaling_constant_lambda_G2E_FF = 0.1 * 0.0001 * original_timestep;  //..0.2ns
        //float biological_conductance_scaling_constant_lambda_E2E_FF = 0.00005 * original_timestep; // 1ns
        //float biological_conductance_scaling_constant_lambda_E2E_FB = 0.1 * 0.0001 * original_timestep; //0.2ns
        //float biological_conductance_scaling_constant_lambda_E2E_L  = 0.000001 * original_timestep; //0.02ns
        //float biological_conductance_scaling_constant_lambda_E2I_L  = 0.001 * original_timestep; // 20ns
        //float biological_conductance_scaling_constant_lambda_I2E_L  = 0.005 * original_timestep; // 100ns

        float biological_conductance_scaling_constant_lambda_G2E_FF;
        float biological_conductance_scaling_constant_lambda_E2E_FF;
        float biological_conductance_scaling_constant_lambda_E2E_FB;
        float biological_conductance_scaling_constant_lambda_E2E_L ;
        float biological_conductance_scaling_constant_lambda_E2I_L ;
        float biological_conductance_scaling_constant_lambda_I2E_L ;


        // is the re-adjust the scaling factor come from optimization? 
        float layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[number_of_layers-1];

        // Aki's model fixed at 40ns for all layers. Use these values for now and change to 40ns if needed
        float layerwise_biological_conductance_scaling_constant_lambda_E2I_L[number_of_layers];

        // Aki's model fixed at 80 ns. Use these values for now and change to 80ns if needed 
        float layerwise_biological_conductance_scaling_constant_lambda_I2E_L[number_of_layers];


        // Tau G = Synaptic Conductance Decay TIME CONSTANT for each synapse type (#1) (Seconds)
        // Most of these values are set to 150ms for trace-like learning. Other than Exc->Inh and Inh->Exc
        float decay_term_tau_g_G2E_FF;
        float decay_term_tau_g_E2E_FF;
        float decay_term_tau_g_E2E_FB;
        float decay_term_tau_g_E2E_L;
        float decay_term_tau_g_E2I_L;
        float decay_term_tau_g_I2E_L;


        /*
         *
         *  Defining the Spiking Model
         *
         */

        // Create the SpikingModel
        SpikingModel* model;

        // Choose an input neuron type
        ImagePoissonInputSpikingNeurons* input_neurons;

        // Choose your neuron type
        LIFSpikingNeurons* lif_spiking_neurons;

        // Choose your synapse type
        ConductanceSpikingSynapses * conductance_spiking_synapses;

        // setup STDP. According to the https://sites.google.com/view/spike-simulator/home/simple-example, add plastiscity rule need to be here
        evans_stdp_plasticity_parameters_struct * STDP_PARAMS;

        EvansSTDPPlasticity* evans_stdp;

        // ADD ANY ACTIVITY MONITORS OR PLASTICITY RULES YOU WISH FOR
        SpikingActivityMonitor* spike_monitor;
        SpikingActivityMonitor* input_spike_monitor;

        /*
           SETUP PROPERTIES AND CREATE NETWORK:
         * Note:
         * All Neuron, Synapse and STDP types have associated parameters structures.
         * These structures are defined in the header file for that class and allow us to set properties.
         */
        // ADD INPUT NEURONS 
        TimerWithMessages * adding_input_neurons_timer;

        // GaborFilter result: Need to include this into a this code

        image_poisson_input_spiking_neuron_parameters_struct * image_poisson_input_spiking_group_params;

        // SETTING UP NEURON GROUPS
        // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
        TimerWithMessages * adding_neurons_timer;

        lif_spiking_neuron_parameters_struct * EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS;
        lif_spiking_neuron_parameters_struct * INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS;

        // Create populations of excitatory and inhibitory neurons for the defined number of layers
        vector<int> EXCITATORY_NEURONS;
        vector<int> INHIBITORY_NEURONS;
        TimerWithMessages * adding_synapses_timer;

        conductance_spiking_synapse_parameters_struct * G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
        conductance_spiking_synapse_parameters_struct * E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
        conductance_spiking_synapse_parameters_struct * E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
        conductance_spiking_synapse_parameters_struct * E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
        conductance_spiking_synapse_parameters_struct * I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
        conductance_spiking_synapse_parameters_struct * E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;


    private:
        void load_weights(SpikingModel* Model, std::string weightloc); 
        void equalize_rates( ImagePoissonInputSpikingNeurons* input_neurons, float target);
        void load_run_config_parameters();
        void set_model_parameters();
        void define_spiking_model();
        void setup_STDP();
        void setup_neuron_groups();
        void setup_synapses();
        void define_synapses_parameters();
        void make_synapses_connections();

    public:
        aki_model();
        ~aki_model();

        void run_spiking_model(bool binary_output_only);
};
