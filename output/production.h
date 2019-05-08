#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iterator>
#include <fstream>

//
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

struct Synapse {
  int pre_ID;
  float weight;
  int delay; // in time_step
};

class production
{
    private: //.. binary output from spike

        string neuron_dir;
        string input_dir;
        string synapse_dir;

        //.. input neurons ....
        string in_SpikeTimes_file;
        string in_SpikeIDs_file;

        ifstream in_SpikeTimes;
        ifstream in_SpikeIDs;

        //.. layer 1 and above 
        string SpikeTimes_file;
        string SpikeIDs_file;

        ifstream SpikeTimes;
        ifstream SpikeIDs;

        //.. synapses 
        string SynapticWeights_file;
        string SynapticDelays_file;
        string PresynapticIDs_file;
        string PostsynapticIDs_file;

        ifstream SynapticWeights;
        ifstream SynapticDelays;
        ifstream PresynapticIDs;
        ifstream PostsynapticIDs;

        //..
        multimap<int, float> map_in_spk_ID_T; //.. map input neuron spike time and ID
        multimap<int, float> map_spk_ID_T; //.. map other neuron spike time and ID
        multimap<int, Synapse> map_Synapse; //.. map post_ID with other information of a synapses 
        Synapse synapse_infor;

    private: //.. for Polychronous groups
        int num_PG;  //.. number of polychronous group
        int max_number_of_connections_per_pair; //.. max number of synapses per pair 
        int min_num_afferent_per_neuron; //.. minimum number of afferent neuron per neuron. ??? not sure yet if 1->1 is OK for PG
        vector<int> neuron_with_synapses; //.. ID of neurons with any number of synapses
        map<int, vector<int>> neuron_with_all_afferent; //.. map each neuron with all of its afferent neurons

        void find_neuron_with_synapses();
        vector<int> get_neuron_with_synapses() {return neuron_with_synapses;}
        void find_all_afferent_neuron_for_a_neuron();

    public:
        production();
        ~production();

        void set_min_num_afferent_per_neuron(int in) {min_num_afferent_per_neuron = in;} 
        void set_max_number_of_connections_per_pair(int in) {max_number_of_connections_per_pair = in;} 

        //..
        void SetIntputBinaryFile(string dir);  //..open spike output binary data. 
        void read_in_SpikeTimes_data(); //.. read input spiking information from Spike binary output
        void read_SpikeTimes_data(); //.. read spiking information from Spike binary output
        void read_Synapses_data(); //.. read synapse information from Spike binary output

        //..
        void analyze_weight_change_after_STDP(string, string); //.. production ntuple of weight before and after STDP

        //..
        int find_PG(); // find polychronous group
};
