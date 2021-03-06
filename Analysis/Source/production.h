#ifndef PRODUCTION_H
#define PRODUCTION_H

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

//
//#include "PolychromousGroup.h"

using namespace std;

struct Synapse {
    int pre_ID;
    float weight;
    int delay; // in time_step
    float pre_spikeTime; //pre_ID spike time in timestep

    Synapse() {
        pre_ID = 0;
        weight = 0;
        delay = 0;
        pre_spikeTime = -1;
    }

    bool operator< ( const Synapse &s ) const { return pre_ID < s.pre_ID; } //..used for equal_range
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
        map<int, vector<Synapse>> neuron_with_all_afferent; //.. map each neuron with all of its afferent neurons

        vector<int> neuron_with_spikes; //.. ID of neurons with any number of spikes
        map<int, vector<float>> neuron_with_all_spikesTime; //.. map each neuron with all of its spike time

        //.. adding spike time into synapses. pair<post_ID, spike_time> is mapped to a vector contains all its synapses 
        //.. with preID spike at time consistent with it's delay. It's a mutimap since there could be more than one synapses
        //.. send spike to fire the post_ID
        multimap<pair<int, float>, vector<Synapse>> synapses_with_spikes; 

        //
        void find_post_neuron_with_synapses();
        void find_all_afferent_neuron_for_a_neuron();
        void find_all_spikeTime_of_a_neuron(multimap<int, float>&);

        template <typename M, typename Vtype>
            Vtype find_a_key_in_a_map(M m, int id);

        void combine_spikeTime_Synapses();
        void clear();

        //
        //PolychromousGroup *PG; 

    public:
        production();
        ~production();

        void set_min_num_afferent_per_neuron(int in) {min_num_afferent_per_neuron = in;} 
        void set_max_number_of_connections_per_pair(int in) {max_number_of_connections_per_pair = in;} 

        //..
        void SetIntputBinaryFile(string &);  //..open spike output binary data. 
        void read_in_SpikeTimes_data(); //.. read input spiking information from Spike binary output
        void read_SpikeTimes_data(); //.. read spiking information from Spike binary output
        void read_Synapses_data(); //.. read synapse information from Spike binary output

        //..
        void analyze_weight_change_after_STDP(string &, string &); //.. production ntuple of weight before and after STDP
        //..
        map<int, vector<pair<Synapse, Synapse>>> find_duplicated_pre_post_pairs(); //.. map each neuron with all its afferent neurons that has two synases per pair
        //..
        void Fig_9(string &, string &); //.. reproduce Fig.9 of the paper 

        //..
        int find_PG(); // find polychronous group

        ClassDef(production, 1)
};

#endif
