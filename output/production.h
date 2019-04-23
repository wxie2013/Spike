#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iterator>
#include <fstream>

using namespace std;

struct Synapse {
  int pre_ID;
  float weight;
  int delay; // in time_step
};

class production
{
    private:
        //.. binary output from spike
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
        int min_num_afferent_per_neuron; //.. minimum number of afferent neuron per neuron. ??? not sure yet if 1->1 is OK for PG
        vector<int> neuron_with_synapses; //.. ID of neurons with any number of synapses
        map<int, vector<int>> neuton_with_all_afferent; //.. map each neuron with all of its afferent neurons

        void find_neuron_with_synapses();
        void find_all_afferent_neuron_for_a_neuron();

    public:
        production(string dir);
        ~production();

        void read_binary_data();  //.. read binary data produced directly from spike

        void set_min_num_afferent_per_neuron(int in) {min_num_afferent_per_neuron = in;} 

        vector<int> get_neuron_with_synapses() {return neuron_with_synapses;}

        int find_PG(); // find polychronous group
};
