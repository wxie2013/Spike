#include <iostream>
#include <map>
#include <string>
#include <iterator>
#include <fstream>

using namespace std;

struct Synapse {
  int post_ID;
  float weight;
  int delay; // in time_step
};

class production
{
    private:
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

        multimap<int, float> map_in_spk_ID_T; //.. map input neuron spike time and ID
        multimap<int, float> map_spk_ID_T; //.. map other neuron spike time and ID
        multimap<int, Synapse> map_Synapse; //.. map pre_ID with other information of a synapses 
        Synapse synapse_infor;

    public:
        production(string dir);
        ~production();
        void read_data();
};
