#include <fstream>
#include <iostream>
#include "TNtuple.h"
using namespace std;

void plot(string dir)
{
    string neuron_dir = dir+"/neuron_dir/";
    string input_dir = dir+"/input_dir/";
    string synapse_dir = dir+"/synapse_dir/";

    //.. input neurons ....
    string in_SpikeTimes_file = input_dir+"SpikeTimes.bin";
    string in_SpikeIDs_file = input_dir+"SpikeIDs.bin";
    ifstream in_SpikeTimes(in_SpikeTimes_file, ios::in | ios::binary);
    ifstream in_SpikeIDs(in_SpikeIDs_file, ios::in | ios::binary);

    //.. layer 1 and above 
    string SpikeTimes_file = neuron_dir+"SpikeTimes.bin";
    string SpikeIDs_file = neuron_dir+"SpikeIDs.bin";
    ifstream SpikeTimes(SpikeTimes_file, ios::in | ios::binary);
    ifstream SpikeIDs(SpikeIDs_file, ios::in | ios::binary);

    //.. synapses 
    string SynapticWeights_file = synapse_dir+"SynapticWeights.bin";
    string SynapticDelays_file = synapse_dir+"SynapticDelays.bin";
    string PresynapticIDs_file = synapse_dir+"PresynapticIDs.bin";
    string PostsynapticIDs_file = synapse_dir+"PostsynapticIDs.bin";

    ifstream SynapticWeights(SynapticWeights_file, ios::in | ios::binary);
    ifstream SynapticDelays(SynapticDelays_file, ios::in | ios::binary);
    ifstream PresynapticIDs(PresynapticIDs_file, ios::in | ios::binary);
    ifstream PostsynapticIDs(PostsynapticIDs_file, ios::in | ios::binary);

    //..
    float in_spkT, spkT, weight;
    int in_spkID, spkID, delay, pre_ID, post_ID;

    //..
    cout<<" .. input neurons .."<<endl;
    TNtuple* in_nt = new TNtuple("in_nt", "", "id:t");
    while (in_SpikeTimes.good()) {
        in_SpikeTimes.read((char*)&in_spkT, sizeof(float));
        in_SpikeIDs.read((char*)&in_spkID, sizeof(int));

        if(in_SpikeTimes.eof()) {
            in_SpikeTimes.close();
            in_SpikeIDs.close();
            break;
        }
        in_nt->Fill(in_spkID, in_spkT);
    }

    //..
    cout<<" .. neurons .."<<endl;
    TNtuple* nt = new TNtuple("nt", "", "id:t");

    while (SpikeTimes.good()) {
        SpikeTimes.read((char*)&spkT, sizeof(float));
        SpikeIDs.read((char*)&spkID, sizeof(int));

        if(SpikeTimes.eof()) {
            SpikeTimes.close();
            SpikeIDs.close();
            break;
        }
        nt->Fill(spkID, spkT);
    }

    //..
    cout<<" .. synapses ..."<<endl;
    TNtuple* synapse = new TNtuple("synapse", "", "preid:postid:w:delay");

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
        synapse->Fill(pre_ID, post_ID, weight, delay);
    }
}
