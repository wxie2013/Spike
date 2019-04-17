#include "production.h"

production::production(string dir)
{
    neuron_dir = dir+"/neuron_dir/";
    input_dir = dir+"/input_dir/";
    synapse_dir = dir+"/synapse_dir/";

    //.. input neurons ....
    in_SpikeTimes_file = input_dir+"SpikeTimes.bin";
    in_SpikeIDs_file = input_dir+"SpikeIDs.bin";

    in_SpikeTimes.open(in_SpikeTimes_file, ios::in | ios::binary);
    in_SpikeIDs.open(in_SpikeIDs_file, ios::in | ios::binary);

    //.. layer 1 and above 
    SpikeTimes_file = neuron_dir+"SpikeTimes.bin";
    SpikeIDs_file = neuron_dir+"SpikeIDs.bin";

    SpikeTimes.open(SpikeTimes_file, ios::in | ios::binary);
    SpikeIDs.open(SpikeIDs_file, ios::in | ios::binary);

    //.. synapses 
    SynapticWeights_file = synapse_dir+"SynapticWeights.bin";
    SynapticDelays_file = synapse_dir+"SynapticDelays.bin";
    PresynapticIDs_file = synapse_dir+"PresynapticIDs.bin";
    PostsynapticIDs_file = synapse_dir+"PostsynapticIDs.bin";

    SynapticWeights.open(SynapticWeights_file, ios::in | ios::binary);
    SynapticDelays.open(SynapticDelays_file, ios::in | ios::binary);
    PresynapticIDs.open(PresynapticIDs_file, ios::in | ios::binary);
    PostsynapticIDs.open(PostsynapticIDs_file, ios::in | ios::binary);

    //..
    if(!in_SpikeTimes.good() || !in_SpikeIDs.good() || !SpikeTimes.good() || !SpikeIDs.good() ||
       !SynapticWeights.good() || !SynapticDelays.good() || !PresynapticIDs.good() || !PostsynapticIDs.good()) {
        cout<<" !!!! one of the input binary data file does not exist,  exit !!!!"<<endl;
        exit(0);
    }
}
//__
production::~production()
{
}
//__
void production::read_data()
{
    //..
    float in_spkT, spkT, weight;
    int in_spkID, spkID, delay, pre_ID, post_ID;

    //.. map input neuron spike time and ID
    multimap<int, float> map_in_spk_ID_T;
    while (in_SpikeTimes.good()) {
        in_SpikeTimes.read((char*)&in_spkT, sizeof(float));
        in_SpikeIDs.read((char*)&in_spkID, sizeof(int));

        if(in_SpikeTimes.eof()) {
            in_SpikeTimes.close();
            in_SpikeIDs.close();
            break;
        }
        map_in_spk_ID_T.insert(make_pair(-in_spkID, in_spkT)); //.. "-" because in synapses collection, input neuron ID is labelled as negative
    }
    cout<<" --input spike neuron size: "<<map_in_spk_ID_T.size()<<endl;


    //.. map other neuron spike time and ID
    multimap<int, float> map_spk_ID_T;
    while (SpikeTimes.good()) {
        SpikeTimes.read((char*)&spkT, sizeof(float));
        SpikeIDs.read((char*)&spkID, sizeof(int));

        if(SpikeTimes.eof()) {
            SpikeTimes.close();
            SpikeIDs.close();
            break;
        }
        map_spk_ID_T.insert(make_pair(spkID, spkT));
    }
    cout<<" --spike neuron size: "<<map_spk_ID_T.size()<<endl;

    //.. map pre_ID with other information of a synapses 
    Synapse synapse_infor;
    multimap<int, Synapse> map_Synapse;
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
        synapse_infor.post_ID = post_ID;
        synapse_infor.weight = weight;
        synapse_infor.delay = delay;

        map_Synapse.insert(make_pair(pre_ID, synapse_infor));
    }
    cout<<" --synapses size: "<<map_Synapse.size()<<endl;
}
