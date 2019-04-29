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

    //..
    Synapse synapse_infor;  //. holder for synapse postID, delay and weight
    multimap<int, float> map_in_spk_ID_T;  //.. map of ID and spiking time of input neuron
    multimap<int, float> map_spk_ID_T;  //.. map of ID and spiking time of neuron
    multimap<int, Synapse> map_Synapse; //.. map of presynaptic ID and (postsynaptic ID, weight, delay)

    //...
    num_PG = 0;
    min_num_afferent_per_neuron = 1; //..default 1.0 for now.
}
//__
production::~production()
{
}
//__ read binary data produced directly from spike
void production::read_binary_data()
{
    //..
    float in_spkT, spkT, weight;
    int in_spkID, spkID, delay, pre_ID, post_ID;

    //.. map input neuron spike time and ID
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
        synapse_infor.pre_ID = pre_ID;
        synapse_infor.weight = weight;
        synapse_infor.delay = delay;

        map_Synapse.insert(make_pair(post_ID, synapse_infor));
    }
    cout<<" --synapses size: "<<map_Synapse.size()<<endl;
}

//__ find unique IDs of neurons with any number of synapse in the network
void production::find_neuron_with_synapses()
{
    for(multimap<int, Synapse>::iterator it = map_Synapse.begin(); it != map_Synapse.end(); it++) 
        neuron_with_synapses.push_back(it->first);

    //.. now remove duplicates..
    sort(neuron_with_synapses.begin(), neuron_with_synapses.end()); //.. sort it first as a general procedure
    vector<int>::iterator ip = unique(neuron_with_synapses.begin(), neuron_with_synapses.end()); //.. call unique function. ip is the address of the last unique element
    neuron_with_synapses.resize(distance(neuron_with_synapses.begin(), ip)); 

    cout<<" --- "<<neuron_with_synapses.size() << " of unique neurons in the network with synapses "<<endl;
}

//__ map each neuron with all of its afferent neurons 
void production::find_all_afferent_neuron_for_a_neuron()
{
    find_neuron_with_synapses();

    pair<multimap<int, Synapse>::iterator, multimap<int, Synapse>::iterator> range;
    vector<int> id;

    for(unsigned int i = 0; i< neuron_with_synapses.size(); i++) {

        id.clear(); //.. clear before handling each neuron

        range = map_Synapse.equal_range(neuron_with_synapses[i]); //.. all afferent neuron of a post synaptic neuron. The map is sorted already

        for(multimap<int, Synapse>::iterator it = range.first; it!=range.second; it++) {
            id.push_back(it->second.pre_ID);
        }
        neuton_with_all_afferent.insert(make_pair(neuron_with_synapses[i], id));

        //.. find duplicated key-value pairs
        std::pair<std::vector<int>::iterator,std::vector<int>::iterator> bounds;
        vector<int> nid = id;
        sort(nid.begin(), nid.end());
        for(unsigned int j = 0; j<nid.size(); j++) {
            bounds=std::equal_range (nid.begin(), nid.end(), nid[j]);
            if((bounds.second - bounds.first) >=2) 
                cout<<neuron_with_synapses[i]<<" "<<nid[j]<<endl;
        }
    }

}

//__ find polychronous group 
int production::find_PG()
{
    find_all_afferent_neuron_for_a_neuron();
    return num_PG;
}
