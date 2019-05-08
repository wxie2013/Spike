#include "production.h"
#include <tuple>


production::production()
{
    num_PG = 0;
    min_num_afferent_per_neuron = 1; //..default 1.0 for now.
    max_number_of_connections_per_pair = 2; //.. default 2 synapses per pair (original Aki's model)
}
//__
production::~production()
{
}

//__ read input spiking nueron data 
void production::read_in_SpikeTimes_data()
{
    float in_spkT;
    int in_spkID;

    //..
    in_SpikeTimes.open(in_SpikeTimes_file, ios::in | ios::binary);
    in_SpikeIDs.open(in_SpikeIDs_file, ios::in | ios::binary);

    if(!in_SpikeTimes.good() || !in_SpikeIDs.good()) { 
        cout<<" !!!! following input binary data file does not exist,  exit !!!!"<<endl;
        cout<<in_SpikeTimes_file<<endl;
        cout<<in_SpikeIDs_file<<endl;
        exit(0);
    }
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
}

//__ read spiking nueron data 
void  production::read_SpikeTimes_data()
{
    float spkT;
    int spkID;

    //..
    SpikeTimes.open(SpikeTimes_file, ios::in | ios::binary);
    SpikeIDs.open(SpikeIDs_file, ios::in | ios::binary);

    if(!SpikeTimes.good() || !SpikeIDs.good()) {
        cout<<" !!!! following input binary data file does not exist,  exit !!!!"<<endl;
        cout<<SpikeTimes_file<<endl;
        cout<<SpikeIDs_file<<endl;
        exit(0);
    }
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
}

//__ read synapses data 
void production::read_Synapses_data()
{
    float in_spkT, spkT, weight;
    int in_spkID, spkID, delay, pre_ID, post_ID;

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
        synapse_infor.pre_ID = pre_ID;
        synapse_infor.weight = weight;
        synapse_infor.delay = delay;

        map_Synapse.insert(make_pair(post_ID, synapse_infor));
    }
    cout<<" --input spike neuron size: "<<map_Synapse.size()<<endl;
}

//__ open data file from Spike output 
void production::SetIntputBinaryFile(string dir)
{
    neuron_dir = dir+"/neuron_dir/";
    input_dir = dir+"/input_dir/";
    synapse_dir = dir+"/synapse_dir/";

    //.. input neurons ....
    in_SpikeTimes_file = input_dir+"SpikeTimes.bin";
    in_SpikeIDs_file = input_dir+"SpikeIDs.bin";

    //.. layer 1 and above 
    SpikeTimes_file = neuron_dir+"SpikeTimes.bin";
    SpikeIDs_file = neuron_dir+"SpikeIDs.bin";

    //.. synapses 
    SynapticWeights_file = synapse_dir+"SynapticWeights.bin";
    SynapticDelays_file = synapse_dir+"SynapticDelays.bin";
    PresynapticIDs_file = synapse_dir+"PresynapticIDs.bin";
    PostsynapticIDs_file = synapse_dir+"PostsynapticIDs.bin";

}

//__ find unique post_IDs of neurons with any number of synapse in the network
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

//__ map each post neuron with all of its afferent neurons 
void production::find_all_afferent_neuron_for_a_neuron()
{
    find_neuron_with_synapses();

    pair<multimap<int, Synapse>::iterator, multimap<int, Synapse>::iterator> range;
    vector<int> pre_ids;

    for(unsigned int i = 0; i< neuron_with_synapses.size(); i++) {

        pre_ids.clear(); //.. clear before handling each neuron

        range = map_Synapse.equal_range(neuron_with_synapses[i]); //.. all afferent neuron of a post synaptic neuron. The map is sorted already

        for(multimap<int, Synapse>::iterator it = range.first; it!=range.second; it++) {
            pre_ids.push_back(it->second.pre_ID);
        }
        neuron_with_all_afferent.insert(make_pair(neuron_with_synapses[i], pre_ids));

        //.. find duplicated key-value pairs
        std::pair<std::vector<int>::iterator,std::vector<int>::iterator> bounds;
        vector<int> nid = pre_ids;
        sort(nid.begin(), nid.end());
        for(unsigned int j = 0; j<nid.size(); j++) {
            bounds=std::equal_range (nid.begin(), nid.end(), nid[j]);
            if((bounds.second - bounds.first) >max_number_of_connections_per_pair) 
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

//__
void production::analyze_weight_change_after_STDP(string dir1, string dir2)
{
    TFile f("analyse_weight_change_after_STDP.root", "RECREATE");

    TNtuple nt("nt", "", "preid:postid:w1:w2:delay1:delay2");

    //.. get synapse information from 1st file 
    SetIntputBinaryFile(dir1);
    read_Synapses_data();
    multimap<int, Synapse> map1 = map_Synapse;

    //.. clear the map_Synapse before reading the 2nd file 
    map_Synapse.clear();

    //.. get synapse information from 1st file 
    SetIntputBinaryFile(dir2);
    read_Synapses_data();
    multimap<int, Synapse> map2 = map_Synapse;

    //.. the two map are syncronized ..
    multimap<int, Synapse>::iterator it1  = map1.begin();
    multimap<int, Synapse>::iterator it2  = map2.begin();

    while(it1 != map1.end()) {

        int postid1 = it1->first;
        int preid1 = it1->second.pre_ID;
        float w1 = it1->second.weight;
        int delay1 = it1->second.delay;

        int postid2 = it2->first;
        int preid2 = it2->second.pre_ID;
        float w2 = it2->second.weight;
        int delay2 = it2->second.delay;

        if(postid1!=postid2 || preid1!=preid2 || delay1!=delay2) {
            cout<<" !!! the two files are not syncronized, exit !!!"<<endl;
            exit(0);
        }

        nt.Fill(preid1, postid1, w1, w2, delay1, delay2);

        it1++;
        it2++;
    }

    nt.Write();
    f.Close();

    cout<<" -- created: analyse_weight_change_after_STDP.root ----"<<endl;
}
