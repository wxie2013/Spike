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

//__
void production::clear()
{
    map_in_spk_ID_T.clear();
    map_spk_ID_T.clear();
    map_Synapse.clear();

    //..
    neuron_with_all_afferent.clear();
    neuron_with_synapses.clear();
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
        cout<<" !!!! following input binary data file does not exist !!!!"<<endl;
        cout<<in_SpikeTimes_file<<endl;
        cout<<in_SpikeIDs_file<<endl;
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
        cout<<" !!!! following input binary data file does not exist!!!!"<<endl;
        cout<<SpikeTimes_file<<endl;
        cout<<SpikeIDs_file<<endl;
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

        //.. when doing EE_L, i.e. lateral connection between excited neurons, 
        //.. there's no protection against self connection . The fraction of 
        //.. this self connection is 0.01%, thus ignored
        if(pre_ID == post_ID) 
            continue;

        //
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
void production::SetIntputBinaryFile(string &dir)
{
    clear(); //.. clear all global maps and vector before reading in new SPIKE inputs

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
void production::find_post_neuron_with_synapses()
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
    pair<multimap<int, Synapse>::iterator, multimap<int, Synapse>::iterator> range;
    vector<Synapse> pre_ids;

    for(unsigned int i = 0; i< neuron_with_synapses.size(); i++) {

        pre_ids.clear(); //.. clear before handling each neuron

        range = map_Synapse.equal_range(neuron_with_synapses[i]); //.. all afferent neuron of a post synaptic neuron. The map is sorted already

        for(multimap<int, Synapse>::iterator it = range.first; it!=range.second; it++) {
            pre_ids.push_back(it->second);
        }
        neuron_with_all_afferent.insert(make_pair(neuron_with_synapses[i], pre_ids));
    }
}

//__  find duplicated key-value pairs
map<int, vector<pair<Synapse, Synapse>>> production::find_duplicated_pre_post_pairs()
{
    //
    find_post_neuron_with_synapses();
    find_all_afferent_neuron_for_a_neuron();

    //
    map<int, vector<pair<Synapse, Synapse>>> neuron_with_all_afferent_with_2_synapses;

    for(map<int, vector<Synapse>>::iterator it = neuron_with_all_afferent.begin(); it !=neuron_with_all_afferent.end(); it++) {

        std::pair<std::vector<Synapse>::iterator,std::vector<Synapse>::iterator> bounds;

        vector<Synapse> nid = it->second;
        sort(nid.begin(), nid.end());

        vector<pair<Synapse, Synapse>> two_synapses;

        for(unsigned int i = 0; i<nid.size(); i++) {
            bounds=std::equal_range (nid.begin(), nid.end(), nid[i]);

            if((bounds.second - bounds.first) == max_number_of_connections_per_pair) {
                unsigned int index = bounds.first - nid.begin();

                two_synapses.push_back(make_pair(nid[index], nid[index+1]));
                i++; //.. no need to go through again the 2nd synapse of the same pair
            }
        }

        //..
        neuron_with_all_afferent_with_2_synapses.insert(make_pair(it->first, two_synapses));
    }

    return neuron_with_all_afferent_with_2_synapses;
}

//__ find polychronous group 
int production::find_PG()
{
    find_duplicated_pre_post_pairs();
    return num_PG;
}

//__
void production::analyze_weight_change_after_STDP(string &dir1, string &dir2)
{
    string file1 = dir1;
    string file2 = dir2;
    replace(file1.begin(), file1.end(), '/', '_'); //.. note: use '' instead of ""
    replace(file2.begin(), file2.end(), '/', '_'); //.. note: use '' instead of ""

    string outfile = "compare_"+file1+"_"+file2+".root";
    TFile f(outfile.c_str(), "RECREATE");

    TNtuple nt("nt", "", "preid:postid:w1:w2:delay1:delay2");

    //.. get synapse information from 1st file 
    SetIntputBinaryFile(dir1);
    read_Synapses_data();
    multimap<int, Synapse> map1 = map_Synapse;

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
            //.. a small bug in the latest version of Spike, only a few synapses affect, thus ignore ..
            //.. in the CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE of Synapse.cpp,  total_probability reduce after each pre_ID is assigned. 
            //.. because of the precision, the last total_probability can be negative, e.g. -1e-7. 
            //.. then the randval > probability_trace can never be satisfied. Then last pre_ID will have the initial value, i.e. 0
            //.. I tried double precision on total_probability and pre_neuron_probabilities and still cannot completely solve the problem. 
            //.. since the number is small and only happens to input neurons, I choose to ignore it. 
            if(!(
                        (postid1==0 && postid2==-1) || \
                        (postid1==-1 && postid2==0) || \
                        (preid1==0 && preid2==-1) || \
                        (preid1==-1 && preid2==0)
                )) { 
                cout<<" !!! the two files are not syncronized, exit !!!"<<endl;
                exit(0);
            }
        }

        nt.Fill(preid1, postid1, w1, w2, delay1, delay2);

        it1++;
        it2++;
    }

    nt.Write();
    f.Close();

    cout<<" -- created: "<<outfile<<" ----"<<endl;
}

//__
void production::Fig_9(string &dir1, string &dir2)
{
    TFile f("Fig_9.root", "RECREATE");

    TNtuple nt("nt", "", "preid:postid:w1:w2:tw1:tw2:delay1:delay2");

    //.. get synapse information from 1st file 
    cout<<" loading synapses in: "<<dir1<<endl;
    SetIntputBinaryFile(dir1);
    read_Synapses_data();
    map<int, vector<pair<Synapse, Synapse>>> m1 = find_duplicated_pre_post_pairs();

    //.. get synapse information from 2nd file 
    cout<<" loading synapses in: "<<dir2<<endl;
    SetIntputBinaryFile(dir2);
    read_Synapses_data();
    map<int, vector<pair<Synapse, Synapse>>> m2 = find_duplicated_pre_post_pairs();

    if(m1.size() != m2.size()) {
        cout<<" !!! the two map should have the same size, exit !!!"<<endl;
        exit(0);
    }
    //..
    map<int, vector<pair<Synapse, Synapse>>>::iterator it1  = m1.begin();
    map<int, vector<pair<Synapse, Synapse>>>::iterator it2  = m2.begin();

    //..
    while(it1 != m1.end()) {

        int postid = it1->first;
        for(unsigned int i = 0; i<it1->second.size(); i++) {
            //.. sanitary check
            if(it1->first != it2->first || 
               it1->second[i].first.pre_ID != it2->second[i].first.pre_ID || 
               it1->second[i].second.pre_ID != it2->second[i].second.pre_ID || 
               it1->second[i].first.delay != it2->second[i].first.delay ||
               it1->second[i].second.delay != it2->second[i].second.delay) {

                cout<<" !! the pre_ID and post_ID should be syncronized in these two maps, exit !!"<<endl;
                exit(0);
            }
            //
            int preid1 = it1->second[i].first.pre_ID;
            int preid2 = it1->second[i].second.pre_ID;

            float w1 = it1->second[i].first.weight;
            float w2 = it1->second[i].second.weight;

            float tw1 = it2->second[i].first.weight;  //.. trained ..
            float tw2 = it2->second[i].second.weight;  //.. trained ..

            float delay1 = it1->second[i].first.delay;
            float delay2 = it1->second[i].second.delay;

            if(preid1 !=preid2) {
                cout<<" !!! the two preids, i.e. id1: "<<preid1<<" and id2: "<<preid2<<" need to be identical, exit(0) "<<endl;
                exit(0);
            }

            nt.Fill(preid1, postid, w1, w2, tw1, tw2, delay1, delay2); 
        }

        it1++;
        it2++;
    }



    nt.Write();
    f.Close();

    cout<<" -- created: Fig_9.root ----"<<endl;
}
