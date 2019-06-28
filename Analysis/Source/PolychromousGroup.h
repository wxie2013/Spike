//..Polychronmous group
class PolychromousGroup
{
    private:
        //.. each neuron and its decendents.  
        //.. pair<post_ID, spike_time> is mapped to a vector contains all its synapses with preID spike at time consistent with it's delay
        //.. the last decendents ends at input neuron.
        vector<map<pair<int, float>, vector<Synapse>>> neuron_id; 

        PolychromousGroup *PG;

    public:

    public:
        PolychromousGroup() {;}
        ~PolychromousGroup() {;}
        
        PolychromousGroup* get_PG() {return PG;}

};

