#include "aki_model.C"

int main(int argc, char** argv)
{
    unsigned int n_epoch = 1; //.. number of epochs
    bool is_ActivityMonitor_on = false; //.. when ON, it consume lots of REM
    bool load_existing_synapses = true; //.. save time
    int which_stimulus = -1; //.. -1 mean all. num > 0, means only that simulus is used as input
    bool binary_output_only=true;  //.. txt format is only needed for debugging

    if(argc!=1) {
        cout<<" ... need 5 inputs:...."<<endl;
        cout<<"   n_epoch: (1 - inf): "<<endl; cin >> n_epoch;
        cout<<"   is_ActivityMonitor_on (0/1): "<<endl; cin >>is_ActivityMonitor_on;
        cout<<"   load_existing_synapses (0/1): "<<endl; cin >>load_existing_synapses;
        cout<<"   which_stimulus (-1:all/number(>0)): "<<endl; cin >> which_stimulus;
        cout<<"   binary_output_only (0/1): "<<endl; cin >> binary_output_only;
    }

    cout<<" ... user inputs ...."<<endl;
    cout<<"   n_epoch: (1 - inf): "<<n_epoch<<endl;  
    cout<<"   is_ActivityMonitor_on (0/1): "<<is_ActivityMonitor_on<<endl; 
    cout<<"   load_existing_synapses (0/1): "<<load_existing_synapses<<endl; 
    cout<<"   which_stimulus (-1:all/number(>0)): "<<which_stimulus<<endl; 
    cout<<"   binary_output_only (0/1): "<<binary_output_only<<endl; 
    cout<<" ..................."<<endl;

    //..
    aki_model model(load_existing_synapses);

    model.activate_ActivityMonitor(is_ActivityMonitor_on);
    model.run_spiking_model(binary_output_only, which_stimulus, n_epoch);

    return 0;
}

