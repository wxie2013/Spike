#include "aki_model.C"

int main(int argc, char** argv)
{
    if(argc!=1 && argc!=4) {
        cout<<" ... need three inputs: load_existing_synapses & binary_output_only=true ...."<<endl;
        cout<<" e.g. ./run_aki_model 0 0 0"<<endl;
        exit(0);
    }
    bool is_ActivityMonitor_on = false; //.. when ON, it consume lots of REM
    bool load_existing_synapses = true; //.. save time
    bool binary_output_only=true;  //.. txt format is only needed for debugging

    if(argc!=1) {
        load_existing_synapses = atoi(argv[1]);
        is_ActivityMonitor_on = atoi(argv[2]);;
        binary_output_only = atoi(argv[3]);
    } 

    //..
    aki_model model(load_existing_synapses);

    model.activate_ActivityMonitor(is_ActivityMonitor_on);
    model.run_spiking_model(binary_output_only);

    return 0;
}

