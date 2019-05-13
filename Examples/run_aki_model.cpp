#include "aki_model.C"

int main(int argc, char** argv)
{
    if(argc!=1 && argc!=3) {
        cout<<" ... need two inputs: load_existing_synapses & binary_output_only=true ...."<<endl;
        cout<<" e.g. ./run_aki_model 0  0"<<endl;
        exit(0);
    }
    bool load_existing_synapses = true;
    bool binary_output_only=true;

    if(argc!=1) {
        load_existing_synapses = atoi(argv[1]);
        binary_output_only = atoi(argv[2]);
    } 

    //..
    aki_model model(load_existing_synapses);
    model.run_spiking_model(binary_output_only);

    return 0;
}

