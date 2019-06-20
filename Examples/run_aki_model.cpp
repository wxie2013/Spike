#include "aki_model.C"

//...
void assign_value(auto &var)
{
    int value = 0;

    cin.ignore(); //.. without this, it will go all the way to the end of a loop without user input. need to be at the beginning

    if (std::cin.peek() == '\n') { //check if next character is newline
        cout<<"     use default value: "<<var<<endl;
    } else if (cin >> value) { 
        var = value;
    } else {
        cout<<" ... invalid inputs, need to be integer ..."<<endl; //error handling
    }

}

int main(int argc, char** argv)
{
    //default value
    unsigned int n_epoch = 1; //.. number of epochs
    bool is_ActivityMonitor_on = false; //.. when ON, it consume lots of REM
    bool load_existing_synapses = true; //.. save time
    int which_stimulus = -1; //.. -1 mean all. num > 0, means only that simulus is used as input
    bool binary_output_only=true;  //.. txt format is only needed for debugging

    string option;
    do {
        cout<<" ___ default value? [Y/N]:"; 
        cin>>option;

        if(option=="N" || option=="n") {
            cout<<"   n_epoch: (1 - inf. Default: 1): "; assign_value(n_epoch); 
            cout<<"   is_ActivityMonitor_on (0/1. Default: 0): "; assign_value(is_ActivityMonitor_on);
            cout<<"   load_existing_synapses (0/1. Default: 1): "; assign_value(load_existing_synapses);
            cout<<"   which_stimulus (-1:all/number(>0). Default: -1): "; assign_value(which_stimulus);
            cout<<"   binary_output_only (0/1. Default: 1): "; assign_value(binary_output_only);
        } else if(option=="Y" || option=="y"){
            cout<<" +++ OK. using default value ++"<<endl;
        } else {
            cout<<"please answer Y or N "<<endl;
        }
    } while(option!="Y" && option!="y" && option!="N" && option!="n");



    cout<<" ... user inputs ...."<<endl;
    cout<<"   n_epoch: (1 - inf. Default: 1): "<<n_epoch<<endl;  
    cout<<"   is_ActivityMonitor_on (0/1. Default: 0): "<<is_ActivityMonitor_on<<endl; 
    cout<<"   load_existing_synapses (0/1. Default: 1): "<<load_existing_synapses<<endl; 
    cout<<"   which_stimulus (-1:all/number(>0). Default: -1): "<<which_stimulus<<endl; 
    cout<<"   binary_output_only (0/1. Default: 1): "<<binary_output_only<<endl; 
    cout<<" ..................."<<endl;

    //..
    aki_model model(load_existing_synapses);

    model.activate_ActivityMonitor(is_ActivityMonitor_on);
    model.run_spiking_model(binary_output_only, which_stimulus, n_epoch);

    return 0;
}

