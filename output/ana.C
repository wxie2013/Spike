R__LOAD_LIBRARY(./lib_production.so)

void ana(string dir="")
{
    production pd;
    pd.SetIntputBinaryFile(dir);
    pd.read_in_SpikeTimes_data();
    pd.read_SpikeTimes_data();
    pd.read_Synapses_data();

    pd.set_max_number_of_connections_per_pair(2); //.. maximum two snapses per pair
    pd.find_PG();
}

