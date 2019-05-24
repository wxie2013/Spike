#include "production.C"

void ana(string dir="")
{
    production pd;
    pd.SetIntputBinaryFile(dir);
    pd.read_in_SpikeTimes_data();
    pd.read_SpikeTimes_data();
    pd.read_Synapses_data();

    pd.set_max_number_of_connections_per_pair(1);
    pd.find_PG();
}

