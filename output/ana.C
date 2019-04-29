#include "production.C"

void ana(string dir="")
{
    production pd(dir);
    pd.read_binary_data();
    pd.find_PG();
}

