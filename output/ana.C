#include "production.C"

void ana(string dir="")
{
    production pd(dir);
    pd.read_data();
}

