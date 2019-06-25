R__LOAD_LIBRARY(./lib_production.so)

#include "production.h"

void ana_Fig_9(const char* in1 = "", const char* in2 = "")
{
    string dir1(in1);
    string dir2(in2);

    production pd;

    pd.Fig_9(dir1, dir2);
}

