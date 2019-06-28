R__LOAD_LIBRARY(../lib/libproduction.so)

#include "production.h"

void analyze_weight_change_after_STDP(const char* in1 = "", const char * in2 = "")
{
    //.. note the "{}" instead of ()
    vector<string> dir1 {"old_Start0.000000End2.000000/"};

    /*
    vector<string> dir1 { 
            "Start0.000000End2.000000/", 
            "Start0.000000End4.000000/", 
            "Start0.000000End6.000000/", 
            "Start0.000000End8.000000/", 
            "Start0.000000End10.000000/", 
            "Start0.000000End12.000000/", 
            "Start0.000000End14.000000/", 
            "Start0.000000End16.000000/", 
            "Start0.000000End18.000000/", 
            "Start0.000000End20.000000/", 
            "Start0.000000End22.000000/", 
            "Start0.000000End24.000000/", 
            "Start0.000000End26.000000/", 
            "Start0.000000End28.000000/", 
            "Start0.000000End30.000000/", 
            "Start0.000000End32.000000/", 
            "Start0.000000End34.000000/", 
            "Start0.000000End36.000000/", 
            "Start0.000000End38.000000/", 
            "Start0.000000End40.000000/", 
            "Start0.000000End42.000000/", 
            "Start0.000000End44.000000/", 
            "Start0.000000End46.000000/", 
            "Start0.000000End48.000000/", 
            "Start0.000000End50.000000/" 
    };
    //.. note the "{}" instead of ()
    vector<string> dir1 { 
            "Start0.000000End2.000000/",
            "Start2.000000End4.000000/",
            "Start4.000000End6.000000/",
            "Start6.000000End8.000000/",
            "Start8.000000End10.000000/",
            "Start10.000000End12.000000/",
            "Start12.000000End14.000000/",
            "Start14.000000End16.000000/",
            "Start16.000000End18.000000/",
            "Start18.000000End20.000000/",
            "Start20.000000End22.000000/",
            "Start22.000000End24.000000/",
            "Start24.000000End26.000000/",
            "Start26.000000End28.000000/",
            "Start28.000000End30.000000/",
            "Start30.000000End32.000000/",
            "Start32.000000End34.000000/",
            "Start34.000000End36.000000/",
            "Start36.000000End38.000000/",
            "Start38.000000End40.000000/",
            "Start40.000000End42.000000/",
            "Start42.000000End44.000000/",
            "Start44.000000End46.000000/",
            "Start46.000000End48.000000/",
            "Start48.000000End50.000000/"
    };
    */

    //..
    //vector<string> dir2(dir1.size(), "./");
    vector<string> dir2{"old_Start0.000000End10.000000/"};

    string data1(in1);
    string data2(in2);

    //..
    production pd;

    if(dir1.size() != dir2.size()) {
        cout<<" ___ the two input arrays need to have the same size. They are different. exit ____"<<endl;
        exit(0);
    }
    for(unsigned short i = 0; i<dir1.size(); i++) {
        string tmp_dir1 = data1+dir1[i];
        string tmp_dir2 = data2+dir2[i];

        cout<<tmp_dir1<<endl;
        cout<<tmp_dir2<<endl;

        pd.analyze_weight_change_after_STDP(tmp_dir1, tmp_dir2);
    }
}

