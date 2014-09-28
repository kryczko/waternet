#include "storage.h"
#include "analysis.h"
#include "helper.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

void OHdistro(Args& args) {
    
    Information info = args.arg_info;
    TimeSteps time_steps = args.arg_time_steps;
    vector<int> bins ( info.OH_bins );
    int counter = 0;
    double binsize = info.max_OH / info.OH_bins;
    for (int i = 0; i < info.OH_bins; i ++) {
        bins[i] = 0;
    }
    double sum = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.local_H_neighbors.size(); k++) {
                Hydrogen& H = Hvec[O.local_H_neighbors[k]];
                double dist = OHdist(O, H, info);
                sum += dist;
                int bin = dist / binsize;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.OH_output.c_str());
    output << "# distance\tprobability\taverage distance\n\n";
    for (int i = 0; i < info.OH_bins; i ++) {
        output << i*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i*binsize + binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    cout << "Outputted OH distribution data file.\n\n";
    
    
}
