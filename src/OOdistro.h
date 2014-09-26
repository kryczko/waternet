#ifndef OODISTRO_H_
#define OODISTRO_H_

#include "storage.h"
#include "analysis.h"
#include "helper.h"

using namespace std;


void  OOdistro(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    vector<int> bins ( info.OO_bins );
    int counter = 0;
    double binsize = info.max_OO / info.OO_bins;
    for (int i = 0; i < info.OO_bins; i ++) {
        bins[i] = 0;
    }
    double sum = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.bonded_O_neighbors.size(); k++) {
                Oxygen& O2 = Ovec[O.bonded_O_neighbors[k]];
                double dist = OOdist(O, O2, info);
                sum += dist;
                int bin = dist / binsize;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.OO_output.c_str());
    output << "# distance\tprobability\taverage distance\n\n";
    for (int i = 0; i < info.OO_bins; i ++) {
        output << i*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i*binsize + binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    
    
}
#endif