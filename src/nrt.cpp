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


void nrt(Args& args) {
    
    Information info = args.arg_info;
    TimeSteps time_steps = args.arg_time_steps;
    
    vector<int> time_counter;
    
    int vector_size = info.time_length;
    if (info.time_length > info.n_frames) {
        vector_size = (int) info.n_frames;
    } 
    double normalizer = 0;
    time_counter.resize(vector_size - 1);
    vector<int> time_counts = set_zero(time_counter);
    cout << "--- Network reorganization time progress ---\n";
    int length = time_steps.size() / info.num_starts;
    
    for (int average_step = 0; average_step < time_steps.size(); average_step += length) {
        int start = average_step + 1;
        int stop = info.time_length + average_step + 1;
        if (stop > time_steps.size()) {
            stop = time_steps.size();
        }
        O_vector& Ovec_before = time_steps[average_step].O_atoms;
        for (int i = 0; i < Ovec_before.size();  i++) {
            normalizer += Ovec_before[i].nearest_neighbors.size();
        }
        for (int step = start; step < stop - 1; step ++) {
            O_vector& Ovec_after = time_steps[step].O_atoms;
            for (int atom = 0; atom < Ovec_before.size(); atom ++) {
                Oxygen& O_before = Ovec_before[atom];
                Oxygen& O_after = Ovec_after[atom];
                for (int n = 0; n < O_before.nearest_neighbors.size(); n ++) {
                    int neighbor_before = O_before.nearest_neighbors[n];
                    for (int m = 0; m < O_after.nearest_neighbors.size(); m ++) {
                        int neighbor_after = O_after.nearest_neighbors[m];
                        if (neighbor_before == neighbor_after) {
                            time_counts[step - start] ++;
                        }
                    }
                }
            }
        }
        if ( average_step < time_steps.size()) {
            cout << " --- |< " << (int) ( 100 * (double) average_step / (double) time_steps.size()) << " % >| ---\r";
            flush(cout);
        } else if (average_step == time_steps.size() - 1) {
            cout << " --- |< 100 % >| ---\n";
            flush(cout);
        }
        
    }
    cout << " --- |< 100 % >| ---\n\n";    
    ofstream output;
    output.open(info.nrt_output.c_str());
    output << "# Averaged " << 1 << " times.\n\n";
    for (int i = 0; i < time_counts.size(); i ++) {
        output << i*info.time_step / 1000.0 << "\t" << 100.0 * time_counts[i] / normalizer << "\n";
    }
    cout << "Outputted network reorganization time data file.\n\n";
    
}
