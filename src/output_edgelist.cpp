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


void output_edgelist(Args& args) {
    Information info = args.arg_info;
    TimeSteps time_steps = args.arg_time_steps;
    if (info.create_edgelist) {
        ofstream output;
        output.open(info.edgelist_output_filename.c_str());
    
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            H_vector& Hvec = time_steps[i].H_atoms;
            output << "Frame:\t" << i << "\t" << num_edges(Ovec) << "\n";
            for (int j = 0; j < Ovec.size(); j ++) {
                for (int k = 0; k < Ovec[j].bonded_O_neighbors.size(); k ++) {
                    output << Ovec[j].ID << "\t" << Ovec[j].bonded_O_neighbors[k] << "\n";
                }
            }
        }
        cout << "Outputted edge list to data file given.\n\n";
        output.close();
    }
    out_count(info, time_steps);
}
