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


void  degree_distro(Args& args) {
    
    Information info = args.arg_info;
    TimeSteps time_steps = args.arg_time_steps;
    
    double incr = (info.cell_block_end - info.cell_block_start) / info.num_cell_blocks;
    for (int cell = 0; cell < info.num_cell_blocks; cell ++) {
        double start_z = info.cell_block_start + cell*incr;
        double end_z = info.cell_block_start + (cell+1)*incr;
        int max_degree = 0;
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                int degree = O.bonded_O_neighbors.size() + O.out_degree;
                if (degree > max_degree) {
                    max_degree = degree;
                }
            }
        }
        vector<int> degree_counts ( max_degree );
        for (int i = 0; i < max_degree; i ++) {
            degree_counts[i] = 0;
        }
        int counter = 0;
        double sum = 0;
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                if (O.z_coords > start_z && O.z_coords < end_z) {
                    int degree = O.bonded_O_neighbors.size() + O.out_degree;
                    degree_counts[degree] ++;
                    sum += degree;
                    counter ++;
                }
            }
        }
        ofstream output;
        string filename = info.degree_output + "region_" + to_string(cell) + ".dat";
        output.open(filename.c_str());
        output << "# degree\tprobability\taverage degree\n# for region " << start_z << " to " << end_z << "[Angstroms]\n\n";
        for (int i = 0; i < max_degree; i ++) {
            output << i << "\t" << degree_counts[i] / (double) counter << "\t" << sum / counter << "\n";
        }
        output.close();
    }
    cout << "Outputted degree distribution data files.\n\n";
    
    
    
}
