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


void output_graphfile(Args& args) {
    Information info = args.arg_info;
    TimeStep time_step = args.arg_time_step;
    
    O_vector& Ovec = time_step.O_atoms;
    int n_bins = info.num_cell_blocks;
    double bin_size = (info.cell_block_end - info.cell_block_start) / n_bins;
    ofstream output;
    output.open(info.gephi_output.c_str());
    
    output << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" xmlns:viz=\"http://www.gexf.net/1.2draft/viz\">\n"
                << "<meta lastmodifieddate=\"2013-11-21\">\n"
                << "<creator> Kevin Ryczko </creator>\n"
                << "<description> Social Network Visualization </description>\n"
                << "</meta>\n"
                << "<graph mode=\"static\" defaultedgetype=\"directed\">\n"
                << "<nodes>\n";
    for (int i = 0; i < info.num_oxygen; i ++) {
        Oxygen& O = Ovec[i];
        vector<int>& bonds = O.bonded_O_neighbors;
        if (O.z_coords < info.starting_z) {
            output << "<node id=\"" << i << "\" label=\"" << (int) ((O.z_coords + info.lattice_z) / bin_size) << "\"/>\n";
        } else {
            output << "<node id=\"" << i << "\" label=\"" << (int) (O.z_coords / bin_size) << "\"/>\n";
        }
    }
    int count = 0;
    output << "</nodes>\n" << "<edges>\n";
    for (int i = 0; i < info.num_oxygen; i ++) {
        Oxygen& O = Ovec[i];
        vector<int>& bonds = O.bonded_O_neighbors;
        for (int j = 0; j < bonds.size(); j ++) {
            output << "<edge id=\"" << count << "\" source=\"" << i
                            << "\" target=\"" << bonds[j] << "\"/>\n";
            count ++;
        }
    }
    output << "</edges>\n" << "</graph>\n" << "</gexf>";
    output.close();
    cout << "Outputted Gephi graph file.\n\n";
    
}