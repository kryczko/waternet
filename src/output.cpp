#include <fstream>
#include <iostream>
#include "storage.h"
using namespace std;

int num_edges(O_vector& Ovec) {
    int count = 0;
    for (int j = 0; j < Ovec.size(); j ++) {
        for (int k = 0; k < Ovec[j].bonded_O_neighbors.size(); k ++) {
            count ++;
        }
    }
    return count;
}

bool output_graphfile(Information& info, TimeStep& time_step) {
    O_vector& Ovec = time_step.O_atoms;
    ofstream output;
    output.open("h_bond_network.gexf");
    
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
        output << "<node id=\"" << i << "\" label=\"" << (int) O.z_coords << "\"/>\n";
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
    return true;
}

bool output_edgelist(Information& info, TimeSteps& time_steps) {
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
    
    if (info.output_gephi) {
        output_graphfile(info, time_steps[time_steps.size() - 1]);
        cout << "Gephi graph file created.\n\n";
    }
    return true;
}