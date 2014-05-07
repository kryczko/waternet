#include <fstream>
#include "storage.h"
using namespace std;

bool output_edgelist(Information& info, TimeSteps& time_steps) {
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        
        ofstream output;
        output.open("output.dat");
        
        for (int j = 0; j < Ovec.size(); j ++) {
            for (int k = 0; k < Ovec[j].nearest_neighbors.size(); k ++) {
                output << Ovec[j].ID << "\t" << Ovec[j].nearest_neighbors[k] << "\n";
            }
        }
        output.close();
    }
    return true;
}