#include <fstream>
#include <iostream>
#include "storage.h"

using namespace std;

bool create_edgelist(Information& info, TimeSteps& time_steps) {
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        cout << i << "\t" << Ovec[0].x_coords[i] << "\t" << Ovec[0].y_coords[i] << "\t" << Ovec[0].z_coords[i] << "\n";
    }
    return true;
}