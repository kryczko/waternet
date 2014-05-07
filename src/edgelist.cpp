#include <fstream>
#include <iostream>
#include <cmath>
#include "storage.h"

using namespace std;

int pbc_round(double input) {
	int i  = input;

	if (abs(input - i) >= 0.5) {
		if (input > 0) {i += 1;}
		if (input < 0) {i -= 1;}
	}
return i;
}

void nearest_neighbors(Information& info, Oxygen& O, O_vector& Ovec, int& timestep) {
    for (int i = 0; i < Ovec.size(); i ++) {
        Oxygen& O2 = Ovec[i];
        if (O.ID != O2.ID) {
            double dx = O.x_coords[timestep] - O2.x_coords[timestep];
            double dy = O.y_coords[timestep] - O2.y_coords[timestep];
            double dz = O.z_coords[timestep] - O2.z_coords[timestep];
            
            dx -= info.lattice_x * pbc_round(dx/info.lattice_x);
            dy -= info.lattice_y * pbc_round(dy/info.lattice_y);
            dz -= info.lattice_z * pbc_round(dz/info.lattice_z);
            
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            if (dist < 3) {
                O.nearest_neighbors.push_back(O2.ID);
            }
        }
    }
}

bool create_edgelist(Information& info, TimeSteps& time_steps) {
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            nearest_neighbors(info, O, Ovec, i);
        }
    }
    cout << "Computed nearest neighbor lists...\n\n";
    return true;
}