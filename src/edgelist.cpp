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

double OOdist(Oxygen& O1, Oxygen& O2, Information& info) {
    double dx = O1.x_coords - O2.x_coords;
    double dy = O1.y_coords - O2.y_coords;
    double dz = O1.z_coords - O2.z_coords;
    
    dx -= info.lattice_x * pbc_round(dx/info.lattice_x);
    dy -= info.lattice_y * pbc_round(dy/info.lattice_y);
    dz -= info.lattice_z * pbc_round(dz/info.lattice_z);
    
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

double OHdist(Oxygen& O, Hydrogen& H, Information& info) {
    double dx = O.x_coords - H.x_coords;
    double dy = O.y_coords - H.y_coords;
    double dz = O.z_coords - H.z_coords;
    
    dx -= info.lattice_x * pbc_round(dx/info.lattice_x);
    dy -= info.lattice_y * pbc_round(dy/info.lattice_y);
    dz -= info.lattice_z * pbc_round(dz/info.lattice_z);
    
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

// radians to degrees
const double rtd = 57.2957795;

double angle_between(Oxygen& O1, Oxygen& O2, Hydrogen& H, Information& info) {
    double odx = O2.x_coords - O1.x_coords;
    double ody = O2.y_coords - O1.y_coords;
    double odz = O2.z_coords - O1.z_coords;
    double hdx = H.x_coords - O1.x_coords;
    double hdy = H.y_coords - O1.y_coords;
    double hdz = H.z_coords - O1.z_coords;
    
    odx -= info.lattice_x * pbc_round(odx/info.lattice_x);
    ody -= info.lattice_y * pbc_round(ody/info.lattice_y);
    odz -= info.lattice_z * pbc_round(odz/info.lattice_z);
    hdx -= info.lattice_x * pbc_round(hdx/info.lattice_x);
    hdy -= info.lattice_y * pbc_round(hdy/info.lattice_y);
    hdz -= info.lattice_z * pbc_round(hdz/info.lattice_z);
    
    double oodist = OOdist(O1, O2, info);
    double ohdist = OHdist(O1, H, info);
    
    double dot = odx*hdx + ody*hdy + odz*hdz;
    double angle = acos( dot / ( oodist*ohdist )) * rtd;
    
    return angle;
}

void nearest_neighbors(Information& info, Oxygen& O, O_vector& Ovec) {
    for (int i = 0; i < Ovec.size(); i ++) {
        Oxygen& O2 = Ovec[i];
        if (O.ID != O2.ID) {
            double dist = OOdist(O, O2, info);
            if (dist < 3.6) {
                O.nearest_neighbors.push_back(O2.ID);
            }
        }
    }
}

void find_local_H(Information& info, Oxygen& O, H_vector& Hvec) {
    for (int i = 0; i < Hvec.size(); i ++) {
        Hydrogen& H = Hvec[i]; 
        double dist = OHdist(O, H, info);
        if (dist < 1.2) {
            O.local_H_neighbors.push_back(H.ID);
        }
    }
}

void find_H_bonds(Oxygen& O, H_vector& Hvec, O_vector& Ovec, Information& info) {
    for (int i = 0; i < O.nearest_neighbors.size(); i ++) {
        Oxygen& O2 = Ovec[O.nearest_neighbors[i]];
        for (int j = 0; j < O2.local_H_neighbors.size(); j ++) {
            Hydrogen& H = Hvec[O2.local_H_neighbors[j]];
            double dist = OHdist(O, H, info);
            if (dist < 2.4) {
                O.bonded_H_neighbors.push_back(H.ID);
            }
        }
    }
}

void check_final_bonds(Oxygen& O, O_vector& Ovec, H_vector& Hvec, Information& info) {
    for (int i = 0; i < O.nearest_neighbors.size(); i ++) {
        Oxygen& O2 = Ovec[O.nearest_neighbors[i]];
        for (int j = 0; j < O.bonded_H_neighbors.size(); j ++) {
            Hydrogen& H = Hvec[O.bonded_H_neighbors[j]];
            double angle = angle_between(O, O2, H, info);
            if (angle < 30.0) {
                O.bonded_O_neighbors.push_back(O2.ID);
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
            nearest_neighbors(info, O, Ovec);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            find_local_H(info, O, Hvec);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            find_H_bonds(O, Hvec, Ovec, info);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            check_final_bonds(O, Ovec, Hvec, info);
        }
        
    }
    cout << "Computed nearest neighbour lists...\n\n";
    return true;
}