#include "storage.h"
#include "analysis.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "helper.h"

using namespace std;

double metals_avg_right(TimeSteps& time_steps, Information& info) {
    double sum = 0;
    double counter = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        double max_m = find_max_m(time_steps[i].M_atoms, info);
        for (auto& M : time_steps[i].M_atoms) {
            if (abs(M.z_coords - max_m) < 1.0) {
                sum += M.z_coords;
                counter += 1;
            }
        }
    }
    
    return sum / counter;
} 

double metals_avg_left(TimeSteps& time_steps, Information& info) {
    double sum = 0;
    double counter = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        double min_m = find_min_m(time_steps[i].M_atoms,info);
        for (auto& M : time_steps[i].M_atoms) {
            if (abs(M.z_coords - min_m) < 1.0) {
                sum += M.z_coords;
                counter += 1;
            }
        }
    }
    
    return sum / counter;
} 

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

double average(vector<double>& vec) {
    double sum = 0;
    for (double& elem : vec) {
        sum += elem;
    }
    return sum / vec.size(); 
}

vector<int> set_zero(vector<int>& vec) {
    for (int i = 0; i < vec.size(); i ++) {
        vec[i] = 0;
    }
    return vec;
}

int max(vector<int> vec) {
    int max = 0;
    for (int i = 0; i < vec.size(); i ++) {
        if (vec[i] > max) {
            max = vec[i];
        }
    }
    return max;
}

int num_edges(O_vector& Ovec) {
    int count = 0;
    for (int j = 0; j < Ovec.size(); j ++) {
        for (int k = 0; k < Ovec[j].bonded_O_neighbors.size(); k ++) {
            count ++;
        }
    }
    return count;
}

void out_count(Information& info, TimeSteps& time_steps) {
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.bonded_O_neighbors.size(); k++) {
                Oxygen& O2 = Ovec[O.bonded_O_neighbors[k]];
                O2.out_degree ++;
            }
        }
    }
}

double wrap(double value, double lattice) {
    if (value > lattice) {
        return value - lattice;
    } else if ( value < 0.0) {
        return value + lattice;
    } 
    return value;
}

void remove_zeros(vector<double>& dens, vector<double>& pos) {
    for (int i = 0; i < dens.size(); i ++) {
        if (dens[i] < 0.1) {
            dens.erase(dens.begin()+i);
            pos.erase(pos.begin()+i);
        }
    }
}

double find_max_m(M_vector& M_atoms, Information& info) {
    double val = 0;
    double dummy;
    for (auto& M : M_atoms) {
        if (M.z_coords < 5.0 ) {
            dummy = M.z_coords + info.lattice_z;
        }
        if (dummy > val) {
            val = dummy;
        }
    }
    return val;
}

double find_min_m(M_vector& M_atoms, Information& info) {
    double val = 500; // some large value
    double dummy;
    for (auto& M : M_atoms) {
        if (M.z_coords < 5.0) {
            dummy = M.z_coords + info.lattice_z;
        } 
        if (dummy < val) {
            val = dummy;
        }
    }
    return val;
}

double min(double x, double y) {
    if (x < y) {
        return x;
    } 
    return y;
}

void unwrap(Information& info, TimeSteps& time_steps) {
    for (int i = 1; i < time_steps.size(); i ++) {
        O_vector& Ovec0 = time_steps[i - 1].O_atoms;
        O_vector& Ovec1 = time_steps[i].O_atoms;
        H_vector& Hvec0 = time_steps[i - 1].H_atoms;
        H_vector& Hvec1 = time_steps[i].H_atoms;
        if (i == 1) {
            for (int j = 0; j < Ovec0.size(); j ++) {
                Ovec0[j].unwrap_x = Ovec0[j].x_coords;
                Ovec0[j].unwrap_y = Ovec0[j].y_coords;
                Ovec0[j].unwrap_z = Ovec0[j].z_coords;
            }
            for (int j = 0; j < Hvec0.size(); j ++) {
                Hvec0[j].unwrap_x = Hvec0[j].x_coords;
                Hvec0[j].unwrap_y = Hvec0[j].y_coords;
                Hvec0[j].unwrap_z = Hvec0[j].z_coords;
            }
        }
        for (int j = 0; j < Ovec0.size() ; j ++) {
            if (Ovec1[j].x_coords - Ovec0[j].unwrap_x > (info.lattice_x / 2.0)) {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords - info.lattice_x;
            } else if (Ovec1[j].x_coords - Ovec0[j].unwrap_x < (-info.lattice_x / 2.0)) {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords + info.lattice_x;
            } else {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords;
            }
            if (Ovec1[j].y_coords - Ovec0[j].unwrap_y > (info.lattice_y / 2.0)) {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords - info.lattice_y;
            } else if (Ovec1[j].y_coords - Ovec0[j].unwrap_y < (-info.lattice_y / 2.0)) {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords + info.lattice_y;
            } else {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords;
            }
            if (Ovec1[j].z_coords - Ovec0[j].unwrap_z > (info.lattice_z / 2.0)) {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords - info.lattice_z;
            } else if (Ovec1[j].z_coords - Ovec0[j].unwrap_z < (-info.lattice_z / 2.0)) {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords + info.lattice_z;
            } else {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords;
            }
        }
        for (int j = 0; j < Hvec0.size() ; j ++) {
            if (Hvec1[j].x_coords - Hvec0[j].unwrap_x > (info.lattice_x / 2.0)) {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords - info.lattice_x;
            } else if (Hvec1[j].x_coords - Hvec0[j].unwrap_x < (-info.lattice_x / 2.0)) {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords + info.lattice_x;
            } else {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords;
            }
            if (Hvec1[j].y_coords - Hvec0[j].unwrap_y > info.lattice_y / 2.0) {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords - info.lattice_y;
            } else if (Hvec1[j].y_coords - Hvec0[j].unwrap_y < -info.lattice_y / 2.0) {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords + info.lattice_y;
            } else {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords;
            }
            if (Hvec1[j].z_coords - Hvec0[j].unwrap_z > info.lattice_z / 2.0) {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords - info.lattice_z;
            } else if (Hvec1[j].z_coords - Hvec0[j].unwrap_z < -info.lattice_z / 2.0) {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords + info.lattice_z;
            } else {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords;
            }
        }
    }
    if (info.write_unwrapped_xyz) {
        ofstream output;
        cout << "Writing unwrapped coordinates to: " << info.unwrapped_coords << "\n\n";
        output.open(info.unwrapped_coords.c_str());
        for (int i = 0; i < time_steps.size(); i++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            H_vector& Hvec = time_steps[i].H_atoms;
            output << info.num_oxygen + info.num_hydrogen << "\n\n";
            for (int j = 0; j < Ovec.size(); j ++) {
                output << "O" << "\t" << Ovec[j].unwrap_x << "\t" << Ovec[j].unwrap_y << "\t" << Ovec[j].unwrap_z << "\n";
            }
            for (int j = 0; j < Hvec.size(); j ++) {
                output << "H" << "\t" << Hvec[j].unwrap_x << "\t" << Hvec[j].unwrap_y << "\t" << Hvec[j].unwrap_z << "\n";
            }
        }
        cout << "Finished writing unwrapped coordinates.\n\n";  
        output.close();
    }
}
