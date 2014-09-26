#include <fstream>
#include <iostream>
#include "storage.h"

using namespace std;

int check_time_steps(Information& info) {
    ifstream input;
    input.open(info.input_filename.c_str());
    string blah;
    int step_counter = 0;
    while(!input.eof()) {
        input >> blah;
        if (blah == "H") {
            step_counter ++;
        }
    }
    cout << step_counter / (double) info.num_hydrogen << " time steps found.\n\n";
    info.n_frames = step_counter / (double) info.num_hydrogen;
    return step_counter / (double) info.num_hydrogen;
}

bool read_datafile(Information& info, TimeSteps& time_steps) {
    
    ifstream input;
    input.open(info.input_filename.c_str());

    int s = check_time_steps(info);
    
    // declare everything up front
    time_steps.resize(s);
    for (int i = 0; i < time_steps.size(); i ++) {
        time_steps[i].H_atoms.resize(info.num_hydrogen);
        time_steps[i].O_atoms.resize(info.num_oxygen);
        time_steps[i].M_atoms.resize(info.num_metals);
    }
    
    int steps = 0;
    int n_hyd = 0;
    int n_oxy = 0;
    int n_met = 0;
    
    string stuff;
    
    while(!input.eof()) {
        input >> stuff;
        H_vector& Hvec = time_steps[steps].H_atoms;
        O_vector& Ovec = time_steps[steps].O_atoms;
        M_vector& Mvec = time_steps[steps].M_atoms;
        
        if (stuff == "H") {
            Hydrogen& H = Hvec[n_hyd];
            H.ID = n_hyd;
            input >> H.x_coords >> H.y_coords >> H.z_coords;
            n_hyd ++;
        } else if (stuff == "O") {
            Oxygen& O = Ovec[n_oxy];
            O.ID = n_oxy;
            input >> O.x_coords >> O.y_coords >> O.z_coords;
            n_oxy ++;
        } else if (stuff == info.metal) {
            Metal& M = Mvec[n_met];
            M.ID = n_met;
            M.name = info.metal;
            input >> M.x_coords >> M.y_coords >> M.z_coords;
            n_met ++;
        } else if (n_hyd == info.num_hydrogen && n_oxy == info.num_oxygen && n_met == info.num_metals ) {
            steps ++;
            n_hyd = 0;
            n_oxy = 0;
            n_met = 0;
        } 
    }
    cout << "Read data file...\n\n";
    return true;
}

