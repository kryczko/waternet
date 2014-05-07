#include <fstream>
#include <iostream>
#include "storage.h"

using namespace std;

int check_time_steps(Information& info) {
    ifstream input;
    input.open(info.filename.c_str());
    string blah;
    int step_counter = 0;
    while(!input.eof()) {
        input >> blah;
        if (blah == "H") {
            step_counter ++;
        }
    }
    return step_counter / (double) info.num_hydrogen;
}

bool read_datafile(Information& info, TimeSteps& time_steps) {
    
    ifstream input;
    input.open(info.filename.c_str());
    int s = check_time_steps(info);
    
    // declare everything up front
    time_steps.resize(s);
    for (int i = 0; i < time_steps.size(); i ++) {
        time_steps[i].H_atoms.resize(info.num_hydrogen);
        time_steps[i].O_atoms.resize(info.num_oxygen);
        for (int j = 0; j < info.num_hydrogen; j ++) {
            time_steps[i].H_atoms[j].x_coords.resize(s);
            time_steps[i].H_atoms[j].y_coords.resize(s);
            time_steps[i].H_atoms[j].z_coords.resize(s);
        }
        for (int j = 0; j < info.num_oxygen; j ++) {
            time_steps[i].O_atoms[j].x_coords.resize(s);
            time_steps[i].O_atoms[j].y_coords.resize(s);
            time_steps[i].O_atoms[j].z_coords.resize(s);
        }
        
    }
    int steps = 0;
    int n_hyd = 0;
    int n_oxy = 0;
    string stuff;
    while(!input.eof()) {
        input >> stuff;
        H_vector& Hvec = time_steps[steps].H_atoms;
        O_vector& Ovec = time_steps[steps].O_atoms;
        
        if (stuff == "H") {
            Hydrogen& H = Hvec[n_hyd];
            H.ID = n_hyd;
            input >> H.x_coords[steps] >> H.y_coords[steps] >> H.z_coords[steps];
            n_hyd ++;
        } else if (stuff == "O") {
            Oxygen& O = Ovec[n_oxy];
            O.ID = n_oxy;
            input >> O.x_coords[steps] >> O.y_coords[steps] >> O.z_coords[steps];
            n_oxy ++;
        } else if (n_hyd == info.num_hydrogen && n_oxy == info.num_oxygen) {
            steps ++;
            n_hyd = 0;
            n_oxy = 0;
        } 
    }
    cout << "Read data file...\n\n";
    return true;
}