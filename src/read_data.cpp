#include <fstream>
#include <iostream>
#include "storage.h"

using namespace std;

bool read_datafile(Information& info, TimeSteps& time_steps) {
    
    ifstream input;
    input.open(info.filename.c_str());
    int steps = 0;
    int n_hyd = 0;
    int n_oxy = 0;
    string stuff;
    while(!input.eof()) {
        input >> stuff;
        time_steps.resize(steps + 1);
        H_vector& Hvec = time_steps[steps].H_atoms;
        O_vector& Ovec = time_steps[steps].O_atoms;
        Hvec.resize(info.num_hydrogen);
        Ovec.resize(info.num_oxygen);
        
        if (stuff == "H") {
            Hydrogen& H = Hvec[n_hyd];
            H.ID = n_hyd;
            double x,y,z;
            input >> x >> y >> z;
            H.x_coords.push_back(x);
            H.y_coords.push_back(y);
            H.z_coords.push_back(z);
            n_hyd ++;
        } else if (stuff == "O") {
            O_vector& Ovec = time_steps[steps].O_atoms;
            Oxygen& O = Ovec[n_oxy];
            O.ID = n_oxy;
            double x,y,z;
            input >> x >> y >> z;
            O.x_coords.push_back(x);
            O.y_coords.push_back(y);
            O.z_coords.push_back(z);
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