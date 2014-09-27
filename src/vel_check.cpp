#include "storage.h"
#include "analysis.h"
#include "helper.h"
#include "vel_check.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;


void  vel_check(Args& args) {
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    VHelp_Vector vhelper(time_steps.size() - 1);
    for (int i = 1; i < time_steps.size(); i ++) {
        O_vector& Ovec_then = time_steps[i].O_atoms;
        O_vector& Ovec_now = time_steps[i-1].O_atoms;
        H_vector& Hvec_then = time_steps[i].H_atoms;
        H_vector& Hvec_now = time_steps[i-1].H_atoms;
        for (int j = 0; j < Ovec_now.size(); j ++) {
            double dx = Ovec_then[j].x_coords - Ovec_now[j].x_coords;
            dx -= info.lattice_x * pbc_round(dx/info.lattice_x);
            double dy = Ovec_then[j].y_coords - Ovec_now[j].y_coords;
            dy -= info.lattice_y * pbc_round(dy/info.lattice_y);
            
            double dz = Ovec_then[j].z_coords - Ovec_now[j].z_coords;
            dz -= info.lattice_z * pbc_round(dz/info.lattice_z);
            
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            double velO = dist / info.time_step;
            vhelper[i].velocities.push_back(velO);
        }
        for (int j = 0; j < Hvec_now.size(); j ++) {
            double dx = Hvec_then[j].x_coords - Hvec_now[j].x_coords;
            dx -= info.lattice_x * pbc_round(dx/info.lattice_x);
            
            double dy = Hvec_then[j].y_coords - Hvec_now[j].y_coords;
            dy -= info.lattice_y * pbc_round(dy/info.lattice_y);
            
            double dz = Hvec_then[j].z_coords - Hvec_now[j].z_coords;
            dz -= info.lattice_z * pbc_round(dz/info.lattice_z);
            
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            double velH = dist / info.time_step;
            
            vhelper[i].velocities.push_back(velH);
        }
    }
    
    ofstream output;
    output.open("output/vel_test.xyz");
    for (int i = 0; i < vhelper.size(); i ++) {
        output << "Timestep: " << i << "\tAvg vel: " << average(vhelper[i].velocities) << "\n";
    }
    output.close();
}