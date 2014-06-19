#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "storage.h"
#include "edgelist.h"
#include "dynamics.h"

using namespace std;

double average(vector<int>& myvec) {
    double sum = 0;
    for (int i = 0; i < myvec.size(); i ++) {
        sum += myvec[i];
    }
    return sum / myvec.size();
}

bool H_group_dynamics(Information& info, TimeSteps& time_steps, H_group_info& hgi) {
    ofstream trajectory;
    trajectory.open(info.H_group_trajectory_filename.c_str());
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            if (Ovec[j].local_H_neighbors.size() == 3) {
                
                if ( hgi.last_oxygen == -1 ) {
                    hgi.last_oxygen = j; 
                    hgi.last_timestep = i;
                } else if ( i > hgi.last_timestep ) {
                    if ( j == hgi.last_oxygen ) {
                        hgi.counter ++;
                    } else {
                        hgi.counts.push_back(hgi.counter);
                        hgi.Os.push_back(hgi.last_oxygen);
                        hgi.counter = 0;
                        hgi.last_oxygen = j;
                        hgi.last_timestep = i;
                    }
                }
                
                
                trajectory << 4 << "\n\n";
                trajectory << "O" << "\t" << Ovec[j].x_coords << "\t" << Ovec[j].y_coords << "\t" << Ovec[j].z_coords << "\n";
                trajectory << "H" << "\t" << Hvec[Ovec[j].local_H_neighbors[0]].x_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[0]].y_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[0]].z_coords << "\n";
                trajectory << "H" << "\t" << Hvec[Ovec[j].local_H_neighbors[1]].x_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[1]].y_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[1]].z_coords << "\n";
                trajectory << "H" << "\t" << Hvec[Ovec[j].local_H_neighbors[2]].x_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[2]].y_coords << "\t" << Hvec[Ovec[j].local_H_neighbors[2]].z_coords << "\n";
            } 
        }
    }
    ofstream test;
    test.open("test_time.dat");
    for (int i = 0; i < hgi.Os.size(); i ++) {
        test << hgi.Os[i] << "\t" << hgi.counts[i] << "\n";
    }
    
    cout << "AVERAGE H-GROUP LIFETIME = " << average(hgi.counts) << "\n\n";
    test.close();
    trajectory.close();
    return true;
}

bool OH_group_dynamics(Information& info, TimeSteps& time_steps) {
    
    
    return true;
}

bool main_dynamics_func(Information& info, TimeSteps& time_steps, H_group_info& hgi) {
    if (info.H_group_dynamics) {
        H_group_dynamics(info, time_steps, hgi);
        cout << "H-group dynamic calculations have been done.\n\n";
    }
    else if (info.OH_group_dynamics) {
        OH_group_dynamics(info, time_steps);
        cout << "OH-group dynamic calculations have been done.\n\n";
    }
    return true;
}