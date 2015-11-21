#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <pthread.h>
#include <vector>
//#include "omp.h"

#include "storage.h"
#include "edgelist.h"
#include "helper.h"
#include "gephifile.h"
#include "degree_z.h"
#include "OOdistro.h"
#include "OHdistro.h"
#include "HOHdistro.h"
#include "degree_distro.h"
#include "density.h"
#include "msd.h"
#include "orientation.h"
#include "sdf.h"
#include "nrt.h"
#include "vel_check.h"
#include "output_edgelist.h"

using namespace std;

void ODownHDown(Args& args) {
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    double average_left = args.avg_left + info.lattice_z;
    double average_right = args.avg_right;
    double metal_dist = average_right - (average_left - info.lattice_z);
    double max_dist = 3.0; // Static for now; in Angstroms
    double distro[2] = {0}; // 0 == O, 1 == H
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int k = 0; k < Ovec.size(); k ++) {
            Oxygen& O = Ovec[k];
            double dist1 = abs(O.z_coords - average_left);
            double dist2 = abs(O.z_coords - average_right);
            double dist_from_closest_metal = min(dist1, dist2);
            if (dist_from_closest_metal < max_dist) {
                vector<int> bonded_H = Ovec[k].bonded_H_neighbors;
                for (int j = 0; j < bonded_H.size(); j ++) {
                    Hydrogen& H = Hvec[bonded_H[j]];
                    if (O.z_coords < H.z_coords) {
                        distro[0] ++;
                    } else {
                        distro[1] ++;
                    }
                }
            }
        }
    }
    cout << distro[0] << "\t" << distro[1] << "\n";
}

bool main_analysis(Information& info, TimeSteps& time_steps) {
    Args args;
    args.arg_info = info;
    args.arg_time_step = time_steps[0]; 
    args.arg_time_steps = time_steps;
   
    void (*function_ptr)(Args&);
    vector<void (*)(Args&)> function_ptrs;
    
    
    
    if (info.create_edgelist) {
        function_ptr = &output_edgelist;
        function_ptrs.push_back(function_ptr);
    }
    
   
    if (info.output_gephi) {
        function_ptr = &output_graphfile;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.degree_z) {
        function_ptr = &degree_respect_z;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.OODistro) {
        function_ptr = &OOdistro;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.OHDistro) {
        function_ptr = &OHdistro;
        function_ptrs.push_back(function_ptr);
    }
    if (info.HOHDistro) {
        function_ptr = &HOHdistro;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.degree_distro) {
        function_ptr = &degree_distro;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.density) {
        function_ptr = &density;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.msd) {
        function_ptr = &msd;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.orientation) {
        function_ptr = &orientation;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.orientation_1D) {
        function_ptr = &orientation_1D;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.sdf) {
        function_ptr = &sdf;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.network_reorganization_time) {
        function_ptr = &nrt;
        function_ptrs.push_back(function_ptr);
    }
    
    if (info.network_reorganization_time) {
        function_ptr = &nrt;
        function_ptrs.push_back(function_ptr);
    }
    
    args.avg_left = metals_avg_left(time_steps, info);
    args.avg_right = metals_avg_right(time_steps, info);
    
    if (info.degree_z_from_metal) {
        function_ptr = &degree_respect_metal;
        function_ptrs.push_back(function_ptr);
    }

    if (info.density_from_metal) {
        function_ptr = &zdens_from_metal;
        function_ptrs.push_back(function_ptr);
    }    
    //#pragma omp parallel for
    for (int i = 0; i < function_ptrs.size(); i ++) {
        function_ptrs[i](args);
    }
    ODownHDown(args);

    return true;
}
