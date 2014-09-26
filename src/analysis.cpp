#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <pthread.h>

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

bool main_analysis(Information& info, TimeSteps& time_steps) {
    Args args;
    args.arg_info = info;
    args.arg_time_step = time_steps[0]; 
    args.arg_time_steps = time_steps;
   
    //THIS IS FOR OPENMP
     // vector of function pointers
    vector<void(*)(Args&)> vec_fp;
    // vector of decisions if functions are called
    vector<bool> dec;
    
    void (*fpointer)(Args&);
    
    fpointer = output_edgelist;
    vec_fp.push_back(fpointer);
    
    fpointer = output_graphfile;
    vec_fp.push_back(fpointer);
    dec.push_back(info.output_gephi);

    fpointer = degree_respect_z;
    vec_fp.push_back(fpointer);
    dec.push_back(info.degree_z);
    
    fpointer = OOdistro;
    vec_fp.push_back(fpointer);
    dec.push_back(info.OODistro);
    
    fpointer = OHdistro;
    vec_fp.push_back(fpointer);
    dec.push_back(info.OHDistro);
    
    fpointer = HOHdistro;
    vec_fp.push_back(fpointer);
    dec.push_back(info.HOHDistro);
    
    fpointer = degree_distro;
    vec_fp.push_back(fpointer);
    dec.push_back(info.degree_distro);
    
    fpointer = density;
    vec_fp.push_back(fpointer);
    dec.push_back(info.density);
    
    fpointer = msd;
    vec_fp.push_back(fpointer);
    dec.push_back(info.msd);
    
    fpointer = orientation;
    vec_fp.push_back(fpointer);
    dec.push_back(info.orientation);
    
    fpointer = orientation_1D;
    vec_fp.push_back(fpointer);
    dec.push_back(info.orientation_1D);
    
    fpointer = sdf;
    vec_fp.push_back(fpointer);
    dec.push_back(info.sdf);
    
    fpointer = nrt;
    vec_fp.push_back(fpointer);
    dec.push_back(info.network_reorganization_time);
   
    fpointer = vel_check;
    vec_fp.push_back(fpointer);
    dec.push_back(false); 
    
    fpointer = zdens_from_metal;
    vec_fp.push_back(fpointer);
    dec.push_back(true);
    
    fpointer = degree_respect_metal;
    vec_fp.push_back(fpointer);
    dec.push_back(true);
    /*if (info.output_gephi) {
        output_graphfile(info, time_steps[time_steps.size() - 1]);
        cout << "Gephi graph file created.\n\n";
    }
    if (info.degree_z) {
        degree_respect_z(info, time_steps);
        cout << "Outputted degrees with respect to the z-axis.\n\n";
    }
    if (info.OODistro) {
        OOdistro(info, time_steps);
        cout << "Outputted neighboring O-O distance distribution.\n\n";
    }
    if (info.OHDistro) {
        OHdistro(info, time_steps);
        cout << "Outputted local O-H distance distribution.\n\n";
    }
    if (info.HOHDistro) {
        HOHdistro(info, time_steps);
        cout << "Outputted local H-O-H angle distribution.\n\n";
    }
    if (info.degree_distro) {
        degree_distro(info, time_steps);
        cout << "Outputted cumulative degree distribution.\n\n";
    }
    if (info.density) {
        density(info, time_steps);
        cout << "Outputted density profile data.\n\n";
    }
    if (info.msd) {
        msd(info, time_steps);
        cout << "Outputted mean square displacement data.\n\n";
    } 
    if (info.orientation) {
        orientation(info, time_steps);
        cout << "Outputted 2D orientation data files.\n\n";
    }
    if (info.orientation_1D) {
        orientation_1D(info, time_steps);
        cout << "Outputted 1D orientation data files.\n\n";
    }
    if (info.sdf) {
        sdf(info, time_steps);
        cout << "Outputted spacial distribution function data files.\n\n";
    }
    if (info.network_reorganization_time) {
        nrt(info, time_steps);
        cout << "Outputted network reorganization time data file.\n\n";
    }*/
    
    // OPENMP parallelization
    //omp_set_dynamic(0);
    //omp_set_num_threads(info.num_threads);
    //#pragma omp parallel for 
    for (int i = 0; i < vec_fp.size(); i ++) {  
        if (dec[i]) {
            vec_fp[i](args);
        }
    }


    // pthread parallelization   
    /*Args args;
    args.arg_info = info;
    args.arg_time_step = time_steps[0]; 
    args.arg_time_steps = time_steps;
    pthread_t threads[vec_fp.size()];
    for (int i = 0; i < vec_fp.size(); i ++) {
        if (dec[i]) {
            int thread = pthread_create(&threads[i], NULL, vec_fp[i], (void *) &args);
            if (thread) {
                cout << "Error with thread.\n";
            }
        }
    }
    for (auto& thread : threads) {
        pthread_join(thread, NULL);
    }*/
    return true;
}
