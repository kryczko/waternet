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
   
    if (info.create_edgelist) {
        output_edgelist(args);
    }
    
    if (info.output_gephi) {
        output_graphfile(args);
    }
    if (info.degree_z) {
        degree_respect_z(args);
    }
    if (info.OODistro) {
        OOdistro(args);
    }
    if (info.OHDistro) {
        OHdistro(args);
    }
    if (info.HOHDistro) {
        HOHdistro(args);
    }
    if (info.degree_distro) {
        degree_distro(args);
    }
    if (info.density) {
        density(args);
    }
    if (info.msd) {
        msd(args);
    } 
    if (info.orientation) {
        orientation(args);
    }
    if (info.orientation_1D) {
        orientation_1D(args);
    }
    if (info.sdf) {
        sdf(args);
    }
    if (info.network_reorganization_time) {
        nrt(args);
    }
    
    zdens_from_metal(args);
    degree_respect_metal(args);
    
    

    return true;
}
