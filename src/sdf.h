#ifndef SDF_H_
#define SDF_H_

#include "storage.h"
#include "analysis.h"
#include "helper.h"

using namespace std;


void sdf(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    double Obins[info.sdf_bins][info.sdf_bins], Hbins[info.sdf_bins][info.sdf_bins];
    for (int i = 0; i < info.sdf_bins; i ++) {
        for (int j = 0; j < info.sdf_bins; j ++) {
            Obins[i][j] = 0;
            Hbins[i][j] = 0;
        }
    }
    int O_count = 0, H_count = 0;
    double xinc = info.lattice_x / info.sdf_bins, yinc = info.lattice_y / info.sdf_bins;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            if (O.z_coords < info.sdf_end && O.z_coords > info.sdf_start) {
                int xbin = O.x_coords / xinc;
                int ybin = O.y_coords / yinc;
                Obins[xbin][ybin] ++;
                O_count ++;
            }
        }
        for (int j = 0; j < Hvec.size(); j ++) {
            Hydrogen& H = Hvec[j];
            if (H.z_coords < info.sdf_end && H.z_coords > info.sdf_start) {
                int xbin = H.x_coords / xinc;
                int ybin = H.y_coords / yinc;
                Hbins[xbin][ybin] ++;
                H_count ++;
            }
        }
    }
    ofstream O_output;
    ofstream H_output;
    string Ofile = info.sdf_output + "_O.dat";
    string Hfile = info.sdf_output + "_H.dat";
    O_output.open(Ofile.c_str());
    H_output.open(Hfile.c_str());
    for (int i = 0; i < info.sdf_bins; i ++) {
        for (int j = 0; j < info.sdf_bins; j ++) {
            O_output << i*xinc << "  " << j*yinc << "  " << "  " <<  100.0 * Obins[i][j] / (double) O_count << "\n";
            H_output << i*xinc << "  " << j*yinc << "  " << "  " <<  100.0 * Hbins[i][j] / (double) H_count << "\n";

        }
    }

    O_output.close();
    H_output.close();
    cout << "Outputted spacial distribution function data files.\n\n";
    
    
}

#endif