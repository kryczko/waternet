#include "storage.h"
#include "analysis.h"
#include "helper.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>


using namespace std;

double HOHangle(Oxygen& O, Hydrogen& H1, Hydrogen& H2, Information& info) {
    double dx1 = H1.x_coords - O.x_coords;
    double dy1 = H1.y_coords - O.y_coords;
    double dz1 = H1.z_coords - O.z_coords;
    double dx2 = H2.x_coords - O.x_coords;
    double dy2 = H2.y_coords - O.y_coords;
    double dz2 = H2.z_coords - O.z_coords;
    dx1 -= info.lattice_x * pbc_round(dx1/info.lattice_x);
    dy1 -= info.lattice_y * pbc_round(dy1/info.lattice_y);
    dz1 -= info.lattice_z * pbc_round(dz1/info.lattice_z);
    dx2 -= info.lattice_x * pbc_round(dx2/info.lattice_x);
    dy2 -= info.lattice_y * pbc_round(dy2/info.lattice_y);
    dz2 -= info.lattice_z * pbc_round(dz2/info.lattice_z);
    
    double dist1 = OHdist(O, H1, info);
    double dist2 = OHdist(O, H2, info);
    
    double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
    double angle = acos ( dot / ( dist1*dist2 ) ) * rtd;
    return angle;
}

void  HOHdistro(Args& args) {
    // 360 degrees, bin for each degree
    Information info = args.arg_info;
    TimeSteps time_steps = args.arg_time_steps;
    
    vector<int> bins ( info.HOH_bins );
    double binsize = 360.0 / info.HOH_bins;
    int counter = 0;
    for (int i = 0; i < info.HOH_bins; i ++) {
        bins[i] = 0;
    }
    double sum = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            if (O.local_H_neighbors.size() == 2) {
                Hydrogen& H1 = Hvec[O.local_H_neighbors[0]];
                Hydrogen& H2 = Hvec[O.local_H_neighbors[1]];
                double angle = HOHangle(O, H1, H2, info);
                sum += angle;
                int bin = angle / binsize;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.HOH_output.c_str());
    output << "# Angle\tProbability\tAverage angle\n\n";
    for (int i = 0; i < info.HOH_bins; i ++) {
        output << i*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << (i + 1)*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    cout << "Outputting HOH angle data file.\n\n";
    
    
    
}