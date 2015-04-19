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


void orientation_1D(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    double incr = (info.cell_block_end - info.cell_block_start) / info.num_cell_blocks;
    for (int cell = 0; cell < info.num_cell_blocks; cell ++) {
        double start_z = info.cell_block_start + cell*incr;
        double end_z = info.cell_block_start + (cell+1)*incr;
        
        vector<double> posx, posy, posz, vecx, vecy, vecz;
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            H_vector& Hvec = time_steps[i].H_atoms;
            for (int j = 0; j < Ovec.size(); j++) {
                Oxygen& O = Ovec[j];
                if (O.local_H_neighbors.size() == 2 && O.z_coords > start_z && O.z_coords < end_z) {
                    Hydrogen& H1 = Hvec[O.local_H_neighbors[0]];    
                    Hydrogen& H2 = Hvec[O.local_H_neighbors[1]];
                    double ohx1 = H1.x_coords - O.x_coords;
                    double ohy1 = H1.y_coords - O.y_coords;
                    double ohz1 = H1.z_coords - O.z_coords;
                    double ohx2 = H2.x_coords - O.x_coords;
                    double ohy2 = H2.y_coords - O.y_coords;
                    double ohz2 = H2.z_coords - O.z_coords;
                    ohx1 -= info.lattice_x * pbc_round(ohx1/info.lattice_x);
                    ohy1 -= info.lattice_y * pbc_round(ohy1/info.lattice_y);
                    ohz1 -= info.lattice_z * pbc_round(ohz1/info.lattice_z);
                    ohx2 -= info.lattice_x * pbc_round(ohx2/info.lattice_x);
                    ohy2 -= info.lattice_y * pbc_round(ohy2/info.lattice_y);
                    ohz2 -= info.lattice_z * pbc_round(ohz2/info.lattice_z);
                    vecx.push_back(ohx1 + ohx2);
                    vecy.push_back(ohy1 + ohy2);
                    vecz.push_back(ohz1 + ohz2);
                    posx.push_back(wrap(O.x_coords + (ohx1 + ohx2) / 2, info.lattice_x));
                    posy.push_back(wrap(O.y_coords + (ohy1 + ohy2) / 2, info.lattice_y));
                    posz.push_back(wrap(O.z_coords + (ohz1 + ohz2) / 2, info.lattice_z));
                }
            }
        }
        int angle_bin[180/3] = {};
        double count = 0;
        for (int i = 0; i < vecx.size(); i ++) {
            double dot = vecz[i];
            double norm = sqrt( vecx[i]*vecx[i] + vecy[i]*vecy[i] + vecz[i]*vecz[i] );
            double angle = acos ( dot / ( norm ) ) * rtd;
            angle_bin[(int) (angle/3.0)] ++;
            count ++;
        
        }
        ofstream output;
        string filename = "output/orientation_1D_region_" + to_string(cell) + ".dat";
        output.open(filename.c_str());
        output << "# 1D orientation for region " << cell << "\n";
        output << "# Region is from " << start_z << " to " << end_z << "\n\n";
        for (int i = 0 ; i < 180/3; i ++) {
            output << cos((double) i*3 / rtd) << "   " << angle_bin[i] / count << "\n";
        }
        output.close();
    }
    cout << "Outputted 1D orientation data files.\n\n";
}

void orientation(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    vector<double> posx, posy, posz, vecx, vecy, vecz;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j++) {
            Oxygen& O = Ovec[j];
            if (O.local_H_neighbors.size() == 2) {
                Hydrogen& H1 = Hvec[O.local_H_neighbors[0]];    
                Hydrogen& H2 = Hvec[O.local_H_neighbors[1]];
                double ohx1 = H1.x_coords - O.x_coords;
                double ohy1 = H1.y_coords - O.y_coords;
                double ohz1 = H1.z_coords - O.z_coords;
                double ohx2 = H2.x_coords - O.x_coords;
                double ohy2 = H2.y_coords - O.y_coords;
                double ohz2 = H2.z_coords - O.z_coords;
                ohx1 -= info.lattice_x * pbc_round(ohx1/info.lattice_x);
                ohy1 -= info.lattice_y * pbc_round(ohy1/info.lattice_y);
                ohz1 -= info.lattice_z * pbc_round(ohz1/info.lattice_z);
                ohx2 -= info.lattice_x * pbc_round(ohx2/info.lattice_x);
                ohy2 -= info.lattice_y * pbc_round(ohy2/info.lattice_y);
                ohz2 -= info.lattice_z * pbc_round(ohz2/info.lattice_z);
                vecx.push_back(ohx1 + ohx2);
                vecy.push_back(ohy1 + ohy2);
                vecz.push_back(ohz1 + ohz2);
                posx.push_back(wrap(O.x_coords + (ohx1 + ohx2) / 2, info.lattice_x));
                posy.push_back(wrap(O.y_coords + (ohy1 + ohy2) / 2, info.lattice_y));
                posz.push_back(wrap(O.z_coords + (ohz1 + ohz2) / 2, info.lattice_z));
            }
        }
    }
    double xzbins[info.orient_x_bins][info.orient_z_bins];
    int xzcounts[info.orient_x_bins][info.orient_z_bins];
    double yzbins[info.orient_y_bins][info.orient_z_bins];
    int yzcounts[info.orient_y_bins][info.orient_z_bins];
    double xybins[info.orient_x_bins][info.orient_y_bins];
    int xycounts[info.orient_x_bins][info.orient_y_bins];
    for (int i = 0; i < info.orient_x_bins; i ++) {
        for (int j = 0; j < info.orient_z_bins; j ++) {
            xzbins[i][j] = 0;
            xzcounts[i][j] = 0;
        }
    }
    for (int i = 0; i < info.orient_y_bins; i ++) {
        for (int j = 0; j < info.orient_z_bins; j ++) {
            yzbins[i][j] = 0;
            yzcounts[i][j] = 0;
        }
    }
    for (int i = 0; i < info.orient_x_bins; i ++) {
        for (int j = 0; j < info.orient_y_bins; j ++) {
            xybins[i][j] = 0;
            xycounts[i][j] = 0;
        }
    }
    double xinc = info.lattice_x / info.orient_x_bins;
    double yinc = info.lattice_y / info.orient_y_bins;
    double zinc = info.lattice_z / info.orient_z_bins;
    for (int i = 0; i < vecx.size(); i ++) {
        double dot = vecz[i];
        double norm = sqrt( vecx[i]*vecx[i] + vecy[i]*vecy[i] + vecz[i]*vecz[i] );
        double angle = acos ( dot / ( norm ) ) * rtd;
        
        int xbin = posx[i] / xinc;
        int ybin = posy[i] / yinc;
        int zbin = posz[i] / zinc;
        if (vecx[i] < 0) {
            angle += 180.0;
            xzbins[xbin][zbin] += angle;
            xzcounts[xbin][zbin] ++;
            angle -= 180.0;
        } else {
            xzbins[xbin][zbin] += angle;
            xzcounts[xbin][zbin] ++;
        } if (vecy[i] < 0) {
            angle += 180.0;
            yzbins[ybin][zbin] += angle;
            yzcounts[ybin][zbin] ++;
            angle -= 180.0;
        } else {
            yzbins[ybin][zbin] += angle;
            yzcounts[ybin][zbin] ++;
        }
        
        xybins[xbin][ybin] += angle;
        xycounts[xbin][ybin] ++;
    }
    ofstream outputxz, outputyz, outputxy;
    string xzfile = info.orientation_filename + "xz.dat";
    string yzfile = info.orientation_filename + "yz.dat";
    string xyfile = info.orientation_filename + "xy.dat";
    outputxz.open(xzfile.c_str());
    outputyz.open(yzfile.c_str());
    outputxy.open(xyfile.c_str());
    for (int i = 0; i < info.orient_x_bins; i ++) {
        for (int j = 0; j < info.orient_z_bins; j ++) {
            if (xzcounts[i][j]) {
                outputxz << i*xinc << "\t" << j*zinc << "\t" << xzbins[i][j] / xzcounts[i][j] << "\n";
            } else {
                outputxz << i*xinc << "\t" << j*zinc << "\t" << 0 << "\n";
            }
        }
    }
    for (int i = 0; i < info.orient_y_bins; i ++) {
        for (int j = 0; j < info.orient_z_bins; j ++) {
            if (yzcounts[i][j]) {
                outputyz << i*yinc << "\t" << j*zinc << "\t" << yzbins[i][j] / yzcounts[i][j] << "\n";
            } else {
                outputyz << i*yinc << "\t" << j*zinc << "\t" << 0 << "\n";
            }
        }
    }
    for (int i = 0; i < info.orient_x_bins; i ++) {
        for (int j = 0; j < info.orient_y_bins; j ++) {
            if (xycounts[i][j]) {
                outputxy << i*xinc << "\t" << j*yinc << "\t" << xybins[i][j] / xycounts[i][j] << "\n";
            } else {
                outputxy << i*xinc << "\t" << j*yinc << "\t" << 0 << "\n";
            }
            
        }
    }
    
    outputxz.close();
    outputyz.close();
    outputxz.close();
    cout << "Outputted 2D orientation data files.\n\n";
    
    
}
