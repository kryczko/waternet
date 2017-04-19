#include "storage.h"
#include "analysis.h"
#include "helper.h"
#include "msd.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

msd_vec msd_vector;

void  msd(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    unwrap(info, time_steps);
    double incr = (info.cell_block_end - info.cell_block_start) / info.num_cell_blocks;
    for (int cell = 0; cell < info.num_cell_blocks; cell ++) {
        double start_z = info.cell_block_start + cell*incr;
        double end_z = info.cell_block_start + (cell+1)*incr;
    
        msd_vector.resize(info.num_blocks);
        for (int i = 0; i < msd_vector.size(); i ++) {
            msd_vector[i].declare(time_steps.size());
        }
    
        ofstream output;
        string filename = info.msd_filename + "region_" + to_string(cell) + ".dat"; 
        output.open(filename.c_str());
        int length = time_steps.size() / info.num_blocks;
    
        vector<double> COMx, COMy, COMz;
        COMx.resize(0);
        COMy.resize(0);
        COMz.resize(0);
    
    
        for (int len = 0; len < time_steps.size(); len ++) {
            int step = len;
            O_vector& Ovec = time_steps[step].O_atoms;
            H_vector& Hvec = time_steps[step].H_atoms;
        
            double dxo = 0, dyo = 0, dzo = 0, dxh = 0, dyh = 0, dzh = 0;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                dxo += O.unwrap_x;
                dyo += O.unwrap_y;
                dzo += O.unwrap_z;
            }
            for (int j = 0; j < Hvec.size(); j ++) {
                Hydrogen& H = Hvec[j];
                dxh += H.unwrap_x;
                dyh += H.unwrap_y;
                dzh += H.unwrap_z;
            }
            if (info.heavy_water) {
                double comx = (16*dxo + 2*dxh) / (16*Ovec.size() + 2*Hvec.size());
                COMx.push_back(comx);
                double comy = (16*dyo + 2*dyh) / (16*Ovec.size() + 2*Hvec.size());
                COMy.push_back(comy);
                double comz = (16*dzo + 2*dzh) / (16*Ovec.size() + 2*Hvec.size());
                COMz.push_back(comz);
            } else {
                double comx = (16*dxo + dxh) / (16*Ovec.size() + Hvec.size());
                COMx.push_back(comx);
                double comy = (16*dyo + dyh) / (16*Ovec.size() + Hvec.size());
                COMy.push_back(comy);
                double comz = (16*dzo + dzh) / (16*Ovec.size() + Hvec.size());
                COMz.push_back(comz);
            }
        }
        cout << "--- MSD progress ---\n";
        int print_count = 0;
        int starting_step, final_step;
        for (int nb = 0; nb < info.num_blocks; nb ++) {
            starting_step = nb*length;
            final_step = starting_step + (int) info.msd_time / info.time_step;
            // if (info.full_msd) {
            //     final_step = time_steps.size();
            // } else {
            //     final_step = (nb + 1)*length;
            // }
            O_vector& Ovec1 = time_steps[starting_step].O_atoms;
            H_vector& Hvec1 = time_steps[starting_step].H_atoms;
            vector<double>& msd_data = msd_vector[nb].msd_step;
            for (int len = starting_step; len < final_step; len ++) {
                int step = len;
                int ocounter = 0, hcounter = 0;
                O_vector& Ovec = time_steps[step].O_atoms;
                H_vector& Hvec = time_steps[step].H_atoms;
                for (int j = 0; j < Ovec.size(); j ++) {
                    Oxygen& O1 = Ovec1[j];
                    Oxygen& O2 = Ovec[j];
                    if (O1.unwrap_z > start_z && O1.unwrap_z < end_z) {
                        double dx1 = O1.unwrap_x - COMx[starting_step];
                        double dy1 = O1.unwrap_y - COMy[starting_step];
                        double dz1 = O1.unwrap_z - COMz[starting_step];
                        double dx2 = O2.unwrap_x - COMx[step];
                        double dy2 = O2.unwrap_y - COMy[step];
                        double dz2 = O2.unwrap_z - COMz[step];
                        double dx = dx2 - dx1;
                        double dy = dy2 - dy1;
                        double dz = dz2 - dz1;
                        ocounter ++;
                        msd_data[step - starting_step] += dx*dx + dy*dy + dz*dz;
                    }
                }
                for (int j = 0; j < Hvec.size(); j ++) {
                    Hydrogen& H1 = Hvec1[j];
                    Hydrogen& H2 = Hvec[j];
                    if (H1.unwrap_z > start_z && H1.unwrap_z < end_z) {
                    
                        double dx1 = H1.unwrap_x - COMx[starting_step];
                        double dy1 = H1.unwrap_y - COMy[starting_step];
                        double dz1 = H1.unwrap_z - COMz[starting_step];   
                        double dx2 = H2.unwrap_x - COMx[step];
                        double dy2 = H2.unwrap_y - COMy[step];
                        double dz2 = H2.unwrap_z - COMz[step];
                        double dx = dx2 - dx1;
                        double dy = dy2 - dy1;
                        double dz = dz2 - dz1;
                        hcounter ++;
                        msd_data[step - starting_step] += dx*dx + dy*dy + dz*dz;
                    }
                }
                msd_data[step - starting_step] /= (double) (ocounter + hcounter) ;
            }
            if ( nb < info.num_blocks) {
                cout << "Cell block: " << cell + 1<<  " of " << info.num_cell_blocks << " --- |< " << (int) ( 100 * (double) nb / (double) info.num_blocks) << "% >| ---\r";
                flush(cout);
            } else if (nb == info.num_blocks - 1) {
                cout << "Cell block: " << cell  + 1<<  " of " << info.num_cell_blocks << " --- |< 100% >| ---\n";
                flush(cout);
            }
            
        }
        cout << "Cell block: " << cell + 1<<  " of " << info.num_cell_blocks << " --- |< 100% >| ---\n\n";
        vector<double> averaged_msd(0);
        vector<int> counts(0);
        // if (info.full_msd) {
        //     for (int i = 0; i < time_steps.size(); i ++) {
        //         averaged_msd.push_back(0);
        //         counts.push_back(0);
        //     }
        // } else {
        for (int i = 0; i < length; i ++) {
            averaged_msd.push_back(0);
            counts.push_back(0);
        }
        // }
        for (int i = 0; i < info.num_blocks; i ++) {
                vector<double>& msd_data = msd_vector[i].msd_step;
                for (int j = 0; j < msd_data.size(); j ++) {
                    averaged_msd[j] += msd_data[j];
                    if (msd_data[j] != 0) {
                        counts[j] ++;
                    }
                }
        }
        averaged_msd[0] = 0;
        double slope_sum = 0;
        double slope_count = 0;
        for (int i = 1; i < averaged_msd.size() - 1; i ++) {
            double time_now = i*info.time_step / 1000.0;
            double time_then = (i-1)*info.time_step / 1000.0;
            averaged_msd[i] /= (double) counts[i];
            double msd_now = averaged_msd[i];
            double msd_then = averaged_msd[i-1];
            double slope = (msd_now - msd_then) / (time_now - time_then);
            slope_sum += slope;
            slope_count ++;
        }
        double dc = slope_sum / (1e8*6.0*slope_count), exp_dc = 2.30e-09;
        output << "# Self diffusion coefficient for region " << cell << ": " << slope_sum / (1e8*6.0*slope_count) << " m^2 / s\n";
        output << "# Experimental diffusion coefficient: 2.30e-09 m^2 /s, Percent difference: " <<  100.0 * abs(dc - exp_dc) / exp_dc << "%\n";
        output << "# Region: " << start_z << " to " << end_z << "[Angstroms] \n\n";
        for (int i = 1; i < averaged_msd.size() - 1; i ++) {
            output << i*info.time_step / 1000.0 << "\t" << averaged_msd[i] << "\n";
        }
        output.close();
        msd_vector.clear();
    }
    cout << "Outputted mean square displacement data.\n\n";
    
    
}