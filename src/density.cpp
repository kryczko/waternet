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

void  zdens_from_metal(Information& info, TimeSteps& time_steps) {
    double average_left = metals_avg_left(time_steps);
    double average_right = metals_avg_right(time_steps, info);
    double metal_dist = average_right - average_left;
    double water_z_dist = info.lattice_z - metal_dist;
    int nbins = info.density_bins/2;
    vector<int> zbins(nbins), O_zbins(nbins), H_zbins(nbins);
    for (int i = 0; i < nbins; i ++) {
        zbins[i] = 0;
        O_zbins[i] = 0;
        H_zbins[i] = 0;
    }
    double zinc = 0.5 * water_z_dist / (nbins);
    double conversion, O_conversion = 1, H_conversion;
    if (info.heavy_water) {
        H_conversion = 2;
        conversion = 20.0e-6/(6.023e23*1.0e-30);
    } else {
        H_conversion = 1;
        conversion = 18.0e-6/(6.023e23*1.0e-30);
    }
    double zvol = info.lattice_x*info.lattice_y*zinc;
    ofstream zoutput;
    string filename = "output/dens_from_metal.dat";
    zoutput.open(filename.c_str());
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            double z = min(abs(O.z_coords - average_left), abs(O.z_coords - average_right));
            int zbin = z / zinc;
            zbins[zbin] ++;
            O_zbins[zbin] ++;
        }
        for (int j = 0; j < Hvec.size(); j ++) {
            Hydrogen& H = Hvec[j];
            double z = min(abs(H.z_coords - average_left), abs(H.z_coords - average_right));
            int zbin = z / zinc;
            H_zbins[zbin] ++;
        }
    }
    for (int i = 0; i < nbins; i ++) {
        zoutput << i*zinc << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\n";  
    }
    
    cout << "got here\n";
    zoutput.close();
    cout << "Outputted density with respect to metal data file.\n\n";
}


void  density(Args& args) {
    
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
    double incr = (info.cell_block_end - info.cell_block_start) / info.num_cell_blocks;
    for (int cell = 0; cell < info.num_cell_blocks; cell ++) {
        double start_z = info.cell_block_start + cell*incr;
        double end_z = info.cell_block_start + (cell+1)*incr;
        vector<int> xbins (info.density_bins), ybins (info.density_bins),
                    O_xbins(info.density_bins), O_ybins(info.density_bins),
                    H_xbins(info.density_bins), H_ybins(info.density_bins);
        for (int i = 0; i < info.density_bins; i ++) {
            xbins[i] = 0;
            ybins[i] = 0;
            O_xbins[i] = 0;
            O_ybins[i] = 0;
            H_xbins[i] = 0;
            H_ybins[i] = 0;
        
        }
        double xbinsize = info.lattice_x / info.density_bins;
        double ybinsize = info.lattice_y / info.density_bins; 
        
        double conversion, O_conversion = 1, H_conversion;
        if (info.heavy_water) {
            H_conversion = 2;
            conversion = 20.0e-6/(6.023e23*1.0e-30);
        } else {
            H_conversion = 1;
            conversion = 18.0e-6/(6.023e23*1.0e-30);
        }
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            H_vector& Hvec = time_steps[i].H_atoms;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                if (O.z_coords > start_z && O.z_coords < end_z) {
                    double x = wrap(O.x_coords, info.lattice_x);
                    double y = wrap(O.y_coords, info.lattice_y);
                    int xbin = x / xbinsize;
                    int ybin = y / ybinsize;
                    xbins[xbin] ++;
                    ybins[ybin] ++;
                    O_xbins[xbin] ++;
                    O_ybins[ybin] ++;
                }
            }
            for (int j = 0; j < Hvec.size(); j ++) {
                Hydrogen& H = Hvec[j];
                if (H.z_coords > start_z && H.z_coords < end_z) {
                    double x = wrap(H.x_coords, info.lattice_x);
                    double y = wrap(H.y_coords, info.lattice_y);
                    int xbin = x / xbinsize;
                    int ybin = y / ybinsize;
                    H_xbins[xbin] ++;
                    H_ybins[ybin] ++;
                }
            }
        }
        
        double xvol = xbinsize*info.lattice_y*(end_z-start_z), yvol = ybinsize*info.lattice_x*(end_z - start_z);
        double xsum = 0, ysum = 0, O_xsum = 0, O_ysum = 0, H_xsum = 0, H_ysum = 0;
        int xcount = 0, ycount = 0, O_xcount = 0, O_ycount = 0, H_xcount = 0, H_ycount = 0;
        for (int i = 0; i < info.density_bins; i ++ ) {
            xsum += xbins[i]*conversion / (xvol*info.n_frames);
            ysum += ybins[i]*conversion / (yvol*info.n_frames);
            O_xsum += O_xbins[i]*O_conversion / (xvol*info.n_frames);
            O_ysum += O_ybins[i]*O_conversion / (yvol*info.n_frames);
            H_xsum += H_xbins[i]*H_conversion / (xvol*info.n_frames);
            H_ysum += H_ybins[i]*H_conversion / (yvol*info.n_frames);
            if (xbins[i]) {
                xcount ++;
            }
            if (ybins[i]) {
                ycount ++;
            }
            if (O_xbins[i]) {
                O_xcount ++;
            }
            if (O_ybins[i]) {
                O_ycount ++;
            }
            if (H_xbins[i]) {
                H_xcount ++;
            }
            if (H_ybins[i]) {
                H_ycount ++;
            }
        }
        ofstream xoutput, youtput;
        string xoutputfile = info.xdens_output + "_region_" + to_string(cell) + ".dat";
        string youtputfile = info.ydens_output + "_region_" + to_string(cell) + ".dat";
        
        xoutput.open(xoutputfile.c_str());
        youtput.open(youtputfile.c_str());
        
        for (int i = 0; i < info.density_bins; i ++) {
            xoutput << i*xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*O_conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*H_conversion / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            xoutput << i*xbinsize + xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*O_conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*H_conversion / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            youtput << i*ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*O_conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*H_conversion / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
            youtput << i*ybinsize + ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*O_conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*H_conversion / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
        }
        xoutput.close();
        youtput.close();
    }
    vector<int> zbins (info.density_bins), O_zbins(info.density_bins), H_zbins(info.density_bins);
    for (int i = 0; i < info.density_bins; i ++) {
        zbins[i] = 0;
        O_zbins[i] = 0; 
        H_zbins[i] = 0;
        
    }
    double zbinsize = info.lattice_z / info.density_bins;
    
    double conversion, O_conversion = 1, H_conversion;
    if (info.heavy_water) {
        H_conversion = 2;
        conversion = 20.0e-6/(6.023e23*1.0e-30);
    } else {
        H_conversion = 1;
        conversion = 18.0e-6/(6.023e23*1.0e-30);
    }
    
    
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            double z = wrap(O.z_coords, info.lattice_z);
            int zbin = z / zbinsize;
            zbins[zbin] ++;
            O_zbins[zbin] ++;
        }
        for (int j = 0; j < Hvec.size(); j ++) {
            Hydrogen& H = Hvec[j];
            double z = wrap(H.z_coords, info.lattice_z);
            int zbin = z / zbinsize;
            H_zbins[zbin] ++;
        }
    }
    double zvol = zbinsize*info.lattice_x*info.lattice_y;
    double zsum = 0, O_zsum = 0, H_zsum = 0;
    int zcount = 0, O_zcount = 0, H_zcount = 0;
    for (int i = 0; i < info.density_bins; i ++ ) {
        zsum += zbins[i]*conversion / (zvol*info.n_frames);
        O_zsum += O_zbins[i]*O_conversion / (zvol*info.n_frames);
        H_zsum += H_zbins[i]*H_conversion / (zvol*info.n_frames);
        if (zbins[i]) {
            zcount ++;
        }
        if (O_zbins[i]) {
            O_zcount ++;
        }
        if (H_zbins[i]) {
            H_zcount ++;
        }
    }
    vector<double> z_dens;
    vector<double> where;
    ofstream xoutput, youtput, zoutput;
    xoutput.open(info.xdens_output.c_str());
    youtput.open(info.ydens_output.c_str());
    zoutput.open(info.zdens_output.c_str());
    if (info.fix_plots) {
        for (int i = 0; i < info.density_bins; i ++) {
            if (i*zbinsize > info.starting_z) {
                if (zbins[i]) {
                    z_dens.push_back(zbins[i]*conversion / (zvol*info.n_frames));
                    where.push_back(i*zbinsize);    
                }
                
                zoutput << i*zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
                zoutput << i*zbinsize + zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";   
            }
        }
        for (int i = 0; i < info.density_bins; i ++) {
            if (i*zbinsize < info.starting_z) {
                if (zbins[i]) {
                
                z_dens.push_back(zbins[i]*conversion / (zvol*info.n_frames));
                where.push_back(i*zbinsize + info.lattice_z);
            }
                zoutput << i*zbinsize + info.lattice_z << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
                zoutput << i*zbinsize + zbinsize + info.lattice_z << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";   
            }
        }
    }
    else {
        for (int i = 0; i < info.density_bins; i ++) {
            zoutput << i*zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
            zoutput << i*zbinsize + zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*O_conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*H_conversion / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
        }
    }
    zoutput.close();    
    cout << "Outputted density data files.\n\n";
    
    
}
