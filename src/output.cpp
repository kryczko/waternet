#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "storage.h"
#include "edgelist.h"

using namespace std;

int num_edges(O_vector& Ovec) {
    int count = 0;
    for (int j = 0; j < Ovec.size(); j ++) {
        for (int k = 0; k < Ovec[j].bonded_O_neighbors.size(); k ++) {
            count ++;
        }
    }
    return count;
}

bool output_graphfile(Information& info, TimeStep& time_step) {
    O_vector& Ovec = time_step.O_atoms;
    int n_bins = info.num_cell_blocks;
    double bin_size = (info.cell_block_end - info.cell_block_start) / n_bins;
    ofstream output;
    output.open(info.gephi_output.c_str());
    
    output << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" xmlns:viz=\"http://www.gexf.net/1.2draft/viz\">\n"
                << "<meta lastmodifieddate=\"2013-11-21\">\n"
                << "<creator> Kevin Ryczko </creator>\n"
                << "<description> Social Network Visualization </description>\n"
                << "</meta>\n"
                << "<graph mode=\"static\" defaultedgetype=\"directed\">\n"
                << "<nodes>\n";
    for (int i = 0; i < info.num_oxygen; i ++) {
        Oxygen& O = Ovec[i];
        vector<int>& bonds = O.bonded_O_neighbors;
        if (O.z_coords < info.starting_z) {
            output << "<node id=\"" << i << "\" label=\"" << (int) ((O.z_coords + info.lattice_z) / bin_size) << "\"/>\n";
        } else {
            output << "<node id=\"" << i << "\" label=\"" << (int) (O.z_coords / bin_size) << "\"/>\n";
        }
    }
    int count = 0;
    output << "</nodes>\n" << "<edges>\n";
    for (int i = 0; i < info.num_oxygen; i ++) {
        Oxygen& O = Ovec[i];
        vector<int>& bonds = O.bonded_O_neighbors;
        for (int j = 0; j < bonds.size(); j ++) {
            output << "<edge id=\"" << count << "\" source=\"" << i
                            << "\" target=\"" << bonds[j] << "\"/>\n";
            count ++;
        }
    }
    output << "</edges>\n" << "</graph>\n" << "</gexf>";
    output.close();
    return true;
}

bool degree_respect_z(Information& info, TimeSteps& time_steps) {
    ofstream output;
    output.open(info.degree_z_output.c_str());
    int n_bins = info.degree_bins;
    double binsize = info.lattice_z / n_bins;
    vector<int> bincounts ( n_bins ), counts ( n_bins );
    vector<double> degrees ( n_bins );
    for (int i = 0; i < n_bins; i ++) {
        bincounts[i] = 0;
        counts[i] = 0;
    }
    int total_counts = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            int degree = O.bonded_O_neighbors.size() + O.out_degree;
            int bin = (int) (O.z_coords / binsize);
            bincounts[bin] += degree;
            counts[bin] ++;
        }
    }
    double sum = 0;
    for (int i = 0; i < n_bins; i ++) {
        degrees[i] = (double) bincounts[i] / counts[i];
        sum += (double) bincounts[i] / counts[i];
    }
    if (info.fix_plots) {
        for (int i = 0; i < n_bins; i ++) {
            degrees[i] = (double) bincounts[i] / counts[i];
            if (counts[i] && i*binsize > info.starting_z) {
                output << i*binsize << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
                output << i*binsize + binsize << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
            } else if (!counts[i] && i*binsize > info.starting_z) {
                    output << i*binsize << "\t" << 0 << "\t" << sum / n_bins << "\n";
                    output << i*binsize + binsize << "\t" << 0 << "\t" << sum / n_bins << "\n";
            }
        }
        for (int i = 0; i < n_bins; i ++) {
            degrees[i] = (double) bincounts[i] / counts[i];
            if (counts[i] && i*binsize < info.starting_z) {
                output << i*binsize + info.lattice_z << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
                output << i*binsize + binsize + info.lattice_z << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
            } else if (!counts[i] && i*binsize < info.starting_z){
                    output << i*binsize + info.lattice_z << "\t" << 0 << "\t" << sum / n_bins << "\n";
                    output << i*binsize + info.lattice_z + binsize << "\t" << 0 << "\t" << sum / n_bins << "\n";
            }
        }
    }
    else {
        for (int i = 0; i < n_bins; i ++) {
            degrees[i] = (double) bincounts[i] / counts[i];
            if (counts[i]) {
                output << i*binsize << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
                output << i*binsize + binsize << "\t" << degrees[i] << "\t" << sum / n_bins << "\n";
            } else {
                    output << i*binsize << "\t" << 0 << "\t" << sum / n_bins << "\n";
                    output << i*binsize + binsize << "\t" << 0 << "\t" << sum / n_bins << "\n"; 
            }
        }
    }
    output.close();
    return true;
}

void out_count(Information& info, TimeSteps& time_steps) {
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.bonded_O_neighbors.size(); k++) {
                Oxygen& O2 = Ovec[O.bonded_O_neighbors[k]];
                O2.out_degree ++;
            }
        }
    }
}

bool OOdistro(Information& info, TimeSteps& time_steps) {
    vector<int> bins ( info.OO_bins );
    int counter = 0;
    double binsize = info.max_OO / info.OO_bins;
    for (int i = 0; i < info.OO_bins; i ++) {
        bins[i] = 0;
    }
    double sum = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.bonded_O_neighbors.size(); k++) {
                Oxygen& O2 = Ovec[O.bonded_O_neighbors[k]];
                double dist = OOdist(O, O2, info);
                sum += dist;
                int bin = dist / binsize;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.OO_output.c_str());
    output << "# distance\tprobability\taverage distance\n\n";
    for (int i = 0; i < info.OO_bins; i ++) {
        output << i*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i*binsize + binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    return true;
}

bool OHdistro(Information& info, TimeSteps& time_steps) {
    vector<int> bins ( info.OH_bins );
    int counter = 0;
    double binsize = info.max_OH / info.OH_bins;
    for (int i = 0; i < info.OH_bins; i ++) {
        bins[i] = 0;
    }
    double sum = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            for (int k = 0; k < O.local_H_neighbors.size(); k++) {
                Hydrogen& H = Hvec[O.local_H_neighbors[k]];
                double dist = OHdist(O, H, info);
                sum += dist;
                int bin = dist / binsize;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.OH_output.c_str());
    output << "# distance\tprobability\taverage distance\n\n";
    for (int i = 0; i < info.OH_bins; i ++) {
        output << i*binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i*binsize + binsize << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    return true;
}

// radians to degrees
const double rtd = 57.2957795;

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

bool HOHdistro(Information& info, TimeSteps& time_steps) {
    // 360 degrees, bin for each degree
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
    return true;
}

bool degree_distro(Information& info, TimeSteps& time_steps) {
    double incr = (info.cell_block_end - info.cell_block_start) / info.num_cell_blocks;
    for (int cell = 0; cell < info.num_cell_blocks; cell ++) {
        double start_z = info.cell_block_start + cell*incr;
        double end_z = info.cell_block_start + (cell+1)*incr;
        int max_degree = 0;
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                int degree = O.bonded_O_neighbors.size() + O.out_degree;
                if (degree > max_degree) {
                    max_degree = degree;
                }
            }
        }
        vector<int> degree_counts ( max_degree );
        for (int i = 0; i < max_degree; i ++) {
            degree_counts[i] = 0;
        }
        int counter = 0;
        double sum = 0;
        for (int i = 0; i < time_steps.size(); i ++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            for (int j = 0; j < Ovec.size(); j ++) {
                Oxygen& O = Ovec[j];
                if (O.z_coords > start_z && O.z_coords < end_z) {
                    int degree = O.bonded_O_neighbors.size() + O.out_degree;
                    degree_counts[degree] ++;
                    sum += degree;
                    counter ++;
                }
            }
        }
        ofstream output;
        string filename = info.degree_output + "region_" + to_string(cell) + ".dat";
        output.open(filename.c_str());
        output << "# degree\tprobability\taverage degree\n# for region " << start_z << " to " << end_z << "[Angstroms]\n\n";
        for (int i = 0; i < max_degree; i ++) {
            output << i << "\t" << degree_counts[i] / (double) counter << "\t" << sum / counter << "\n";
        }
        output.close();
    }
    return true;
}

double wrap(double value, double lattice) {
    if (value > lattice) {
        return value - lattice;
    } else if ( value < 0.0) {
        return value + lattice;
    } 
    return value;
}

bool density(Information& info, TimeSteps& time_steps) {
    vector<int> xbins (info.density_bins), ybins (info.density_bins), zbins (info.density_bins),
                O_xbins(info.density_bins), O_ybins(info.density_bins), O_zbins(info.density_bins),
                H_xbins(info.density_bins), H_ybins(info.density_bins), H_zbins(info.density_bins);
    for (int i = 0; i < info.density_bins; i ++) {
        xbins[i] = 0;
        ybins[i] = 0;
        zbins[i] = 0;
        O_xbins[i] = 0;
        O_ybins[i] = 0;
        O_zbins[i] = 0;
        H_xbins[i] = 0;
        H_ybins[i] = 0;
        H_zbins[i] = 0;
        
    }
    double xbinsize = info.lattice_x / info.density_bins;
    double ybinsize = info.lattice_y / info.density_bins;
    double zbinsize = info.lattice_z / info.density_bins;
    
    double conversion, O_conversion = 16.0e-6/(6.023e23*1.0e-30), H_conversion;
    if (info.heavy_water) {
        H_conversion = 2.0e-6/(6.023e23*1.0e-30);
        conversion = 20.0e-6/(6.023e23*1.0e-30);
    } else {
        H_conversion = 1.0e-6/(6.023e23*1.0e-30);
        conversion = 18.0e-6/(6.023e23*1.0e-30);
    }
    
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            double x = wrap(O.x_coords, info.lattice_x);
            double y = wrap(O.y_coords, info.lattice_y);
            double z = wrap(O.z_coords, info.lattice_z);
            int xbin = x / xbinsize;
            int ybin = y / ybinsize;
            int zbin = z / zbinsize;
            xbins[xbin] ++;
            ybins[ybin] ++;
            zbins[zbin] ++;
            O_xbins[xbin] ++;
            O_ybins[ybin] ++;
            O_zbins[zbin] ++;
        }
        for (int j = 0; j < Hvec.size(); j ++) {
            Hydrogen& H = Hvec[j];
            double x = wrap(H.x_coords, info.lattice_x);
            double y = wrap(H.y_coords, info.lattice_y);
            double z = wrap(H.z_coords, info.lattice_z);
            int xbin = x / xbinsize;
            int ybin = y / ybinsize;
            int zbin = z / zbinsize;
            H_xbins[xbin] ++;
            H_ybins[ybin] ++;
            H_zbins[zbin] ++;
        }
    }
    double xvol = xbinsize*info.lattice_y*info.lattice_z, yvol = ybinsize*info.lattice_x*info.lattice_z, zvol = zbinsize*info.lattice_x*info.lattice_y;
    double xsum = 0, ysum = 0, zsum = 0, O_xsum = 0, O_ysum = 0, O_zsum = 0, H_xsum = 0, H_ysum = 0, H_zsum = 0;
    int xcount = 0, ycount = 0, zcount = 0, O_xcount = 0, O_ycount = 0, O_zcount = 0, H_xcount = 0, H_ycount = 0, H_zcount = 0;
    for (int i = 0; i < info.density_bins; i ++ ) {
        xsum += xbins[i]*conversion / (xvol*info.n_frames);
        ysum += ybins[i]*conversion / (yvol*info.n_frames);
        zsum += zbins[i]*conversion / (zvol*info.n_frames);
        O_xsum += O_xbins[i]*conversion / (xvol*info.n_frames);
        O_ysum += O_ybins[i]*conversion / (yvol*info.n_frames);
        O_zsum += O_zbins[i]*conversion / (zvol*info.n_frames);
        H_xsum += H_xbins[i]*0.5*conversion / (xvol*info.n_frames);
        H_ysum += H_ybins[i]*0.5*conversion / (yvol*info.n_frames);
        H_zsum += H_zbins[i]*0.5*conversion / (zvol*info.n_frames);
        if (xbins[i]) {
            xcount ++;
        }
        if (ybins[i]) {
            ycount ++;
        }
        if (zbins[i]) {
            zcount ++;
        }
        if (O_xbins[i]) {
            O_xcount ++;
        }
        if (O_ybins[i]) {
            O_ycount ++;
        }
        if (O_zbins[i]) {
            O_zcount ++;
        }
        if (H_xbins[i]) {
            H_xcount ++;
        }
        if (H_ybins[i]) {
            H_ycount ++;
        }
        if (H_zbins[i]) {
            H_zcount ++;
        }
    }
    ofstream xoutput, youtput, zoutput;
    xoutput.open(info.xdens_output.c_str());
    youtput.open(info.ydens_output.c_str());
    zoutput.open(info.zdens_output.c_str());
    if (info.fix_plots) {
        for (int i = 0; i < info.density_bins; i ++) {
            xoutput << i*xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*conversion*0.5 / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            xoutput << i*xbinsize + xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*conversion*0.5 / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            youtput << i*ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*conversion*0.5 / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
            youtput << i*ybinsize + ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*conversion*0.5 / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
            if (i*zbinsize > info.starting_z) {
                zoutput << i*zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
                zoutput << i*zbinsize + zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";   
            }
        }
        for (int i = 0; i < info.density_bins; i ++) {
            if (i*zbinsize < info.starting_z) {
                zoutput << i*zbinsize + info.lattice_z << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
                zoutput << i*zbinsize + zbinsize + info.lattice_z << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";   
            }
        }
    }
    else {
        for (int i = 0; i < info.density_bins; i ++) {
            xoutput << i*xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*conversion*0.5 / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            xoutput << i*xbinsize + xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << O_xbins[i]*conversion / (xvol*info.n_frames) << "\t" << H_xbins[i]*conversion*0.5 / (xvol*info.n_frames) << "\t" << xsum / xcount << "\t" << O_xsum / O_xcount << "\t" << H_xsum / H_xcount << "\n";
            youtput << i*ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*conversion*0.5 / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
            youtput << i*ybinsize + ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << O_ybins[i]*conversion / (yvol*info.n_frames) << "\t" << H_ybins[i]*conversion*0.5 / (yvol*info.n_frames) << "\t" << ysum / ycount << "\t" << O_ysum / O_ycount << "\t" << H_ysum / H_ycount << "\n";
            zoutput << i*zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
            zoutput << i*zbinsize + zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << O_zbins[i]*conversion / (zvol*info.n_frames) << "\t" << H_zbins[i]*conversion*0.5 / (zvol*info.n_frames) << "\t" << zsum / zcount << "\t" << O_zsum / O_zcount << "\t" << H_zsum / H_zcount << "\n";
        }
    }
    xoutput.close();
    youtput.close();
    zoutput.close();
    return true;
}

void unwrap(Information& info, TimeSteps& time_steps) {
    for (int i = 1; i < time_steps.size(); i ++) {
        O_vector& Ovec0 = time_steps[i - 1].O_atoms;
        O_vector& Ovec1 = time_steps[i].O_atoms;
        H_vector& Hvec0 = time_steps[i - 1].H_atoms;
        H_vector& Hvec1 = time_steps[i].H_atoms;
        if (i == 1) {
            for (int j = 0; j < Ovec0.size(); j ++) {
                Ovec0[j].unwrap_x = Ovec0[j].x_coords;
                Ovec0[j].unwrap_y = Ovec0[j].y_coords;
                Ovec0[j].unwrap_z = Ovec0[j].z_coords;
            }
            for (int j = 0; j < Hvec0.size(); j ++) {
                Hvec0[j].unwrap_x = Hvec0[j].x_coords;
                Hvec0[j].unwrap_y = Hvec0[j].y_coords;
                Hvec0[j].unwrap_z = Hvec0[j].z_coords;
            }
        }
        for (int j = 0; j < Ovec0.size() ; j ++) {
            if (Ovec1[j].x_coords - Ovec0[j].unwrap_x > (info.lattice_x / 2.0)) {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords - info.lattice_x;
            } else if (Ovec1[j].x_coords - Ovec0[j].unwrap_x < (-info.lattice_x / 2.0)) {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords + info.lattice_x;
            } else {
                Ovec1[j].unwrap_x = Ovec1[j].x_coords;
            }
            if (Ovec1[j].y_coords - Ovec0[j].unwrap_y > (info.lattice_y / 2.0)) {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords - info.lattice_y;
            } else if (Ovec1[j].y_coords - Ovec0[j].unwrap_y < (-info.lattice_y / 2.0)) {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords + info.lattice_y;
            } else {
                Ovec1[j].unwrap_y = Ovec1[j].y_coords;
            }
            if (Ovec1[j].z_coords - Ovec0[j].unwrap_z > (info.lattice_z / 2.0)) {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords - info.lattice_z;
            } else if (Ovec1[j].z_coords - Ovec0[j].unwrap_z < (-info.lattice_z / 2.0)) {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords + info.lattice_z;
            } else {
                Ovec1[j].unwrap_z = Ovec1[j].z_coords;
            }
        }
        for (int j = 0; j < Hvec0.size() ; j ++) {
            if (Hvec1[j].x_coords - Hvec0[j].unwrap_x > (info.lattice_x / 2.0)) {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords - info.lattice_x;
            } else if (Hvec1[j].x_coords - Hvec0[j].unwrap_x < (-info.lattice_x / 2.0)) {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords + info.lattice_x;
            } else {
                Hvec1[j].unwrap_x = Hvec1[j].x_coords;
            }
            if (Hvec1[j].y_coords - Hvec0[j].unwrap_y > info.lattice_y / 2.0) {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords - info.lattice_y;
            } else if (Hvec1[j].y_coords - Hvec0[j].unwrap_y < -info.lattice_y / 2.0) {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords + info.lattice_y;
            } else {
                Hvec1[j].unwrap_y = Hvec1[j].y_coords;
            }
            if (Hvec1[j].z_coords - Hvec0[j].unwrap_z > info.lattice_z / 2.0) {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords - info.lattice_z;
            } else if (Hvec1[j].z_coords - Hvec0[j].unwrap_z < -info.lattice_z / 2.0) {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords + info.lattice_z;
            } else {
                Hvec1[j].unwrap_z = Hvec1[j].z_coords;
            }
        }
    }
    if (info.write_unwrapped_xyz) {
        ofstream output;
        cout << "Writing unwrapped coordinates to: " << info.unwrapped_coords << "\n\n";
        output.open(info.unwrapped_coords.c_str());
        for (int i = 0; i < time_steps.size(); i++) {
            O_vector& Ovec = time_steps[i].O_atoms;
            H_vector& Hvec = time_steps[i].H_atoms;
            output << info.num_oxygen + info.num_hydrogen << "\n\n";
            for (int j = 0; j < Ovec.size(); j ++) {
                output << "O" << "\t" << Ovec[j].unwrap_x << "\t" << Ovec[j].unwrap_y << "\t" << Ovec[j].unwrap_z << "\n";
            }
            for (int j = 0; j < Hvec.size(); j ++) {
                output << "H" << "\t" << Hvec[j].unwrap_x << "\t" << Hvec[j].unwrap_y << "\t" << Hvec[j].unwrap_z << "\n";
            }
        }
        cout << "Finished writing unwrapped coordinates.\n\n";  
        output.close();
    }
}

struct MsdStep {
    vector<double> msd_step;
    void declare(int len) {
            msd_step.clear();
            for (int i = 0; i < len; i ++) {
                msd_step.push_back(0);
            }
    }
};
    
typedef std::vector<MsdStep> msd_vec;

msd_vec msd_vector;

bool msd(Information& info, TimeSteps& time_steps) {
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
                double comx = (8*dxo + 2*dxh) / (8*Ovec.size() + 2*Hvec.size());
                COMx.push_back(comx);
                double comy = (8*dyo + 2*dyh) / (8*Ovec.size() + 2*Hvec.size());
                COMy.push_back(comy);
                double comz = (8*dzo + 2*dzh) / (8*Ovec.size() + 2*Hvec.size());
                COMz.push_back(comz);
            } else {
                double comx = (8*dxo + dxh) / (8*Ovec.size() + Hvec.size());
                COMx.push_back(comx);
                double comy = (8*dyo + dyh) / (8*Ovec.size() + Hvec.size());
                COMy.push_back(comy);
                double comz = (8*dzo + dzh) / (8*Ovec.size() + Hvec.size());
                COMz.push_back(comz);
            }
        }
        cout << "--- MSD progress ---\n";
        int print_count = 0;
        int starting_step, final_step;
        for (int nb = 0; nb < info.num_blocks; nb ++) {
            starting_step = nb*length;
            if (info.full_msd) {
                final_step = time_steps.size();
            } else {
                final_step = (nb + 1)*length;
            }
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
        if (info.full_msd) {
            for (int i = 0; i < time_steps.size(); i ++) {
                averaged_msd.push_back(0);
                counts.push_back(0);
            }
        } else {
            for (int i = 0; i < length; i ++) {
                averaged_msd.push_back(0);
                counts.push_back(0);
            }
        }
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
    return true;
}

bool orientation(Information& info, TimeSteps& time_steps) {
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
    return true;
}
bool sdf(Information& info, TimeSteps& time_steps) {
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
    return true;
}

bool orientation_1D(Information& info, TimeSteps& time_steps) {
    return true;
}

bool output_edgelist(Information& info, TimeSteps& time_steps) {
    ofstream output;
    output.open(info.edgelist_output_filename.c_str());
    
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        output << "Frame:\t" << i << "\t" << num_edges(Ovec) << "\n";
        for (int j = 0; j < Ovec.size(); j ++) {
            for (int k = 0; k < Ovec[j].bonded_O_neighbors.size(); k ++) {
                output << Ovec[j].ID << "\t" << Ovec[j].bonded_O_neighbors[k] << "\n";
            }
        }
    }
    cout << "Outputted edge list to data file given.\n\n";
    output.close();
    
    out_count(info, time_steps);
    
    if (info.output_gephi) {
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
    return true;
}