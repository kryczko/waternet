#include <fstream>
#include <iostream>
#include <cmath>
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
    int n_bins = info.label_bins;
    double bin_size = info.lattice_z / n_bins;
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
        output << "<node id=\"" << i << "\" label=\"" << (int) (O.z_coords / bin_size) << "\"/>\n";
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
    vector<int> bins ( 360 );
    int counter = 0;
    for (int i = 0; i < 360; i ++) {
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
                int bin = angle;
                bins[bin] ++;
                counter ++;
            }
        }
    }
    ofstream output;
    output.open(info.HOH_output.c_str());
    output << "# Angle\tProbability\tAverage angle\n\n";
    for (int i = 0; i < 360; i ++) {
        output << i << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i + 1 << "\t" << bins[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
    return true;
}

bool degree_distro(Information& info, TimeSteps& time_steps) {
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
            int degree = O.bonded_O_neighbors.size() + O.out_degree;
            degree_counts[degree] ++;
            sum += degree;
            counter ++;
        }
    }
    ofstream output;
    output.open(info.degree_output.c_str());
    output << "# degree\tprobability\taverage degree\n\n";
    for (int i = 0; i < max_degree; i ++) {
        output << i - 0.5 << "\t" << degree_counts[i] / (double) counter << "\t" << sum / counter << "\n";
        output << i + 0.5 << "\t" << degree_counts[i] / (double) counter << "\t" << sum / counter << "\n";
    }
    output.close();
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
    vector<int> xbins (info.density_bins), ybins (info.density_bins), zbins (info.density_bins);
    for (int i = 0; i < info.density_bins; i ++) {
        xbins[i] = 0;
        ybins[i] = 0;
        zbins[i] = 0;
    }
    double xbinsize = info.lattice_x / info.density_bins;
    double ybinsize = info.lattice_y / info.density_bins;
    double zbinsize = info.lattice_z / info.density_bins;
    
    double conversion;
    if (info.heavy_water) {
        conversion = 20.0e-6/(6.023e23*1.0e-30);
    } else {
        conversion = 18.0e-6/(6.023e23*1.0e-30);
    }
    
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
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
        }
    }
    double xvol = xbinsize*info.lattice_y*info.lattice_z, yvol = ybinsize*info.lattice_x*info.lattice_z, zvol = zbinsize*info.lattice_x*info.lattice_y;
    double xsum = 0, ysum = 0, zsum = 0;
    for (int i = 0; i < info.density_bins; i ++ ) {
        xsum += xbins[i]*conversion / (xvol*info.n_frames);
        ysum += ybins[i]*conversion / (yvol*info.n_frames);
        zsum += zbins[i]*conversion / (zvol*info.n_frames);
    }
    ofstream xoutput, youtput, zoutput;
    xoutput.open(info.xdens_output.c_str());
    youtput.open(info.ydens_output.c_str());
    zoutput.open(info.zdens_output.c_str());
    xoutput << "# xpos\tdensity (g/cc)\taverage\n\n";
    youtput << "# ypos\tdensity (g/cc)\taverage\n\n";
    zoutput << "# zpos\tdensity (g/cc)\taverage\n\n";
    for (int i = 0; i < info.density_bins; i ++) {
        xoutput << i*xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << xsum / info.density_bins << "\n";
        xoutput << i*xbinsize + xbinsize << "\t" << xbins[i]*conversion / (xvol*info.n_frames) << "\t" << xsum / info.density_bins << "\n";
        youtput << i*ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << ysum / info.density_bins << "\n";
        youtput << i*ybinsize + ybinsize << "\t" << ybins[i]*conversion / (yvol*info.n_frames) << "\t" << ysum / info.density_bins << "\n";
        zoutput << i*zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << zsum / info.density_bins << "\n";
        zoutput << i*zbinsize + zbinsize << "\t" << zbins[i]*conversion / (zvol*info.n_frames) << "\t" << zsum / info.density_bins << "\n";
    }
    xoutput.close();
    youtput.close();
    zoutput.close();
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
    return true;
}