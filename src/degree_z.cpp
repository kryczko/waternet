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

void  degree_respect_metal(Args& args) {
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    double average_l = args.avg_left + info.lattice_z;
    double average_r = args.avg_right;
    double metal_d = average_r - (average_l - info.lattice_z);
    double water_z_d = abs(info.lattice_z - metal_d);
    ofstream output;
    string filename = "output/degree_wrt_metal.dat";
    output.open(filename.c_str());
    
    int n_bins = info.degree_bins / 2;
    double zinc = 0.5 * water_z_d / (n_bins);
    vector<int> bincounts ( n_bins ), counts ( n_bins );
    vector<double> degrees ( n_bins );
    for (int i = 0; i < n_bins; i ++) {
        bincounts[i] = 0;
        counts[i] = 0;
        degrees[i] = 0;
    }
    int total_counts = 0;
    for (int i = 0; i < time_steps.size(); i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            int degree = O.bonded_O_neighbors.size() + O.out_degree;
            double dist1 = abs(O.z_coords - average_l);
            double dist2 = abs(O.z_coords - average_r);
            double z = min(dist1, dist2);
            int bin = z / zinc;
            degrees[bin] += degree;
            counts[bin] ++;
        }
    }
    for (int i = 0; i < n_bins; i ++) {
        if (counts[i]) {
            output << i*zinc << "\t" << 2.0 * degrees[i] / counts[i] << "\n";
        } else {
            output << i*zinc << "\t" << 0 << "\n";
        }
    }
    cout << "Outputted degree with respect to metal data file.\n\n";
    
    output.close();
}

void  degree_respect_z(Args& args) {
    Information& info = args.arg_info;
    TimeSteps& time_steps = args.arg_time_steps;
    
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
    cout << "Outputted degree with respect to z data file.\n\n";
    
}