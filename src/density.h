#ifndef DENSITY_H_
#define DENSITY_H_

#include "helper.h"

class DensityDistro {
public:

    vector<double> distro;
    void normalize() {
        double sum = 0;
        for (auto& val : this->distro)
            sum += val;
        for(auto& val : this->distro)
            val /= sum;
    }

    DensityDistro(int nbins) {
        for (int i = 0; i < nbins; i ++) {
            distro.push_back(0.0);
        }
    }
};

class AllDistros {
public:
    vector<double> z;
    vector<DensityDistro> all_distros;

    vector<double> avg_distro() {
        DensityDistro return_distro(this->all_distros[0].distro.size());
        for (auto& d : this->all_distros) {
            for (int i = 0; i < d.distro.size(); i ++) {
                return_distro.distro[i] += d.distro[i];
            }
        }
        for (auto& elem : return_distro.distro) {
            elem /= (double) this->all_distros.size();
        }
        return return_distro.distro;
    }

    vector<double> std_err() {
        vector<double> returnDistro;
        vector<double> averages = this->avg_distro();
        vector<vector<double>> distros;
        for (int i = 0; i < this->all_distros[0].distro.size(); i ++) {
        	vector<double> std_err_distro;
            for (auto& d : this->all_distros) {
                std_err_distro.push_back(d.distro[i]);
            }
            distros.push_back(std_err_distro);
        }
        for (int i = 0; i < distros.size(); i ++) {
        	returnDistro.push_back(standardError(distros[i]));
        }
        return returnDistro;
    }
};

void  density(Args& args);
void  zdens_from_metal(Args& args);

#endif