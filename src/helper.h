#ifndef HELPER_H_
#define HELPER_H_

#include <vector>
#include "storage.h"
#include "analysis.h"
#include <cmath> 
#include <iostream>
using namespace std;

// radians to degrees
const double rtd = 57.2957795;

int pbc_round(double input);
double OOdist(Oxygen& O1, Oxygen& O2, Information& info);
double OHdist(Oxygen& O, Hydrogen& H, Information& info);
double angle_between(Oxygen& O1, Oxygen& O2, Hydrogen& H, Information& info);
int n_times(int);

template <class T> double average(vector<T>& vec) {
    double sum = 0;
    for (double& elem : vec) {
        sum += elem;
    }
    return sum / vec.size(); 
}

template <class T> vector<double> halfAveraged(vector<T>&v) { 
    vector<double> half_avg_vec;
    for (int i = 0; i < v.size() / 2; i ++) {
        half_avg_vec.push_back(0.5 * (v[2*i] + v[2*i+1]));
    }
    return half_avg_vec;
}

template <class T> double standardDeviation(vector<T>& v) {
    double avg = average(v);
    double std_sum = 0;
    for (auto& elem : v) {
        std_sum += pow(elem - avg, 2);
    }
    return sqrt(std_sum / v.size());
}

template<class T> vector<T> makeUncorrelated(vector<T>& vec) {
    int iter = n_times(vec.size()) - 1;
    int num_of_blocks = pow(2, iter);
    int size_of_block = vec.size() / num_of_blocks;

    vector<T> avg_vecs = vec;
    for (int i = 0; i < iter; i ++) {
       avg_vecs = halfAveraged(avg_vecs);
    }
    return avg_vecs;
}
 
template <class T> double standardError(vector<T>& vec) {
    int iter = n_times(vec.size()) - 1;
    int num_of_blocks = pow(2, iter);
    int size_of_block = vec.size() / num_of_blocks;

    vector<T> avg_vecs = vec;
    for (int i = 0; i < iter; i ++) {
       avg_vecs = halfAveraged(avg_vecs);
    }
    double err = standardDeviation(avg_vecs) / sqrt(size_of_block - 1);
    err *= sqrt((size_of_block * num_of_blocks) / (double)vec.size() );
    return err;

}

template <class T> double sum(vector<T>& v) {
    double sum = 0;
    for (auto& elem : v) {
        sum += elem;
    }
    return sum;
}

template<class T> T max(vector<T>& vec) {
    T max = 0;
    for (int i = 0; i < vec.size(); i ++) {
        if (vec[i] > max) {
            max = vec[i];
        }
    }
    return max;
}

vector<int> set_zero(vector<int>& vec);
int max(vector<int> vec);
int num_edges(O_vector& Ovec);
void out_count(Information& info, TimeSteps& time_steps);
double wrap(double value, double lattice);
void remove_zeros(vector<double>& dens, vector<double>& pos);
double find_max_m(M_vector& M_atoms, Information& info);
double find_min_m(M_vector& M_atoms, Information& info);
double min(double x, double y);
void unwrap(Information& info, TimeSteps& time_steps);
double metals_avg_right(TimeSteps& time_steps, Information& info);
double metals_avg_left(TimeSteps& time_steps, Information& info);
#endif