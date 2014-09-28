#ifndef HELPER_H_
#define HELPER_H_

#include <vector>
#include "storage.h"
#include "analysis.h"

using namespace std;

// radians to degrees
const double rtd = 57.2957795;

int pbc_round(double input);
double OOdist(Oxygen& O1, Oxygen& O2, Information& info);
double OHdist(Oxygen& O, Hydrogen& H, Information& info);
double angle_between(Oxygen& O1, Oxygen& O2, Hydrogen& H, Information& info);
double average(vector<double>& vec);
vector<int> set_zero(vector<int>& vec);
int max(vector<int> vec);
int num_edges(O_vector& Ovec);
void out_count(Information& info, TimeSteps& time_steps);
double wrap(double value, double lattice);
void remove_zeros(vector<double>& dens, vector<double>& pos);
double find_max_m(M_vector& M_atoms, Information& info);
double find_min_m(M_vector& M_atoms);
double min(double x, double y);
void unwrap(Information& info, TimeSteps& time_steps);
double metals_avg_right(TimeSteps& time_steps, Information& info);
double metals_avg_left(TimeSteps& time_steps);
#endif