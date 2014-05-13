#ifndef EDGELIST_H_
#define EDGELIST_H_

bool create_edgelist(Information& info, TimeSteps& time_steps);
double OHdist(Oxygen& O, Hydrogen& H, Information& info);
double OOdist(Oxygen& O, Oxygen& O2, Information& info);
int pbc_round(double input);

#endif