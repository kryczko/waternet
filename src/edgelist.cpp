#include <fstream>
#include <iostream>
#include <cmath>
#include <thread>
#include <pthread.h>
#include "storage.h"
#include "helper.h"
//#include "omp.h"

using namespace std;


void nearest_neighbors(Information& info, Oxygen& O, O_vector& Ovec) {
    for (int i = 0; i < Ovec.size(); i ++) {
        Oxygen& O2 = Ovec[i];
        if (O.ID != O2.ID) {
            double dist = OOdist(O, O2, info);
            if (dist < 3.6) {
                O.nearest_neighbors.push_back(O2.ID);
            }
        }
    }
}

void find_local_H(Information& info, Oxygen& O, H_vector& Hvec) {
    for (int i = 0; i < Hvec.size(); i ++) {
        Hydrogen& H = Hvec[i]; 
        double dist = OHdist(O, H, info);
        if (dist < 1.2) {
            O.local_H_neighbors.push_back(H.ID);
        }
    }
}

void find_H_bonds(Oxygen& O, H_vector& Hvec, O_vector& Ovec, Information& info) {
    for (int i = 0; i < O.nearest_neighbors.size(); i ++) {
        Oxygen& O2 = Ovec[O.nearest_neighbors[i]];
        for (int j = 0; j < O2.local_H_neighbors.size(); j ++) {
            Hydrogen& H = Hvec[O2.local_H_neighbors[j]];
            double dist = OHdist(O, H, info);
            if (dist < 2.4) {
                O.bonded_H_neighbors.push_back(H.ID);
            }
        }
    }
}

void check_final_bonds(Oxygen& O, O_vector& Ovec, H_vector& Hvec, Information& info) {
    for (int i = 0; i < O.nearest_neighbors.size(); i ++) {
        Oxygen& O2 = Ovec[O.nearest_neighbors[i]];
        for (int j = 0; j < O.bonded_H_neighbors.size(); j ++) {
            Hydrogen& H = Hvec[O.bonded_H_neighbors[j]];
            double angle = angle_between(O, O2, H, info);
            if (angle < 30.0) {
                O.bonded_O_neighbors.push_back(O2.ID);
            }
        }
    }
    
}

void perform(Information& info, TimeSteps& time_steps, int start, int end) {
    for (int i = start; i < end; i ++) {
        O_vector& Ovec = time_steps[i].O_atoms;
        H_vector& Hvec = time_steps[i].H_atoms;
        
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            nearest_neighbors(info, O, Ovec);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            find_local_H(info, O, Hvec);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            find_H_bonds(O, Hvec, Ovec, info);
        }
        for (int j = 0; j < Ovec.size(); j ++) {
            Oxygen& O = Ovec[j];
            check_final_bonds(O, Ovec, Hvec, info);
        }
    }
        
}
    


bool create_edgelist(Information& info, TimeSteps& time_steps) {
    if (!info.create_edgelist) {
        return false;
    }

    int chunk = time_steps.size() / info.num_threads;
    //omp_set_num_threads(info.num_threads);
    int begin, end;
    //#pragma omp parallel for 
    for (int i = 0; i < info.num_threads; i ++) {
        begin = i*chunk;
        if (i == info.num_threads - 1) {
            end = time_steps.size();
        } else {
            end = (i+1)*chunk;
        }
        perform(info, time_steps, begin, end);
    }
    cout << "Computed nearest neighbour lists...\n\n";    
    return true;
}

