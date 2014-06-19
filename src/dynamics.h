#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include<vector>

struct H_group_info {
    int last_oxygen;
    int counter;
    int last_timestep;
    std::vector<int> counts, Os;
    H_group_info() {
        last_oxygen = -1;
        counter = 0;
    }

};

struct OH_group_info {

    OH_group_info() {
        
    }

};

bool main_dynamics_func(Information& info, TimeSteps& time_steps, H_group_info& hgi);


#endif