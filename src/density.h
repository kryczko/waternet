#ifndef DENSITY_H_
#define DENSITY_H_

void  density(Args& args);
void  zdens_from_metal(Information& info, TimeSteps& time_steps);
double metals_avg_left(TimeSteps& time_steps);
double metals_avg_right(TimeSteps& time_steps, Information& info);

#endif