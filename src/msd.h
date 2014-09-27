#ifndef MSD_H_
#define MSD_H_

void  msd(Args& args);

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


#endif