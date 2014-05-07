#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <vector>


struct Information {
    std::string filename;
    int num_oxygen, num_hydrogen;
    double lattice_x, lattice_y, lattice_z;
    int n_frames;
    
    Information() {
        // declare incorrect values, so we need to read things in to obtain information
        filename = "not_a_file";
        n_frames = num_oxygen = num_hydrogen = -1;
        lattice_x = lattice_y = lattice_z = -1.0; 
    }
};

struct Oxygen {
    int ID;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
    
    Oxygen() {
        ID = -1;
    }
};

struct Hydrogen {
    int ID;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
    
    Hydrogen() {
        ID = -1;
    }
};

typedef std::vector<Oxygen> O_vector;
typedef std::vector<Hydrogen> H_vector;

struct TimeStep {
    O_vector O_atoms;
    H_vector H_atoms;
};

typedef std::vector<TimeStep> TimeSteps;

#endif

