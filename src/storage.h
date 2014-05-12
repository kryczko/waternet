#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <vector>


struct Information {
    std::string input_filename;
    std::string edgelist_output_filename, in_degree, out_degree, cumulative_degree;
    bool output_gephi;
    int num_oxygen, num_hydrogen, timesteps;
    double lattice_x, lattice_y, lattice_z;
    int n_frames;
    
    Information() {
        // declare incorrect values, so we need to read things in to obtain information
        edgelist_output_filename = input_filename = in_degree = out_degree = cumulative_degree = "not_a_file";
        n_frames = num_oxygen = num_hydrogen = -1;
        lattice_x = lattice_y = lattice_z = -1.0; 
        output_gephi = false;
    }
};

struct Hydrogen {
    int ID;
    double x_coords;
    double y_coords;
    double z_coords;
    
    Hydrogen() {
        ID = -1;
        x_coords = y_coords = z_coords = -1.0;
    }
};

typedef std::vector<Hydrogen> H_vector;

struct Oxygen {
    int ID;
    double x_coords;
    double y_coords;
    double z_coords;
    std::vector<int> nearest_neighbors;
    std::vector<int> bonded_O_neighbors;
    std::vector<int> bonded_H_neighbors;
    std::vector<int> local_H_neighbors;
    Oxygen() {
        ID = -1;
        x_coords = y_coords = z_coords = -1.0;
    }
};

typedef std::vector<Oxygen> O_vector;

struct TimeStep {
    O_vector O_atoms;
    H_vector H_atoms;
};

typedef std::vector<TimeStep> TimeSteps;

#endif

