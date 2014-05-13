#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <vector>


struct Information {
    std::string input_filename;
    std::string edgelist_output_filename, in_degree, out_degree, cumulative_degree, gephi_output, degree_z_output;
    std::string OO_output, OH_output, HOH_output, degree_output, xdens_output, ydens_output, zdens_output;
    bool output_gephi, degree_z, OODistro, OHDistro, HOHDistro, degree_distro, density, heavy_water;
    int num_oxygen, num_hydrogen, timesteps, label_bins, degree_bins, OO_bins, OH_bins, HOH_bins, density_bins;
    double lattice_x, lattice_y, lattice_z, max_OO, max_OH;
    int n_frames;
    
    Information() {
        // declare incorrect values, so we need to read things in to obtain information
        edgelist_output_filename = input_filename = in_degree = out_degree = cumulative_degree = gephi_output = degree_z_output = "not_a_file";
        OO_output = OH_output = HOH_output = degree_output = xdens_output = ydens_output = zdens_output = "notafile";
        n_frames = num_oxygen = num_hydrogen = -1;
        lattice_x = lattice_y = lattice_z = max_OO = max_OH = -1.0; 
        output_gephi = degree_z = OODistro = OHDistro = HOHDistro = degree_distro = density = heavy_water = false;
        label_bins = degree_bins = OO_bins = OH_bins = HOH_bins = density_bins = 0;
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
    int out_degree;
    double x_coords;
    double y_coords;
    double z_coords;
    std::vector<int> nearest_neighbors;
    std::vector<int> bonded_O_neighbors;
    std::vector<int> bonded_H_neighbors;
    std::vector<int> local_H_neighbors;
    Oxygen() {
        ID = -1;
        out_degree = 0;
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

