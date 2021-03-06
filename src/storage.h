#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <vector>


struct Information {
    std::string input_filename, metal;
    std::string edgelist_output_filename, in_degree, out_degree, cumulative_degree, gephi_output, degree_z_output, msd_filename, orientation_filename;
    std::string OO_output, OH_output, HOH_output, degree_output, xdens_output, ydens_output, zdens_output;
    bool create_edgelist, output_gephi, degree_z, OODistro, OHDistro, HOHDistro, degree_distro, density, heavy_water, msd, orientation;
    bool H_group_dynamics, OH_group_dynamics, write_unwrapped_xyz, fix_plots;
    int num_oxygen, num_hydrogen, timesteps, num_blocks, label_bins, degree_bins, OO_bins, OH_bins, HOH_bins, density_bins, orient_x_bins, orient_y_bins, orient_z_bins;
    double lattice_x, lattice_y, lattice_z, max_OO, max_OH, time_step, cell_block_start, cell_block_end;
    int n_frames, num_threads;
    std::string H_group_trajectory_filename, OH_group_trajectory_filename;
    std::string unwrapped_coords;
    bool sdf, network_reorganization_time, orientation_1D, full_msd, degree_z_from_metal, density_from_metal;
    double sdf_start, sdf_end, starting_z, time_length, msd_time;
    int sdf_bins, num_cell_blocks, num_starts, num_metals;
    std::string sdf_output, nrt_output;
    
    Information() {
        // declare incorrect values, so we need to read things in to obtain information
        edgelist_output_filename = input_filename = in_degree = out_degree = cumulative_degree = gephi_output = degree_z_output = unwrapped_coords = "not_a_file";
        OO_output = OH_output = HOH_output = degree_output = xdens_output = ydens_output = zdens_output = msd_filename = orientation_filename = sdf_output = nrt_output = "notafile";
        n_frames = num_oxygen = num_hydrogen = sdf_bins = -1;
        H_group_dynamics = OH_group_dynamics = write_unwrapped_xyz = sdf = network_reorganization_time = orientation_1D = create_edgelist = degree_z_from_metal = density_from_metal = false;
        lattice_x = lattice_y = lattice_z = max_OO = max_OH = time_step = sdf_start = sdf_end = cell_block_start = cell_block_end = -1.0; 
        output_gephi = degree_z = OODistro = OHDistro = HOHDistro = degree_distro = density = heavy_water = msd = orientation = fix_plots = full_msd = network_reorganization_time = false;
        label_bins = degree_bins = OO_bins = num_blocks = OH_bins = HOH_bins = density_bins = orient_x_bins = orient_y_bins = orient_z_bins = 0;
        starting_z = time_length = msd_time -1.0;
        num_cell_blocks = 1;
        num_starts = num_metals = -1;
    }
};

struct Metal {
    std::string name; 
    int ID;
    double x_coords;
    double y_coords;
    double z_coords;
    double unwrap_x;
    double unwrap_y;
    double unwrap_z;
    
    Metal() {
        ID = -1;
        x_coords = y_coords = z_coords = -1.0;
    }
};
typedef std::vector<Metal> M_vector;

struct Hydrogen {
    int ID;
    double x_coords;
    double y_coords;
    double z_coords;
    double unwrap_x;
    double unwrap_y;
    double unwrap_z;
    
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
    double unwrap_x;
    double unwrap_y;
    double unwrap_z;
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
    M_vector M_atoms;
};

typedef std::vector<TimeStep> TimeSteps;

struct Args {
    Information arg_info;
    TimeStep arg_time_step;
    TimeSteps arg_time_steps;
    double avg_left, avg_right;
    Args() {
        avg_left = avg_right = 0;
    }
};

#endif

