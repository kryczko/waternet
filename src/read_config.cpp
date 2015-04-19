#include <fstream>
#include <string>
#include <yaml-cpp/yaml.h>
#include <iostream>

#include "storage.h"

using namespace std;
using namespace YAML;

template<class T>
static inline void parse(const Node& node, const char* key, T& value) {
    if (node.FindValue(key)) {
        node[key] >> value;
    } else {
        printf("'%s' was not found!", key);
    }
}


void parse_inputfile(Information& info, const Node& node) {
    const Node& modules = node["modules"];
    parse(modules, "metal", info.metal);
    parse(modules, "num_threads", info.num_threads);
    parse(modules, "create_edgelist", info.create_edgelist);
    parse(modules, "num_cell_blocks", info.num_cell_blocks);
    parse(modules, "cell_block_start", info.cell_block_start);
    parse(modules, "cell_block_end", info.cell_block_end);
    parse(modules, "output_gephi", info.output_gephi);
    parse(modules, "degree_z", info.degree_z);
    parse(modules, "degree_z_from_metal", info.degree_z_from_metal);
    parse(modules, "degree_bins", info.degree_bins);
    parse(modules, "fix_plots", info.fix_plots);
    parse(modules, "starting_z", info.starting_z);
    parse(modules, "OOdistro", info.OODistro);
    parse(modules, "OO_bins", info.OO_bins);
    parse(modules, "OO_max_dist", info.max_OO);
    parse(modules, "OHdistro", info.OHDistro);
    parse(modules, "OH_bins", info.OH_bins);
    parse(modules, "OH_max_dist", info.max_OH);
    parse(modules, "HOHdistro", info.HOHDistro);
    parse(modules, "HOH_bins", info.HOH_bins);
    parse(modules, "degree_distro", info.degree_distro);
    parse(modules, "num_cell_blocks", info.num_cell_blocks);
    parse(modules, "density", info.density);
    parse(modules, "density_from_metal", info.density_from_metal);
    parse(modules, "density_bins", info.density_bins);
    parse(modules, "heavy_water", info.heavy_water);
    parse(modules, "mean_square_displacement", info.msd);
    parse(modules, "time_step", info.time_step);
    parse(modules, "num_blocks", info.num_blocks);
    parse(modules, "full_msd", info.full_msd);
    parse(modules, "write_unwrapped_xyz", info.write_unwrapped_xyz);
    parse(modules, "orientation_2D", info.orientation);
    parse(modules, "orient_x_bins", info.orient_x_bins);
    parse(modules, "orient_y_bins", info.orient_y_bins);
    parse(modules, "orient_z_bins", info.orient_z_bins);
    parse(modules, "H_group_dynamics", info.H_group_dynamics);
    parse(modules, "OH_group_dynamics", info.OH_group_dynamics);
    parse(modules, "spacial_distribution_function", info.sdf);
    parse(modules, "sdf_z_start", info.sdf_start);
    parse(modules, "sdf_z_end", info.sdf_end);
    parse(modules, "sdf_bins", info.sdf_bins);
    parse(modules, "network_reorganization_time", info.network_reorganization_time);
    parse(modules, "orientation_1D", info.orientation_1D);
    parse(modules, "time_length", info.time_length);
    parse(modules, "num_starts", info.num_starts);
    
    const Node& files = node["filenames"];
    parse(files, "input", info.input_filename);
    parse(files, "edgelist_output", info.edgelist_output_filename);
    parse(files, "gephi_output", info.gephi_output);
    parse(files, "degree_z_output", info.degree_z_output);
    parse(files, "OOdistro_output", info.OO_output);
    parse(files, "OHdistro_output", info.OH_output);
    parse(files, "HOHdistro_output", info.HOH_output);
    parse(files, "degree_distro_output", info.degree_output);
    parse(files, "x_density_output", info.xdens_output);
    parse(files, "y_density_output", info.ydens_output);
    parse(files, "z_density_output", info.zdens_output);
    parse(files, "mean_square_displacement_output", info.msd_filename);
    parse(files, "orientation_2D_output", info.orientation_filename);
    parse(files, "H_group_trajectory", info.H_group_trajectory_filename);
    parse(files, "OH_group_trajectory", info.OH_group_trajectory_filename);
    parse(files, "unwrapped_coords", info.unwrapped_coords);
    parse(files, "spacial_distribution_output", info.sdf_output);
    parse(files, "network_reorganization_time_distro_output", info.nrt_output);
    
    const Node& lattice_consts = node["lattice_constants"];
    parse(lattice_consts, "x", info.lattice_x);
    parse(lattice_consts, "y", info.lattice_y);
    parse(lattice_consts, "z", info.lattice_z);
    const Node& atoms = node["n_atoms"];
    parse(atoms, "O", info.num_oxygen);
    parse(atoms, "H", info.num_hydrogen);
    parse(atoms, "metal", info.num_metals);
}

bool read_configfile(Information& info) {
    fstream file("input.yaml", fstream::in);
    if (!file.is_open()) {
        return false;
    }
    Parser parser(file);
    Node root;
    parser.GetNextDocument(root);
    parse_inputfile(info, root);
    cout << "Configuration file \"input.yaml\" read.\n\n";
    return true;
}

