#include <fstream>
#include <string>
#include <yaml-cpp/yaml.h>

#include "storage.h"

using namespace std;
using namespace YAML;

template<class T>
static inline void parse(const Node& node, const char* key, T& value, bool optional = false) {
    if (node.FindValue(key)) {
        node[key] >> value;
    } else {
        printf("'%s' was not found!", key);
    }
}


void parse_inputfile(Information& info, const Node& node) {
    const Node& modules = node["modules"];
    parse(modules, "output_gephi", info.output_gephi);
    parse(modules, "label_bins", info.label_bins);
    parse(modules, "degree_z", info.degree_z);
    parse(modules, "degree_bins", info.degree_bins);
    parse(modules, "OOdistro", info.OODistro);
    parse(modules, "OO_bins", info.OO_bins);
    parse(modules, "OO_max_dist", info.max_OO);
    parse(modules, "OHdistro", info.OHDistro);
    parse(modules, "OH_bins", info.OH_bins);
    parse(modules, "OH_max_dist", info.max_OH);
    parse(modules, "HOHdistro", info.HOHDistro);
    parse(modules, "HOH_bins", info.HOH_bins);
    parse(modules, "degree_distro", info.degree_distro);
    parse(modules, "density", info.density);
    parse(modules, "density_bins", info.density_bins);
    parse(modules, "heavy_water", info.heavy_water);
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
    const Node& lattice_consts = node["lattice_constants"];
    parse(lattice_consts, "x", info.lattice_x);
    parse(lattice_consts, "y", info.lattice_y);
    parse(lattice_consts, "z", info.lattice_z);
    const Node& atoms = node["n_atoms"];
    parse(atoms, "O", info.num_oxygen);
    parse(atoms, "H", info.num_hydrogen);
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