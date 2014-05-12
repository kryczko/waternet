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
    const Node& files = node["filenames"];
    parse(files, "input", info.input_filename);
    parse(files, "edgelist_output", info.edgelist_output_filename);
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