// global includes
#include <iostream>
#include <fstream>
#include <string>

// local includes
#include "storage.h"
#include "read_config.h"
#include "read_data.h"
#include "edgelist.h"
#include "output.h"
#include "dynamics.h"

using namespace std;

// declaration of data structures
Information info;
TimeSteps time_steps;
H_group_info hgi;

void print_welcome_message() {
    cout << "\n\t\tWelcome to --waternet--\n\n";
}

int main() {
    print_welcome_message();
    if (!read_configfile(info)) {
        cout << "Error reading configuration file, make sure it is in your working directory!\n\n";
        return 0;
    }
    if (!read_datafile(info, time_steps)) {
        cout << "Error reading/storing datafile, exiting...\n\n";
        return 0;
    }
    if (!create_edgelist(info, time_steps)) {
        cout << "Error creating edgelists, exiting...\n\n";
        return 0;
    }
    if (!output_edgelist(info, time_steps)) {
        cout << "Error outputting edgelists, exiting...\n\n";
        return 0;
    } 
    if (!main_dynamics_func(info, time_steps, hgi)) {
        cout << "Error with dynamics, exiting...\n\n";
        return 0;
    }
    return 0;
}