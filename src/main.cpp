// global includes
#include <iostream>
#include <fstream>
#include <string>

// local includes
#include "storage.h"
#include "read_config.h"
#include "read_data.h"
#include "edgelist.h"

using namespace std;

// declaration of data structures
Information info;
TimeSteps time_steps;

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
    return 0;
}