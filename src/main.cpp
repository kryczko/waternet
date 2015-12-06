// global includes
#include <iostream>
#include <fstream>
#include <string>

// local includes
#include "storage.h"
#include "read_config.h"
#include "read_data.h"
#include "edgelist.h"
#include "analysis.h"

using namespace std;

// declaration of data structures
Information info;
TimeSteps time_steps;

void print_welcome_message() {
    cout << "\x1b[1;44;40m\n\n";
cout << "                _    _  ___ _____ ___________ _   _  _____ _____\n";
cout << "               | |  | |/ _ \\_   _|  ___| ___ \\ \\ | ||  ___|_   _|\n";
cout << "               | |  | / /_\\ \\| | | |__ | |_/ /  \\| || |__   | |\n";  
cout << "               | |/\\| |  _  || | |  __||    /| . ` ||  __|  | |\n";  
cout << "               \\  /\\  / | | || | | |___| |\\ \\| |\\  || |___  | |\n";  
cout << "                \\/  \\/\\_| |_/\\_/ \\____/\\_| \\_\\_| \\_/\\____/  \\_/\n\n";  
cout << "              WRITTEN BY K. RYCZKO, WITH HELP FROM I. TAMBLYN             \n";
cout << "              UNIVERSITY OF ONTARIO INSTITUTE OF TECHNOLOGY               \n";
cout << "            COMPUTATIONAL LABORATORY FOR ENERGY AND NANOSCIENCE           \n\n\n";
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
        cout << "Error creating edgelists, some modules will fail...\n\n";
    }
    if (!main_analysis(info, time_steps)) {
        cout << "Error running analysis, exiting...\n\n";
        return 0;
    } 
    cout << "\033[0m";
    return 0;
}