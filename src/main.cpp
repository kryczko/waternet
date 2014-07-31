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
    cout << "\033[1;46m\033[1;37m\n\n";
    cout << " ------------------------------------------------------------------------ \n";      
    cout << "|  W         W    A  TTTTTTTTT EEEEEEE R R R  N      N EEEEEEE TTTTTTTTT |\n";
    cout << "| W           W  A A     T     E       R   R  N N    N E           T     |\n";
    cout << "| W     W     W A   A    T     E       R  R   N  N   N E           T     |\n";
    cout << "|  W   W W   W AAAAAAA   T     EEEEEEE R R    N   N  N EEEEEEE     T     |\n";
    cout << "|   W W   W W A       A  T     E       R   R  N    N N E           T     |\n";
    cout << "|    W     W A         A T     EEEEEEE R    R N      N EEEEEEE     T     |\n";
    cout << " ------------------------------------------------------------------------ \n\n";
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
    if (!output_edgelist(info, time_steps)) {
        cout << "Error outputting edgelists, exiting...\n\n";
        return 0;
    } 
    if (!main_dynamics_func(info, time_steps, hgi)) {
        cout << "Error with dynamics, exiting...\n\n";
        return 0;
    }
    cout << "\033[0m";
    return 0;
}