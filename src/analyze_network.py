#!/usr/bin/env python

import matplotlib.pyplot as plt
import networkx as nx
from networkx import *
import time
import sys

t0 = time.clock()


def compute_clustering_coefficient(Directed_G, Undirected_G):
    return nx.average_clustering(Undirected_G)
    
def average_shortest_path(Directed_G, Undirected_G):
    t1 = time.clock()
    path_sum = 0
    num_paths = 0
    no_paths = 0
    counter = 0
    for node0 in Undirected_G:
        counter += 1
        for node1 in Undirected_G:
            if (nx.has_path(Undirected_G, node0, node1) and node0 != node1):
                num_paths += 1
                sp = nx.shortest_path_length(Undirected_G, node0, node1)
                path_sum += sp
            else:
                no_paths += 1
    t2 = time.clock()
    if (num_paths != 0):
        return (path_sum / float(num_paths))

def report_final_stats():
    t6 = time.clock()
    print '%f seconds elapsed in total\n' % (t6-t0)

def run_main():
    file = str(sys.argv[1])
    f = open(file, 'r')
    print "\nReading inputfile:", file, "..."
    
    total_edgelist = []
    for line in f.readlines():
        if ("Frame:" in str(line)):
            num_edges = int(line.split()[2])
            edgelist = []
        else:
            edgelist.append((int(line.split()[0]), int(line.split()[1])))
            if (num_edges == len(edgelist)):
                total_edgelist.append(edgelist)
    
    average_clustering = 0
    average_path_length = 0
    average_diameter = 0
    average_number_of_edges = 0
    average_number_of_nodes = 0
    timesteps = len(total_edgelist)
    counter = 0
    #plt.ion()
    #plt.show()
    ac = []
    apl = []
    dia = []
    ae = []
    an = []
    ad = []
    second_counter = 0
    for edgelist in total_edgelist: 
        if (second_counter % 10 == 0):
            Directed_G = nx.DiGraph(edgelist)
            Undirected_G = Directed_G.to_undirected()
            #plt.figure(figsize=(8,8))
            #nx.draw(Directed_G,pos=nx.spring_layout(Directed_G))
            #plt.draw()
            #time.sleep(0.1)

            # compute other things
            average_clustering += compute_clustering_coefficient(Directed_G, Undirected_G)
            average_path_length += average_shortest_path(Directed_G, Undirected_G)
            #average_diameter += nx.diameter(Undirected_G);
            average_number_of_edges += nx.number_of_edges(Undirected_G);
            average_number_of_nodes += nx.number_of_nodes(Undirected_G);
            average_degree = float(average_number_of_edges) / average_number_of_nodes
            print "Timestep: ", second_counter
            #print "Average diameter:", average_diameter / (counter + 1)
            #dia.append(average_diameter / (counter + 1))
            print "Average clustering coefficient:", average_clustering / (counter + 1)
            ac.append(average_clustering / (counter + 1))
            print "Average path length:", average_path_length / (counter + 1)
            apl.append(average_path_length / (counter + 1))
            print "Average number of edges:", float(average_number_of_edges) / (counter + 1)
            ae.append(average_number_of_edges / (counter + 1))
            print "Average number of nodes involved in network:", float(average_number_of_nodes) / (counter + 1)
            an.append(average_number_of_nodes / (counter + 1))
            print "Degrees:\n[0, 1, 2, 3, 4 ,5]"
            ad.append(average_degree)
            print nx.degree_histogram(Undirected_G)
           
            report_final_stats()
            counter += 1
        second_counter += 1
    opl = open("clustering_time.dat", 'w')
    i = 0
    for elem in ac:
        opl.write("%i\t%f\n" % (i, elem))
        i += 1
    avpl = open("path_length_time.dat", 'w')
    j = 0
    for elem in apl:
        avpl.write("%i\t%f\n" % (j, elem))
        j += 1

def main():
    run_main()
if __name__ == "__main__":
    main()
