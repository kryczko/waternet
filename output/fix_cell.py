#!/usr/bin/python

import sys, string

f = open(str(sys.argv[1]), 'r')
dim = int(sys.argv[2])
push = float(sys.argv[3])

g = open('new_data.dat', 'w')

x1 = []
y1 = []
val1 = []
val3 = []
val4 = []
val5 = []
val6 = []

x2 = []
y2 = []
val2 = []

count1 = 0
count2 = 0

if (dim == 1):
	for line in f:
		i = float(line.split()[0])
		if (i < 5.0):
			i += push
			x1.append(i)
			val1.append(float(line.split()[1]))
			val3.append(float(line.split()[2]))
			val4.append(float(line.split()[3]))
		else:
			x2.append(i)
			val2.append(float(line.split()[1]))
			val5.append(float(line.split()[2]))
                        val6.append(float(line.split()[3]))
    	for elem in x2:
		g.write("%f\t%f\t%f\t%f\n" % (elem, val2[count1], val5[count1], val6[count1]))
	        count1 += 1
    	for elem in x1:
	    	g.write("%f\t%f\t%f\t%f\n" % (elem, val1[count2], val3[count2], val4[count2]))
	        count2 += 1	

elif (dim == 2):
    for line in f:
    	i = float(line.split()[1])
    	if (i < 5.0):
    		i += push
    		y1.append(i)
	        x1.append(float(line.split()[0]))
    		val1.append(float(line.split()[2]))
    	else:
    		y2.append(i)
		x2.append(float(line.split()[0]))
    		val2.append(float(line.split()[2]))

    for elem in x2:
        g.write("%f\t%f\t%f\n" % (elem, y2[count1], val2[count1]))
        count1 += 1
    for elem in x1:
        g.write("%f\t%f\t%f\n" % (elem, y1[count2],val1[count2]))
        count2 += 1
	
