#!/usr/bin/env python

import os, sys
import matplotlib
from matplotlib import rc
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm
import string
from matplotlib import colors

nooa = 100

circles = []

g = open(str(sys.argv[2]), 'r')
for line in g:
    circle = plt.Circle((float(line.split()[0]), float(line.split()[1])), radius=0.6, fc='g', alpha=0.5)
    circles.append(circle)
    

f = open(str(sys.argv[1]), 'r')


xcoords = []
ycoords = []
zcoords = []
hbonds = []

for line in f.readlines():
    x = float(line.split()[0])
    y = float(line.split()[1])
    #z = float(line.split()[2])
    hb = float(line.split()[2])

    xcoords.append(x)
    ycoords.append(y)
    #zcoords.append(z)
    hbonds.append(hb)	

plt.figure()

fig = plt.gcf()
for elem in circles:
    fig.gca().add_artist(elem)
xlist = xcoords
ylist = ycoords               #np.linspace(-2., 1., 100)
#zlist = zcoords               #np.linspace(-1., 1., 100)
hblist = hbonds

xi = np.linspace(min(xcoords), max(xcoords), 100)
yi = np.linspace(min(ycoords), max(ycoords), 100)

hbi = griddata (xlist, ylist, hblist, xi, yi)

colormap = colors.ListedColormap(['white', 'blue', 'green', 'yellow', 'orange', 'red'])

v = np.linspace(0, max(hblist), 6, endpoint=True)
CS = plt.contour(xi,yi,hbi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,hbi,v,cmap=colormap,alpha=1.0)
plt.colorbar(CS)
#plt.hlines(13.38, min(xcoords), max(xcoords), linewidth = 4)
#plt.hlines(34.2, min(xcoords), max(xcoords), linewidth = 4)

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
##plt.title('Density Contour', size = 'x-large')
#plt.xlabel(r'x-axis [\AA]', size = 'x-large')
#plt.ylabel(r'z-axis [\AA]', size = 'x-large')

#plt.clabel(CS, inline=1, fontsize=10)

plt.show()

