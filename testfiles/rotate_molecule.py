#!/usr/bin/env python

import numpy as np
import random
from ase import Atoms
from ase.io import write

# starting config
# O  5.44223  5.61731  4.95808
# H  6.39922  5.49052  5.04453
#H  5.35756  6.27071  4.24718

def pbc_round(input):
     i = int(input)
     if (abs(input - i) >= 0.5):
         if (input > 0):
             i += 1
         if (input < 0):
             i -= 1
     return i

vec_pbc_round = np.vectorize(pbc_round)

def giveRandRot(coords):
    theta = random.uniform(0, 2 * np.pi)
    R_x = np.array([[1., 0., 0.], [0., np.cos(theta), -np.sin(theta)], [0., np.sin(theta), np.cos(theta)]])
    R_y = np.array([[np.cos(theta), 0., np.sin(theta)], [0., 1., 0.], [-np.sin(theta), 0., np.cos(theta)]])
    R_z = np.array([[np.cos(theta), -np.sin(theta), 0.], [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
    for i in range(3):
        coords[i] = np.matmul(R_x, coords[i])
        # coords[i] = np.matmul(R_y, coords[i])
        # coords[i] = np.matmul(R_z, coords[i])
    return coords

Os_coord = np.array([0.0, 0.0, 0.0])
H1s_coord = np.array([1.0, 0.0, 0.0])
H2s_coord = np.array([-0.25, 0.96,0.0])
coords = np.array([Os_coord, H1s_coord, H2s_coord])

f = open('orientation_test.xyz', 'w')
distrox = np.zeros(180)
distroy = np.zeros(180)
distroz = np.zeros(180)

for i in range(50000):
    f.write('3\n\n')
    new_coords = giveRandRot(coords)
    atoms = Atoms('OH2', coords, cell=[20., 20., 20.], pbc=[1,1,1])
    atoms.wrap()
    new_coords = atoms.get_positions()
    f.write('O\t%f\t%f\t%f\n' % (new_coords[0][0], new_coords[0][1], new_coords[0][2]))
    f.write('H\t%f\t%f\t%f\n' % (new_coords[1][0], new_coords[1][1], new_coords[1][2]))
    f.write('H\t%f\t%f\t%f\n' % (new_coords[2][0], new_coords[2][1], new_coords[2][2]))
    dr1 = new_coords[1] - new_coords[0]
    dr2 = new_coords[2] - new_coords[0]
    dr1 -= vec_pbc_round(dr1 / 20.) * 20.
    dr2 -= vec_pbc_round(dr2 / 20.) * 20.
    totr = dr1 + dr2
    norm = np.sqrt(totr[0]*totr[0] + totr[1]*totr[1] + totr[2]*totr[2])
    angle = np.arccos(totr[0] / norm) * 57.2958
    distrox[int(np.arccos(totr[0] / norm) * 57.2958)] += 1
    distroy[int(np.arccos(totr[1] / norm) * 57.2958)] += 1
    distroz[int(np.arccos(totr[2] / norm) * 57.2958)] += 1

for i in range(len(distrox)):
    print i, distrox[i] / np.sum(distrox), distroy[i] / np.sum(distroy), distroz[i] / np.sum(distroz)
