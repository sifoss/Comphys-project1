#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

mat = np.loadtxt("outfile.txt")

"""
infile = open("outfile.txt", "r")

x = []
u = []
d = []

for line in infile: 
    x.append(float(line.split()[0]))
    u.append(float(line.split()[1]))
    d.append(float(line.split()[2]))
""" 

import matplotlib.pyplot as plt
plt.plot(mat[:,0],mat[:,1],mat[:,0],mat[:,2])
plt.show()