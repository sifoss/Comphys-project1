#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:15:59 2018

@author: gunnar
"""
import numpy as np

n = 10

x = np.linspace(0,1,n+2)

h = 1./(n+1)

print(h,x[2]-x[1])

a = np.zeros(n)-1
b = np.zeros(n)+2
c = np.zeros(n)-1

a[0]  = 0
c[-1] = 0

f = 100*np.exp(-10*x)
u = 1 - (1-np.exp(-10))*x - np.exp(-10*x)




import matplotlib.pyplot as plt
plt.plot(x,u)
plt.show()
