#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import time
from numba import jit

def Exact(x):
    return 1- (1-np.exp(-10))*x - np.exp(-10*x)

def F(x):
    return 100*np.exp(-10*x)

#@jit
def general_elimination_find_x(a,b,c,d):
    start = time.time()
    for i in range(1,n):    #foreward substitution
	    cof   = a[i]/b[i-1]
	    b[i] -= cof*c[i-1]
	    d[i+1] -= cof*d[i]
    
    for i in range(n-1,0,-1):           #backward substitution
	    d[i] -= (c[i-1]/b[i])*d[i+1]
    
    end = time.time()
    
    d[1:-1] = d[1:-1]/b
    
    print(end-start)
    
    return d

n = 10
x = np.linspace(0,1,n+2)
h = 1./(n+1)


# de tre diagonalenei dentridiagonale matrisen, fra venstre til høyre
a = np.zeros(n)-1
b = np.zeros(n)+2
c = np.zeros(n)-1

a[0]  = 0 #randbetingelser
c[-1] = 0

analytical = Exact(x)

b_vec = h**2*F(x)  #b-vector i ligningen Ax=b
b_vec[0]  = 0
b_vec[-1] = 0


x_vec = general_elimination_find_x(a,b,c,b_vec) # x i likningen Ax=b

import matplotlib.pyplot as plt
plt.plot(x,analytical,x,x_vec)
plt.show()