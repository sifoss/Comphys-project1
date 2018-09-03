import numpy as np
import time

n = 10

x = np.linspace(0,1,n+2)
h = 1./(n+1)


# diagonalene i matrisen
a = np.zeros(n)-1
b = np.zeros(n)+2
c = np.zeros(n)-1

a[0]  = 0
c[-1] = 0

f = 100*np.exp(-10*x)
u = 1- (1-np.exp(-10))*x - np.exp(-10*x)

d = h**2*f
d[0]  = 0
d[-1] = 0


start = time.time()

for i in range(1,n):
	cof   = a[i]/b[i-1]
	b[i] -= cof*c[i-1]
	d[i+1] -= cof*d[i]




"""for i in range(1,n):
	#cof      = c[-i-1]/b[-i]
	#c[-i-1] -= cof*b[-i]
	#print d[-i-2]
	d[-i-2] -= (c[-i-1]/b[-i])*d[-i-1]
	#d[-i-2] = (d[-i-2] - (c[-i-1]/b[-i])*d[-i-1])/b[-i]
	#d[n-1] = d[n-1]-(c[n-2]/b[n-1])*d[n-1]
	print(i)"""

for i in range(n-1,0,-1):
	d[i] -= (c[i-1]/b[i])*d[i+1]

print d

d[1:-1] = d[1:-1]/b
end = time.time()

print(end-start)

import matplotlib.pyplot as plt
plt.plot(x,u,x,d)
plt.show()