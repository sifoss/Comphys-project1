import numpy as np
import time

n = 100000

x = np.linspace(0,1,n+2)
h = 1./(n+1)


# diagonalene i matrisen
a = -1
b = np.zeros(n)+2
c = -1


f = 100*np.exp(-10*x)
u = 1- (1-np.exp(-10))*x - np.exp(-10*x)

d = h**2*f
d[0]  = 0
d[-1] = 0


start = time.time()

for i in range(1,n):
	cof   = -1/b[i-1]
	b[i] += cof
	d[i+1] -= cof*d[i]

d[-2] = d[-2]/b[-1]

for i in range(1,n):
	d[-i-2] += d[-i-1]/b[-i]

d[1:-1] = d[1:-1]/b
end = time.time()

print(end-start)

#import matplotlib.pyplot as plt
#plt.plot(x,u,x,d)
#plt.show()