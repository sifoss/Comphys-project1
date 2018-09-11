import numpy as np, matplotlib.pyplot as plt
import pandas as pd

def exercise_b():
    for i in range(3):
        infile = np.loadtxt('u_exact_approx_%d.txt'% 10**(i+1))

        x = np.asarray(infile[:, 0])
        approx = np.asarray(infile[:, 2])
    
        plt.plot(x, approx, label = 'n = %d' % 10**(i+1))

    exact = np.asarray(infile[:, -1])

    plt.plot(x, exact, 'r--', label = 'Analytical solution')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('comparison of numerical solution with analytical solution')
    plt.legend()
    plt.show()
    
exercise_b()



