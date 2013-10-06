import matplotlib.pyplot as plt
import numpy as np
import os


f = open('../Project2-build/output_eigenstates.txt', 'r')

i = 0
for line in f:
    if i == 0:
        column = line.split()
        m = int(column[0])
        n = int(column[1])
        rho_max = float(column[2])
        omega_r = float(column[3])

        h = rho_max/(n+1)

        v = np.zeros((m,n))
        i = 1
    else:
        column = line.split()
        for j in range(m):
            v[j,i-1] = float(column[j])
        
        i += 1

x = np.linspace(0, rho_max, n)

plt.plot(x, v[0])
plt.hold('on')
for j in range(1,m):
    plt.plot(x, v[j])

plt.title(r'$\omega_r$ = %.2f' % omega_r)
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\psi(\rho)$')
plt.legend([(r'$\psi_{%d}(\rho)$' % n) for n in range(m)])
plt.show()
        
