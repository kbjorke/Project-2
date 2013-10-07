import matplotlib.pyplot as plt
import numpy as np
import os


f = open('../Project2-build/output_stability.txt', 'r')

eigenvalues = [3.0, 7.0, 11.0]

i = 0
for line in f:
    column = line.split()
    if i == 0:
        resolution = int(column[0])
        rho_max = float(column[1])
        m = int(column[2])

        n_list = np.zeros(resolution)
        error = np.zeros((m,resolution))
        counter_transformations = np.zeros(resolution)
        i = 1
    else:
        n_list[i-1] = int(column[0])
        counter_transformations[i-1] = int(column[m+1])
        for j in range(m):
            error[j,i-1] = np.abs(float(column[j+1]) -
                                  eigenvalues[j])/eigenvalues[j]
        
        i += 1

error = np.log10(error)

plt.figure(1)
plt.plot(n_list, error[0])
plt.hold('on')
for j in range(1,m):
    plt.plot(n_list, error[j])

#plt.title(r'$\rho_{max} = %.2f$' % rho_max)
plt.xlabel(r'n',fontsize=16)
plt.ylabel(r'$log(|\epsilon_{rel}|)$', fontsize=20)
plt.legend([(r'$\lambda_{%d}$' % n) for n in range(m)])
plt.savefig('stability_6_householder.eps', format='eps')

fit = np.polyfit(n_list, counter_transformations, 2)
print fit
counter_fit = np.poly1d(fit)
print counter_fit

plt.figure(2)
plt.plot(n_list, counter_transformations,
        n_list, counter_fit(n_list))
plt.xlabel('n', fontsize=16)
plt.ylabel('Transformations', fontsize=16)
plt.savefig('transform.eps', format='eps')

plt.show()
        
