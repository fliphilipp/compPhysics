# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# approximate results
r = np.linspace(0,5,100);
Z = 2
cfunscreened = Z**3 * 4 * r**2 * np.exp(- 2 * Z * r)
Z = 27/16
cfoptimized = Z**3 * 4 * r**2 * np.exp(- 2 * Z * r)

# input file
filename = 'file-dist-nuc.dat'

# import data
data = np.loadtxt(filename)


# initial size of plot 
fig = plt.figure(figsize=(10,8))

plt.hist(data, bins=200, normed=1, label='Monte Carlo')
plt.plot(r, cfunscreened, 'r-', lw=2, label='central field approximation, unscreened')
plt.plot(r, cfoptimized, 'k--', lw=2, label='central field approximation, optimized')

plt.xlim([0,5])

# labels
plt.xlabel(r'$r$ / $[a.u.]$', fontsize=20)
plt.ylabel(r'probability $\rho(r)$', fontsize=20)

# legend
plt.legend(fontsize=16)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)


# display the plot
fig.tight_layout()
plt.savefig('fig-dist-nuc.pdf')
plt.show()
