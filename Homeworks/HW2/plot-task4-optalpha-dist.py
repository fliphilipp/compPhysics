# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# import data
optalpha = np.loadtxt('file-optimized-alphas-10e8-98-100.dat')
# optalpha = np.loadtxt('file-optimized-alphas-5x10e6-95-1000.dat')
optalpha = optalpha[optalpha > 0.13]
optalpha = optalpha[optalpha < 0.15]
std = np.std(optalpha)
mean = np.mean(optalpha)

# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.hist(optalpha, bins=25, label=r'distribution of optimized $\alpha$')
lab = r"optimal $\alpha = %.4f \pm %.4f$" % (mean, std)
[y1, y2] = plt.gca().get_ylim()
plt.plot([mean, mean], [y1, y2], 'r--', lw=2, label=lab)

# labels
plt.xlabel(r'optimized $\alpha$', fontsize=20)
plt.ylabel(r'frequency', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.legend(fontsize=14)


# display the plot
fig.tight_layout()
plt.savefig('fig-opt-alpha-distribution.pdf')
plt.savefig('fig-opt-alpha-distribution.pdf')
plt.show()
