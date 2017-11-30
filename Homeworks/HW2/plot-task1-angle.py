# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-cosAngle.dat'

# import data
data = np.loadtxt(filename)


# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.hist(data, bins=100, normed=True, label='Monte Carlo')

# labels
plt.xlabel(r'$cos(\theta)$', fontsize=20)
plt.ylabel(r'probability', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)


# display the plot
fig.tight_layout()
plt.savefig('fig-dist-angle.pdf')
plt.show()
