# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-alpha-eloc.dat'

# import data
data = np.loadtxt(filename)
alphaval = data[:,0]
locen = data[:,1]
lower = data[:,2]
upper = data[:,3]


# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.plot(alphaval, lower, 'b:', label='lower bound standard error')
plt.plot(alphaval, upper, 'b:', label='upper bound standard error')
plt.plot(alphaval, locen, 'k-', label='energy')

# labels
plt.xlabel(r'alpha value', fontsize=20)
plt.ylabel(r'average local energy', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.legend(fontsize=16)


# display the plot
fig.tight_layout()
plt.savefig('fig-different-alpha.pdf')
plt.show()
