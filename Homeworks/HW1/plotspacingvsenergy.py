# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
#import matplotlib.pylab as plt
import numpy as np
import matplotlib.pylab as plt

# input file
filename = 'latticeparameter-vs-energy.dat'

# import data
data = np.loadtxt(filename)

# initial size of plot window
fig = plt.figure(figsize=(12,8))

# plot
plt.plot(data[:,0], data[:,1])

# labels
plt.xlabel(r'lattice spacing [Å]', fontsize=18)
plt.ylabel(r'potential energy [eV]', fontsize=18)

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.savefig('fig-latticeparameter-vs-potEnergy.pdf',bbox_inches='tight')
plt.show()
