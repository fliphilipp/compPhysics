# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
#import matplotlib.pylab as plt
import numpy as np
import matplotlib.pylab as plt

# input file
filename1 = 'latticeparameter-vs-energy-widerange.dat'
filename2 = 'latticeparameter-vs-energy-narrowrange.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)

# initial size of plot window
fig = plt.figure(figsize=(12,5))

# plot
plt.subplot(1,2,1)
plt.plot(data1[:,0], data1[:,1])
# labels
plt.xlabel(r'lattice spacing [Å]', fontsize=18)
plt.ylabel(r'potential energy [eV]', fontsize=18)
# tick fontsize
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.subplot(1,2,2)
plt.plot(data2[:,0], data2[:,1])
plt.scatter(4.0323, -860.2721,s=150, color='red')
plt.arrow(4.5, 0, 0.85*(4.0323-4.5), 0.85 * (-860.2721),head_width=0.1, head_length=60, fc='k', ec='k')
plt.text(3.6, 0, 'lattice parameter = 4.0323 Å \n     potential = -860.27 eV',fontsize=15)
# labels
plt.xlabel(r'lattice spacing [Å]', fontsize=18)
plt.ylabel(r'potential energy [eV]', fontsize=18)
# tick fontsize
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)


# display the plot
plt.tight_layout()
plt.savefig('fig-latticeparameter-vs-potEnergy.pdf',bbox_inches='tight')
plt.show()
