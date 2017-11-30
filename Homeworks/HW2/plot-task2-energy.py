# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-locEn.dat'

# import data
data = np.loadtxt(filename)


# initial size of plot 
fig = plt.figure(figsize=(8,3))

plt.plot(data,'k-')

# labels
plt.xlabel(r'time step', fontsize=20)
plt.ylabel(r'local Energy', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)


# display the plot
fig.tight_layout()
plt.savefig('fig-locEn.pdf')
plt.show()
