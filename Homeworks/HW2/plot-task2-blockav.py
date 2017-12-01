# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-block_s.dat'

# import data
data = np.loadtxt(filename)


# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.plot(data[:,0],data[:,1],'k-')

# labels
plt.xlabel(r'block size', fontsize=20)
plt.ylabel(r'statistical inefficiency', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)


# display the plot
fig.tight_layout()
plt.savefig('fig-blockav.pdf')
plt.show()
