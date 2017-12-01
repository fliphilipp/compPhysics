# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-block_s.dat'

# import data
data = np.loadtxt(filename)
blockav = data[0:len(data)-2,:]
stat_ineff = data[len(data)-1,1]


# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.plot(blockav[:,0],blockav[:,1], 'k-', label='block statistical inefficiency')
lab = "limit for large block size = %.2f" % stat_ineff
plt.plot([min(blockav[:,0]), max(blockav[:,0])], [stat_ineff, stat_ineff], 'r:', lw=2, label=lab)

# labels
plt.xlabel(r'block size', fontsize=20)
plt.ylabel(r'statistical inefficiency', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.legend(fontsize=16)


# display the plot
fig.tight_layout()
plt.savefig('fig-blockav.pdf')
plt.show()
