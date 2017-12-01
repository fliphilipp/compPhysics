# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'file-corrfunc.dat'

# import data
data = np.loadtxt(filename)
corrfunc = data[0:len(data)-2]
stat_ineff = int(data[len(data)-1])


# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.plot(corrfunc, 'k-', label='correlation function')
lab = "statistical inefficiency = %d" % stat_ineff
plt.plot([stat_ineff, stat_ineff], [-0.05, 1.05], 'r:', label=lab)

# labels
plt.xlabel(r'sequence length', fontsize=20)
plt.ylabel(r'correlation', fontsize=20)

plt.ylim(-0.05, 1.05)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.legend(fontsize=16)


# display the plot
fig.tight_layout()
plt.savefig('fig-autocorr.pdf')
plt.show()
