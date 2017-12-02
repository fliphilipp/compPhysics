# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# import data
data50 = np.loadtxt('file-alphaval50.dat')
data75 = np.loadtxt('file-alphaval75.dat')
data90 = np.loadtxt('file-alphaval90.dat')
data99 = np.loadtxt('file-alphaval99.dat')

iteration = data50[:,0];
alpha50 = data50[:,1];
alpha75 = data75[:,1];
alpha90 = data90[:,1];
alpha99 = data99[:,1];

# initial size of plot 
fig = plt.figure(figsize=(8,6))

plt.plot(iteration, alpha50,'k-',label=r'$\beta = 0.5$')
plt.plot(iteration, alpha75,'g-.',label=r'$\beta = 0.75$')
plt.plot(iteration, alpha90,'r--',label=r'$\beta = 0.9$')
plt.plot(iteration, alpha99,'b:',label=r'$\beta = 0.99$')

plt.ylim(0.05, 0.3)

# labels
plt.xlabel(r'time step', fontsize=20)
plt.ylabel(r'alpha', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.legend(fontsize=14)


# display the plot
fig.tight_layout()
plt.savefig('fig-alpha-optimization.pdf')
plt.savefig('fig-alpha-optimization.png')
plt.show()
