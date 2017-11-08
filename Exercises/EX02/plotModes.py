# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'mode_tavg.dat'

# import data
data1 = np.loadtxt(filename1)

# initial size of plot window
plt.figure(figsize=(12,8))

##########################################################
# plot modes energy
for i in range(1,data1.shape[1]):
    plt.loglog(data1[:,0], data1[:,i],'-',alpha=0.5)
# plt.loglog(data1[:,0],np.sum(data1[:,1:], axis=1),':')

# labels
plt.xlabel(r'Time $t$', fontsize=18)
plt.ylabel(r'Time average of energy in modes, $\langle E_k \rangle_t$', fontsize=18)
plt.title(r'$\alpha = 0.1 \rightarrow$ equipartition of energy!', fontsize=24)
#plt.title(r'$\alpha = 0.01 \rightarrow$ no equipartition of energy!', fontsize=24)

# axis limits
# plt.ylim([10^-30,10^2])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.tight_layout()
plt.savefig('task4-time-averages.pdf')
plt.show()
