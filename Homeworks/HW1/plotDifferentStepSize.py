# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'energy.dat'
filename2 = 'energy10e-2.dat'
filename3 = 'energy10e-3.dat'
filename4 = 'pos1.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)

# initial size of plot window
plt.figure(figsize=(12,8))

##########################################################
# plot for dt = 0.1
plt.subplot(2,1,1)

plt.plot(data1[:,0], data1[:,1],'r:', label='potential energy')
plt.plot(data1[:,0], data1[:,2],'b:', label='kinetic energy')
plt.plot(data1[:,0], data1[:,3],'k-', label='total energy')

plt.legend()

# labels
plt.xlabel(r'time $t$', fontsize=12)
plt.ylabel(r'energy', fontsize=12)
plt.title(r'energies')

#axis limits
#plt.xlim([0.9,1])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


##########################################################
# plot for dt = 0.0001
plt.subplot(2,1,2)

plt.plot(data4[:,0], data4[:,1],'r-', label='atom1, dim1')
plt.plot(data4[:,0], data4[:,2],'b-', label='atom1, dim2')
plt.plot(data4[:,0], data4[:,3],'k-', label='atom1, dim3')

plt.plot(data4[:,0], data4[:,4],'r:', label='atom2, dim1')
plt.plot(data4[:,0], data4[:,5],'b:', label='atom2, dim2')
plt.plot(data4[:,0], data4[:,6],'k:', label='atom2, dim3')

plt.legend()

# labels
plt.xlabel(r'time $t$', fontsize=12)
plt.ylabel(r'position atom1', fontsize=12)
plt.title(r'position of first atom')

#axis limits
# plt.xlim([0,500])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.tight_layout()
plt.savefig('fig-different-dt.pdf')
plt.show()
