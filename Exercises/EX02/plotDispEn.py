# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'disp.dat'
filename2 = 'energy.dat'
filename3 = 'modes.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)

# initial size of plot window
plt.figure(figsize=(12,8))

##########################################################
# plot trajectories
plt.subplot(3,1,1)
for i in range(1,data1.shape[1]):
    plt.plot(data1[:,0], data1[:,i],'-')

# labels
plt.xlabel(r'time $t$', fontsize=12)
plt.ylabel(r'displacement $q_i$', fontsize=12)

#axis limits
plt.xlim([0,500])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

##########################################################
# plot energies
plt.subplot(3,1,2)
plt.plot(data2[:,0], data2[:,1],':',label='potential energy',alpha=0.5)
plt.plot(data2[:,0], data2[:,2],':',label='kinetic energy',alpha=0.5)
plt.plot(data2[:,0], data2[:,3],'-',label='total energy')
#plt.plot(data2[:,0], np.ones((data2.shape[0],1)) * (data1.shape[1] - 1),'-',label='total energy')

# labels
plt.xlabel(r'time $t$', fontsize=12)
plt.ylabel(r'total energies', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
#plt.xlim([0,10000])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

##########################################################
# plot energy modes
plt.subplot(3,1,3)
for i in range(1,data3.shape[1]):
    plt.plot(data3[:,0], data3[:,i],'-',alpha=0.5,label='mode ' + str(i))

# labels
plt.xlabel(r'time $t$', fontsize=12)
plt.ylabel(r'energy in modes $E_k(t)$', fontsize=12)
plt.title(r'$\alpha = 0.01 \longrightarrow$ energy transfer among normal modes!', fontsize=16)

# legend
plt.legend(loc='upper right')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.tight_layout()
plt.savefig('task3-disp-energy-modes-001.pdf')
plt.show()
