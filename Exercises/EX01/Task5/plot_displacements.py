# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'disp.dat'
filename2 = 'pe.dat'
filename3 = 'ke.dat'
filename4 = 'pspec.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)

# initial size of plot window
plt.figure(figsize=(12,8))

##########################################################
# plot trajectories
plt.subplot(3,1,1)
plt.plot(data1[:,0], data1[:,1],'-',label='Atom 1')
plt.plot(data1[:,0], data1[:,2],'-',label='Atom 2')
plt.plot(data1[:,0], data1[:,3],'-',label='Atom 3')

# labels
plt.xlabel('Time / [dim. unit]', fontsize=12)
plt.ylabel('Displacement / [dim. unit]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
plt.xlim([0,50])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

##########################################################
# plot energies
plt.subplot(3,1,2)
plt.plot(data2[:,0], data2[:,1],'--',label='potential energy')
plt.plot(data3[:,0], data3[:,1],':',label='kinetic energy')
plt.plot(data2[:,0], data2[:,1] + data3[:,1],'-',label='total energy')

# labels
plt.xlabel('Time / [dim. unit]', fontsize=12)
plt.ylabel('Energy / [dim. unit]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
plt.xlim([0,50])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

##########################################################
# plot power spectrum
plt.subplot(3,1,3)
plt.plot(data4[:,0], data4[:,1],'-',label='Atom 1')
plt.plot(data4[:,0], data4[:,2],'--',label='Atom 2')
plt.plot(data4[:,0], data4[:,3],':',label='Atom 3')

# labels
plt.xlabel('Frequency / [dim. unit]', fontsize=12)
plt.ylabel('Power spectrum / [dim. unit]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
plt.xlim([-0.4,0.4])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.tight_layout()
plt.savefig('task5-energy-pspec.pdf')
plt.show()
