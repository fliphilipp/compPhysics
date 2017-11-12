# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'energy-dt-1.dat'
filename2 = 'energy-dt-2.dat'
filename3 = 'energy-dt-3.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)

# initial size of plot 
fig = plt.figure(figsize=(10,8))

##########################################################
# plot for dt = 0.1
ax11left = plt.subplot2grid((3, 2), (0, 0), colspan=1)

ax11left.plot(data1[:,0], data1[:,1],'b-')
ax11left.set_xlabel(r'time $[ps]$', fontsize=12)
ax11left.set_ylabel(r'$E_{pot}$ $[eV]$', fontsize=12, color='b')
ax11left.tick_params('y', colors='b')

ax11right = ax11left.twinx()
ax11right.plot(data1[:,0], data1[:,2],'r-')
ax11right.set_ylabel(r'$E_{kin}$ $[eV]$', fontsize=12, color='r')
ax11right.tick_params('y', colors='r')
ax11right.set_title(r'potential and kinetic energies, $dt = 0.1$')

ax12 = plt.subplot2grid((3, 2), (0, 1), colspan=1)
ax12.plot(data1[:,0], data1[:,3],'k-')
ax12.set_xlabel(r'time $[ps]$', fontsize=12)
ax12.set_ylabel(r'$E_{tot}$ $[eV]$', fontsize=12, color='k')
ax12.set_title(r'total energy, $dt = 0.1$')

##########################################################
# plot for dt = 0.01
ax21left = plt.subplot2grid((3, 2), (1, 0), colspan=1)

ax21left.plot(data2[:,0], data2[:,1],'b-')
ax21left.set_xlabel(r'time $[ps]$', fontsize=12)
ax21left.set_ylabel(r'$E_{pot}$ $[eV]$', fontsize=12, color='b')
ax21left.tick_params('y', colors='b')

ax21right = ax21left.twinx()
ax21right.plot(data2[:,0], data2[:,2],'r-')
ax21right.set_ylabel(r'$E_{kin}$ $[eV]$', fontsize=12, color='r')
ax21right.tick_params('y', colors='r')
ax21right.set_title(r'potential and kinetic energies, $dt = 0.01$')

ax22 = plt.subplot2grid((3, 2), (1, 1), colspan=1)
ax22.plot(data2[:,0], data2[:,3],'k-')
ax22.set_xlabel(r'time $[ps]$', fontsize=12)
ax22.set_ylabel(r'$E_{tot}$ $[eV]$', fontsize=12, color='k')
ax22.set_ylim(-835, -833.5)
ax22.set_title(r'total energy, $dt = 0.01$')

##########################################################
# plot for dt = 0.01
ax31left = plt.subplot2grid((3, 2), (2, 0), colspan=1)

ax31left.plot(data3[:,0], data3[:,1],'b-')
ax31left.set_xlabel(r'time $[ps]$', fontsize=12)
ax31left.set_ylabel(r'$E_{pot}$ $[eV]$', fontsize=12, color='b')
ax31left.tick_params('y', colors='b')

ax31right = ax31left.twinx()
ax31right.plot(data3[:,0], data3[:,2],'r-')
ax31right.set_ylabel(r'$E_{kin}$ $[eV]$', fontsize=12, color='r')
ax31right.tick_params('y', colors='r')
ax31right.set_title(r'potential and kinetic energies, $dt = 0.001$')

ax32 = plt.subplot2grid((3, 2), (2, 1), colspan=1)
ax32.plot(data3[:,0], data3[:,1] + data3[:,2],'k-')
ax32.set_xlabel(r'time $[ps]$', fontsize=12)
ax32.set_ylabel(r'$E_{tot}$ $[eV]$', fontsize=12, color='k')
ax32.set_ylim(-832.5, -832.3)
ax32.set_title(r'total energy, $dt = 0.001$')


# display the plot
fig.tight_layout()
plt.savefig('fig-different-dt.pdf')
plt.show()
