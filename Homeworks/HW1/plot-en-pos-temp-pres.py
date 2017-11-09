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

# initial size of plot 
fig = plt.figure(figsize=(12,8))

##########################################################
# plot energies
ax11left = plt.subplot2grid((3, 2), (0, 0), colspan=1)

ax11left.plot(data1[:,0], data1[:,1],'b-')
ax11left.set_xlabel(r'time $[ps]$', fontsize=12)
ax11left.set_ylabel(r'$E_{pot}$ $[eV]$', fontsize=12, color='b')
ax11left.tick_params('y', colors='b')

ax11right = ax11left.twinx()
ax11right.plot(data1[:,0], data1[:,2],'r-')
ax11right.set_ylabel(r'$E_{kin}$ $[eV]$', fontsize=12, color='r')
ax11right.tick_params('y', colors='r')
ax11right.set_title(r'$E_{pot}$ and $E_{kin}$')

ax12 = plt.subplot2grid((3, 2), (0, 1), colspan=1)
ax12.plot(data1[:,0], data1[:,3],'k-')
ax12.set_xlabel(r'time $[ps]$', fontsize=12)
ax12.set_ylabel(r'$E_{pot}$ $[eV]$', fontsize=12, color='k')
ax12.set_title(r'$E_{tot}$')

##########################################################
# plot particle positons for atoms 1 and 2
ax2 = fig.add_subplot(312)

ax2.plot(data4[:,0], data4[:,1],'-', color='green', label='atom1, dim1')
ax2.plot(data4[:,0], data4[:,2],'-', color='blue', label='atom1, dim2')
ax2.plot(data4[:,0], data4[:,3],'-', color='orange', label='atom1, dim3')

ax2.plot(data4[:,0], data4[:,4],':', color='green', label='atom2, dim1')
ax2.plot(data4[:,0], data4[:,5],':', color='blue', label='atom2, dim2')
ax2.plot(data4[:,0], data4[:,6],':', color='orange', label='atom2, dim3')

ax2.legend()

# labels
ax2.set_xlabel(r'time $t$', fontsize=12)
ax2.set_ylabel(r'atoms 1 and 2 $[Å]$', fontsize=12)
ax2.set_title(r'particle positions')


##########################################################
# plot temperature in Kelvin and Pressure in kPa
ax31 = fig.add_subplot(313)

ax31.plot(data1[:,0], data1[:,4],'b-', label=r'temp')
ax31.set_xlabel(r'time $t$', fontsize=12)
ax31.set_ylabel(r'temperature $[K]$', fontsize=12, color='b')
ax31.tick_params('y', colors='b')

ax32 = ax31.twinx()
ax32.plot(data1[:,0], data1[:,5],'r:', label=r'pressure')
ax32.set_ylabel(r'pressure $[kPa]$', fontsize=12, color='r')
ax32.tick_params('y', colors='r')
ax32.set_title(r'temperature and pressure')


# display the plot
fig.tight_layout()
plt.savefig('fig-energy-positions-temp.pdf')
plt.show()
