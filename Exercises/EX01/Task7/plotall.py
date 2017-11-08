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
plt.plot(data1[:,0], data1[:,1],'r--',label='Oxygen atom 1')
plt.plot(data1[:,0], data1[:,2],'k-',label='Carbon atom')
plt.plot(data1[:,0], data1[:,3],'r:',label='Oxygen atom 2')

# labels
plt.xlabel(r'Time / $[ps]$', fontsize=12)
plt.ylabel(r'Displacement / [Å]', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
plt.xlim([0,0.3])

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
plt.xlabel(r'Time / $[ps]$', fontsize=12)
plt.ylabel(r'Energy / $[eV]$', fontsize=12)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
plt.xlim([0,0.3])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

##########################################################
# plot power spectrum
plt.subplot(3,1,3)
plt.plot(data4[:,0], data4[:,1],'r--',label='Oxygen atom 1')
plt.plot(data4[:,0], data4[:,2],'k-',label='Carbon atom')
plt.plot(data4[:,0], data4[:,3],'r:',label='Oxygen atom 2')
freqs = data4[:,0]
freqs1 = freqs[(freqs > 25) & (freqs < 50)]
freqs2 = freqs[freqs > 50]
arr = data4[:,1] + data4[:,2] + data4[:,3]
arr1 = np.array(arr[(freqs > 25) & (freqs < 50)])
arr2 = np.array(arr[freqs > 50])
maxfreq1 = freqs1[np.argmax(arr1)]
maxfreq2 = freqs2[np.argmax(arr2)]
print('Max intensities at frequencies ', maxfreq1, ' and ', maxfreq2)

# labels
plt.xlabel(r'Frequency / [$ps^{-1}$]', fontsize=12)
plt.ylabel('Power spectrum', fontsize=12)

plt.text(maxfreq1-7, max(arr1)*0.6, r'$39.063 \frac{1}{ps}$',fontsize=14)
plt.text(maxfreq2-7, max(arr2), r'$74.769 \frac{1}{ps}$',fontsize=14)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

# axis limits
plt.xlim([-100,100])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot
plt.tight_layout()
plt.savefig('task7-disp-energy-pspec.pdf')
plt.show()
