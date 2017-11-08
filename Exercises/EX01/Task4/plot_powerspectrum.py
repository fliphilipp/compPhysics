# plot the powerspectrum
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input files function
filename1 = 'function.dat'
filename2 = 'function005.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)

# initial size of plot window
plt.figure(figsize=(12,8))

# function subplot
plt.subplot(211)

# plot
plt.plot(data1[:,0], data1[:,1],'k-',label=r'$\Delta t = 0.1$')
plt.plot(data2[:,0], data2[:,1],'r:',label=r'$\Delta t = 0.05$')

# labels
plt.xlabel('t / [arb. unit]', fontsize=16)
plt.ylabel('h(t) / [arb. unit]', fontsize=16)

# legend
plt.legend()

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

######################################################
# input files powerspectrum
filename1 = 'powerspectrum.dat'
filename2 = 'powerspectrum005.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)

# powerspectrum subplot
plt.subplot(212)

# plot
plt.plot(data1[:,0], data1[:,1],'k-',label=r'$\Delta t = 0.1$')
plt.plot(data2[:,0], data2[:,1],'r:',label=r'$\Delta t = 0.05$')

plt.xlim([-8,8])

# labels
plt.xlabel('Frequency / [arb. unit]', fontsize=16)
plt.ylabel('Power spectrum / [arb. unit]', fontsize=16)

# legend
plt.legend()

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

# display the plot
plt.tight_layout()

plt.savefig('plot.pdf')
plt.show()
