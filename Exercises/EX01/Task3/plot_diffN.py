# plot the powerspectrum
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input files function
filename1 = 'functionf2.dat'
filename2 = 'functionN258.dat'
filename3 = 'functionN58.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)

# initial size of plot window
plt.figure(figsize=(12,8))

# function subplot
plt.subplot(211)

# plot
plt.plot(data1[:,0], data1[:,1],'k-',label=r'$N=250$')
plt.plot(data2[:,0], data2[:,1],'b--',label=r'$N=258$')
plt.plot(data3[:,0], data3[:,1],'r:',label=r'$N=58$')

# labels
plt.xlabel('t / [arb. unit]', fontsize=16)
plt.ylabel('h(t) / [arb. unit]', fontsize=16)

# legend
plt.legend(loc='center left')

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

######################################################
# input files powerspectrum
filename1 = 'powerspectrumf2.dat'
filename2 = 'powerspectrumN258.dat'
filename3 = 'powerspectrumN58.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)

# powerspectrum subplot
plt.subplot(212)

# plot
plt.plot(data1[:,0], data1[:,1],'k-',label=r'$N=250$')
plt.plot(data2[:,0], data2[:,1],'b--',label=r'$N=258$')
plt.plot(data3[:,0], data3[:,1],'r:',label=r'$N=58$')

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

plt.savefig('diffN.pdf')
plt.show()
