# plot the powerspectrum
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename1 = 'powerspectrumf2.dat'
filename2 = 'powerspectrumf1.dat'
filename3 = 'powerspectrumphi.dat'

# import data
data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
plt.plot(data1[:,0], data1[:,1],'k-',label=r'$f=2,\phi=0$')
plt.plot(data2[:,0], data2[:,1],'b-',label=r'$f=1,\phi=0$')
plt.plot(data3[:,0], data3[:,1],'r--',label=r'$f=2,\phi=\pi/2$')

# labels
plt.xlabel('Frequency / [arb. unit]', fontsize=20)
plt.ylabel('Power spectrum / [arb. unit]', fontsize=20)

# legend
plt.legend()

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

# display the plot

plt.savefig('powerspectrum.pdf')
plt.show()
