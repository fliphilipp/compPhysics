# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
#import matplotlib.pylab as plt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# input file
filename1 = 'largeBox.dat'

# import data
data = np.loadtxt(filename1)

# initial size of plot window
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d') 

# plot atoms
ax.scatter(data[:,0], data[:,1], data[:,2], c='red', marker='.')

# labels
ax.set_xlabel(r'$x$ - position', fontsize=18)
ax.set_ylabel(r'$y$ - position', fontsize=18)
ax.set_zlabel(r'$z$ - position', fontsize=18)

# remove background gridlines
ax.set_axis_off()

# display the plot
plt.savefig('largeBox.png',bbox_inches='tight')
plt.show()
