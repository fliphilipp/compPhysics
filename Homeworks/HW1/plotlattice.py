# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
#import matplotlib.pylab as plt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# input file
filename1 = 'pos.dat'

# import data
data = np.loadtxt(filename1)

spacing = data[1,0]
ndim = max(data[:,0]) / spacing

# initial size of plot window
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d') 

##########################################################

# for i in range(0,data.shape[0]):
# 	x = data[i,0]
# 	y = data[i,1]
# 	z = data[i,2]
# 	if (x % spacing < 0.1) and (y % spacing < 0.1) and (z % spacing < 0.1):
# 		ax.plot([x, x + spacing], [y, y], [z, z], 'k-', lw=0.7)
# 		ax.plot([x, x], [y, y + spacing], [z, z], 'k-', lw=0.7)
# 		ax.plot([x, x], [y, y], [z, z + spacing], 'k-', lw=0.7)
# 		ax.plot([x, x + spacing], [y, y + spacing], [z, z], 'k:', lw=0.3)
# 		ax.plot([x, x + spacing], [y, y], [z, z + spacing], 'k:', lw=0.3)
# 		ax.plot([x, x], [y, y + spacing], [z, z + spacing], 'k:', lw=0.3)
# 		ax.plot([x + spacing, x], [y, y + spacing], [z, z], 'k:', lw=0.3)
# 		ax.plot([x + spacing, x], [y, y], [z, z + spacing], 'k:', lw=0.3)
# 		ax.plot([x, x], [y + spacing, y], [z, z + spacing], 'k:', lw=0.3)


for i in range(0,data.shape[0]):
	x = data[i,0]
	y = data[i,1]
	z = data[i,2]
	iscorner = (x % spacing < 0.1) and (y % spacing < 0.1) and (z % spacing < 0.1)
	a = 0.3
	if x == 0 and iscorner: 
		ax.plot([x, x - spacing], [y, y], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y], [z, z + spacing], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x - spacing], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x - spacing], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x - spacing, x], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x - spacing, x], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y + spacing, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)

	if y == 0 and iscorner: 
		ax.plot([x, x + spacing], [y, y], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y - spacing], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y], [z, z + spacing], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x + spacing], [y, y - spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x + spacing], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y, y - spacing], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y - spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y - spacing, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)

	if z == 0 and iscorner: 
		ax.plot([x, x + spacing], [y, y], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y], [z, z - spacing], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x + spacing], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x + spacing], [y, y], [z, z - spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z - spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y], [z, z - spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y + spacing, y], [z, z - spacing], 'k:', lw=0.3, alpha=a)

	if iscorner:
		ax.plot([x, x + spacing], [y, y], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x], [y, y], [z, z + spacing], 'k-', lw=0.7, alpha=a)
		ax.plot([x, x + spacing], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x + spacing], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y, y + spacing], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y + spacing], [z, z], 'k:', lw=0.3, alpha=a)
		ax.plot([x + spacing, x], [y, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)
		ax.plot([x, x], [y + spacing, y], [z, z + spacing], 'k:', lw=0.3, alpha=a)

# plot atoms
ax.scatter(data[:,0], data[:,1], data[:,2], c='red', marker='.')

# labels
ax.set_xlabel(r'$x$ - position', fontsize=18)
ax.set_ylabel(r'$y$ - position', fontsize=18)
ax.set_zlabel(r'$z$ - position', fontsize=18)

# remove background gridlines
ax.set_axis_off()

# display the plot
plt.savefig('lattice-positions.png',bbox_inches='tight')
plt.show()
