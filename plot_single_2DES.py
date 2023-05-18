import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
from os.path import *

AU_2_EV = 27.211396132

#   global plot properties
rcParams = {
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'font.weight': 'bold',
    'axes.labelsize': 14,
    'axes.labelweight': 'extra bold'
}
plt.rcParams.update(rcParams)


data_file = sys.argv[1]

fig, ax = plt.subplots(figsize=np.array((4,4))*1.2)

min_pos = None
ideal_max_pos = (2.1, 2.1)
ref_bottom_corner = None
colormap_name = 'jet'

#   import data, first try to load an uncompressed compressed version
data = np.loadtxt(data_file)
dim = int(np.sqrt(data.shape[0]))
data[:, 0:2] *= AU_2_EV

print(np.min(data[:, 2]))
print(np.max(data[:, 2]))

#   recenter data at experimental maximum
#   further data are then centered by it's bottom corner
if min_pos is None:
    min_loc = np.argmax(data[:, 2])
    min_pos = data[min_loc, 0:2].copy()
    data[:, 0:2] += (ideal_max_pos- min_pos)
    ref_bottom_corner = np.min(data, axis=0)[0:2]
bottom_corner = np.min(data, axis=0)[0:2]
data[:, 0:2] += (ref_bottom_corner - bottom_corner)

X = data[:,0].reshape(dim,dim)
Y = data[:,1].reshape(dim,dim)
Z = data[:,2].reshape(dim,dim)
Z /= np.max(Z)
Z[np.where(Z < 0)] = 0.0

contour = ax.contourf(X,Y,Z, dim,cmap=colormap_name)

ax.set_xticks([2.0, 2.1, 2.2, 2.3])
ax.set_yticks([2.0, 2.1, 2.2, 2.3])

# ax.axis('equal')
for spine in ax.spines.values():
    spine.set_edgecolor('white')

# fig.colorbar(contour)
fig.tight_layout()
plt.subplots_adjust(wspace=0.05,
                    hspace=0.05)

fig.savefig('plot.png', dpi=300)
plt.show()
