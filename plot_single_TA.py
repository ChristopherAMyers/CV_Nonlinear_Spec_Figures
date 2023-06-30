import matplotlib.pyplot as plt
import numpy as np
from os.path import *
plt.style.use('style.mplstyle')

AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

fig, ax = plt.subplots(figsize=np.array((3.0,2.9))*1.5)

min_pos = None
ideal_max_pos = 2.1
ref_bottom_corner = None
colormap_name = 'jet'

#   import data
data_file = '/Users/cmyers/Library/CloudStorage/OneDrive-UniversityofCaliforniaMerced/crysel_violet/2DES/shao-yu_replication/gs-aimd/stripped/results/vee__2nd_order_cumulant_transient_absorption_spec.txt'

data = np.loadtxt(data_file)
data[:, 1] *= AU_2_EV
data[:, 0] *= AU_2_FS

#   shift the data to only include time >= 0
#   this is mostely for the experimental data
zero_idx = 0
for n, time in enumerate(data[:, 0]):
    if time >= 0:
        zero_idx = n
        break
data = data[zero_idx:]
dim_x = len(set(data[:, 0]))
dim_y = len(set(data[:, 1]))

#   recenter data at experimental maximum
min_loc = np.argmax(data[0:dim_y, 2])
min_pos = data[min_loc, 1]
print(min_pos, zero_idx)
data[:, 1] += (ideal_max_pos- min_pos)
data[:, 0]

#   create data grids
X = data[:,0].reshape(dim_x,dim_y)
Y = data[:,1].reshape(dim_x,dim_y)
Z = data[:,2].reshape(dim_x,dim_y)
Z /= np.max(Z)
Z[np.where(Z < 0)] = 0.0

contour = ax.contourf(X,Y,Z, 100,cmap=colormap_name)

#   labels and colorbar
fig.colorbar(contour, ax=ax, ticks=np.arange(0, 1.01, 0.2))
ax.set_ylabel('$\omega_3$ (eV)')
ax.set_xlabel('Time (ps)')
fig.tight_layout()

# fig.savefig('png/transient_abs.png', dpi=300)
plt.show()
