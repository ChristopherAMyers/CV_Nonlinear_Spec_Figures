import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np
import global_settings as GS
from os.path import *

AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

#   global plot properties
fontsize = plt.rcParams['axes.labelsize']

data_files = (
    GS.data_root_dir + 'gs-aimd/stripped/results/vee__2nd_order_cumulant_transient_absorption_spec.txt', 
    GS.data_root_dir + 'gs-aimd/mm_4hb/results/vee__2nd_order_cumulant_transient_absorption_spec.txt', 
    GS.data_root_dir + 'gs-aimd/mm_C4/results/vee__2nd_order_cumulant_transient_absorption_spec.txt', 
    GS.data_root_dir + 'gs-aimd/mm_qm1/results/vee__2nd_order_cumulant_transient_absorption_spec.txt',
    GS.data_root_dir + 'gs-aimd/qm2/results/vee__2nd_order_cumulant_transient_absorption_spec.txt',
    GS.data_root_dir + 'experimental_data/transient_abs.txt'
    )
dt = 4 #    in femtoseconds
titles = ['Stripped', 'H-bonded', 'Closest 4', 'QM1', 'QM2', 'EXP']

length = 1.8*len(data_files) + 0.5
fig, ax_grid = plt.subplots(1, len(data_files), figsize=np.array((length,2.9))*1.2)
ax_grid = np.array([ax_grid])



min_pos = None
ideal_max_pos = 2.1
ref_bottom_corner = None
colormap_name = 'jet'

idx = -1
for i, ax_row in enumerate(ax_grid):
    for j, ax in enumerate(ax_row):
        idx += 1
        #   import data, first try to load an uncompressed compressed version
        data_file = data_files[j]
        if isfile(data_file):
            data = np.loadtxt(data_file)
        else:
            data = np.loadtxt(data_file + '.gz')
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

        X = data[:,0].reshape(dim_x,dim_y)
        Y = data[:,1].reshape(dim_x,dim_y)
        Z = data[:,2].reshape(dim_x,dim_y)
        Z /= np.max(Z)
        Z[np.where(Z < 0)] = 0.0
        # exit()

        contour = ax.contourf(X,Y,Z, 100,cmap=colormap_name)

        if i+1 != len(ax_grid):
            ax.set_xticks([])
        else:
            pass
            # ax.set_xlabel('Time (fs)')
            # ax.set_xticks([0, 250, 500])
        if j != 0:
            ax.set_yticks([])
        else:
            ax.set_yticks([1.9, 2.0, 2.1, 2.2, 2.3])
            # ax.set_ylabel('$\omega_3$ (eV)')
        
        #   MM environment label
        ax.text(0.05, 0.95, titles[idx],
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes,
            color='white', fontsize=16),

        ax.set_xlim(0, 500)
        ax.set_ylim(1.88, 2.32)

        #   custom tick labels are needed to prevent overlaps between subplots
        ax.set_xticks([0.0, 250, 500])
        ax.set_xticklabels([' ', ' ', ' '])
        ax.text(0.0, 1.84, '0', ha='left', fontsize=fontsize)
        ax.text(250, 1.84, '250', ha='center', fontsize=fontsize)
        ax.text(500, 1.84, '500', ha='right', fontsize=fontsize)
        # txt = ax.text(0.05, 0.15, titles[i],
        #             verticalalignment='top', horizontalalignment='left',
        #             transform=ax.transAxes,
        #             color='white', fontsize=fontsize+4)
        # txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])

for ax in ax_grid.flatten():
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

# fig.colorbar(contour, use_gridspec=False)
fig.supylabel('$\omega_3$ (eV)')
fig.supxlabel('Time (ps)')
fig.tight_layout()
# plt.subplots_adjust(wspace=0.1,
                    # hspace=0.5)

fig.savefig('png/transient_abs.png', dpi=300)
plt.show()
