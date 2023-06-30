import matplotlib.pyplot as plt
import numpy as np
from os.path import *
import global_settings as GS
'''
    This script creates 2DES plot grid.
    Each column is a fixed t2 fime and each row
    is a fixed computation (H-Bonds, QM1, QM2, etc.)
'''

AU_2_EV = 27.211396132

#   global plot properties
rcParams = {
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'font.weight': 'bold',
    'axes.labelsize': 16,
    'axes.labelweight': 'extra bold'
}
plt.rcParams.update(rcParams)

#   each key is the text that will be placed on each plot.
#   within each key are the lists of data files to use.
plot_info = {
    'Peripheral': (
    GS.data_root_dir + 'gs-aimd/mm_4hb/results/vee__2DES_12.dat',
    GS.data_root_dir + 'gs-aimd/mm_4hb/results/vee__2DES_50.dat',
    GS.data_root_dir + 'gs-aimd/mm_4hb/results/vee__2DES_100.dat'),

    'Axial': (
    GS.data_root_dir + 'gs-aimd/mm_C4/results/vee__2DES_12.dat',
    GS.data_root_dir + 'gs-aimd/mm_C4/results/vee__2DES_50.dat',
    GS.data_root_dir + 'gs-aimd/mm_C4/results/vee__2DES_100.dat'),

    'Full MM': (
    GS.data_root_dir + 'gs-aimd/mm_qm1/results/vee__2DES_12.dat',
    # GS.data_root_dir + 'gs-aimd/mm_qm1/results/vee__2DES_1.dat',
    GS.data_root_dir + 'gs-aimd/mm_qm1/results/vee__2DES_50.dat',
    GS.data_root_dir + 'gs-aimd/mm_qm1/results/vee__2DES_100.dat'),

    'Experiment': (
    GS.data_root_dir + 'experimental_data/twoDSpec_Converted/cropped_data/twoDSpecMat_7.dat',
    GS.data_root_dir + 'experimental_data/twoDSpec_Converted/cropped_data/twoDSpecMat_37.dat',
    GS.data_root_dir + 'experimental_data/twoDSpec_Converted/cropped_data/twoDSpecMat_67.dat',),
}
t2_times = [0, 200, 400] # t2 times in femtoseconds

fig_height = 1.8*len(plot_info) + 0.5
n_rows = len(plot_info)
n_columns = np.max([len(x) for x in plot_info.values()])
fig, ax_grid = plt.subplots(
    n_rows, 
    n_columns, 
    figsize=np.array((6.2,fig_height))*1.2,
    )

# ideal_max_pos = (2.088511, 2.088511)
ideal_max_pos = (2.05149689, 2.03809733)
ref_bottom_corner = None
colormap_name = 'jet'

x_lims = []
y_lims = []

for i, (title, files) in enumerate(plot_info.items()):
    min_pos = None
    print("Plotting: ", title)
    for j, data_file in enumerate(files):
        ax = ax_grid[i, j]

        #   import data, first try to load an uncompressed compressed version
        if isfile(data_file):
            data = np.loadtxt(data_file)
        else:
            data = np.loadtxt(data_file + '.gz')
        dim = int(np.sqrt(data.shape[0]))
        data[:, 0:2] *= AU_2_EV
        # exit()

        unique_x = set(data[:, 0])
        unique_y = set(data[:, 1])
        dim_x = len(unique_x)
        dim_y = len(unique_y)

        if min_pos is None:
            min_loc = np.argmax(data[:, 2])
            min_pos = data[min_loc, 0:2].copy()
            shift = ideal_max_pos- min_pos
        if i != len(plot_info) - 1:
            data[:, 0:2] += shift

        print(dim_x, dim_y)
        X = data[:,0].reshape(dim_x,dim_y)
        Y = data[:,1].reshape(dim_x,dim_y)
        Z = data[:,2].reshape(dim_x,dim_y)
        # exit()
        Z /= np.max(Z)
        Z[np.where(Z < 0)] = 0.0

        contour = ax.contourf(X,Y,Z, dim,cmap=colormap_name)

        #   x-ticks and labels are only on the bottom plots
        if i == len(plot_info) - 1:
            ax.set_xticks([1.9, 2.0, 2.1, 2.2, 2.3])
            ax.set_xlabel('$\omega_1$ (eV)')
        else:
            ax.set_xticks([])

        #   y-ticks and labels are only on the left-most plots
        if j == 0:
            ax.set_yticks([1.9, 2.0, 2.1, 2.2, 2.3])
            ax.set_ylabel('$\omega_3$ (eV)')
        else:
            ax.set_yticks([])

        #   T2 time label
        ax.text(0.95, 0.01, '{:d} fs'.format(t2_times[j]),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes,
            color='white', fontsize=16)
        
        #   MM environment label
        if j == 0:
            ax.text(0.05, 0.95, title,
                verticalalignment='top', horizontalalignment='left',
                transform=ax.transAxes,
                color='white', fontsize=16)

        x_lims.append(ax.get_xlim())
        y_lims.append(ax.get_ylim())
        ax.set_xlim(1.85, 2.32)
        ax.set_ylim(1.85, 2.32)

        ax.plot(ax.get_ylim(), ax.get_ylim(), color='white', linestyle='--', linewidth=1.5)

min_min_freq = min(np.min(x_lims), np.min(y_lims))
min_max_freq = min(np.max(x_lims), np.max(y_lims))

for ax in ax_grid.flatten():
    pass
    # ax.set_xlim(xmin=min_min_freq, xmax=min_max_freq)
    # ax.set_xbound(lower=min_min_freq, upper=min_max_freq)
    # ax.set_ylim(ymin=min_min_freq, ymax=min_max_freq)
    # ax.axis('equal')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')


# fig.colorbar(contour)
fig.tight_layout()
plt.subplots_adjust(wspace=0.05,
                    hspace=0.05)

fig.savefig('png/2DES_solvent_poster.png', dpi=400)
plt.show()
