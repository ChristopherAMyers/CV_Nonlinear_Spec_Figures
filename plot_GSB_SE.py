import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import global_settings as GS
from os.path import *

AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

data_folders = (
    GS.data_root_dir + 'gs-aimd/stripped_SE/results', 
    GS.data_root_dir + 'gs-aimd/mm_4hb_SE/results', 
    GS.data_root_dir + 'gs-aimd/mm_C4_SE/results', 
    GS.data_root_dir + 'gs-aimd/mm_qm1_SE/results', 
    GS.data_root_dir + 'gs-aimd/qm2_SE/results')

dt = 4 #    in femtoseconds
titles = ['Stripped', 'H-bonded', 'Closest 4', 'QM1', 'QM2']
colors = ['black', '#21ADEF', '#D321FF', 'red', 'blue']
time_number = [97, 97, 100, 100, 100]

fig, ax = plt.subplots(figsize=np.array((7.5, 4))*1.0)

ideal_max_pos = 2.1

data_regular = None
data_SE = None
data_GB = None
for idx in range(len(data_folders)):
    #   import data, first try to load an uncompressed compressed version
    data_file = join(data_folders[idx], 'vee__2nd_order_cumulant_transient_absorption_spec.txt')
    print(data_file)
    if isfile(data_file):
        data = np.loadtxt(data_file)
    else:
        data = np.loadtxt(data_file + '.gz')

    dim_x = len(set(data[:, 0]))
    dim_y = len(set(data[:, 1]))
    n_times = int(len(data)/dim_y)
    data[:, 1] *= AU_2_EV
    print(dim_x, dim_y)

    data_GB = data[0:dim_y, 1:]
    data_GB[:, 1] /= np.max(data_GB[:, 1])

    start = (time_number[idx]-1)*dim_y
    end = time_number[idx]*dim_y
    data_SE = data[start:end, 1:]
    print(data.shape, data_SE.shape)
    data_SE[:, 1] /= np.max(data_SE[:, 1])

    max_freq = data_GB[np.argmax(data_GB[:, 1]), 0]
    shift = (ideal_max_pos - max_freq)
    data_GB[:, 0] += shift
    data_SE[:, 0] += shift

    ax.plot(data_GB[:, 0], data_GB[:, 1], color=colors[idx], label=titles[idx])
    ax.plot(data_SE[:, 0], data_SE[:, 1], color=colors[idx], linestyle='--')

# fig.colorbar(contour)
ax.legend(fontsize=plt.rcParams['axes.labelsize']*0.9)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel("Intensity (Arb. Units)")
ax.set_ylim(0, 1.1)
ax.set_xlim(1.6, 2.5)
ax.set_xticks(np.arange(1.6, 2.5, 0.2))
ax.set_yticks([])
fig.tight_layout()

fig.savefig('png/aimd_GB_SE.png', dpi=300)
plt.show()
