import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import global_settings as GS
from os.path import *
from scipy.interpolate import interp1d
from scipy.optimize import root

AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

#   global plot properties

data_folders = (
    GS.data_root_dir + 'stripped/results', 
    GS.data_root_dir + 'mm_4hb/results', 
    GS.data_root_dir + 'mm_C4/results', 
    GS.data_root_dir + 'mm_qm1/results', 
    GS.data_root_dir + 'qm2/results')
dt = 4 #    in femtoseconds
titles = ['Stripped', 'H-bonded', 'Pi Solvent', 'QM1', 'QM2']
colors = ['black', '#21ADEF', '#D321FF', 'red', 'blue']

fig, ax = plt.subplots(figsize=np.array((5, 7))*1.1)

min_pos = None
ideal_max_pos = 2.10
ref_bottom_corner = None
colormap_name = 'jet'
interpolate = True

for i in range(len(data_folders)):
    #   import data, first try to load an uncompressed compressed version
    data_file = join(data_folders[i], 'vee__2nd_order_cumulant_transient_absorption_spec.txt')
    print(data_file)
    if isfile(data_file):
        data = np.loadtxt(data_file)
    else:
        data = np.loadtxt(data_file + '.gz')
    data[:, 1] *= AU_2_EV
    
    #   reshape data
    n_pts_per_time = len(set(data[:, 1]))
    times = np.array(sorted(set(data[:, 0])))*AU_2_FS
    n_times = len(times)
    data = data[:, 1:3].reshape((n_times, n_pts_per_time, 2))

    emissions = []
    absorptions = []
    emissions_interps = []
    absorptions_interps = []
    shift = None
    for j, spectra in enumerate(data):
        spectra[:, 1] /= np.max(spectra[:, 1])
        # spectra[:, 0] += shift_freq
        #   first split the spectrum by from the location of the maximum
        idx_of_maximum = np.argmax(spectra[:, 1])
        spectra_1 = spectra[0:idx_of_maximum]
        spectra_2 = spectra[idx_of_maximum:]

        #   seach for the first half (emission)
        freq_1, intensity_1 = spectra_1.T
        intensity_1 -= 0.5
        tmp_intensity = intensity_1.copy()
        tmp_intensity[tmp_intensity < 0] = 1E10  #  negative values will not pass the minimum search (next line)
        idx_pos = np.argmin(tmp_intensity)  #   closest to zero, but positive
        if intensity_1[idx_pos - 1] < 0:
            idx_neg = idx_pos - 1
        else:
            idx_neg = idx_pos + 1
        emiss_f1, emiss_f2 = freq_1[idx_pos], freq_1[idx_neg]
        emiss_i1, emiss_i2 = intensity_1[idx_pos], intensity_1[idx_neg]
        emiss_avg = 0.5*(emiss_f1 + emiss_f2)

        #   now search the seond half (absorption)
        freq_2, intensity_2 = spectra_2.T
        intensity_2 -= 0.5
        tmp_intensity = intensity_2.copy()
        tmp_intensity[tmp_intensity < 0] = 1E10  #  negative values will not pass the minimum search (next line)
        idx_pos = np.argmin(tmp_intensity)  #   closest to zero, but positive
        if intensity_2[idx_pos - 1] < 0:
            idx_neg = idx_pos - 1
        else:
            idx_neg = idx_pos + 1
        abs_f1, abs_f2 = freq_2[idx_pos], freq_2[idx_neg]
        abs_i1, abs_i2 = intensity_2[idx_pos], intensity_2[idx_neg]
        abs_avg = 0.5*(abs_f1 + abs_f2)
        
        #   finalize by solving for the line that connects the two points
        slope_e = (emiss_i1 - emiss_i2)/(emiss_f1 - emiss_f2)
        solution_e = emiss_f1 - emiss_i1/slope_e
        slope_a = (abs_i1 - abs_i2)/(abs_f1 - abs_f2)
        solution_a = abs_f1 - abs_i1/slope_a
        if shift is None:
            shift = (ideal_max_pos - solution_a)
            
        emissions_interps.append(solution_e + shift)
        absorptions_interps.append(solution_a + shift)

        emissions.append(emiss_avg + shift)
        absorptions.append(abs_avg + shift)

    if interpolate:
        ax.plot(times, absorptions_interps, color=colors[i], linestyle='--')
        ax.plot(times, emissions_interps, color=colors[i], label=titles[i])
    else:
        ax.plot(times, absorptions, color=colors[i], linestyle='--')
        ax.plot(times, emissions, color=colors[i], label=titles[i])

    


# fig.colorbar(contour)
ax.legend(ncol=2)
ax.set_xlabel('Time (fs)')
ax.set_ylabel("$E_\mathrm{GSB}$ / $E_\mathrm{SE} $ (eV)")
ax.set_ylim(1.8, 2.13)
ax.set_xlim(0, 500)
fig.tight_layout()

fig.savefig('SS_vs_Time.png', dpi=300)
plt.show()
