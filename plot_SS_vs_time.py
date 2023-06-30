import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import global_settings as GS
from os.path import *
from os import makedirs
from scipy.interpolate import interp1d
from scipy.optimize import root

AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

def run_analysis(ax=None):
    data_folders = (
        # GS.data_root_dir + 'gs-aimd/stripped/results', 
        # GS.data_root_dir + 'gs-aimd/mm_4hb/results', 
        # GS.data_root_dir + 'gs-aimd/mm_C4/results', 
        # GS.data_root_dir + 'gs-aimd/mm_qm1/results', 
        # GS.data_root_dir + 'gs-aimd/qm2/results',
    
        GS.data_root_dir + 'gs-aimd-longer/stripped/', 
        GS.data_root_dir + 'gs-aimd-longer/mm_4hb/', 
        GS.data_root_dir + 'gs-aimd-longer/mm_C4/', 
        GS.data_root_dir + 'gs-aimd-longer/mm_qm1/', 
        GS.data_root_dir + 'gs-aimd-longer/qm2/',
    )
    dt = 4 #    in femtoseconds
    titles = ['Stripped', 'H-bonded', 'Pi Solvent', 'QM1', 'QM2']
    save_names = ['stripped', 'mm_4hb', 'mm_C4', 'mm_qm1', 'qm2']
    colors = ['black', '#21ADEF', '#D321FF', 'red', 'blue']


    min_pos = None
    ideal_max_pos = 2.10
    ref_bottom_corner = None
    colormap_name = 'jet'
    interpolate = True

    #   results are stored in directory
    makedirs('stokes_shift', exist_ok=True)

    for i in range(len(data_folders)):
        #   import data, first try to load an uncompressed compressed version
        data_file = join(data_folders[i], 'vee__2nd_order_cumulant_transient_absorption_spec.txt')
        print(data_file)
        if isfile(data_file):
            data = np.loadtxt(data_file)
        else:
            data = np.loadtxt(data_file + '.gz')
        data[:, 1] *= AU_2_EV
        print(data.shape)
        
        #   reshape data
        n_pts_per_time = len(set(data[:, 1]))
        times = np.array(sorted(set(data[:, 0])))*AU_2_FS
        n_times = len(times)
        data = data[:, 1:3].reshape((n_times, n_pts_per_time, 2))

        emissions = []
        absorptions = []
        emissions_interps = []
        absorptions_interps = []
        keep_times = []
        shift = None
        for j, spectra in enumerate(data):
            
            if np.sum(np.isnan(spectra)) > 0:
                print("Skipping ", j)
                continue

            spectra[:, 1] /= np.max(spectra[:, 1])
            # spectra[:, 0] += shift_freq

            #   First split the spectrum by from the location of the maximum.
            #   We also need to make sure that the maximum is not at the
            #   end points. This comes up when the transient-abs is periodic.
            #   Also, we keep trimming until the last and first 10% of points
            #   are also no greater than 1/2.

            dim = spectra.shape[0]
            idx_of_maximum = np.argmax(spectra[:, 1])
            max_first_10 = np.max(spectra[0:int(dim*0.1)])
            max_last_10 = np.max(spectra[int(dim*0.9):])
            while idx_of_maximum/dim < 0.1 or idx_of_maximum/dim > 0.9 or max_last_10 > 0.5 or max_first_10 > 0.5:
                spectra = spectra[1:]
                spectra = spectra[0:-1]
                dim = spectra.shape[0]
                idx_of_maximum = np.argmax(spectra[:, 1])
                max_first_10 = np.max(spectra[0:int(dim*0.1), 1])
                max_last_10 = np.max(spectra[int(dim*0.9):, 1])
                
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

            try:
                abs_f1, abs_f2 = freq_2[idx_pos], freq_2[idx_neg]
            except:
                fig2, ax2 = plt.subplots()
                ax2.plot(spectra[:, 1])
                plt.show()
                exit()
            abs_i1, abs_i2 = intensity_2[idx_pos], intensity_2[idx_neg]
            abs_avg = 0.5*(abs_f1 + abs_f2)
            
            #   finalize by solving for the line that connects the two points
            slope_e = (emiss_i1 - emiss_i2)/(emiss_f1 - emiss_f2)
            solution_e = emiss_f1 - emiss_i1/slope_e
            slope_a = (abs_i1 - abs_i2)/(abs_f1 - abs_f2)
            solution_a = abs_f1 - abs_i1/slope_a
            if shift is None:
                shift = (ideal_max_pos - solution_a)

            if abs_avg + shift < 1.5 and i == 3:
                print("PLOT: ", j, idx_neg, idx_pos)
                fig2, ax2 = plt.subplots()
                ax2.plot(*spectra.T, marker='.')
                ax2.hlines(0, 2.0, 4.0)
                plt.show()
                exit()
                
            emissions_interps.append(solution_e + shift)
            absorptions_interps.append(solution_a + shift)

            emissions.append(emiss_avg + shift)
            absorptions.append(abs_avg + shift)

            keep_times.append(times[j])


        if ax is not None:
            if interpolate:
                ax.plot(keep_times, absorptions_interps, color=colors[i], linestyle='--')
                ax.plot(keep_times, emissions_interps, color=colors[i], label=titles[i])
            else:
                ax.plot(keep_times, absorptions, color=colors[i], linestyle='--')
                ax.plot(keep_times, emissions, color=colors[i], label=titles[i])

        out_data = np.array([keep_times, absorptions, emissions, absorptions_interps, emissions_interps]).T
        np.savetxt('stokes_shift/{:s}.txt'.format(save_names[i]), out_data, 
                header='time (fs) absorption(nearest) emission(nearest) absorption(interp) emission(interp)')


if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=np.array((5, 7))*1.1)
    run_analysis(ax)

    ax.legend(ncol=2)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel("$E_\mathrm{GSB}$ / $E_\mathrm{SE} $ (eV)")
    ax.set_ylim(1.78, 2.13)
    ax.set_xlim(-5, 1300)

    fig.tight_layout()
    fig.savefig('png/SS_vs_Time.png', dpi=300)
    plt.show()
