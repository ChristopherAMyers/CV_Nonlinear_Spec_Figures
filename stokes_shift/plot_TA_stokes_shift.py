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

def crop_data(spectra):
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
    return spectra

def run_analysis(ax=None):

    ideal_max_pos = 2.10
    interpolate = True

    #   results are stored in directory
    makedirs('SS_results', exist_ok=True)

    # for i in range(len(data_folders)):
    for i, (job_name, job) in enumerate(GS.job_info.items()):
        if i < 0:
            continue
        #   import data, first try to load an uncompressed compressed version
        data_file = join(job['job_dir'], 'vee__2nd_order_cumulant_transient_absorption_spec.txt')
        print(f'\nAnalyzing TA for file:\n{data_file}')
        if isfile(data_file):
            data = np.loadtxt(data_file)
        else:
            data = np.loadtxt(data_file + '.gz')
        data[:, 1] *= AU_2_EV

        data_file_SE = join(job['job_dir_SE'], 'vee__2nd_order_cumulant_transient_absorption_spec.txt')
        if isfile(data_file_SE):
            data_SE = np.loadtxt(data_file_SE)
        else:
            data_SE = np.loadtxt(data_file_SE + '.gz')
        data_SE[:, 1] *= AU_2_EV
        data_GB = np.copy(data)
        #   ground state bleach is the full spectra minus the SE spectra
        data_GB[:, 2] -= data_SE[:, 2]

        
        #   reshape data
        n_pts_per_time = len(set(data[:, 1]))
        times = np.array(sorted(set(data[:, 0])))*AU_2_FS
        n_times = len(times)
        data = data[:, 1:3].reshape((n_times, n_pts_per_time, 2))
        data_SE = data_SE[:, 1:3].reshape((n_times, n_pts_per_time, 2))
        data_GB = data_GB[:, 1:3].reshape((n_times, n_pts_per_time, 2))


        #   these will be our collected results
        SE_maximums = []
        GB_maximums = []
        SE_maximums_interp = []
        GB_maximums_interp = []
        keep_times = []

        #   loop over each row in the TA
        shift = None
        for j in range(len(data_SE)):
            spectra_SE = data_SE[j].T
            spectra_GB = data_GB[j].T
            
            if np.sum(np.isnan(spectra_SE)) > 0 or np.sum(np.isnan(spectra_GB)) > 0:
                print(f"\tSkipping TA time-step {j}, contains NaN")
                continue

            spectra_SE[1] /= np.max(spectra_SE[1])
            spectra_GB[1] /= np.max(spectra_GB[1])

            spectra_SE = crop_data(spectra_SE.T).T
            spectra_GB = crop_data(spectra_GB.T).T

            max_freqs = []
            max_freqs_interp = []
            for spectra in (spectra_SE, spectra_GB):
                #   find the neareast frequency as the maximum
                max_idx = np.argmax(spectra[1])
                max_freq = spectra[0][max_idx]
                max_freqs.append(max_freq)

                #   use a cubic spline to get more accurate maximum peak
                start, end = max_idx - 10, max_idx + 10
                sub_freq = spectra[0][start:end]
                interp_f = interp1d(sub_freq, spectra[1][start:end], kind='cubic')
                interp_x = np.linspace(sub_freq[0], sub_freq[-1], 1000)
                interp_y = interp_f(interp_x)
                max_freqs_interp.append(interp_x[interp_y.argmax()])

                # fig2, ax2 = plt.subplots()
                # ax2.scatter(sub_freq, spectra[1][start:end])
                # ax2.plot(interp_x, interp_y)
                # ax2.plot(spectra[0], spectra[1])
                # plt.show()
                # exit()

            if shift is None:
                shift = (ideal_max_pos - max_freqs[0])

            #   record results
            SE_maximums.append(max_freqs[0] + shift)
            GB_maximums.append(max_freqs[1] + shift)
            SE_maximums_interp.append(max_freqs_interp[0] + shift)
            GB_maximums_interp.append(max_freqs_interp[1] + shift)
            keep_times.append(times[j])

        GB_maximums_interp = np.array(GB_maximums_interp)
        SE_maximums_interp = np.array(SE_maximums_interp)

        delta_omega = SE_maximums_interp[0] - SE_maximums_interp
        SS = (delta_omega[-1] - delta_omega)/delta_omega[-1]

        if ax is not None:
            # ax.plot(keep_times, GB_maximums_interp, color=colors[i], linestyle='--')
            # ax.plot(keep_times, SE_maximums_interp, color=colors[i], label=titles[i])

            ax.plot(keep_times, SS, **job['plt_args'])

        out_data = np.array([keep_times, GB_maximums, SE_maximums, GB_maximums_interp, SE_maximums_interp]).T
        np.savetxt(join('SS_results', f'SE_method_{job_name:s}.txt'), out_data, 
                header='time (fs) absorption(nearest) emission(nearest) absorption(interp) emission(interp)')



if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=np.array((5, 7))*1.1)
    run_analysis(ax)

    ax.legend(ncol=2)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel("$E_\mathrm{GSB}$ / $E_\mathrm{SE} $ (eV)")
    # ax.set_ylim(1.78, 2.13)
    ax.set_xlim(-5, 1300)

    fig.tight_layout()
    fig.savefig('png/SS_vs_Time.png', dpi=300)
    plt.show()
