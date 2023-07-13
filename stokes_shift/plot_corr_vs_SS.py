import matplotlib.pyplot as plt
import numpy as np
import global_settings as GS
from os.path import join
from scipy.interpolate import interp1d
from scipy.integrate import simpson
from correlation_function import compute_corr_func


AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254
AU_2_CM = 219474.63

#   use explicit separation of SE and GSB, as opposed to Jessica Anna's method using total Trans. Abs.
use_SE_GSB = True
#   use only the SE signal as the decay signal, as opposed to using the diff between SE and GSB
use_SE_only = True

if use_SE_GSB:
    method = 'SE_method'    # 'SE_method' or 'JA_method'
else:
    method = 'JA_method'


if __name__ == '__main__':

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=np.array((8, 9))*0.9, height_ratios=[1, 0.5])

    unshifted_diff = []
    shifted_diff = []
    # for i in range(len(plot_info)):
    for i, (job_name, job) in enumerate(GS.job_info.items()):
        if i not in [0, 4]:
            continue
        
        data_file = join('SS_results', f'{method}_{job_name}.txt')
        plt_args = job['plt_args']
    

        # data_file, plt_args = plot_info[i]
    
        data = np.loadtxt(data_file)
        times = data[:, 0]
        keep = (times <= 1500)
        data = data[keep]
        times = times[keep]

        GSB = data[:, 3]
        SE = data[:, 4]
        freq = GSB - SE
        if use_SE_only:
            ss = (SE - SE[-1]) / (SE[0] - SE[-1])
        else:
            ss = (freq - freq[-1])/(freq[0] - freq[-1])
        ax1.plot(times, ss, **plt_args)

        vee_file = join(job['job_dir'], 'vee_traj1.dat')
        vee_data = np.loadtxt(vee_file, usecols=0)
        corr_func, corr_times = compute_corr_func(vee_data, 2000)
        corr_func /= corr_func[0]
        ax1.plot(corr_times, corr_func, **plt_args, linestyle=(0, (1, 1)))

        interp_ss = interp1d(times, ss, bounds_error=False, fill_value=0.0)
        interp_cf = interp1d(corr_times, corr_func, bounds_error=False, fill_value=0.0)

        #   account for the phase difference by minimizing the std. dev. of the
        #   difference between the two signals
        sub_times = np.linspace(400, 800, 500)
        interpolated_cf = interp_cf(sub_times)
        shift_range = [0] + np.arange(-20, 20, 1).tolist()
        std_diff = []
        integral_diff = []
        for s in shift_range:
            diff = np.abs(interp_ss(sub_times + s) - interpolated_cf)
            std = np.std(diff)
            integral = simpson(diff, sub_times)/(sub_times[-1] - sub_times[0])
            std_diff.append(std)
            integral_diff.append(integral)
        idx_min_std = np.argmin(std_diff)
        unshifted_diff.append(integral_diff[0])
        shifted_diff.append(integral_diff[idx_min_std])
        print(f'Unshifted: {unshifted_diff[-1]:.4f} shifted: {shifted_diff[-1]:.4f}')


        w = 0.5
        ax2.bar(i, unshifted_diff[-1], width=w, **plt_args, alpha=0.5)
        ax2.bar(i, shifted_diff[-1], width=w, **plt_args, alpha=1.0)

        # combo_times = np.linspace(0, 1500, 1000)
        # diff = np.abs(interp_ss(combo_times + shift) - interp_cf(combo_times))
        # ax2.plot(combo_times, diff, **plt_args)

    ax1.set_xlim(0, 800)
    ax1.set_xlabel('Time (fs)')
    ax1.set_ylabel('SS(t), C(t)')

    ax2.set_xticks(range(len(GS.job_info)))
    ax2.set_xticklabels([x for x in GS.job_info])
    ax2.set_ylabel('Difference')

    
    # ax2.set_xlim(0, 800)
    # ax2.set_xlabel('Time (fs)')
    # ax2.set_ylabel('Difference')

    ax1.legend(ncol=2)
    fig.tight_layout()
    plt.show()