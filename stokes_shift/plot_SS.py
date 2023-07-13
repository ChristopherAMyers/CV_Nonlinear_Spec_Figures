import matplotlib.pyplot as plt
import numpy as np
from os.path import join, isfile
from FourierTransform import explicitFourierTransform
from exp_fit import fit_decay, fit_sinusoidal_decay
from plot_JA_method import run_analysis as run_JA
from plot_TA_stokes_shift import run_analysis as run_TA
import global_settings as GS


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

DEBUG = False

def plot_on_ax(ax, times, GSB, SE, 
               ax0 : plt.Axes = None,
               ax1 : plt.Axes = None, 
               ax2 : plt.Axes = None, 
               ax3 : plt.Axes = None, 
               plt_args={}):
    freq = GSB - SE

    if use_SE_only:
        ss = (SE - SE[-1]) / (SE[0] - SE[-1])
    else:
        ss = (freq - freq[-1])/(freq[0] - freq[-1])
    

    if ax0 is not None:
        ax0.plot(times, GSB, **plt_args, linestyle='--')
        ax0.plot(times, SE, **plt_args)
        ax0.set_xlabel('Time (fs)')
        ax0.set_ylabel('SE / GSB (eV)')


    if DEBUG:
        t_exp = np.linspace(10, 2000, 50)
        y = []
        for x in t_exp:
            y.append(fit_decay(times, ss, x)[1])
            print(y[-1])
        fig2, ax2 = plt.subplots()
        ax2.plot(t_exp, y, marker='.')
        plt.show()
        exit()

    #   fit exponential to stokes shift
    # popt, rrmsd, exp_fit = fit_decay(times, ss)
    popt, rrmsd, exp_fit = fit_sinusoidal_decay(times, ss)
    print(f'RRMSD for {plt_args["label"]}: {rrmsd:.5f}')
    for p in popt:
        print(f'    {p:.5f}')
    ss_residule = (ss - exp_fit)

    if ax1 is not None:
        ax1.plot(times, ss, **plt_args)
        ax1.plot(times, exp_fit, color='black')
        # ax1.plot(times, exp_fit, color='black', linestyle='--')
        # ax1.plot(times, ss_residule, color=colors[i], linestyle='-')
        ax1.set_xlabel('Time (fs)')
        ax1.set_ylabel('Stokes Shift (eV)')

    if ax2 is not None:
        ax2.plot(times, ss_residule, **plt_args)
        ax2.set_xlabel('Time (fs)')
        ax2.set_ylabel('Residual')

    if ax3 is not None:
        fft_res = explicitFourierTransform(times/AU_2_FS, ss_residule, np.linspace(300/AU_2_CM, 800/AU_2_CM, 1000))
        omega = fft_res['freq']*AU_2_CM
        ft = np.abs(fft_res['ft'])
        ft /= ft.max()
        ax3.plot(omega*0.96, ft, **plt_args)
        ax3.set_ylim(0, 1.1)
        ax3.set_xlim(omega.min(), omega.max())
        ax3.set_ylabel('Intensity')
        ax3.set_xlabel('Frequency (cm$^{-1}$)')

if __name__ == '__main__':

    fig, ax = plt.subplots(2, 2, figsize=np.array([16, 10])*0.8)
    ax = ax.flatten()
    
    for i, (job_name, job) in enumerate(GS.job_info.items()):
        if i < 0: continue
        if i == 4:
            DEBUG = False

        data_file = join('SS_results', f'{method}_{job_name}.txt')
        plt_args = job['plt_args']

        if not isfile(data_file):
            if use_SE_GSB:
                run_TA()
            else:
                run_JA()

        data = np.loadtxt(data_file)
        times = data[:, 0]
        keep = (times <= 1500)
        data = data[keep]
        times = times[keep]

        plot_on_ax(ax, times, data[:, 3], data[:, 4], *ax, plt_args=plt_args)

    exp_ft = np.loadtxt('exp_FT.csv', delimiter=',')
    # ax[3].plot(exp_ft[:, 0], exp_ft[:, 1], color='gray', zorder=-1)
    ax[3].fill_between(exp_ft[:, 0], exp_ft[:, 1], color='lightgray', zorder=-2)

    ax[-1].legend()
    fig.tight_layout()
    plt.show()