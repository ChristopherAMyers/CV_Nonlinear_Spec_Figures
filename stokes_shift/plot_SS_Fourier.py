import matplotlib.pyplot as plt
import numpy as np
from plot_JA_method import run_analysis
import global_settings as GS
from scipy.optimize import curve_fit
from os.path import join
from FourierTransform import explicitFourierTransform

def exp_func(t, a, b, c):
    return a * np.exp(-t/b) + c


AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254
AU_2_CM = 219474.63

fig, ax = plt.subplots(figsize=np.array((7, 5))*1.1)

data_files = (
    join('SS_results', 'JA_method_stripped.txt'),
    join('SS_results', 'JA_method_mm_4hb.txt'),
    join('SS_results', 'JA_method_mm_C4.txt'),
    join('SS_results', 'JA_method_mm_qm1.txt'),
    join('SS_results', 'JA_method_qm2.txt'),
)
dt = 4 #    in femtoseconds
titles = ['Stripped', 'H-bonded', 'Pi Solvent', 'QM1', 'QM2']
save_names = ['stripped', 'mm_4hb', 'mm_C4', 'mm_qm1', 'qm2']
colors = ['black', '#21ADEF', '#D321FF', 'red', 'blue']

#   first run stokes shift plot
# run_analysis()


for i in range(len(data_files)):
    if i < 2: continue
    data = np.loadtxt(data_files[i])
    times = data[:, 0]
    keep = (times <= 1000)
    data = data[keep]
    times = times[keep]/AU_2_FS
    # ss_original = data[:, 3] - data[:, 4]
    ss_original = data[:, 4]

    #   fit exponential to stokes shift
    popt, pcov = curve_fit(exp_func, times, ss_original, maxfev=10000, p0=(0.1, 10000, 2))
    print(popt)
    exp_fit = exp_func(times, *popt)
    # ss_residule = (ss_original - exp_fit)*np.exp(-times*0.005)
    ss_residule = (ss_original - exp_fit)

    # ax.plot(times*AU_2_FS, ss_original, label=titles[i], color=colors[i])
    # ax.plot(times*AU_2_FS, exp_fit, color='black', linestyle='--')
    # ax.plot(times*AU_2_FS, ss_residule, color=colors[i], linestyle='-')
    # ax.set_xlabel('Time (fs)')
    # ax.set_ylabel('Stokes Shift (eV)')

    # fft_res = fourierTransform(times, ss_residule, len(times), interp1d, {'bounds_error': False, 'fill_value': 0.0})

    fft_res = explicitFourierTransform(times, ss_residule, np.linspace(300/AU_2_CM, 800/AU_2_CM, 1000))
    omega = fft_res['freq']*AU_2_CM
    ft = np.abs(fft_res['ft'])
    ft /= ft.max()
    ax.plot(omega, ft)
    ax.set_ylim(0, 1.1)
    ax.set_xlim(omega.min(), omega.max())


    # break

ax.legend()
fig.tight_layout()
plt.show()