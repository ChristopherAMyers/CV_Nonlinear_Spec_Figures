import matplotlib.pyplot as plt
import numpy as np
from plot_SS_vs_time import run_analysis
import global_settings as GS
from scipy.fft import fft, fftfreq

def fit_to_exponential():
    pass


AU_2_EV = 27.211396132
AU_2_FS = 0.02418884254

fig, ax = plt.subplots(figsize=np.array((7, 5))*1.1)

data_files = (
    'stokes_shift/stripped.txt', 
    'stokes_shift/mm_4hb.txt', 
    'stokes_shift/mm_C4.txt', 
    'stokes_shift/mm_qm1.txt', 
    'stokes_shift/qm2.txt')
dt = 4 #    in femtoseconds
titles = ['Stripped', 'H-bonded', 'Pi Solvent', 'QM1', 'QM2']
save_names = ['stripped', 'mm_4hb', 'mm_C4', 'mm_qm1', 'qm2']
colors = ['black', '#21ADEF', '#D321FF', 'red', 'blue']

#   first run stokes shift plot
run_analysis()


for i in range(len(data_files)):
    if i != 4:
        continue
    data = np.loadtxt(data_files[i])
    times = data[:, 0]
    diff_interp = data[:, 3] - data[:, 4]
    diff_interp -= np.mean(diff_interp)
    # ax.plot(times, diff_interp, label=titles[i], color=colors[i])

    N = times.shape[0]
    T = times[1] - times[0]
    fourier_transform = fft(diff_interp)[0:N//2]
    fourier_transform_freq = fftfreq(N, T)[:N//2]
    # ax.plot(fourier_transform_freq, 2.0/N * np.abs(fourier_transform), label=titles[i])
    

    break

ax.legend()
plt.show()