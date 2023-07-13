import numpy as np
import matplotlib.pyplot as plt
from os.path import join
import global_settings as GS
from exp_fit import fit_decay, fit_sinusoidal_decay

EV_2_AU = 0.0367493
FS_2_HA = 2.418884326505*10.0**(-2.0)
MD_STEP = 4.0*FS_2_HA
DECAY_LENGTH = 500.0*FS_2_HA

def compute_corr_func(data: np.ndarray, max_time: float = None):
    '''
        Computes correlation function from a veritical excitation 
        energy (VEE) trajectory

        Parameters
        ----------
        data : np.ndarray
            VEE trajectory (in EV)
        max_time : float
            crop all results so that the maximum time does not exceed 
            this value (in femtoseconds)

        Returns
        -------
        corr_data : np.ndarray
            Correlation function array with first row being the 
            times (in fs) and the second being the correlation function
            itself (in eV^2)

        corr_data_pos : np.ndarray
            Same as `corr_data`, but only for time >= 0
    '''

    #   import VEE data
    dim = data.shape[0]

    #   shift by the mean and convert to HA
    delta_U = (data - np.mean(data))

    #   form (symmetric) correlation function
    correlation = np.correlate(delta_U, delta_U, 'full')/dim
    times = np.arange(-dim+1, dim)*MD_STEP
    correlation *= np.exp(-np.abs(times)/DECAY_LENGTH)

    #   the part of the correlation functions that has positive times
    correlation_pos = correlation[dim-1:]
    times_pos = times[dim-1:]/FS_2_HA

    if max_time is not None:
        keep = times_pos < 2000
        correlation_pos, times_pos = correlation_pos[keep], times_pos[keep]

    return correlation_pos, times_pos

if __name__ == '__main__':

    data_folders = (
        join(GS.data_root_dir, 'gs-aimd-longer', 'stripped'), 
        join(GS.data_root_dir, 'gs-aimd-longer', 'mm_4hb'), 
        join(GS.data_root_dir, 'gs-aimd-longer', 'mm_C4'), 
        join(GS.data_root_dir, 'gs-aimd-longer', 'mm_qm1'), 
        join(GS.data_root_dir, 'gs-aimd-longer', 'qm2'),
    )

    data = np.loadtxt(join(data_folders[3], 'vee_traj1.dat'), usecols=0)
    correlation, times = compute_corr_func(data, 2000)

    n_smooth = 1
    if n_smooth == 1:
        smoothed_corr = correlation
        smoothed_time = times
    else:
        smooth_dim = len(correlation) - n_smooth
        smoothed_corr = np.zeros(smooth_dim)
        smoothed_time = times[0:smooth_dim]
        for n in range(smooth_dim):
            smoothed_corr[n] = np.mean(correlation[n:n+n_smooth])

    smoothed_corr -= smoothed_corr[-1]
    smoothed_corr /= smoothed_corr[0]

    fig, ax = plt.subplots()

    # ax.plot(times, correlation, marker='.')
    ax.plot(smoothed_time, smoothed_corr)


    # params, residual, fit = fit_decay(smoothed_time, smoothed_corr)
    params, rrmsd, fit = fit_sinusoidal_decay(smoothed_time, smoothed_corr)
    ax.plot(smoothed_time, fit)

    exp_names = ['A_gauss', 't_gauss', 'A_exp', 't_exp']
    print(f"\nFinal RRMSD = {rrmsd:.8f}")
    print("\nExponential terms")
    for i in range(4):
        print(f'{exp_names[i]:10s} {params[i]:.4f}')

    sin_params = params[4:].reshape((4, -1))
    order = np.argsort(sin_params[:, 0])
    sin_params = sin_params[order]
    sin_param_names = ['Freq', 'Amp', 'Phase', 'Decay']

    print("\nSinusoidal terns")
    print(f'         {"Freq":>7s}  {"Amp":>7s}  {"Phase":>7s}  {"Decay":>7s}')
    for i, sp in enumerate(sin_params):
        print(f'Wave {i+1:2d}  {sp[0]:7.1f}  {sp[1]:7.3f}  {sp[2]*180/np.pi:7.1f}  {sp[3]:7.2f}')

    # for i in range(0, len(sin_params), 4):
    #     for j in range(4):
    #         print(f'{sin_param_names[j]:6s} {sin_params[i + j]:10.3f}')
    #     print()



    # ref_corr = np.loadtxt('vee_MD_correlation_function_cl.dat')
    # plt.plot(times, ref_corr/EV_2_AU/EV_2_AU)
    ax.set_xlabel('Time (fs)')
    ax.set_xlim(0)
    fig.tight_layout()
    plt.show()