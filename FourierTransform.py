import numpy as np
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.integrate import simpson


def fourierTransform(t, y, num_points=100, interp_func=None, interp_args=None):

    times = np.copy(t)
    amps = np.copy(y)

    #   if interpolation function is supplied, generate new times and amplitudes
    if interp_func:
        if interp_args is None:
            interp_args = {}
        f = interp_func(t, y, **interp_args)
        times = np.linspace(times[0], times[-1], num_points)
        amps = f(times)
        print(amps)
        
    #   compute fft and it's corresponding frequencies
    dT = times[1] - times[0]
    N = times.shape[0]
    fft_amps = fft(amps, N)[0:N//2]
    fft_freq = fftfreq(N, dT)[0:N//2]

    result = {
        'times': times,
        'amps': amps,
        'fft': fft_amps,
        'freq': fft_freq
    }

    return result

def explicitFourierTransform(t, y, freq, num_points=100, interp_func=None, interp_args=None):

    times = np.copy(t)
    amps = np.copy(y)

    #   if interpolation function is supplied, generate new times and amplitudes
    if interp_func:
        if interp_args is None:
            interp_args = {}
        f = interp_func(t, y, **interp_args)
        times = np.linspace(times[0], times[-1], num_points)
        amps = f(times)
        print(amps)
        
    #   compute fft and it's corresponding frequencies
    dT = times[1] - times[0]
    N = freq.shape[0]
    ft = np.zeros(N)
    for j in range(N):
        integrand = np.exp(1j*freq[j]*times)*y
        ft[j] = simpson(integrand, times)

    result = {
        'times': times,
        'amps': amps,
        'ft': ft,
        'freq': np.copy(freq)
    }

    return result