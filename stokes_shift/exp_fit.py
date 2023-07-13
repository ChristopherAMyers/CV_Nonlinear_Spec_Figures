import numpy as np
from scipy.optimize import curve_fit, minimize

CM_2_INV_FS = 3.0E-5

def gauss_exp_func(t, A_gauss, t_gauss, A_exp, t_exp):
    return A_gauss*np.exp(-t*t/t_gauss**2) + A_exp*np.exp(-t/t_exp)

def fit_decay(t, y, t_exp=None):

    bounds_lower = [0, 0, 0]
    bounds_upper = [np.inf, np.inf, np.inf]
    random_limits_l = [0.0, 50, 0.0]
    random_limits_u = [1.0, 1000, 1.0]
    if t_exp is None:
        func_to_fit = lambda t, A1, t1, A2, t2 : gauss_exp_func(t, A1, t1, A2, t2)
        guess = [0.5, 20, 0.5, 2]
        bounds_lower.append(0.5)
        bounds_upper.append(np.inf)
        random_limits_l.append(50)
        random_limits_u.append(1000)
    else:
        func_to_fit = lambda t, A1, t1, A2 : gauss_exp_func(t, A1, t1, A2, t_exp=t_exp)
        guess = [0.5, 20, 0.5]

    best_params = []
    best_residual = np.inf
    best_fit_func = None
    for n in range(20):
        guess = []
        for l, u in zip(random_limits_l, random_limits_u):
            guess.append(np.random.random()*u + l)
        params = curve_fit(func_to_fit, t, y, guess, bounds=(bounds_lower, bounds_upper))[0]
        fit_func = func_to_fit(t, *params)
        residual = np.sum((y - fit_func)**2)

        if residual < best_residual:
            best_residual = residual
            best_params = np.copy(params)
            best_fit_func = np.copy(fit_func)

    return best_params, best_residual, best_fit_func

def decay_sin(t, *params):
    sin_params = params[4:]
    A_gauss, t_gauss, A_exp, t_exp = params[0:4]
    total_func = gauss_exp_func(t, *params[0:4])
    for i in range(0, len(sin_params), 4):
        freq, amp, phase, decay = sin_params[i:i+4]
        freq *= CM_2_INV_FS
        total_func += amp*np.sin(2*np.pi*freq*t + phase)*np.exp(-t/decay)

    return total_func

def decay_sin_jac(t, *params):
    jac = np.zeros((len(t), len(params)))
    sin_params = params[4:]
    A_gauss, t_gauss, A_exp, t_exp = params[0:4]

    #   exponential function part
    exp_1 = np.exp(-t*t/t_gauss**2)
    exp_2 = np.exp(-t/t_exp)
    jac[:, 0] = exp_1
    jac[:, 1] = A_gauss * exp_1 * (2*t*t/t_gauss**3)
    jac[:, 2] = exp_2
    jac[:, 3] = A_exp * exp_2 * (t/t_exp**2)

    #   sinusoidal function part
    for i in range(0, len(sin_params), 4):
        freq, amp, phase, decay = sin_params[i:i+4]
        freq *= CM_2_INV_FS

        sin_array = np.sin(2*np.pi*freq*t + phase)*np.exp(-t/decay)
        cos_array = np.cos(2*np.pi*freq*t + phase)*np.exp(-t/decay)
        jac[:, 4+i+0] = amp * cos_array * (2*np.pi*t) * CM_2_INV_FS
        jac[:, 4+i+1] =       sin_array
        jac[:, 4+i+2] = amp * cos_array
        jac[:, 4+i+3] = amp * sin_array * (t/decay**2)

    return jac

def numerical_jac_test(t, guess):
    ''' Tests to make sure the derivatives are being taken correctly'''
    act_deriv = decay_sin_jac(t, *guess)[20]
    num_deriv = np.zeros(len(guess))
    for j in range(len(guess)):
        params_p = np.copy(guess)
        params_n = np.copy(guess)
        params_p[j] += 0.001
        params_n[j] -= 0.001
        func_p = decay_sin(t[20], *params_p)
        func_n = decay_sin(t[20], *params_n)
        num_deriv[j] = (func_p - func_n)/0.002
        print("DERIV: ", act_deriv[j], num_deriv[j])
    exit()

def residual_func(t, y, func_params):
    return np.sum((y - decay_sin(t, *func_params))**2) / np.sum(y**2)

def fit_sinusoidal_decay(t, y, t_exp=None):

    guess_wavenumbers = np.array([500, 550, 600, 650])
    # guess_wavenumbers = np.array([350, 530, 580, 600])
    # guess_wavenumbers = np.linspace(300, 700, 7)
    n_waves = len(guess_wavenumbers)

    min_func = lambda args: residual_func(t, y, args)
    def callback(params):
        pass
        print("Callback: ", min_func(params))

    best_params = []
    best_rrmsd = np.inf
    best_fit_func = None
    for n in range(20):
        guess = [0.5, 20, 0.5, 200]
        phase_guesses = np.random.random(n_waves)*2*np.pi
        bounds_u = [np.inf, np.inf, np.inf, np.inf]
        bounds_l = [0, 0, 0, 0]
        for i in range(n_waves):
            bounds_l += [300, 0, 0, 0]
            bounds_u += [800, 2.0, 2*np.pi, np.inf]
            guess += [guess_wavenumbers[i], 0.3/n_waves, phase_guesses[0], 1000.0]
            
        # numerical_jac_test(t, guess)

        try:
            params = curve_fit(decay_sin, t, y, guess, bounds=(bounds_l, bounds_u), maxfev=1E4, jac=decay_sin_jac)[0]
            # res = minimize(min_func, guess,
            #                     method='Nelder-Mead',
            #                     bounds=np.transpose((bounds_l, bounds_u)),
            #                     callback=callback
            #                 )
            # params = res.x
        
            fit_func = decay_sin(t, *params)
            rrmsd = np.sqrt(np.sum((y - fit_func)**2) / np.sum(y**2))
            if rrmsd < best_rrmsd:
                best_rrmsd = rrmsd
                best_params = np.copy(params)
                best_fit_func = np.copy(fit_func)
            print(f"Trial {n+1:3d}, residual = {best_rrmsd}")
        except:
            print("Failed")


    return best_params, best_rrmsd, best_fit_func
