# Trap data - sinusoidal curve fitting

from scipy.optimize import curve_fit

# curve fitting
g_f1 = 1.1; g_A1 = 0.10; g_ph1 = 0.0; g_off1 = 0.0
p1=[g_f1, g_A1, g_ph1, g_off1]

# create the function we want to fit
def my_sin(t, freq, amplitude, phase, offset):
    return np.sin(2*np.pi*t * freq + phase) * amplitude + offset

# now do the fit
fit = curve_fit(my_sin, t, x, p0=p1)
x_fit = my_sin(t, *fit[0])
x_res = x- x_fit



