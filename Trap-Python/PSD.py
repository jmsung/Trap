# Power Spectral Density of Trap data

import numpy as np
import matplotlib.pyplot as plt
from ReadTdms import x, fs, t

def runningMean(x, N):
    return np.convolve(x, np.ones((N,))/N, mode='valid')

N = 100

y = x[1]
my = runningMean(y, N)
mt = runningMean(t, N)

# PSD
ps = np.abs(np.fft.fft(y))**2
dt = 1/fs
f = np.fft.fftfreq(y.size, dt)
idx = np.argsort(f)

mps = runningMean(ps[idx], N)
mf = runningMean(f[idx], N)


# Plot
plt.close('all')
plt.figure(1)
plt.loglog(f[idx], ps[idx], 'k')
#plt.loglog(mf[1000:], mps[1000:], 'y')
plt.xlabel('Frequency (Hz)'), plt.ylabel('Power (V^2/Hz)')

plt.show()

