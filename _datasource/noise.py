"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np

def sim_gaussian_segment(freq, psd, Nstride):
    Ntilde = len(freq)
    df = freq[1] - freq[0]
    sigma = 2*np.sqrt(psd / df)
    stilde = np.random.randn(Ntilde) * sigma + 1.j*np.random.randn(Ntilde) * sigma
    data = np.fft.irfft(stilde)
    Ndata = len(data)
    if Nstride == 0:
        return data
    elif Nstride == Ndata:
        Nstride = 0
    Nolp = Ndata - Nstride
    olp_x = np.arange(Nolp) * np.pi / (2*Nolp)
    overlap = data[Nstride:]
    data[:Nolp] = overlap * np.cos(olp_x) + data[:Nolp] * np.sin(olp_x)
    return data

def sim_gaussian_from_psd(freq, psd, fs, length):
    ret = np.zeros(length)
    df = freq[1] - freq[0]
    Ndata = int(fs / df)
    Ntilde = Ndata // 2 + 1
    Nstride = Ndata // 2
    if Ntilde > len(freq):
        raise ValueError("PSD not compatible with requested sample rate")
    
    length_generated = 0
    data = sim_gaussian_segment(freq, psd, 0)
    while(length_generated < length):
        if length_generated + Nstride < length:
            ret[length_generated:length_generated+Nstride] = data[:Nstride]
        else:
            ret[length_generated:length] = data[:length - length_generated]
        length_generated += Nstride
        data = sim_gaussian_segment(freq, psd, Nstride)
    return ret
