"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from scipy.signal import resample as scipy_resample
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d


def resample(data, fs_old, fs_new):
    if fs_new != fs_old:
        nsamp = int(len(data) * fs_new / fs_old)
        new = scipy_resample(data, nsamp)
        return new
    else:
        return data

def whiten(strain, interp_psd, fs):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, 1./fs)
    df = freqs[1] - freqs[0]
    # freqs1 = np.linspace(0,2048.,Nt/2+1)
    
    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    #norm = 1./np.sqrt(fs/2)
    white_hf = hf / np.sqrt(interp_psd(np.abs(freqs)))# * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    sigmasq = 4 * (hf * hf.conjugate() / interp_psd(freqs)).sum() * df
    return white_ht, sigmasq

def get_psdfun(data, fs, NFFT = None, NOVL = None, window = False):
    data_psd, freqs = get_psd(data, fs = fs, NFFT = NFFT, window=window, NOVL=NOVL)
    return interp1d(freqs, data_psd)

def get_psd(data, fs, NFFT = None, NOVL = None, window = False):
    if NFFT is None:
        NFFT = 4*fs
    if window:
        psd_window = np.blackman(NFFT)
    else:
        psd_window = None
    # and a 50% overlap:
    if NOVL is None:
        NOVL = NFFT/2
    data_psd, freqs = mlab.psd(data, Fs = fs, NFFT = NFFT, window=psd_window, noverlap=NOVL)
    return data_psd, freqs

def padinsert(a,b,length = None):
    if length is None:
        length = max(len(a), len(b))
    if length < max(len(a), len(b)):
        raise ValueError('pad length should not less than max(len(a), len(b))')
    ax = np.pad(a, (0, length - len(a)), mode = 'constant')
    bx = np.pad(b, (0, length - len(b)), mode = 'constant')
    return ax, bx

def cutinsert(a,b,cutpct = None):
    if cutpct is None or cutpct > 1 or cutpct < 0:
        cutpct = 0.5
    b_cut = b[:int(len(b)*cutpct)]
    return padinsert(a, b_cut)        

def correlate_real(stilde, htilde, fs, power_vec, NFFT):
    sigmasq = 1 * (htilde * htilde.conjugate() / power_vec).sum() * fs / NFFT
    corr = 1 * stilde * htilde.conjugate() / power_vec
    corr_time = np.fft.irfft(corr) / np.sqrt(np.abs(sigmasq)) * fs
    return corr_time
