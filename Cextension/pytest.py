#! /Library/Frameworks/Python.framework/Versions/3.6/bin/python3
"""
    
   This is a test script for pyGpc
    
"""

import PyGWCOH as pg
import numpy as np
from pygwcoh._core.skymap import nside2npix, pix2ang
import time

epoch = 1186624818
FS = 4096

nside = 32
npix=nside2npix(nside)
tta, phi = pix2ang(nside, np.arange(npix))
ra = phi-np.pi
de = -tta+np.pi/2

#SNRtime = np.arange(epoch, epoch + 10, 1./FS)
#SNR = np.abs(np.random.randn(SNRtime.size))

SNR_times = np.linspace(epoch - 5, epoch + 5, 10*FS)
SNR = (np.sin(SNR_times) + 1) * (np.cos(SNR_times) + 2) * np.abs(np.random.randn(len(SNR_times)))

times = np.linspace(epoch - 1, epoch + 1, 100)

ifo = 'H1'
t_start = time.time()
Gpc_sngl, SNR_sngl = pg.Gpc_time_pix(SNR, SNR_times, ra, de, times, 'H1', epoch);
print('time cost: {}'.format(time.time() - t_start))
print(Gpc_sngl.shape)
print(SNR_sngl.shape)
