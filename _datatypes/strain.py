"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from . import TimeSeries, TimeFreqSpectrum
from .detector import Detector
from .._core.filter import padinsert, cutinsert, correlate_real, get_psdfun
from .._core.qplane import DEFAULT_FRANGE, DEFAULT_MISMATCH
from .._core import resample
from scipy import signal


# H1-118 L1-150 V1-53
def get_sigma2(ifo):
    if ifo == 'L1':
        return (150 * 2)**2
    if ifo == 'H1':
        return (118 * 2)**2
    if ifo == 'V1':
        return (53 * 2)**2

#-----------------------My Strain Series Class--------------------#
# Must be real series
class gwStrain(TimeSeries):
    def __init__(self, strain, epoch, ifo, fs, info = ''):
        super(gwStrain, self).__init__(value = strain, epoch = epoch, fs = fs, info = info)
        self._ifo = ifo
        Det = Detector(self._ifo)
        self._ifo_latitude = Det.latitude
        self._ifo_longtitude = Det.longtitude
        self._ifo_location = Det.location
        self._ifo_response = Det.response
        self._ifo_antenna_pattern = Det.antenna_pattern
        self._ifo_delay = Det.time_delay_from_earth_center
        self._psdfun_set = lambda x : 1
        self._sigma2 = get_sigma2(ifo)
    
    @property
    def ifo(self):
        return self._ifo

    @property
    def ifo_latitude(self):
        return self._ifo_latitude

    @property
    def sigma2(self):
        return self._sigma2
    
    @property
    def ifo_longtitude(self):
        return self._ifo_longtitude

    @property
    def ifo_response(self):
        return self._ifo_response

    @property
    def ifo_antenna_pattern(self):
        return self._ifo_antenna_pattern

    @property
    def ifo_delay(self):
        return self._ifo_delay
    
    def set_psd(self, psd):
        self._psdfun_set = psd

    @property
    def psdfun_set(self):
        return self._psdfun_set

    def psdfun(self):
        return get_psdfun(self.value, self.fs)

    @property
    def duration(self):
        return self.length

    def resample(self, fs):
        if fs != self.fs:
            new = resample(self.value, self.fs, fs)
            return gwStrain(new, self.epoch, self.ifo, fs, info = self._info)
        else:
            return self

    def matched_filter(self,
                       tmpl,
                       cut = None,
                       window = True,
                       psd = None):
        if not isinstance(tmpl, Template):
            raise TypeError(f'Incomatible type: {type(tmpl)}')
        shift = tmpl.dtpeak
        h = tmpl.value
        if len(h) > self.size:
            s, h = cutinsert(self.value, h)
        elif len(h) < self.size:
            s, h = padinsert(self.value, h)
        else:
            s = self.value
        if psd in ('self',):
            psdfun = self.psdfun()
        if psd in ('set',):
            psdfun = self._psdfun_set
        
        if window:
            try:   
                dwindow = signal.tukey(h.size, alpha=1./8)  # Tukey window preferred, but requires recent scipy version 
            except: 
                dwindow = signal.blackman(h.size)     
        else:
            dwindow = 1
            
        fs = self.fs
        stilde = np.fft.rfft(s * dwindow) / fs
        hrtilde = np.fft.rfft(h.real * dwindow) / fs
        hitilde = np.fft.rfft(h.imag * dwindow) / fs
        datafreq = np.fft.rfftfreq(h.size, 1./fs)
        df = abs(datafreq[1] - datafreq[0])

        power_vec = psdfun(np.abs(datafreq))

        if cut is not None:
            fmin, fmax = cut
            if fmin < min(abs(datafreq)):
                fmin = min(abs(datafreq))
            if fmax > max(abs(datafreq)):
                fmax = max(abs(datafreq))
            kmin = np.where( np.abs(datafreq - fmin) < df)[0][0]
            kmax = np.where( np.abs(datafreq - fmax) < df)[0][0]
            stilde[:kmin] = 0
            stilde[kmax:] = 0
            hrtilde[:kmin] = 0
            hrtilde[kmax:] = 0
            hitilde[:kmin] = 0
            hitilde[kmax:] = 0

        snr_r = correlate_real(stilde, hrtilde, self.fs, power_vec)
        snr_i = correlate_real(stilde, hitilde, self.fs, power_vec)
        SNR = snr_r + 1.j*snr_i
        return gwStrain(SNR, self.epoch + shift, self.ifo, self.fs, info = f'{self.ifo}_SNR')


    def qfilter(self,
                tmpl,
                q,
                psd = None,
                cut = None,
                frange = DEFAULT_FRANGE,
                mismatch = DEFAULT_MISMATCH,
                window = True):
        shift = tmpl.dtpeak
        h = tmpl.value
        if len(h) > self.size:
            s, h = cutinsert(self.value, h)
        elif len(h) < self.size:
            s, h = padinsert(self.value, h)
        else:
            s = self.value
        if psd in ('self',):
            psdfun = self.psdfun()
        if psd in ('set',):
            psdfun = self._psdfun_set
        
        if window:
            try:   
                dwindow = signal.tukey(h.size, alpha=1./8)  # Tukey window preferred, but requires recent scipy version 
            except: 
                dwindow = signal.blackman(h.size)     
        else:
            dwindow = 1
            
        fs = self.fs
        stilde = np.fft.rfft(s * dwindow) / fs
        hrtilde = np.fft.rfft(h.real * dwindow) / fs
        hitilde = np.fft.rfft(h.imag * dwindow) / fs
        datafreq = np.fft.rfftfreq(h.size, 1./fs)
        df = abs(datafreq[1] - datafreq[0])
        NFFT = stilde.size
        power_vec = psdfun(np.abs(datafreq))

        if cut is not None:
            fmin, fmax = cut
            if fmin < min(abs(datafreq)):
                fmin = min(abs(datafreq))
            if fmax > max(abs(datafreq)):
                fmax = max(abs(datafreq))
            kmin = np.where( np.abs(datafreq - fmin) < df)[0][0]
            kmax = np.where( np.abs(datafreq - fmax) < df)[0][0]
            stilde[:kmin] = 0
            stilde[kmax:] = 0
            hrtilde[:kmin] = 0
            hrtilde[kmax:] = 0
            hitilde[:kmin] = 0
            hitilde[kmax:] = 0
        
        outspec = CreateEmptySpectrum(self.ifo)
        for (shift, qtile) in tmpl.iter_fftQPlane(q = q, 
                                                  duration = self.duration,
                                                  fs = fs,
                                                  frange = frange,
                                                  mismatch = mismatch):
            qwindow = qtile.get_window()
            hrwindowed = hrtilde * qwindow
            hiwindowed = hitilde * qwindow
            snr_r = correlate_real(stilde, hrwindowed, fs, power_vec)
            snr_i = correlate_real(stilde, hiwindowed, fs, power_vec)
            snr = snr_r + 1.j*snr_i
            outspec.append(snr, freq, epoch=self.epoch+shift, fs=fs)
        return outspec


class gwStrainSpectrum(TimeFreqSpectrum):
    def __init__(self, ifo, array, epoch, fs, freqs, info = 'StrainSpectrum'):
        super(gwStrainSpectrum, self).__init__(array, epoch, fs, freqs, info)
        self._ifo = ifo
        Det = Detector(self._ifo)
        self._ifo_latitude = Det.latitude
        self._ifo_longtitude = Det.longtitude
        self._ifo_location = Det.location
        self._ifo_response = Det.response
        self._ifo_antenna_pattern = Det.antenna_pattern
        self._ifo_delay = Det.time_delay_from_earth_center
        self._psdfun_set = lambda x : 1
        self._sigma2 = get_sigma2(ifo)

    @property
    def ifo_latitude(self):
        return self._ifo_latitude

    @property
    def sigma2(self):
        return self._sigma2
    
    @property
    def ifo_longtitude(self):
        return self._ifo_longtitude

    @property
    def ifo_response(self):
        return self._ifo_response

    @property
    def ifo_antenna_pattern(self):
        return self._ifo_antenna_pattern

    @property
    def ifo_delay(self):
        return self._ifo_delay

    def set_psd(self, psd):
        self._psdfun_set = psd

    @property
    def psdfun_set(self):
        return self._psdfun_set
        

def CreateEmptySpectrum(ifo):
    array = np.array([])
    freqs = None
    epoch = None
    fs = 1
    empty = gwStrainSpectrum(ifo, array, epoch, fs, freqs, info = f'{ifo}_StrainSpectrum')
    return empty
    
