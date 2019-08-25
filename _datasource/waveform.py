"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from .._datatypes.series import TimeSeries
from .._utils import CEV, CallCommand_With_Output
from .._core import resample, QPlane
from .._core.qplane import DEFAULT_FRANGE, DEFAULT_MISMATCH

def CMD_lalsim_inspiral(exe,  
                        m1,
                        m2,
                        s1z,
                        s2z,
                        D,
                        srate,
                        f_ini,
                        approx):
    CMD = f'{exe} --m1={m1} --m2={m2} \
            --spin1z={s1z} --spin2z={s2z} \
            --distance={D} --sample-rate={srate} \
            --f-min={f_ini} --approx={approx} --inclination=0'
    return CMD


class Template(TimeSeries):
    def __init__(self, m1, m2, s1z, s2z,
                    fini = 20,
                    approx = 'SEOBNRv4',
                    srate = 4096,
                    D = 100,
                    duration = None,
                    info = 'template'):
        self._info = 'template'
        self._m1 = m1
        self._m2 = m2
        self._s1z = s1z
        self._s2z = s2z
        self._fini = fini
        self._srate = srate
        self._approx = approx
        self._D = D
        STATE, data = CallCommand_With_Output(self._CMD, name_out = '_template')
        if len(data) == 0:
            STATE = CEV.PROCESS_FAIL
        self._STATE = STATE
        if STATE is not CEV.SUCCESS:
            # regeneration
            fail = True
            for i in range(1,5):
                fs_retry = int(self._srate * i)
                STATE, data = CallCommand_With_Output(self._fCMD(fs_retry), name_out = '_template')
                if len(data) != 0:
                    fail = False
                    self._STATE = STATE
                    break
            if fail:
                self._STATE = STATE
                raise Exception(f'Cannot generate waveform.')
            else:
                # do resample
                hp = np.asarray(data[:,1])
                hp_new = resample(hp, fs_retry, self._srate)
                hc = np.asarray(data[:,2])
                hc_new = resample(hc, fs_retry, self._srate)
                ht = hp_new + 1.j*hc_new
        else:
            ht = np.asarray(data[:,1]) + 1.j*np.asarray(data[:,2])
        super(Template, self).__init__(ht, 0, self._srate, info = self._info)
        if duration is not None:
            self._check_duration(duration)

    def _check_duration(self, duration):
        if self.dtpeak > abs(duration):
            cut_idx = self.argpeak - int(abs(duration)*self.fs)
            self._value = self._value[cut_idx:]
    
    def _fCMD(self, fs):
        return CMD_lalsim_inspiral(exe = 'lalsim-inspiral',
                                   m1 = self._m1,
                                   m2 = self._m2,
                                   s1z = self._s1z,
                                   s2z = self._s2z,
                                   D = self._D,
                                   srate = fs,
                                   f_ini = self._fini,
                                   approx = self._approx)

    @property
    def distance(self):
        return self._D

    @property
    def _CMD(self):
        return self._fCMD(self._srate)

    @property
    def template(self):
        return self.value

    @property
    def STATE(self):
        return self._STATE
            
    @property
    def m1(self):
        return self._m1
    
    @property
    def m2(self):
        return self._m2
    
    @property
    def s1z(self):
        return self._s1z
    
    @property
    def s2z(self):
        return self._s2z
    
    @property
    def approx(self):
        return self._approx

    @property
    def duration(self):
        return len(self) / self._srate
    
    @property
    def argpeak(self):
        return np.argmax(np.abs(self.value))
    
    @property
    def dtpeak(self):
        return self.time[self.argpeak] - self.time[0]
    
    @property
    def time(self):
        return np.arange(0, len(self) / self.fs, 1./self.fs)

    @property
    def phase(self):
        return np.unwrap(np.angle(self.value))
    
    @property
    def phasedot(self):
        return np.gradient(self.phase) / np.gradient(self.time)

    @property
    def track(self):
        return self.get_track(0,0)

    def get_track(self, epoch_peak, extra_index = 5):
        idx_peak = self.argpeak
        idx_end = min(len(self), idx_peak + extra_index)
        time = self.time - self.time[idx_peak] + epoch_peak
        return self.time[:idx_end], self.phasedot[:idx_end] / (2*np.pi)
    
    def iter_fftQPlane(self, q, duration, fs, frange = None, mismatch = None):
        plane = QPlane(q, duration, fs, frange = frange, mismatch = mismatch)
        # For time shift
        if frange is None:
            frange = DEFAULT_FRANGE
        if mismatch is None:
            mismatch = DEFAULT_MISMATCH
        track_x, track_y = self.track
        track_x -= track_x[0]
        for qtile in plane:
            freq = qtile.frequency
            idx = np.where(np.abs(track_y - freq) == np.min(np.abs(track_y - freq)))[0][0]
            shift = track_x[idx]
            yield (shift, qtile)
