"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from ._datasource import load_data_from_ifo
from ._core.filter import correlate_real
from ._core.skymap import nside2npix, pix2ang, Skymap
from ._core.utdk import calc_sngl_Gpc_and_shift, calc_sngl_shift
from ._datatypes.strain import CreateEmptySpectrum
from ._datatypes.series import TimeFreqSpectrum
from ._utils import interp2d_complex, LOGGER

DEFAULT_SBACK = 0.5
DEFAULT_SFWD = 0.5

class gwStrainCoherent(object):
    def __init__(self, epoch, duration, fs, verbose = False):
        self._epoch = epoch
        self._duration = duration
        self._data = []
        self._verbose = verbose
        self._fs = fs
    
    def __len__(self):
        return len(self._data)

    def __iter__(self):
        for strain in self._data:
            yield strain
        
    def load_data(self, cache = None, ifos = None, channel = 'GATED'):
        if cache is not None:
            pass
        else:
            datadict = load_data_from_ifo(self._epoch, 
                    self._epoch + self._duration, 
                    ifos = ifos, channel = channel, 
                    fs = self._fs)
            for key in datadict:
                self._data.append(datadict[key])
                
    def set_psd(self, refpsd):
        for strain in self:
            if strain.ifo in refpsd:
                strain.set_psd(refpsd[strain.ifo])
            else:
                LOGGER.warning(f'Cannot set psd for strain {strain.ifo}\n')

    def iter_matched_filter(self, tmpl, **kwargs):
        """
        Will cover SNR in gwStrain
        """
        for strain in self:
            yield strain.matched_filter(tmpl, **kwargs)

    def matched_filter(self, tmpl, **kwargs):
        for snr in self.iter_matched_filter(tmpl, **kwargs):
            self._data.appendSNR(snr)
        
    def iter_qfilter(self, tmpl, q, **kwargs):
        for strain in self:
            yield strain.qfilter(tmpl, q, **kwargs)

    def calc_coherent_snr_skymap(self, tmpl, 
                                 nside, gps_trigger, 
                                 trange = None, **kwargs):
        """
        This method will do matched filtering and calculate SNR for each channel.
        Then combine all SNR series and calculate coherent skymap series.

        return
            (SNR_ifo1, SNR_ifo2, ...), coh_SNR[trange, ra & de] 
        """
        psd = kwargs.pop('psd', 'set')
        cut = kwargs.pop('cut', None)
        window = kwargs.pop('window', True)
        npix = nside2npix(nside)
        theta,phi = pix2ang(nside,np.arange(npix))
        ra_pix = phi-np.pi
        de_pix = -theta+np. pi / 2
        ndet = len(self)

        if trange is None:
            sback = min(gps_trigger - self._epoch, DEFAULT_SBACK)
            sfwd = min(self._epoch + self._duration - gps_trigger, DEFAULT_SFWD)
        elif isinstance(trange, np.float) or isinstance(trange, np.int):
            sback = sfwd = trange
        else:
            sback, sfwd = trange
        if sback < 0 or sfwd < 0:
            raise Exception(f'Invalid trange: ({sback}, {sfwd})')
        geocent_times = np.arange(gps_trigger - sback, gps_trigger + sfwd, 1./self._fs)
        ntime = len(geocent_times)
        Gpc_matrix = np.zeros([npix, ndet, 2], np.float)
        snr_matrix = np.zeros([ntime, npix, ndet], np.complex)
        retSNR = []

        for i, strain in enumerate(self):
            #stilde, hrtilde, hitilde, power_vec = \
                #strain.rfft_utils(tmpl, psd, cut, window)
            gwSNR = strain.matched_filter(tmpl, cut = cut, window = window, psd = psd)
            retSNR.append(gwSNR)
            Gpc_sngl, snr_sngl = \
                calc_sngl_Gpc_and_shift(gwSNR, geocent_times, ra_pix, de_pix, gps_trigger)
            Gpc_matrix[:,i,:] = Gpc_sngl * np.sqrt(gwSNR.sigma2)
            snr_matrix[:,:,i] = snr_sngl
        u,s,v = np.linalg.svd(Gpc_matrix)
        u_matrix = np.zeros([ntime, npix, ndet, ndet], u.dtype)
        u_matrix[:] = u
        ndet = u_matrix.shape[-1]
        utdk = np.zeros([ntime, npix, ndet, ndet], np.complex)
        for i in range(u.shape[-1]):
            utdk[:,:,:,i]  = np.multiply(u_matrix[:,:,:,i], snr_matrix)
        utdk = np.sum(utdk, axis = 2)
        utdk2 = np.multiply(utdk, utdk.conjugate()).real
        smap = Skymap(utdk2)
        return retSNR, smap


    def calc_coherent_snr_qspectrum(self, tmpl, 
                                    q, gps_trigger, ra, de,
                                    trange = None, **kwargs):
        """
        This method will apply Qwindow on template and do matched filtering for each qtile.
        Then combin all spectrum and calculate coherent spectrum.

        return
            (spec_ifo1, spec_ifo2, ...), coh_spec[trange]
        """
        psd = kwargs.pop('psd', 'set')
        cut = kwargs.pop('cut', None)
        window = kwargs.pop('window', True)
        frange = kwargs.pop('frange', None)
        mismatch = kwargs.pop('mismatch', None)

        if trange is None:
            sback = np.min(gps_trigger - self._epoch, DEFAULT_SBACK)
            sfwd = np.min(self._epoch + self._duration - gps_trigger, DEFAULT_SFWD)
        elif isinstance(trange, np.float) or isinstance(trange, np.int):
            sback = sfwd = trange
        else:
            sback, sfwd = trange
        if sback < 0 or sfwd < 0:
            raise Exception(f'Invalid trange: ({sback}, {sfwd})')
        geocent_times = np.arange(gps_trigger - sback, gps_trigger + sfwd, 1./self._fs)

        retSPEC = [CreateEmptySpectrum(strain.ifo) for strain in self]
        frequencies = []
        for idxf, (shift, qtile) in tmpl.iter_fftQPlane(q = q, 
                                                duration = self._duration,
                                                fs = self._fs,
                                                frange = frange,
                                                mismatch = mismatch):
            qwindow = qtile.get_window()
            frequencies.append(freq)
            for i, strain in enumerate(self):
                stilde, hrtilde, hitilde, power_vec = \
                    strain.rfft_utils(tmpl, psd, cut, window)
                hrwindowed = hrtilde * qwindow
                hiwindowed = hitilde * qwindow
                snr_r = correlate_real(stilde, hrwindowed, fs, power_vec)
                snr_i = correlate_real(stilde, hiwindowed, fs, power_vec)
                snr = snr_r + 1.j*snr_i
                retSPEC[i].append(snr, freq, epoch=strain.epoch+shift, fs=self._fs)
        frequencies = np.asarray(frequencies)
        ndet = len(self)
        ntime = len(geocent_times)
        nfreq = len(frequencies)
        Gpc_matrix = np.zeros([ndet, 2], np.float)
        spec_matrix = np.zeros([nfreq, ntime, ndet], np.complex)
        for i, strain in enumerate(self):
            ar, delay = strain.ifo_get_at_and_delay(ra, de, 0, gps_trigger)
            Gpc_matrix[i,:] = ar * np.sqrt(strain.sigma2)
            spec_matrix[:, :, i] = retSPEC[i].interp(geocent_times + delay)
        u,s,v = np.linalg.svd(Gpc_matrix)
        u_matrix = np.zeros([nfreq, ntime, ndet, ndet], u.dtype)
        u_matrix[:,:] = u
        coh = np.zeros([nfreq, ntime, ndet, ndet], np.complex)
        for i in range(ndet):
            coh[:,:,:,i] = np.multiply(u[:,:,:,i], spec_matrix)
        coh = np.sum(coh, axis=2)
        coh2 = np.multiply(coh, coh.conjugate()).real
        cohspec_array = np.sqrt(np.sum(coh2[:,:,:2], axis=2))
        coh_SPEC = TimeFreqSpectrum(cohspec_array, geocent_times[0], self._fs, frequencies, info = 'Coherent Q Spectrum')
        if coh2.shape[2] > 2:
            nullspec_array = np.sqrt(np.sum(coh2[:,:,2:], axis=2))
            null_SPEC = TimeFreqSpectrum(nullspec_array, geocent_times[0], self._fs, frequencies, info = 'Null Q Spectrum')
        else:
            null_SPEC = None
        return retSPEC, coh_SPEC, null_SPEC