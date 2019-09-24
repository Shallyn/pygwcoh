"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from ._datasource import load_data_from_ifo, load_data_from_cache
from ._datasource.noise import sim_gaussian_from_psd
from ._core.filter import correlate_real, padinsert, cutinsert
from ._core.skymap import nside2npix, pix2ang, Skymap
from ._core.utdk import calc_sngl_Gpc_and_shift
from ._datatypes.strain import CreateEmptySpectrum, gwStrain
from ._datatypes.series import TimeFreqSpectrum
from ._utils import interp2d_complex, LOGGER
from scipy import signal as scipysignal
import matplotlib.pyplot as plt

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
            datatict = load_data_from_cache(self._epoch, self._epoch + self._duration, 
                    ifos = ifos, fcache = cache, channel = channel, fs = self._fs)
        else:
            datadict = load_data_from_ifo(self._epoch, 
                    self._epoch + self._duration, 
                    ifos = ifos, channel = channel, 
                    fs = self._fs)
            for key in datadict:
                self._data.append(datadict[key])
        # Check duration
        strain = self._data[0]
        if strain.duration != self._duration:
            self._duration = strain.duration
        if strain.epoch != self._epoch:
            self._epoch = strain.epoch

    @property
    def broken(self):
        if len(self) < 2:
            return True
        else:
            return False

    def make_noise_from_psd(self, ifos, psddict):
        for ifo in ifos:
            funcpsd = psddict[ifo]
            length = int(30 * self._fs)
            freq = np.fft.rfftfreq(length, 1./self._fs)
            psd = funcpsd(freq)
            value = sim_gaussian_from_psd(freq, psd, self._fs, int(self._duration * self._fs))
            self._data.append(gwStrain(value, self._epoch, ifo, self._fs, info = f'Fake_{ifo} strain'))

    def calc_expected_snr(self, tmpl_inj, tmpl, gps, ra_inj, de_inj, snr_expected, psi = 0, phic = 0):
        SNR2 = 0
        hinj = tmpl_inj.template * np.exp(1.j*phic)
        hmatch = tmpl.template
        if len(hmatch) > len(hinj):
            hinj, hmatch = cutinsert(hinj, h)
        elif len(hmatch) < len(hinj):
            hinj, hmatch = padinsert(hinj, h)
        else:
            hinj = hinj
        ret = {}
        hrtilde = np.fft.rfft(hmatch.real)
        hitilde = np.fft.rfft(hmatch.imag)
        hfreq = np.fft.rfftfreq(hmatch.size, 1./self._fs)
        df = hfreq[1] - hfreq[0]
        for strain in self:
            at = strain.ifo_antenna_pattern(ra_inj, de_inj, psi, gps)
            signal = at[0]*hinj.real + at[1]*hinj.imag
            stilde = np.fft.rfft(signal)
            power_vec = strain.psdfun_set(hfreq)

            sigmasq_r = 1 * (hrtilde * hrtilde.conjugate() / power_vec).sum() * df
            snr_r = 1 * (stilde * hrtilde.conjugate() / power_vec).sum() * df / np.sqrt(np.abs(sigmasq_r))

            sigmasq_i = 1 * (hitilde * hitilde.conjugate() / power_vec).sum() * df
            snr_i = 1 * (stilde * hitilde.conjugate() / power_vec).sum() * df / np.sqrt(np.sqrt(sigmasq_i))

            snr = snr_r**2 + snr_i**2
            ret[strain.ifo] = np.sqrt(snr)
            SNR2 += snr
        rescaled =  snr_expected / np.sqrt(SNR2)
        for ifo in ret:
            ret[ifo] *= rescaled
        return ret, rescaled

    def calc_expected_track_SNR(self, q, 
                                tmpl_inj, tmpl, gps, 
                                ra_inj, de_inj, 
                                rescaled, 
                                psi = 0, phic = 0,
                                **kwargs):
        frange = kwargs.pop('frange', None)
        mismatch = kwargs.pop('mismatch', None)

        hinj = rescaled * tmpl_inj.template.copy() * np.exp(1.j*phic)
        hmatch = tmpl.template.copy()
        if len(hmatch) > len(hinj):
            hinj, hmatch = cutinsert(hinj, hmatch)
        elif len(hmatch) < len(hinj):
            hinj, hmatch = padinsert(hinj, hmatch)
        else:
            hinj = hinj
        ret = {}
        hrtilde = np.fft.rfft(hmatch.real)
        hitilde = np.fft.rfft(hmatch.imag)
        hfreq = np.fft.rfftfreq(hmatch.size, 1./self._fs)
        df = hfreq[1] - hfreq[0]
        stilde_dict = {}
        power_vec_dict = {}
        for strain in self:
            at = strain.ifo_antenna_pattern(ra_inj, de_inj, psi, gps)
            signal = at[0]*hinj.real + at[1]*hinj.imag
            stilde = np.fft.rfft(signal)
            power_vec = strain.psdfun_set(hfreq)

            stilde_dict[strain.ifo] = stilde
            power_vec_dict[strain.ifo] = power_vec

        frequencies = []
        ret_trackSNR = []
        for shift, qtile in tmpl.iter_fftQPlane(q = q, 
                                                duration = tmpl_inj.duration,
                                                fs = tmpl_inj.fs,
                                                frange = frange,
                                                mismatch = mismatch):
            freq = qtile.frequency
            qwindow = qtile.get_window()
            frequencies.append(freq)
            hrwindowed = hrtilde * qwindow
            hiwindowed = hitilde * qwindow
            SNR2_total = 0
            for strain in self:
                power_vec = power_vec_dict[strain.ifo]
                stilde = stilde_dict[strain.ifo]
                sigmasq_r = 1 * (hrwindowed * hrwindowed.conjugate() / power_vec).sum() * df
                snr_r = 1 * (stilde * hrwindowed.conjugate() / power_vec).sum() * df / np.sqrt(np.abs(sigmasq_r))

                sigmasq_i = 1 * (hiwindowed * hiwindowed.conjugate() / power_vec).sum() * df
                snr_i = 1 * (stilde * hiwindowed.conjugate() / power_vec).sum() * df / np.sqrt(np.abs(sigmasq_i))

                SNR2_total += snr_r.real**2 + snr_i.real**2
            ret_trackSNR.append(np.sqrt(SNR2_total))
        return frequencies, ret_trackSNR



    def make_injection(self, tmpl_inj, tmpl, gps, ra_inj, de_inj, snr_expected,
                        psi = 0, phic = 0):
        SNR2 = 0
        hinj = tmpl_inj.template.copy() * np.exp(1.j*phic)
        hmatch = tmpl.template.copy()
        if len(hmatch) > len(hinj):
            hinj, hmatch = cutinsert(hinj, hmatch)
        elif len(hmatch) < len(hinj):
            hinj, hmatch = padinsert(hinj, hmatch)
        else:
            hinj = hinj
        hrtilde = np.fft.rfft(hmatch.real)
        hitilde = np.fft.rfft(hmatch.imag)
        hfreq = np.fft.rfftfreq(hmatch.size, 1./self._fs)
        df = hfreq[1] - hfreq[0]
        ret = {}
        for strain in self:
            at = strain.ifo_antenna_pattern(ra_inj, de_inj, psi, gps)
            signal = at[0]*hinj.real + at[1]*hinj.imag
            stilde = np.fft.rfft(signal)
            power_vec = strain.psdfun_set(hfreq)

            sigmasq_r = 1 * (hrtilde * hrtilde.conjugate() / power_vec).sum() * df
            snr_r = 1 * (stilde * hrtilde.conjugate() / power_vec).sum() * df / np.sqrt(np.abs(sigmasq_r))

            sigmasq_i = 1 * (hitilde * hitilde.conjugate() / power_vec).sum() * df
            snr_i = 1 * (stilde * hitilde.conjugate() / power_vec).sum() * df / np.sqrt(np.abs(sigmasq_i))
            snr2 = snr_r.real**2 + snr_i.real**2
            ret[strain.ifo] = np.sqrt(snr2)
            SNR2 += snr2
        rescaled =  snr_expected / np.sqrt(SNR2)
        LOGGER.info(f'rescaled distance = {tmpl.distance / rescaled} Mpc\n')
        for strain in self:
            ret[strain.ifo] *= rescaled
            strain.make_injection(tmpl, gps, ra_inj, de_inj, rescaled, 
                        psi = psi, phic = phic)
        return ret, rescaled
                
    def set_psd(self, refpsd):
        for strain in self:
            if strain.ifo in refpsd:
                strain.set_psd(refpsd[strain.ifo])
            else:
                LOGGER.warning(f'Cannot set psd for strain {strain.ifo}\n')

    def plot_psd(self, fsave):
        freqs = np.fft.rfftfreq( int(self._duration * self._fs), 1./self._fs)
        plt.figure(figsize = (10,8))
        for strain in self:
            plt.loglog(freqs, power_vec, label = strain.ifo)
        plt.xlabel('freq [Hz]')
        plt.ylabel('$S_h \[Hz^{-1}\]$')
        plt.legend()
        plt.title('Power Spectral Density')
        plt.savefig(fsave, dpi = 200)
        plt.close()

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
        smap = Skymap(utdk2, geocent_times)
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

        retSPEC = [CreateEmptySpectrum(strain.ifo, info = None) for strain in self]
        frequencies = []
        for shift, qtile in tmpl.iter_fftQPlane(q = q, 
                                                duration = self._duration,
                                                fs = self._fs,
                                                frange = frange,
                                                mismatch = mismatch):
            freq = qtile.frequency
            qwindow = qtile.get_window()
            frequencies.append(freq)
            for i, strain in enumerate(self):
                stilde, hrtilde, hitilde, power_vec = \
                    strain.rfft_utils(tmpl, psd, cut, window)
                hrwindowed = hrtilde * qwindow
                hiwindowed = hitilde * qwindow
                df = self._fs / strain.size
                snr_r = correlate_real(stilde, hrwindowed, power_vec, df)
                snr_i = correlate_real(stilde, hiwindowed, power_vec, df)
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
            Gpc_matrix[i,0] = ar[0] * np.sqrt(strain.sigma2)
            Gpc_matrix[i,1] = ar[1] * np.sqrt(strain.sigma2)
            spec_matrix[:, :, i] = retSPEC[i].interpolate(geocent_times + delay)
        u,s,v = np.linalg.svd(Gpc_matrix)
        u_matrix = np.zeros([nfreq, ntime, ndet, ndet], u.dtype)
        u_matrix[:,:] = u
        coh = np.zeros([nfreq, ntime, ndet, ndet], np.complex)
        for i in range(ndet):
            coh[:,:,:,i] = np.multiply(u_matrix[:,:,:,i], spec_matrix)
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