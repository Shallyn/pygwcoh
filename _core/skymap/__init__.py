"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from .pix import nside2npix, pix2ang, npix2nside
from .skymap import mollview, graticule, MollweideProj
from pathlib import Path
import matplotlib.pyplot as plt


__all__ = ['nside2npix', 
           'pix2ang', 
           'mollview', 
           'graticule', 
           'MollweideProj',
           'Skymap']


class Skymap(object):
    def __init__(self, utdk2, geocent_times):
        # utdk2 [ntimd, npix, ndet]
        self._geocent_times = geocent_times
        coh_snr_all = np.sqrt(np.sum(utdk2[:,:,:2], axis = 2)
        self._coh_snr = coh_snr_all.max(axis=1)
        self._projector  = MollweideProj()
        self._coh_snr_pix = coh_snr_all.max(axis=0))
        self._nside = npix2nside(len(self._coh_snr_pix))
        self._NULL = False
        if utdk2.shape[2] > 2:
            self._NULL = True
            null = np.sum(utdk2[:,:,2:], axis=2)
            null_snr2 = null.sum(axis=0) / utdk2.shape[0]
            self._null_snr_pix = np.sqrt(null_snr2)
    @property
    def max_gps_time(self):
        idx = self._coh_snr.argmax()
        return self._geocent_times[idx]

    @property
    def max_ra_de(self):
        max_de,max_ra = pix2ang(self._nside,np.argmax(self._coh_snr_pix))
        max_de = max_de[0]
        max_ra = max_ra[0]
        return max_ra - np.pi, np.pi/2 - max_de


    @property
    def NULL(self):
        return self._NULL

    def plot_skymap(self, prefix, 
                    plot_peak = True, 
                    ra_inj = None,
                    de_inj = None):
        prefix = Path(prefix)
        if plot_peak:
            max_ra,max_de = self.max_ra_de
            x1,y1 = self._projector.ang2xy(np.array([np.pi/2 - max_de, max_ra + np.pi]))
        mollview(self._coh_snr_pix,title=f'Coherent SNR, max = ({max_ra},{max_de})')
        graticule(coord='G',local=True)
        if plot_peak:
            plt.plot(x1,y1,'rx')

        if ra_inj is not None and de_inj is not None:
            x2,y2 = self._projector.ang2xy(np.array([np.pi/2 - de_inj, np.pi + ra_inj]))
            plt.plot(x2,y2,'r+')
        plt.savefig(prefix/'Skymap_Coherent.png', dpi = 200)
        plt.close()

        if self.NULL:
            mollview(self._null_snr_pix, title = 'null SNR')
            graticule(coord='G',local = True)
            if plot_peak:
                plt.plot(x1, y1, 'rx')
            plt.savefig(prefix/'Skymap_Null.png', dpi = 200)
            plt.close()

