"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from .pix import nside2npix, pix2ang
from .skymap import mollview, graticule, MollweideProj
from pathlib import Path

__all__ = ['nside2npix', 
           'pix2ang', 
           'mollview', 
           'graticule', 
           'MollweideProj',
           'Skymap']


class Skymap(object):
    def __init__(self, utdk2):
        self._projector  = MollweideProj()
        self._coh_snr_pix = np.sqrt(np.sum(utdk2[:,:,:2], axis = 2).max(axis=0))
        self._NULL = False
        if utdk2.shape[2] > 2:
            self._NULL = True
            null = np.sum(utdk2[:,:,2:], axis=2)
            null_snr2 = null.sum(axis=0) / utdk2.shape[0]
            self._null_snr_pix = np.sqrt(null_snr2)

    @property
    def NULL(self):
        return self._NULL

    def plot_skymap(self, prefix, 
                    plot_peak = True, 
                    ra_inj = None,
                    de_inj = None):
        prefix = Path(prefix)
        mollview(self._coh_snr_pix,title=f'Coherent SNR, max = ({max_ra},{max_de})')
        graticule(coord='G',local=True)
        if plot_peak:
            max_de,max_ra = pix2ang(nside,np.argmax(self._coh_snr_pix))
            max_de = max_de[0]
            max_ra = max_ra[0]
            x1,y1 = self._projector.ang2xy(np.array([max_de, max_ra]))
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

