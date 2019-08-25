"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from .pix import nside2npix, pix2ang
from .skymap import mollview, graticule, MollweideProj

__all__ = ['nside2npix', 
           'pix2ang', 
           'mollview', 
           'graticule', 
           'MollweideProj',
           'Skymap']


class Skymap(object):
    def __init__(self, utdk2):
        self._coh_snr_pix = np.sqrt(np.sum(utdk2[:,:,:2], axis = 2).max(axis=0))
        if utdk2.shape[2] > 2:
            null = np.sum(utdk2[:,:,2:], axis=2)
            null_snr2 = null.sum(axis=0) / utdk2.shape[0]
            self._null_snr_pix = np.sqrt(null_snr2)
        

