"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from .._datatypes.detector import Detector
from .._utils import LOGGER, interp1d_complex

try:
    from . import _PyGWCOH as pg
    CEXT = True
except:
    CEXT = False
    LOGGER.warning('Cannot import PyGWCOH\n')

def calc_sngl_Gpc_and_shift(gwSNR, times, ra_pix, de_pix, gps_geocent):
    if CEXT:
        return calc_sngl_Gpc_and_shift_cextension(gwSNR, times, ra_pix, de_pix, gps_geocent)
    else:
        return calc_sngl_Gpc_and_shift_python(gwSNR, times, ra_pix, de_pix, gps_geocent)

def calc_sngl_Gpc_and_shift_python(gwSNR, times, ra_pix, de_pix, gps_geocent):
    Det = Detector(gwSNR.ifo)
    ntime = len(times)
    npix = len(ra_pix)
    Gpc_sngl = np.zeros([npix, 2], np.float)
    snr_sngl = np.zeros([ntime, npix], gwSNR.value.dtype)
    fitp = interp1d_complex(gwSNR.time, gwSNR.value)

    for k, (ra, de) in enumerate(zip(ra_pix, de_pix)):
        ar, delay = gwSNR.ifo_get_at_and_delay(ra, de, 0, gps_geocent)
        Gpc_matrix[k, 0] = ar[0]
        Gpc_matrix[k, 1] = ar[1]
        snr_sngl[:, k] = fitp(times + delay)
    return Gpc_sngl, snr_sngl

def calc_sngl_shift(ifo, snr, times, ra, de, gps_trigger):
    pass

"""
For Cextension interface
"""

def calc_sngl_Gpc_and_shift_cextension(gwSNR, times, ra_pix, de_pix, gps_geocent):
    pass
