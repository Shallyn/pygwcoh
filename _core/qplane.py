"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from math import (log, ceil, pi, isinf, exp)
import warnings
from scipy.interpolate import InterpolatedUnivariateSpline, interp2d, interp1d

DEFAULT_FRANGE = (10, 1200)
DEFAULT_MISMATCH = 0.2

class QObject(object):
    """Base class for Q-transform objects

    This object exists just to provide basic methods for all other
    Q-transform objects.
    """
    # pylint: disable=too-few-public-methods

    def __init__(self, duration, sampling, mismatch=DEFAULT_MISMATCH):
        self.duration = float(duration)
        self.sampling = float(sampling)
        self.mismatch = float(mismatch)

    @property
    def deltam(self):
        """Fractional mismatch between neighbouring tiles

        :type: `float`
        """
        return 2 * (self.mismatch / 3.) ** (1/2.)

class QBase(QObject):
    """Base class for Q-transform objects with fixed Q

    This class just provides a property for Q-prime = Q / sqrt(11)
    """
    def __init__(self, q, duration, sampling, mismatch=DEFAULT_MISMATCH):
        super(QBase, self).__init__(duration, sampling, mismatch=mismatch)
        self.q = float(q)

    @property
    def qprime(self):
        """Normalized Q `(q/sqrt(11))`
        """
        return self.q / 11**(1/2.)


class QPlane(QBase):
    def __init__(self, q, duration, sampling,
                 frange = DEFAULT_FRANGE,
                 mismatch=DEFAULT_MISMATCH):
        super(QPlane, self).__init__(q, duration, sampling, mismatch=mismatch)
        self.frange = [float(frange[0]), float(frange[1])]

        if self.frange[0] == 0:  # set non-zero lower frequency
            self.frange[0] = 50 * self.q / (2 * pi * self.duration)
        if isinf(self.frange[1]) or self.frange[1] > sampling/2:  # set non-infinite upper frequency
            self.frange[1] = self.sampling / 2 / (1 + 1/self.qprime)
        
    def __iter__(self):
        """Iterate over this `QPlane`

        Yields a `QTile` at each frequency
        """
        # for each frequency, yield a QTile
        for freq in self._iter_frequencies():
            yield QTile(self.q, freq, self.duration, self.sampling,
                        mismatch=self.mismatch)
    
    def _iter_frequencies(self):
        minf, maxf = self.frange
        fcum_mismatch = log(maxf / minf) * (2 + self.q**2)**(1/2.) / 2.
        nfreq = int(max(1, ceil(fcum_mismatch / self.deltam)))
        fstep = fcum_mismatch / nfreq
        fstepmin = 1 / self.duration
        # for each frequency, yield a QTile
        for i in range(nfreq):
            yield (minf *
                   exp(2 / (2 + self.q**2)**(1/2.) * (i + .5) * fstep) //
                   fstepmin * fstepmin)

    @property
    def frequencies(self):
        """Array of central frequencies for this `QPlane`

        :type: `numpy.ndarray`
        """
        return np.array(list(self._iter_frequencies()))


class QTile(QBase):
    def __init__(self, q, frequency, duration, sampling,
                 mismatch=DEFAULT_MISMATCH):
        super(QTile, self).__init__(q, duration, sampling, mismatch=mismatch)
        self.frequency = frequency
        
    @property
    def bandwidth(self):
        return 2 * pi ** (1/2.) * self.frequency / self.q
    
    @property
    def ntiles(self):
        tcum_mismatch = self.duration * 2 * pi * self.frequency / self.q
        return next_power_of_two(tcum_mismatch / self.deltam)
    
    @property
    def windowsize(self):
        return 2 * int(self.frequency / self.qprime * self.duration) + 1
    
    def _get_indices(self):
        half = int((self.windowsize - 1) / 2)
        return np.arange(-half, half + 1)

    def get_frequency_indices(self):
        """Returns the index array of interesting frequencies for this row
        """
        return np.round(self._get_indices() + 1 +
                           self.frequency * self.duration).astype(int)

    
    def _get_bisquare_window(self):
        """Generate the bi-square window for this row

        Returns
        -------
        window : `numpy.ndarray`
        """
        # real frequencies
        wfrequencies = self._get_indices() / self.duration
        # dimensionless frequencies
        xfrequencies = wfrequencies * self.qprime / self.frequency
        # normalize and generate bi-square window
        norm = self.ntiles / (self.duration * self.sampling) * (
            315 * self.qprime / (128 * self.frequency)) ** (1/2.)
        return (1 - xfrequencies ** 2) ** 2 * norm

    @property
    def padding(self):
        """The `(left, right)` padding required for the IFFT

        :type: `tuple` of `int`
        """
        pad = self.ntiles - self.windowsize
        return (int((pad - 1)/2.), int((pad + 1)/2.))

    def get_window(self):
        size = int(self.duration * self.sampling/2) + 1
        window = np.zeros(size)
        sngl_window = self._get_bisquare_window()
        window[self.get_frequency_indices()] = sngl_window
        return window


def next_power_of_two(x):
    return 2**(ceil(log(x, 2)))

