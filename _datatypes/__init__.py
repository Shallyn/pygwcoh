"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from .series import TimeSeries, TimeFreqSpectrum
from .strain import gwStrain
from .detector import Detector

__all__ = ['TimeSeries',
           'gwStrain',
           'Detector',
           'TimeFreqSpectrum']
