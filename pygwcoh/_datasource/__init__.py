"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from .waveform import Template
from .datasrc import load_data_from_ifo

__all__ = ['Template', 'load_data_from_ifo']
    