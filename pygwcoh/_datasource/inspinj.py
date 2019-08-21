"""
This is the module for gravitational wave coherent search.
Construct inspinj fake frame data
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from .._utils import infoPrinter, CallCommand, Commander
from pathlib import Path

_INSPINJ_OPTIONS = \
    ['--seed', '--f-lower', '--waveform', 
    '--gps-start-time', '--gps-end-time', '--t-distr', '--time-step', '--time-interval',
    '--l-distr', '--longitude', '--latitude', 
    '--d-distr', '--min-distance', '--max-distance',
    '--snr-distr', '--min-snr', '--max-snr', '--min-coinc-snr',
    '--ligo-psd', '--ligo-fake-psd', '--ligo-start-freq', '--virgo-psd', '--virgo-fake-psd', '--ifos',
    '--i-distr', '--polarization', '--incl-std', '--fixed-inc', '--max-inc',
    '--coa-phase-distr', '--fixed-coa-phase',
    '--m-distr', '--min-mass1', '--max-mass1', '--min-mass2',
    '--min-mtotal', '--max-mtotal', '--fixed-mass1', '--fixed-mass2',
    '--mean-mass1', '--stdev-mass1', '--mean-mass2', '--stdev-mass2',
    '--min-mratio', '--max-mratio', '--mass1-points', '--mass2-points',
    '--diable-spin', '--enable-spin', '--spin-gaussian', '--aligned',
    '--min-spin1', '--max-spin1', '--mean-spin1', '--stdev-spin1',
    '--min-spin2', '--max-spin2', '--mean-spin2', '--stdev-spin2']

_DEFAULT_INSPINJ = 'lalapps_inspinj'

class _gwInspinj(object):
    """
    Generate lalapps_inspinj xml file
    """
    def __init__(self, fname, executable = _DEFAULT_INSPINJ):
        self._file = Path(fname)
        self._optparser = Commander(_INSPINJ_OPTIONS)
        self._exe = executable
    
    @property
    def exists(self):
        return self._file.exists()

    def get_inj(self):
        pass

    def make_inj_file(self, *args, **kwargs):
        """
        Construct Shell Command to call lalapps_inspinj
        """
        options = self._optparser(*args, **kwargs)
        CMD = f'{self._exe} {options} --output {self._file}'
        return CallCommand(CMD)

