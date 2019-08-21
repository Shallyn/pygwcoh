"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from ._datasource import load_data_from_ifo

class gwStrainCoherent(object):
    def __init__(self, epoch, duration, verbose = False):
        self._epoch = epoch
        self._duration = duration
        self._core = []
        self._verbose = verbose
    
    def __len__(self):
        return len(self._core)

    def __iter__(self):
        for strain in self._core:
            yield strain
        
    def load_data(self, cache = None, ifos = None, channel = 'GATED'):
        if cache is None:
            pass
        else:
            datadict = load_data_from_ifo(self._epoch, 
                    self._epoch + self._duration, 
                    ifos = ifos, channel = channel)
            for key in datadict:
                self._core.append(datadict[key])


    def iter_matched_filter(self, tmpl, **kwargs):
        """
        Will cover SNR in gwStrain
        """
        for strain in self:
            yield strain.matched_filter(tmpl, **kwargs)
        
    def iter_qfilter(self, tmpl, q, **kwargs):
        for strain in self:
            yield strain.qfilter(tmpl, q, **kwargs)


