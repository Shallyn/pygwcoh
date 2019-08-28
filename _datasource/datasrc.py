"""
This is the module for gravitational wave coherent search.
Get LIGO data...
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from gwpy.timeseries import TimeSeries
from pathlib import Path
import numpy as np
from .._utils import CEV, LOGGER
import sys
from .._datatypes import gwStrain
from .._core import resample

channel_dict = {'H1_CALIB':'H1:GDS-CALIB_STRAIN', 'H1_GATED':'H1:GDS-GATED_STRAIN', 'H1_IDQ': 'H1:IDQ-PGLITCH_OVL_16_4096' ,\
                'L1_CALIB':'L1:GDS-CALIB_STRAIN', 'L1_GATED':'L1:GDS-GATED_STRAIN', 'L1_IDQ': 'L1:IDQ-PGLITCH_OVL_16_4096' ,\
                'V1_CALIB':'V1:Hrec_hoft_16384Hz', 'V1_GATED':'V1:Hrec_hoft_16384Hz_Gated', 'V1_IDQ': None}
BROKEN_THRE = 0.3

class gwStrainSRC(object):
    def __init__(self, ifo, gpsstart, gpsend, channel, tag = 'kafka'):
        self._shmpath = Path(f'/dev/shm/{tag}')
        self._gpsstart = gpsstart
        self._gpsend = gpsend
        self._fprefix = '{0}-{1}_llhoft'.format(ifo[0], ifo)
        if channel in channel_dict:
            self._channel = channel_dict[channel]
        else:
            self._channel = channel
        self._ifo = ifo
        if self._ifo == 'H1' or self._ifo == 'L1':
            self._frame = '{}_HOFT_C00'.format(self._ifo)
        else:
            self._frame = 'V1Online'
    
    def load_data(self, fs = 4096):
        t_start = int(self._gpsstart)
        t_end = int(self._gpsend)
        sys.stderr.write(f'--Loading {self._channel} {t_start}...{t_end}\n')
        ret = load_data_from_shm(t_start, t_end, 
                                 self._ifo, self._channel, self._shmpath, fs)
        if not isinstance(ret, CEV):
            return ret
        LOGGER.warning('Failed to find data in shm, try using gwpy...\n')
        ret = load_data_from_gwpy(t_start, t_end, 
                                  self._ifo, self._channel, self._frame, fs)
        if not isinstance(ret, CEV):
            return ret
        LOGGER.warning('Failed to find data in gwpy, try manual...\n')
        ret = load_data_manual(self._gpsstart, self._gpsend, self._ifo, self._channel, fs)
        if not isinstance(ret, CEV):
            return ret
        LOGGER.warning('Failed to find data...\n')
        return CEV.PROCESS_FAIL

def load_data_from_ifo(tstart, tend, ifos, channel = 'GATED', fs = 4096):
    sys.stderr.write(f'--:Load data {tstart} .. {tend}\n')
    ret = dict()
    for ifo in ifos:
        gws = gwStrainSRC(ifo, tstart, tend, channel = f'{ifo}_{channel}')
        data = gws.load_data(fs)
        if not isinstance(data, CEV):
            ret[f'{ifo}'] = data
    return ret
    
            
class shmseg(object):
    def __init__(self, name, channel, epoch, fs = 16384):
        self._name = name
        self._ifo = name[:2]
        self._channel = channel
        self._value = np.array([])
        self._epoch = epoch
        self._fs = fs
        self._duration = 0
        self._broken_time = 0

    
    @property
    def deltat(self):
        return 1./self._fs
        
    def append(self, filename, duration):
        try:
            data = TimeSeries.read(filename, self._channel)
            self._value = np.concatenate([self._value, data])
        except:
            data = np.zeros(self._fs * duration)
            self._value = np.concatenate([self._value, data])
            self._broken_time += duration
        self._duration += duration
    
    @property
    def broken(self):
        if self._duration == 0:
            return False
        if self._broken_time / self._duration > BROKEN_THRE:
            return True
    
    def togwStrain(self, fs = 4096):
        if self.broken:
            return CEV.PROCESS_FAIL
        gs = gwStrain(self._value, self._epoch, self._ifo, self._fs, info = f'{self._ifo}_strain')
        if fs != self._fs:
            gs = gs.resample(fs)
        return gs
    



def load_data_from_gwpy(gpsstart, gpsend, ifo, channel, frame, fs = 4096):
    try:
        data = TimeSeries.find(channel, 
                                gpsstart, gpsend, 
                                frametype = frame, 
                                allow_tape=False)
        value = data.value
        srate = data.sample_rate.value
        epoch = data.epoch.value
        #duration = data.duration.value
        ret = gwStrain(value, epoch, ifo, srate, info = f'{ifo}_strain')
        if srate != fs:
            return ret.resample(fs)
        return ret
    except:
        return CEV.PROCESS_FAIL

def load_data_from_shm(gpsstart, gpsend, ifo, channel, shmpath, fs = 4096):
    
    t_start = int(gpsstart)
    t_end = int(gpsend)
    fprefix = '{0}-{1}_llhoft'.format(ifo[0], ifo)
    if 'Gated' in channel or 'GATED' in channel:
        name = f'{ifo}_GATED'
    else:
        name = f'{ifo}_CALIB'
    seg = shmseg(name, channel, gpsstart)
    srcpath = shmpath / ifo
    for t in range(t_start, t_end):
        filename = '{0}/{1}-{2}-{3}'.format(srcpath, fprefix, t, '1.gwf')
        seg.append(filename, 1)
    return seg.togwStrain(fs = fs)
    
def load_data_manual(gpsstart, gpsend, ifo, channel, fs = 4096):
    # This method will reset self.value & self.epoch & self.duration
    filelist, glist = find_data_path(ifo, gpsstart, gpsend)
    if len(filelist) == 0:
            return CEV.PROCESS_FAIL
    argsrt = np.argsort(np.asarray(glist))
    value = np.array([])
    epoch = int(gpsstart)
    for idx in argsrt:
        fname = filelist[idx]
        data = TimeSeries.read(fname, channel)
        srate = data.sample_rate.value
        gf0, gf1 = parse_datafile(fname.name)
        if epoch < gf0:
            dataidx0 = 0
            if len(value) == 0:
                epoch = gf0
        else:
            dataidx0 = int( (epoch - gf0) * srate)
        if gpsend > gf1:
            dataidx1 = data.size
        else:
            dataidx1 = int( (gpsend - gf0) * srate ) + 1
        value = np.concatenate([value, data.value[dataidx0:dataidx1]])
    ret = gwStrain(value, epoch, ifo, srate, info = f'{ifo}_strain')
    if fs != srate:
        return ret.resample(fs)
    return ret
        

def find_data_path(ifo, gpsstart, gpsend):
    if ifo == 'H1' or ifo == 'L1':
        shmpath = Path('/hdfs/frames/O3/hoft/{}'.format(ifo))
        prefix = '{0}-{1}_HOFT_C00-'.format(ifo[0], ifo)
    elif ifo == 'V1':
        shmpath = Path('/hdfs/frames/O3/V1Online')
        prefix = 'V-V1Online-'
    else:
        raise ValueError('Incorrect instrument name {}'.format(ifo))
    g5 = int(gpsstart / 1e5)
    srcpath = None
    for src in shmpath.iterdir():
        gsrc = int(src.name.split(prefix)[-1])
        if g5 == gsrc:
            srcpath = src
            break
    if srcpath is None:
        return [], None
    gpsstart = int(gpsstart)
    gpsend = int(gpsend)
    filelist = []
    glist = []
    datapath = shmpath / srcpath
    for src in datapath.iterdir():
        gf0, gf1 = parse_datafile(src)
        if gpsstart < gf1 and gpsend > gf0:
            filelist.append(src)
            glist.append(gf0)
    return filelist, glist

def parse_datafile(filename):
    tmp = filename.name.split('/')[-1].split('.gwf')[0].split('-')[-2:]
    gf0 = int(tmp[0])
    gf1 = gf0 + int(tmp[1])
    return gf0, gf1

def load_data_from_cache(tstart, tend, ifos, fcache, channel = 'GATED', fs = 4096):
    LOGGER.info(f'load data from {fcache}...')
    cache = LIGOCache(fcache)
    return cache.load_data(tstart, tend, ifos, fs, channel = channel)

# Ligo Cache
class LIGOCache(object):
    def __init__(self, filename):
        self._fcache = filename
        with open(self._fcache, 'r') as f:
            contents = f.read()
        self._contents = contents
                
    def _iter_segment(self, ifo):
        for line in self._contents.split('\n'):
            spt = line.split(' ')
            if len(spt) != 5:
                continue
            name, frame, gps_start_str, gps_duration_str, fname = spt
            ifo_name = name + '1'
            if ifo_name == ifo:
                gps = float(gps_start_str)
                timeseg = TimeSegment(gps, gps + float(gps_duration_str))
                yield frame, timeseg, fname[16:]
            else:
                continue
                
    def load_data(self, gps_start, gps_end, ifos, fs, channel = 'GATED'):
        """
        return dict(ifo : gwStrain)
        """
        #print(self._contents)
        gps_seg = TimeSegment(gps_start, gps_end)
        ret = dict()
        for ifo in ifos:
            channel_name = f'{ifo}_{channel}'
            if channel in channel_dict:
                channel_ifo = channel_dict[channel_name]
            else:
                channel_ifo = channel_name
            value = np.array([])
            epoch = gps_start
            check_epoch = True
            for frame, timeseg, fname in self._iter_segment(ifo):
                seg = gps_seg.check_overlap(timeseg)
                if seg is None:
                    continue
                if check_epoch:
                    epoch = seg.start
                    check_epoch = False
                data = TimeSeries.read(fname, channel_ifo)
                srate = data.sample_rate.value
                idx = seg.get_idx(timeseg.start, srate)
                data = data[idx]
                if fs != srate:
                    data = resample(data, srate, fs)
                value = np.insert(value, len(value), data)
            ret[ifo] = gwStrain(value, epoch, ifo, fs, info = f'{ifo}_strain')
        return ret



class TimeSegment(object):
    def __init__(self, start, end):
        self._start = start
        self._end = end
    
    @property
    def start(self):
        return self._start
    
    @property
    def end(self):
        return self._end
    
    def __str__(self):
        return f'({self.start}, {self.end})'
    
    def check_overlap(self, seg):
        if seg.end < self.start or seg.start > seg.end:
            return None
        if seg.start < self.start:
            start = self.start
        else:
            start = seg.start
        
        if seg.end > self.end:
            end = self.end
        else:
            end = seg.end
        return TimeSegment(start, end)

    def get_slice(self, epoch, fs):
        idx_start = int( fs*(self.start - epoch))
        idx_end = int( fs*(self.end - epoch) )
        return slice(idx_start, idx_end)
    
    def overlap_slice(self, seg, fs):
        """
        return overlap part of the two TimeSegment
        """
        if seg.end < self.start or seg.start > self.end:
            return None
        if seg.start < self.start:
            idx_start = int((self.start - seg.start) * fs)
        else:
            idx_start = 0
        
        if seg.end > self.end:
            idx_end = int((self.end - seg.start) * fs)
        else:
            idx_end = int((seg.end - seg.start) * fs)
        return slice(idx_start, idx_end)

