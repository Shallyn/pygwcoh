"""
This is the module for gravitational wave coherent search.
Get LIGO data...
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from ligo.gracedb.rest import GraceDb
from glue.ligolw import ligolw, lsctables
import sys
from .datasrc import gwStrainSRC
from .._utils import CEV, LOGGER
from glue import gpstime
from astropy.time import Time
import time

def get_nowtime():
    t = time.time()
    return gpstime.GpsSecondsFromPyUTC(t)

def GPS2ISO(gps):
    scTime = Time(int(gps), format = 'gps')
    return scTime.iso

class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass

client = GraceDb()
pipelist = ['gstlal', 'MBTAOnline', 'spiir', 'pycbc']

def get_events_from_time(gstart = None, gend = None, evttag = None, url = None, getid = True):
    if gend is None:
            gend = get_nowtime()
    if gstart is None:
            gstart = gend - 86400
    if evttag is None:
            evttag = ''
    event_list = []
    for pipeline in pipelist:
        evttag = pipeline + ' ' + str(gstart) + ' .. ' + str(gend)
        events = client.events(evttag)
        for evt in events:
            event_list.append(GraceEvent(event = evt))
    return event_list

def get_Sevents_from_time(gstart = None, gend = None):
    if gend is None:
            gend = get_nowtime()
    if gstart is None:
            gstart = gend - 86400
    Sevent_list = []
    Sevttag = f'{gstart} .. {gend}'
    Sevents = client.superevents(Sevttag)
    for Sevt in Sevents:
        ret = GraceSuperEvent(Sevent = Sevt)
        if ret.STAT:
            Sevent_list.append(ret)
    return Sevent_list

class GraceEvent(object):
    def __init__(self, event = None, GraceID = None, verbose = False):
        if event is None and GraceID is None:
            raise Exception('One of event and GraceID must be specified.')
        if event is not None and isinstance(event, dict):
            datatable = event
            self._GraceID = event['graceid']
        else:
            self._GraceID = GraceID
            response = client.event(GraceID)
            datatable = response.json()
        self._verbose = verbose
        try:
            self._load_coinc(datatable)
            self._load_sngl(datatable)
            self._STAT = True
        except:
            self._STAT = False
    
    @property
    def STAT(self):
        return self._STAT
        
    def _load_coinc(self, table):
        data = table['extra_attributes']['CoincInspiral']
        self._snr = data.pop('snr', 0)
        self._end_time = data.pop('end_time', 0) + 1e-9 * data.pop('end_time_ns', 0)
        self._combined_far = data.pop('combined_far', 0)
        ifos = data.pop('ifos', '')
        if len(ifos) == 0:
            self._ifos = []
        else:
            self._ifos = ifos.split(',')
        
    def _load_sngl(self, table):
        data = table['extra_attributes']['SingleInspiral']
        names = self.__dict__
        for ele in data:
            names[ele['ifo']] = DetTable(ele)
    
    def get_sngl(self, name):
        try:
            ret = self.__dict__[name]
            if isinstance(ret, DetTable):
                return ret
            LOGGER.warning(f'Invalid name {name}\n')
            return None
        except:
            LOGGER.warning(f'Invalid name {name}\n')
            return None
    @property
    def end_time(self):
        return self._end_time
    
    @property
    def ifos(self):
        return self._ifos
    
    @property
    def combined_far(self):
        return self._combined_far
    
    @property
    def GraceID(self):
        return self._GraceID
    
    @property
    def snr(self):
        return self._snr
    
    def load_data(self, stepback = 15, stepforward = 15, channel = 'GATED', fs = 4096):
        tstart = self.end_time - abs(stepback)
        tend = self.end_time + abs(stepforward)
        if self._verbose:
            sys.stderr.write(f'--Load data {tstart} .. {tend}\n')
        ret = dict()
        for ifo in self._ifos:
            gws = gwStrainSRC(ifo, tstart, tend, channel = f'{ifo}_{channel}')
            data = gws.load_data(fs)
            if not isinstance(data, CEV):
                ret[f'{ifo}'] = data
        return ret
    
def find_strain_all(gps_start, gps_end, channel = 'GATED', fs = 4096):
    ret = dict()
    for ifo in ['H1','L1','V1']:
        gws = gwStrainSRC(ifo, gps_start, gps_end, channel = f'{ifo}_{channel}')
        ret[f'{ifo}'] = gws.load_data(fs)
        LOGGER.warning(f'Failed to load {ifo} data {gps_start} .. {gps_end}\n')
    return ret

class DetTable(object):
    def __init__(self, table):
        names = self.__dict__
        for key in table:
            names[key] = table[key]
        self.gps = self.end_time + 1e-9 * self.end_time_ns


class GraceSuperEvent(object):
    def __init__(self, Sevent = None, SGraceID = None, verbose = False):
        if Sevent is None and SGraceID is None:
            raise Exception('One of Sevent and SGraceID must be specified.')
        if Sevent is not None and isinstance(Sevent, dict):
            datatable = Sevent
            self._SGraceID = datatable['superevent_id']
        else:
            self._SGraceID = SGraceID
            response = client.superevent(SGraceID)
            datatable = response.json()
        self._verbose = verbose
        self._load_table(datatable)
        
    def _load_table(self, datatable):
        self._GraceID_list = datatable['gw_events']
        self._GraceID_preferred = datatable['preferred_event']
        self._gps_start = datatable['t_start']
        self._gps = datatable['t_0']
        self._gps_end = datatable['t_end']
        self._far = datatable['far']
        self._GraceEvent = GraceEvent(GraceID = self._GraceID_preferred, verbose = self._verbose)
    
    def __iter__(self):
        for gid in self._GraceID_list:
            yield gid
            
    def __len__(self):
        return len(self._GraceID_list)
    
    @property
    def SGraceID(self):
        return self._SGraceID
    
    @property
    def STAT(self):
        return self.Preferred_GraceEvent.STAT
    
    @property
    def Preferred_GraceEvent(self):
        return self._GraceEvent
    
    @property
    def gps_start(self):
        return self._gps_start
    
    @property
    def gps_trigger(self):
        return self._gps
    
    @property
    def gps_end(self):
        return self._gps_end
    
    @property
    def far(self):
        return self._far
        

