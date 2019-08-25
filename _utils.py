"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""
import numpy as np
import time, os, sys
from pathlib import Path
import subprocess
from enum import Enum, auto
import logging
from scipy.interpolate import InterpolatedUnivariateSpline

#------Color Str------#
class COLOR(object):
    LOG = '\033[1;34;47m'
    WARNING = '\033[4;31m'
    DEBUG = '\033[4;31;46m'
    END = '\033[0m'

def logWrapper(log):
    return f'{COLOR.LOG}{log}{COLOR.END}'
LOG = logWrapper('LOG')
INFO = logWrapper('INFO')

def warningWrapper(warning):
    return f'{COLOR.WARNING}{warning}{COLOR.END}'
WARNING = warningWrapper('Warning')
ERROR = warningWrapper('Error')

def debugWrapper(debug):
    return f'{COLOR.DEBUG}{debug}{COLOR.END}'
MESSAGE = debugWrapper('Message')
DEBUG = debugWrapper('DEBUG')

class infoPrinter(object):
    def __init__(self, name = 'logger'):
        self._name = name

    def __call__(self, msg = 'Blank'):
        return sys.stderr.write(f'{msg}')

    def debug(self, msg = 'Blank'):
        return sys.stderr.write(f'{DEBUG}: {msg}')
    
    def warning(self, msg = 'Blank'):
        return sys.stderr.write(f'{WARNING}: {msg}')

    def log(self, msg = 'Blank'):
        return sys.stderr.write(f'{LOG}: {msg}')
    
    def info(self, msg = 'Blank'):
        return sys.stderr.write(f'{INFO}: {msg}')

    def error(self, msg = 'Blank'):
        return sys.stderr.write(f'{ERROR}: {msg}')
LOGGER = infoPrinter()

#------Exception Type------#
class CEV(Enum):
    SUCCESS = auto()
    PROCESS_FAIL = auto()
    VALUE_ERROR = auto()
    NORMAL = auto()
    TIMEOUT = auto()
    UNKNOWN = auto()
    
def CEV_parse_value(val):
    if isinstance(val, int):
        try:
            return CEV(val).name
        except:
            return CEV.UNKNOWN
    if isinstance(val, Iterable):
        ret = []
        for value in val:
            try:
                rst = CEV(value)
            except:
                rst = CEV.UNKNOWN
            ret.append(rst)
        return CEV_Array(ret)
    
class CEV_Array(object):
    def __init__(self, array):
        self._array = array
    
    def __iter__(self):
        for x in self._array:
            yield x
            
    @property
    def name(self):
        return np.array([itr.name for itr in self])
    
    @property
    def value(self):
        return np.array([itr.value for itr in self])
    
    def __str__(self):
        return '{}'.format(self._array)
    
    def __len__(self):
        return len(self._array)
    
    def __repr__(self):
        return self.__str__()
    
    def __format__(self):
        return self.__str__()

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.integer):
            return self._array[key]
        return self._getslice(key)

    def __setitem__(self, key, value):
        if isinstance(value, int):
            value = CEV_parse_value(value)
        if not isinstance(value, CEV):
            raise ValueError('The value to be set must be CEV or int.')
        self._array[key] = value

    def _getslice(self, index):
        if index.start is not None and index.start < 0:
            raise ValueError(('Negative start index ({}) is not supported').format(index.start))        
        return CEV_Array(self._array[index])

#-----Progress-----#
def Progress_with_bar(itr,N):
    arrow = '|'
    pcg_str = '%.2f'%min( 100, float(101*itr/N)) 
    pcg = float(pcg_str)
    for i in range(50):
        if 2 * i < pcg:
            arrow += '>'
        else:
            arrow += ' '
    arrow += '|'
    sys.stderr.write('\r')
    sys.stderr.write(pcg_str+ '%|'+arrow)
    sys.stderr.flush()
    time.sleep(0.02)
    
def Progress(itr,N, remarks = ''):
    pcg_str = '%.2f'%min( 100, float(101*itr/N)) 
    sys.stderr.write('\r')
    sys.stderr.write(f'{remarks}|{pcg_str}%|')
    sys.stderr.flush()
    time.sleep(0.02)
    
def Progress_time(dt, itr, N, remarks = None):
    tr_str = '%.1f'%( (dt+0.02) * (N-itr) / 60)
    pcg_str = '%.2f'%min(100, float(101*itr/N))
    sys.stderr.write('\r')
    if remarks is None:
        printout = pcg_str+ '%|time remain: '+tr_str+' min'
    else:
        printout = pcg_str+ '%|time remain: '+tr_str+' min-' + remarks
    sys.stderr.write(printout)
    sys.stderr.flush()
    time.sleep(0.02)


#------command line tools-------#
def CallCommand(CMD, out = None, err = None, timeout = 60):
    if isinstance(out, str):
        fileout = Path(out)
        pout = open(fileout, "w")
    else:
        pout = subprocess.PIPE
    if isinstance(err, str):
        fileerr = Path(err)
        perr = open(fileerr, "w")
    else:
        perr = subprocess.PIPE
    obj = subprocess.Popen(CMD,stdout=pout,stderr=perr,shell=True)
    t_bgn = time.time()
    while(True):
        if obj.poll() is not None:
            if hasattr(pout, 'close'):
                pout.close()
            if hasattr(perr, 'close'):
                perr.close()
            return True
        time_cost = time.time() - t_bgn
        if time_cost > timeout:
            obj.terminate()
            if hasattr(pout, 'close'):
                pout.close()
            if hasattr(perr, 'close'):
                perr.close()
            return False
        time.sleep(0.5)
        
def CallCommand_With_Output(CMD, name_out = '_out', timeout = 60):
    fname = Path(name_out)
    pout = open(fname, "w")
    obj = subprocess.Popen(CMD, stdout = pout, shell = True)
    t_start = time.time()
    while(True):
        if obj.poll() is not None:
            data = np.loadtxt(name_out)
            pout.close()
            if fname.exists():
                os.remove(fname)
            return CEV.SUCCESS, data
        time_cost = time.time() - t_start
        if time_cost > timeout:
            obj.terminate()
            pout.close()
            if fname.exists():
                os.remove(fname)
            return CEV.TIMEOUT, None
        time.sleep(0.5)

class Commander(object):
    """
    Convert Python function parameters to command options
    options format: --option-format, corresponding to option_format
    """
    def __init__(self, options):
        self._options = options
        self._funcpms_dict = {}
        for opt in options:
            pms = self._opt2pms(opt)
            if pms not in self._funcpms_dict:
                self._funcpms[pms] = opt
            else:
                sys.stderr.write(f'{WARNING}:Duplicate options: {opt}.\n')

    def _opt2pms(self, option):
        if option[:2] != '--':
            raise ValueError(f'Invalid option: {option}')
        pms = '_'.join(option[2:].split('-'))
        return pms

    def __call__(self, *args, **kwargs):
        cmd_options = []
        for pms in kwargs:
            opt = self._funcpms_dict[pms]
            val = kwargs[pms]
            cmd_option.append(opt)
            if val is not None:
                cmd_option.append(f'{val}')
        return ' '.join(cmd_options)

#-----switch method-----#
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

#------Other Method------#
def myresample(t, h, fs):
    dt_new = 1./fs
    t -= t[0]
    h_itp = interp1d(t, h)
    lth = t[-1] - t[0]
    t_new = np.arange(0, lth * 1.1, dt_new)
    idx_end = get_idx(t_new, lth)
    lth_new = len(t_new)
    h_new = np.zeros(lth_new, dtype = h.dtype)
    h_new[:idx_end] = h_itp(t_new[:idx_end])
    return t_new, h_new

class interp1d_complex(object):
    def __init__(self, x, y, w=None, bbox = [None]*2, k=3,ext=0, check_finite = False):
        yreal = y.real
        yimag = y.imag
        self._func_real = InterpolatedUnivariateSpline(x,yreal,w=w,bbox=bbox,k=k,ext=ext,check_finite=check_finite)
        self._func_imag = InterpolatedUnivariateSpline(x,yimag,w=w,bbox=bbox,k=k,ext=ext,check_finite=check_finite)

    def __call__(self, x):
        return self._func_real(x) + 1.j*self._func_imag(x)

class interp2d_complex(object):
    def __init__(self, x, y, z, kind = 'cubic'):
        zreal = z.real
        zimag = z.imag
        self._func_real = interp2d(x,y,zreal, kind = kind)
        self._func_imag = interp2d(x,y,zimag, kind = kind)
    
    def __call__(self, x, y):
        return self._func_real(x,y) + 1.j * self._func_imag(x,y)

