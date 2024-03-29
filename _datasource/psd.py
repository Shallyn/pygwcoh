"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from glue.ligolw import ligolw
from glue.ligolw.array import get_array
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from scipy.interpolate import interp1d
from pathlib import Path

@ligolw_array.use_in
@ligolw_param.use_in
class PSDContentHandler(ligolw.LIGOLWContentHandler):
        pass

def get_refpsd(refpsd):
    retdict = dict()
    xmldoc = ligolw_utils.load_filename(refpsd, contenthandler = PSDContentHandler, verbose = True)
    root_name = u"psd"
    xmldoc, = (elem for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == root_name)
    for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
        if elem.hasAttribute(u"Name") and elem.Name == u"REAL8FrequencySeries":
            ifo = ligolw_param.get_pyvalue(elem, u"instrument")
            # t, = elem.getElementsByTagName(ligolw.Time.tagName)
            a, = elem.getElementsByTagName(ligolw.Array.tagName)
            # dims = a.getElementsByTagName(ligolw.Dim.tagName)
            # f0 = ligolw_param.get_param(elem, u"f0")
            fpsd = interp1d(a.array[0,:], a.array[1,:])
            retdict[str(ifo)] = fpsd
    return retdict

def get_refpsd_xml(xml):
    retdict = dict()
    xmldoc = ligolw_utils.load_filename(refpsd, contenthandler = PSDContentHandler, verbose = True)
    root_name = u"psd"
    xmldoc, = (elem for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == root_name)
    for elem in xmldoc.getElementsByTagName(ligolw.LIGO_LW.tagName):
        if elem.hasAttribute(u"Name") and elem.Name == u"REAL8FrequencySeries":
            ifo = ligolw_param.get_pyvalue(elem, u"instrument")
            # t, = elem.getElementsByTagName(ligolw.Time.tagName)
            a, = elem.getElementsByTagName(ligolw.Array.tagName)
            # dims = a.getElementsByTagName(ligolw.Dim.tagName)
            # f0 = ligolw_param.get_param(elem, u"f0")
            retdict[str(ifo)] = a.array
    return retdict

def get_refpsd_from_dir(fdir, channel = 'CALIB', ret_func = True):
    fdir = Path(fdir)
    fout = []
    for file in fdir.iterdir():
        if channel in file.name.split('_'):
            fout.append(file)
    if len(fout) == 0:
        raise ValueError('No such channel data file')
    out = dict()
    for file in fout:
        for ifo in ['H1', 'L1', 'V1']:
            if ifo in file.name.split('_'):
                datapsd = np.loadtxt(file)
                psd = interp1d(datapsd[0,:], datapsd[1,:])
                out[ifo] = psd
    return out
