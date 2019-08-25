"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

from optparse import OptionParser
from pathlib import Path

import numpy as np
import time, sys, os
import logging
from ._coherent import gwStrainCoherent
from ._utils import LOGGER
from ._datasource import Template, get_refpsd

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

logging.basicConfig(format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                    datefmt="%d-%M-%Y %H:%M:%S", 
                    level=logging.DEBUG)

DEFAULT_NSIDE = 16
DEFAULT_FRANGE = (20, 1200)
DEFAULT_PCOLORBINS = 100
DEFAULT_CMAP = 'jet'
DEFAULT_MISMATCH = 0.1
FIGSIZE_QSCAN = (14,6)

def get_proper_approx(m1, m2):
    if m1 + m2 > 6:
        return 'SEOBNRv4'
    else:
        return 'SpinTaylorT4'


def parseargs(argv):
    parser = OptionParser(description='Waveform Comparation With SXS')
    # Sample rate
    parser.add_option('--sample-rate', type = 'int', default = 4096, help = 'sample rate used.')
    
    # Source parameters
    parser.add_option('--m1', type = 'float', help = 'mass1 of this event, for template generation.')
    parser.add_option('--m2', type = 'float', help = 'mass2 of this event, for template generation.')
    parser.add_option('--s1z', type = 'float', help = 'spin1z of this event, for template generation.')
    parser.add_option('--s2z', type = 'float', help = 'spin2z of this event, for template generation.')
    parser.add_option('--fini', type = 'float', default = 20, help = 'Initial frequency for template generation(natual dimension).')
    parser.add_option('--approx', type = 'str', help = 'approx for template generation.')

    # Localize 
    parser.add_option('--ra', type = 'float', help = 'ra of this event, if added, will use this value.')
    parser.add_option('--de', type = 'float', help = 'dec of this event, if added, will use this value.')

    # Skymap resolution
    parser.add_option('--nside', type = 'int', default = DEFAULT_NSIDE, help = 'Nside for skymap pix.')

    # Qtransform setting
    parser.add_option('--Q', type = 'float', default = 10, help = 'Q value for qtransform')
    parser.add_option('--minf', type = 'float', default = DEFAULT_FRANGE[0], help = 'min Frequency.')
    parser.add_option('--maxf', type = 'float', default = DEFAULT_FRANGE[1], help = 'max Frequency.')
    parser.add_option('--mismatch', type = 'float', default = DEFAULT_MISMATCH, help = 'mismatch for qscan.')

    # Plot setting
    parser.add_option('--cmap', type = 'str', default = DEFAULT_CMAP, help = 'Plot color map type.')
    parser.add_option('--pcolorbins', type = 'int', default = DEFAULT_PCOLORBINS, help = 'color bins for pcolor mesh plot.')
    
    # Data source & Data saving
    parser.add_option('--gps', type = 'float', help = 'gps trigger time for this event.')

    # GraceDB id
    parser.add_option('--Sgraceid', type = 'str', help = 'GraceDB Super Event ID, if added, will load the Preferred event parameters.')
    parser.add_option('--graceid', type = 'str', help = 'GraceDB Event ID, if added, will load such Grace event parameters(Prefer).')

    # For data find
    parser.add_option('--stepback', type = 'int', default = 15, help = 'Used for GraceDB data load.')
    parser.add_option('--stepforward', type = 'int', default = 15, help = 'Used for GraceDB data load.')
    
    # Specify detector
    parser.add_option('--H1', action = 'store_true', help = 'Add detector')
    parser.add_option('--L1', action = 'store_true', help = 'Add detector')
    parser.add_option('--V1', action = 'store_true', help = 'Add detector')

    # Psd, cache
    parser.add_option('--ref-psd', type = 'str', help = 'prefix for reference psd, preferred.')
    parser.add_option('--channel', type = 'str', default = 'GATED', help = 'channel type, if local data used.')
    parser.add_option('--cache', type = 'str', help = 'Data cache for data loading')

    # Output
    parser.add_option('--prefix', type = 'str', default = '.', help = 'prefix for results saving.')

    # Log setting
    # FIXME

    # Other options
    parser.add_option('--injection', action = 'store_true', help = 'If added, will make an injection')
    parser.add_option('--track', action = 'store_true', help = 'If added, will plot track.')
    parser.add_option('--dimless', action = 'store_true', help = 'If added, will use nature dimemsion unit for initial frequency.')

    args = parser.parse_args(argv)
    return args

def main(argv = None):
    # Step.0 set logger...

    # Step.1 parse args...
    logging.info('Parsing args...')
    args, empty = parseargs(argv)
    ra = args.ra
    de = args.de
    
    gps = args.gps
    Sgraceid = args.Sgraceid
    graceid = args.graceid
    sback = args.stepback
    sfwd = args.stepforward
    fs = args.sample_rate
    
    m1 = args.m1
    m2 = args.m2
    s1z = args.s1z
    s2z = args.s2z
    fini = args.fini
    approx = args.approx

    if not args.H1 and not args.L1 and not args.V1:
        ifos = None
    else:
        ifos = []
        if args.H1:
            ifos.append('H1')
        if args.L1:
            ifos.append('L1')
        if args.V1:
            ifos.append('V1')
    
    nside = args.nside
    Q = args.Q
    frange = (args.minf, args.maxf)
    mismatch = args.mismatch
    cmaptype = args.cmap
    pcolorbins = args.pcolorbins
    
    prefix = args.prefix
    refpsd = args.ref_psd
    channel = args.channel
    cache = args.cache

    track = args.track
    injection = args.injection
    dimless = args.dimless

    # Step.2 load data...
    """
    Three ways to load data & source parameters: m1 m2 s1z s2z gps
        1. If graceid & Sgraceid(preferred) was specified,
            will try loading source parameters from GraceDB.
                if failed, will skip to method 2.
                if loaded, will use event end time of this event.
            if source parameters was not specified, will use parameters of this event in GraceDB.
            else, will use input parameters to generate waveform.
        2. If cache was specified, will read cache file to load data.
            Source parameters should be specified
        3. If gps was specified, will trying load file directly(C00).
        4. If injection was set, will make injection.
    """
    if graceid is not None or Sgraceid is not None:
        logging.info('Fetch GraceDB event...')
        from ._datasource.gracedb import GraceEvent, GraceSuperEvent
        # Parsing GraceID
        try:
            if graceid is not None:
                Gevt = GraceEvent(GraceID = graceid, verbose = True)
            else:
                Gevt = GraceSuperEvent(SGraceID = Sgraceid, verbose = True).Preferred_GraceEvent
            sngl = Gevt.get_sngl(Gevt.ifos[0])
        except:
            LOGGER.error('Cannot fetch GraceDB event...exit.\n')
            return -1
        # Setting parameters
        if m1 is None:
            m1 = sngl.mass1
        if m2 is None:
            m2 = sngl.mass2
        if s1z is None:
            s1z = sngl.spin1z
        if s2z is None:
            s2z = sngl.spin2z
        if gps is None:
            gps = Gevt.end_time
        if ifos is None:
            ifos = Gevt.ifos
    logging.info(f'Parameters:\n\t\
                    m1 = {m1}\n\t\
                    m2 = {m2}\n\t\
                    s1z = {s1z}\n\t\
                    s2z = {s2z}\n\t\
                    gps end time: {gps}')

    # Setting fini_SI
    import astropy.constants as cst
    if dimless:
        fini_SI = fini * cst.c.value**3 / ((m1 + m2) * cst.M_sun.value * cst.G.value)
    else:
        fini_SI = fini

    # Setting approx
    if approx is None:
        approx = get_proper_approx(m1, m2)

    # Now let's try loading data
    tstart = gps - sback
    tend = gps + sfwd
    
    # Create Coherent Object
    logging.info(f'Builing coherent strain, {tstart} ... {tend}')
    Strains = gwStrainCoherent(tstart, tend-tstart, fs = fs, verbose = True)
    if ifos is None:
        ifos = ['H1', 'L1', 'V1']

    # Loading data
    logging.info(f'Loading data {ifos}')
    Strains.load_data(cache = cache, ifos = ifos, channel = channel)

    # Setting psd
    if refpsd is None:
        LOGGER.error('reference-psd was not specified\n')
        return -1
    logging.info('Setting psd')
    psddict = get_refpsd(refpsd)
    Strains.set_psd(psddict)

    # Step.3 Making Template
    logging.info('Generating template...')
    tmpl = Template(m1 = m1, m2 = m2, s1z = s1z, s2z = s2z, 
                    fini = fini_SI, approx = approx, srate = fs, 
                    duration = 0.8 * sback)
    logging.info(f'Get template, duration = {tmpl.duration}')

    """
    Preparatory work complete
    Now let's start data analysis.
    """

    # Step.0 Pre-setting...
    # fsave prefix setting
    fsave = Path(prefix)
    if not fsave.exists():
        fsave.mkdir(parents=True)
        logging.info(f'mkdir for output prefix {fsave}')
    else:
        logging.info(f'Output prefix {fsave} exists.')

    # plot setting
    cmap = plt.get_cmap(cmaptype)

    """
    Step.1 Matched filtering & Skymap
    """
    # Call...
    logging.info('Matched filtering & Coherent...')
    SNRs, skymap = Strains.calc_coherent_snr_skymap(tmpl, nside, gps)

    """
    Step.2 Plot SNR & Skymap...
    """
    logging.info('Ploting coherent SNR skymap & SNR time series...')
    skymap.plot_skymap(prefix, plot_peak = True)
    for snr in SNRs:
        snr.plot(epoch = gps, fsave = prefix / f'SNR_{snr.ifo}.png', 
            pset = 'abs')
    return 0

    """
    Step.3 Qscan...
    """

    """
    Step. Coherent Qscan
    """

    return 0
