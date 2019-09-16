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
from ._datasource import Template, get_refpsd, get_refpsd_from_dir

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

logging.basicConfig(format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                    datefmt="%d-%M-%Y %H:%M:%S", 
                    level=logging.INFO)

DEFAULT_Q = 10.0
DEFAULT_NSIDE = 16
DEFAULT_FRANGE = (20, 1200)
DEFAULT_PCOLORBINS = 100
DEFAULT_CMAP = 'jet'
DEFAULT_MISMATCH = 0.1
FIGSIZE_QSCAN = (14,6)
DEFAULT_SBACK = 30.0
DEFAULT_SFWD = 5.0
DEFAULT_SNR = 9.0
DEFAULT_MAX_DURATION = 30.0

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

    # Localize & injection
    parser.add_option('--injection', action = 'store_true', help = 'If added, will make an injection')
    parser.add_option('--ra', type = 'float', help = 'ra of this event, if added, will use this value.')
    parser.add_option('--de', type = 'float', help = 'dec of this event, if added, will use this value.')
    parser.add_option('--gps-inj', type = 'float', help = 'Used for injection.')
    parser.add_option('--m1-inj', type = 'float', help = 'Used for injection')
    parser.add_option('--m2-inj', type = 'float', help = 'Used for injection')
    parser.add_option('--s1z-inj', type = 'float', help = 'Used for injection')
    parser.add_option('--s2z-inj', type = 'float', help = 'Used for injection')
    parser.add_option('--approx-inj', type = 'str', help = 'Used for injection')
    parser.add_option('--psi', type = 'float', help = 'Polarisation, used for injection')
    parser.add_option('--phic', type = 'float', help = 'Initial orbit phase, used for injection')
    parser.add_option('--snr', type = 'float', default = DEFAULT_SNR, help = 'Expected coherent SNR, used for injection.')

    # Skymap resolution
    parser.add_option('--nside', type = 'int', default = DEFAULT_NSIDE, help = 'Nside for skymap pix.')

    # Qtransform setting
    parser.add_option('--Q', type = 'float', default = DEFAULT_Q, help = 'Q value for qtransform')
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
    parser.add_option('--stepback', type = 'float', default = DEFAULT_SBACK, help = 'Used for GraceDB data load.')
    parser.add_option('--stepforward', type = 'float', default = DEFAULT_SFWD, help = 'Used for GraceDB data load.')
    parser.add_option('--gps-start', type = 'float', help = 'You can set this gps start for data load.')
    parser.add_option('--gps-end', type = 'float', help = 'You can set this gps end for data load.')
    
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

    # Template Generation
    parser.add_option('--m1-tmpl', type = 'float', help = 'Parameters for template generation.')
    parser.add_option('--m2-tmpl', type = 'float', help = 'Parameters for template generation.')
    parser.add_option('--s1z-tmpl', type = 'float', help = 'Parameters for template generation.')
    parser.add_option('--s2z-tmpl', type = 'float', help = 'Parameters for template generation.')
    parser.add_option('--fini', type = 'float', default = 20.0, help = 'Initial frequency for template generation(natual dimension).')
    parser.add_option('--approx', type = 'str', help = 'approx for template generation.')

    # Other options
    parser.add_option('--track', action = 'store_true', help = 'If added, will plot track.')
    parser.add_option('--dimless', action = 'store_true', help = 'If added, will use nature dimemsion unit for initial frequency.')
    parser.add_option('--gaussian', action = 'store_true', help = 'If added, will generate gaussian noise using given psd.')
    parser.add_option('--plot-mode', type = 'str', default = 'all', help = 'Option for plot setting [all]')

    args = parser.parse_args(argv)
    return args

def main(argv = None):
    # Step.0 set logger...

    # Step.1 parse args...
    logging.info('Parsing args...')
    args, empty = parseargs(argv)

    injection = args.injection
    gps_inj = args.gps_inj
    ra = args.ra
    de = args.de
    m1_inj = args.m1_inj
    m2_inj = args.m2_inj
    s1z_inj = args.s1z_inj
    s2z_inj = args.s2z_inj
    psi = args.psi
    phic = args.phic
    approx_inj = args.approx_inj
    snr_expected = args.snr

    gps = args.gps
    Sgraceid = args.Sgraceid
    graceid = args.graceid
    sback = args.stepback
    sfwd = args.stepforward
    gps_start = args.gps_start
    gps_end = args.gps_end
    fs = args.sample_rate
    
    m1 = args.m1
    m2 = args.m2
    s1z = args.s1z
    s2z = args.s2z

    gaussian = args.gaussian
    fini = args.fini
    approx = args.approx
    m1_t = args.m1_tmpl
    m2_t = args.m2_tmpl
    s1z_t = args.s1z_tmpl
    s2z_t = args.s2z_tmpl

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
    dimless = args.dimless
    plot_mode = plot_mode_parser(args.plot_mode)

    # Handlings...
    """
    1. Getting source parameters: m1, m2, s1z, s2z, gps
        1.1 If graceid & Sgraceid(preferred) was specified,
            will try loading source parameters from GraceDB.
        1.2 Using parameters of input m1, m2, s1z, s2z, gps(will cover GraceDB parameters)
    """
    # Parsing GraceDB
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
            LOGGER.error('Cannot fetch GraceDB event....\n')
    
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

    """
    2. Checking parameters
    """
    if m1 is None or \
        m2 is None or \
        s1z is None or \
        s2z is None or \
        gps is None:
        LOGGER.error('Parameters is not enough, exit...')
        return -1

    """
    3. Generating Template
        3.1 Checking m1_t, m2_t, s1z_t, s2z_t
        3.2 set fini
        3.3 set approx
        3.4 Generation
    """
    if m1_t is None:
        m1_t = m1
    if m2_t is None:
        m2_t = m2
    if s1z_t is None:
        s1z_t = s1z
    if s2z_t is None:
        s2z_t = s2z
    
    # Setting fini_SI
    import astropy.constants as cst
    if dimless:
        fini_SI = fini * cst.c.value**3 / ((m1 + m2) * cst.M_sun.value * cst.G.value)
    else:
        fini_SI = fini

    # Setting approx
    if approx is None:
        approx = get_proper_approx(m1, m2)

    logging.info(f'Parameters:\n\t\
                    m1 = {m1}, template m1 = {m1_t}\n\t\
                    m2 = {m2}, template m2 = {m2_t}\n\t\
                    s1z = {s1z}, template s1z = {s1z_t}\n\t\
                    s2z = {s2z}, template s2z = {s2z_t}\n\t\
                    template approx = {approx}\n\t\
                    template initial frequency = {fini_SI}\n\t\
                    gps end time: {gps}')
    
    # Making Template
    logging.info('Generating template...')
    tmpl = Template(m1 = m1, m2 = m2, s1z = s1z, s2z = s2z, 
                    fini = fini_SI, approx = approx, srate = fs,
                    duration = DEFAULT_MAX_DURATION)
    logging.info(f'Get template, duration = {tmpl.dtpeak}')
    track_x, track_y = tmpl.track
    trange_peak = [min(20/max(track_y), sback), min(16/max(track_y), sfwd)]

    """
    4. Checking data segment (gps - sfwd, gps + sback)
    """
    if gps_start is not None:
        sfwd = gps - gps_start
    if gps_end is not None:
        sback = gps_end - gps
    if sback < 0:
        LOGGER.warning(f'Invalid stepback, reset stepback = {DEFAULT_SBACK}\n')
        sback = DEFAULT_SBACK
    if sfwd < 0:
        LOGGER.warning(f'Invalid stepforward, reset stepforward = {DEFAULT_SFWD}\n')
        sfwd = DEFAULT_SFWD
    
    # Data segment should be as long as 2 times of the template~
    if sback < 2.1 * tmpl.dtpeak:
        LOGGER.warning(f'Time duration of template is too long.\n')
        sback = 2.1 * tmpl.dtpeak

    # Now setting data segment
    tstart = gps - sback
    tend = gps + sfwd

    # Create Coherent Object
    logging.info(f'Builing coherent strain, {tstart} ... {tend}')
    Strains = gwStrainCoherent(tstart, tend-tstart, fs = fs, verbose = True)
    if ifos is None:
        ifos = ['H1', 'L1', 'V1']

    """
    5. Loading PSD. If fail, exit
    """
    if refpsd is None:
        if Sgraceid is not None:
            pass
        LOGGER.error('reference-psd was not specified\n')
        return -1
    else:
        refpsd = Path(refpsd)
        if refpsd.is_dir():
            psddict = get_refpsd_from_dir(refpsd, channel = channel)
        else:
            psddict = get_refpsd(refpsd)
    

    """
    6. Loading data.
        6.1 If gaussian was set, will generate gaussian noise.
            else, loading data from cache or manually.
        6.2 If injection was set, will make an injection.
            6.2.1 Checking ra, de, psi, phic
            6.2.2 Calculating distance, using expected snr
            6.2.3 Calculating expected SNR and trackSNR for each detector
    """
    # Loading data & making noise
    if gaussian:
        logging.info(f'Making gaussian noise {ifos}')
        Strains.make_noise_from_psd(ifos, psddict)
    else:
        logging.info(f'Loading data {ifos}')
        Strains.load_data(cache = cache, ifos = ifos, channel = channel)
    logging.info('Setting psd')
    Strains.set_psd(psddict)
    # Shoule we make injection?
    if injection:
        # Checking ra, de
        if ra is None:
            LOGGER.warning('You have not set ra, now we set it to 0\n')
            ra = 0
        if de is None:
            LOGGER.warning('You have not set de, now we set it to 0\n')
            de = 0
        if gps_inj is None:
            LOGGER.warning(f'You have not set gps-injection, now we set it to {gps}\n')
            gps_inj = gps
        if psi is None:
            psi = 0
        if phic is None:
            phic = 0
        regeneration = False
        if m1_inj is None:
            m1_inj = m1_t
        else:
            regeneration = True
        if m2_inj is None:
            m2_inj = m2_t
        else:
            regeneration = True
        if s1z_inj is None:
            s1z_inj = s1z_t
        else:
            regeneration = True
        if s2z_inj is None:
            s2z_inj = s2z_t
        else:
            regeneration = True

        if regeneration:
            logging.info('Regenerating template for injection...')
            if approx_inj is None:
                approx_inj = get_proper_approx(m1_inj, m2_inj)
            tmpl_inj = Template(m1 = m1_inj, m2 = m2_inj, s1z = s1z_inj, s2z = s2z_inj, 
                                fini = fini_SI, approx = approx_inj, srate = fs,
                                duration = None)
        else:
            tmpl_inj = tmpl
        logging.info(f'Making injection... expected coherent snr = {snr_expected}')
        expected_snr_dict, rescaled = Strains.make_injection(tmpl_inj, tmpl, gps_inj, ra, de, snr_expected, psi = psi, phic = phic)
        logging.info(f'Injection done, expected SNR')
        for ifo in expected_snr_dict:
            sys.stderr.write(f'\t{ifo}: {expected_snr_dict[ifo]}\n')
        logging.info('Calculating expected track SNR')
        exp_freqs, exp_trackSNR = Strains.calc_expected_track_SNR(q = Q, tmpl_inj = tmpl_inj,
                                        tmpl = tmpl, gps = gps, ra_inj = ra, de_inj = de,
                                        rescaled = rescaled,
                                        psi = psi, phic = phic,
                                        frange = frange, mismatch = mismatch)

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

    """
    Step.1 Matched filtering & Skymap
    """
    # Call...
    logging.info('Matched filtering & Coherent...')
    SNRs, skymap = \
        Strains.calc_coherent_snr_skymap(tmpl, nside, gps)

    """
    Step.2 Plot SNR & Skymap...
    """
    logging.info('Ploting coherent SNR skymap & SNR time series...')
    skymap.plot_skymap(fsave, plot_peak = True)
    skymap.plot_coherent_snr(fsave / f'Coherent_SNR.png')
    for snr in SNRs:
        snr.plot(fsave = fsave / f'SNR_{snr.ifo}.png', 
            pset = 'abs')
        snr.plot(fsave = fsave / f'SNR_{snr.ifo}_zoom.png',
            pset = 'abs', xrange = [gps - trange_peak[0], gps + trange_peak[1]])

    gps_max = skymap.max_gps_time
    logging.info(f'Coherent SNR peak at geocent time {gps_max}, while gps input is {gps}')
    return -1
    """
    Step.3 Qscan...
    """
    logging.info('Q matched filtering & Coherent...')
    trange_duration = [min(max(tmpl.duration*1.2, trange_peak[0]), sback*0.8), min(1, sfwd*0.8)]
    max_ra,max_de = skymap.max_ra_de
    SPECs, cohSPEC, nullSPEC = \
        Strains.calc_coherent_snr_qspectrum(tmpl, q = Q, 
            gps_trigger = gps, ra = max_ra, de = max_de, 
            trange = trange_duration,
            frange = frange)

    """
    Step. Ploting coherent Q spectrum...
    """
    logging.info('Ploting coherent Q spectrum...')
    tspec_plot = np.arange(gps - trange_duration[0], gps + trange_duration[1], 1./fs)
    tdur_range_plot = [tspec_plot[0], tspec_plot[-1]]
    tpeak_range_plot = [gps - trange_peak[0], gps + trange_peak[1]]
    fspec_plot = np.logspace(np.log10(frange[0]), np.log10(frange[1]), 500)
    flabel = f'frequency [Hz]{frange}'
    logging.info('Coherent...')
    cohSPEC.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                          figsize = FIGSIZE_QSCAN, fsave = fsave/'snrQscan_coh.png',
                          cmaptype = cmaptype, pcolorbins = pcolorbins,
                          ylabel = flabel,
                          xlim = tdur_range_plot, ylim = frange)

    cohSPEC.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                          figsize = FIGSIZE_QSCAN, fsave = fsave/'snrQscan_coh_zoom.png',
                          cmaptype = cmaptype, pcolorbins = pcolorbins,
                          ylabel = flabel,
                          xlim = tpeak_range_plot, ylim = frange)

    cohSPEC.plot_spectrum_with_track(tmpl, gps_max, 
                            figsize = FIGSIZE_QSCAN, fsave = fsave/'snrQscan_coh_track.png',
                            cmaptype = cmaptype, pcolorbins = pcolorbins,
                            ylabel = flabel)
    if nullSPEC is not None:
        logging.info('Null...')
        nullSPEC.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                            figsize = FIGSIZE_QSCAN, fsave = fsave/'nullQscan_coh.png',
                            cmaptype = cmaptype, pcolorbins = pcolorbins,
                            ylabel = flabel,
                            xlim = tdur_range_plot, ylim = frange)

        nullSPEC.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                            figsize = FIGSIZE_QSCAN, fsave = fsave/'nullQscan_coh_zoom.png',
                            cmaptype = cmaptype, pcolorbins = pcolorbins,
                            ylabel = flabel,
                            xlim = tpeak_range_plot, ylim = frange)
    
    for spec in SPECs:
        logging.info(f'{spec.ifo}...')
        spec.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                            figsize = FIGSIZE_QSCAN, fsave = fsave/f'Qscan_{spec.ifo}.png',
                            cmaptype = cmaptype, pcolorbins = pcolorbins,
                            ylabel = flabel,
                            xlim = tdur_range_plot, ylim = frange)

        spec.plot_spectrum(times = tspec_plot, freqs = fspec_plot,
                            figsize = FIGSIZE_QSCAN, fsave = fsave/f'Qscan_{spec.ifo}_zoom.png',
                            cmaptype = cmaptype, pcolorbins = pcolorbins,
                            ylabel = flabel,
                            xlim = tpeak_range_plot, ylim = frange)
        delay = spec.ifo_delay(max_ra, max_de, gps)
        spec.plot_spectrum_with_track(tmpl, gps_max + delay, 
                                figsize = FIGSIZE_QSCAN, fsave = fsave/f'Qscan_{spec.ifo}_track.png',
                                cmaptype = cmaptype, pcolorbins = pcolorbins,
                                ylabel = flabel)
    logging.info('Calculating track significance...')
    traceSNR, freqs, backSNR = cohSPEC.calc_trace(tmpl, gps_max)
    traceSNR_int = np.average(traceSNR)

    iter_back = 10
    backSNR_int = []
    for i in range(iter_back):
        logging.info('Calculating background...')
        gps_back = gps_max - 100*i - np.random.random() * 50
        backStrains = gwStrainCoherent(gps_back - sback, sback+sfwd, fs = fs, verbose = False)
        backStrains.load_data(cache = cache, ifos = ifos, channel = channel)
        backSNRs, back_skymap = \
            backStrains.calc_coherent_snr_skymap(tmpl, nside, gps_back)
        max_ra_back, max_de_back = back_skymap.max_ra_de
        back_SPECs, back_cohSPEC, back_nullSPEC = \
            backStrains.calc_coherent_snr_qspectrum(tmpl, q = Q, 
                    gps_trigger = gps_back, ra = max_ra_back, de = max_de_back, 
                    trange = trange_duration,
                    frange = frange)
        backSNR_int += back_cohSPEC.calc_background_track(tmpl)

    backtraceSNR = np.zeros(len(traceSNR))
    count = 0
    for back in backSNR:
        if len(back) == len(traceSNR):
            backtraceSNR = backtraceSNR + back
            count += 1
        backSNR_int.append(np.average(back))
    backtraceSNR = backtraceSNR / count
    plt.plot(freqs, traceSNR, label = 'trace')
    if injection:
        plt.plot(exp_freqs, exp_trackSNR, label = 'expected')
        plt.legend()
    plt.xlabel('frequency [Hz]')
    plt.ylabel('SNR')
    plt.savefig(fsave/'traceSNR.png', dpi = 200)
    plt.close()

    backSNR_int = np.asarray(backSNR_int)
    count_x = np.arange(len(backSNR_int))
    plt.figure(figsize = (8,4))
    plt.scatter(count_x, backSNR_int, marker = 'x', color = 'gray', label = 'background')
    plt.scatter([len(backSNR_int)/2], [traceSNR_int], marker = '.', color = 'red',label = 'foreground')
    plt.legend()
    plt.xticks([])
    plt.ylabel('SNR')
    plt.savefig(fsave/'significance.png', dpi = 200)
    plt.close()
    LOGGER.info(f'Trace SNR = {traceSNR_int}\n')
    LOGGER.info(f'Average background SNR = {np.average(backSNR_int)}\n')

    return 0



class plot_mode_parser(object):
    def __init__(self, plot_mode):
        plot_mode = str(plot_mode).lower()
        if plot_mode in ('all', 'default','every','plot_all', 'plotall',):
            self._all = True
            self._off = True
        elif plot_mode in ('no','noplot','no_plot','nothing','off','mute',):
            self._all = False
            self._off = True
        else:
            self._all = False
            self._off = False
        self._mode = plot_mode
        self._split = plot_mode.split('_')
        if self._split[0] in ('no', 'dont', 'off', 'kill'):
            self._judge = False
        else:
            self._judge = True
    
    @property
    def mode(self):
        return self._mode

    def __call__(self, option):
        if self._all:
            return True
        if self._off:
            return False
        if option in self._split:
            return self._judge
        else:
            return not self._judge