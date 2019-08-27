"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
from .._core import resample
from scipy.signal import resample as scipy_resample
from scipy.interpolate import interp1d, interp2d
from .._utils import interp2d_complex

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec


class Series(object):
    def __init__(self, value, deltax, info = 'Series'):
        value = np.asarray(value)
        if (len(value.shape) != 1 and value.shape[0] != 1):
            raise Exception(f'Shape error: {value.shape}')
        self._value = value
        self._deltax = deltax
        self._info = info
    
    @property
    def value(self):
        return self._value
    
    @property
    def deltax(self):
        return self._deltax

    def __len__(self):
        return len(self._value)

    def __abs__(self):
        return abs(self._value)

    @property
    def size(self):
        return self._value.size
    
    @property
    def x(self):
        return np.arange(0, self.length, self._deltax)

    @property
    def length(self):
        return self._deltax * len(self)

    @property
    def real(self):
        return Series(self._value.real, self._deltax, info = f'Re_{self._info}')

    @property
    def imag(self):
        return Series(self._value.imag, self._deltax, info = f'Im_{self._info}')

    def conjugate(self):
        return Series(self._value.conjugate(), self._deltax, info = f'Conj_{self._info}')

    def __str__(self):
        return f'{self._info}: {self._value}'
    
    def __repr__(self):
        return self.__str__()
    
    def __format__(self):
        return self.__str__()
    
    def __iter__(self):
        for x in self._value:
            yield x

    def __getitem__(self, key):
        if isinstance(key, np.int):
            return self._value[key]
        return self._getslice(key)

    def _getslice(self, index):
        if isinstance(index, slice):        
            if index.start is not None and index.start < 0:
                raise ValueError(('Negative start index ({}) is not supported').format(index.start))
            
            if index.step is not None:
                new_deltax = self.deltax * index.step
            else:
                new_deltax = self.deltax
            return Series(self._value[index], deltax = new_deltax)
        if isinstance(index, np.ndarray):
            if len(index) == 1:
                return self._value[index]
            # Check uniform
            grad = np.gradient(index)
            if max(grad) != min(grad):
                raise ValueError(f'Invalid index for Series: {index}')
            step = index[1] - index[0]
            new_deltax = self._deltax * step
            return Series(self._value[index], deltax = new_deltax)

    
    def __setitem__(self, key, val):
        self._value[key] = val

    def resample(self, new_deltax):
        if new_deltax != self.deltax:
            new = resample(self.value, 1./self.deltax, 1./new_deltax)
            return Series(new, new_deltax, info = self._info)
        else:
            return self


class TimeSeries(Series):
    def __init__(self, value, epoch, fs, info = 'TimeSeries'):
        super(TimeSeries, self).__init__(value, 1./fs, info = info)
        self._epoch = epoch
    
    @property
    def fs(self):
        return int(1./self._deltax)
    
    @property
    def time(self):
        return self._epoch + self.x
    
    @property
    def duration(self):
        return self.length
    
    @property
    def epoch(self):
        return self._epoch
    
    def resample(self, fs_new):
        if fs_new != self.fs:
            new = resample(self.value, self.fs, fs_new)
            return TimeSeries(new, epoch=self.epoch, fs=self.fs, info=self.info)
        else:
            return self

    def plot(self, fsave,
             xrange = None, yrange = None,
             xlabel = None, ylabel = None,
             figsize = None, pset = None,
             title = None):
        if figsize is None:
            figsize = (10, 5)
        if pset in (None, 'origin',):
            val = self.value
        if pset in ('abs', 'snr'):
            val = np.abs(self.value)
        if title is None:
            title = self._info
        plt.figure(figsize = figsize)
        plt.plot(self.time, val)
        plt.xlim(xrange)
        plt.ylim(yrange)
        plt.title(title)
        plt.savefig(fsave, dpi = 200)
        plt.close()

    

class MultiSeries(object):
    def __init__(self, array, deltax, y):
        array = np.asarray(array)
        if len(array.shape) == 1:
            if array.shape[0] > 0:
                array = array.reshape(1, array.size)
                self._isempty = False
            else:
                array = np.array([])
                self._isempty = True
        else:
            self._isempty = False
        self._array = array
        if not self._isempty:
            self._deltax = deltax
            if isinstance(y, np.int) or isinstance(y, np.float):
                y = [y]
            y = np.asarray(y)
            if len(y) > 0:
                if (len(y.shape) != 1 and y.shape[0] != 1):
                    raise Exception(f'Shape error for y: {y.shape}')
                if y.size != self._array.shape[0]:
                    raise Exception(f'Incompatible size for y: {y.size}')
                self._y = y.reshape(y.size)
            else:
                raise Exception(f'Invalid variable: {y}')
        else:
            self._deltax = None
            self._y = None

    @property
    def y(self):
        return self._y

    @property
    def deltax(self):
        return self._deltax

    @property
    def x(self):
        return np.arange(0, self.length, self._deltax)

    def __len__(self):
        return self.shape[1]

    @property
    def ysize(self):
        return self.shape[0]

    @property
    def length(self):
        return self.xsize * self._deltax
    
    @property
    def height(self):
        return self._y[-1] - self._y[0]

    @property
    def xsize(self):
        return self.shape[1]

    @property
    def shape(self):
        return self._array.shape

    def __iter__(self):
        for i in range(self.ysize):
            yield (self.y[i], Series(self._array[i,:], self.deltax))

    def append(self, series, y):
        if not isinstance(series, Series):
            series = Series(series, self.deltax)
        if not self._isempty:
            if len(series) != self.xsize:
                raise Exception(f'Incompatible size: {series.size} != {self.xsize}')
            if series.deltax != self.deltax:
                raise Exception(f'Incompatible deltax: {series.deltax} != {self.deltax}')
            if y > self._y[-1]:
                idx_insert = self.ysize
                self._array = np.insert(self._array, idx_insert, series.value, axis=0)
                self._y = np.insert(self._y, idx_insert, y)
            else:
                idx_insert = np.where(self._y - y >= 0)[0][0]
                if self._y[idx_insert] == y:
                    self._array[idx_insert,:] = series.value
                else:
                    self._array = np.insert(self._array, idx_insert, series.value, axis=0)
                    self._y = np.insert(self._y, idx_insert, y)
        else:
            size = series.size
            self._deltax = series.deltax
            self._array = series.value.reshape(1, size)
            if not isinstance(y, np.int) or not isinstance(y, np.float):
                raise TypeError(f'Invalid type: {type(y)}')
            self._y = np.array([y])



    
class TimeFreqSpectrum(MultiSeries):
    def __init__(self, array, epoch, fs, freqs, info = 'TimeFreqSpectrum'):
        super(TimeFreqSpectrum, self).__init__(array, 1./fs, freqs)
        self._info = info
        if not self._isempty:
            if isinstance(epoch, np.int) or isinstance(epoch, np.float):
                epoch = [epoch]
            epoch = np.asarray(epoch)
            if len(epoch) == 1:
                self._epoch = np.ones(self.ysize) * epoch[0]
            elif len(epoch) == self.ysize:
                self._epoch = epoch
            else:
                raise Exception(f'Incompatible shape for epoch: {epoch.shape}')
        else:
            epoch = None
    @property
    def epoch(self):
        return self._epoch

    @property
    def trange(self):
        epoch_min = min(self.epoch)
        epoch_max = max(self.epoch)
        return epoch_max, epoch_min + self.length

    @property
    def times(self):
        return np.arange(self.trange[0], self.trange[1], self._deltax)
    
    @property
    def fs(self):
        return 1./self.deltax
    
    @property
    def frequencies(self):
        return self.y

    def __iter__(self):
        for i in range(self.ysize):
            yield (self.frequencies[i], TimeSeries(self._array[i,:], self.epoch[i], self.fs, info = self._info))
        
    def append(self, timeseries, freq, epoch = None, fs = None):
        if not isinstance(timeseries, TimeSeries) and epoch is None:
            raise TypeError(f'Invalid type: {timeseries}')
        elif epoch is not None and isinstance(timeseries, np.ndarray):
            value = timeseries
            if fs is None:
                deltax = self.deltax
            else:
                deltax = 1./fs
            size = value.size
        else:
            value = timeseries.value
            epoch = timeseries.epoch
            deltax = timeseries.deltax
            size = timeseries.size
            

        if not self._isempty:
            if size != self.xsize:
                raise Exception(f'Incompatible size: {timeseries.size} != {self.xsize}')
            if deltax != self.deltax:
                raise Exception(f'Incompatible deltax: {timeseries.deltax} != {self.deltax}')
            if freq > self._y[-1]:
                idx_insert = self.ysize
                self._array = np.insert(self._array, idx_insert, value, axis=0)
                self._epoch = np.insert(self._epoch, idx_insert, epoch)
                self._y = np.insert(self._y, idx_insert, freq)
            else:
                idx_insert = np.where(self._y - freq >= 0)[0][0]
                if self._y[idx_insert] == freq:
                    self._array[idx_insert, :] = value
                    self._epoch[idx_insert] = epoch
                else:
                    self._array = np.insert(self._array, idx_insert, value, axis=0)
                    self._epoch = np.insert(self._epoch, idx_insert, epoch)
                    self._y = np.insert(self._y, idx_insert, freq)
        else:
            self._array = value.reshape(1, size)
            self._epoch = np.array([epoch])
            self._y = np.array([freq])
            self._deltax = deltax
            self._isempty = False

    def interpolate(self, t_interp):
        ret = np.zeros([self.ysize, len(t_interp)], self._array.dtype)
        for i, epoch in enumerate(self.epoch):
            ret[i,:] = np.interp(t_interp, self.x + epoch, self._array[i,:])
        return ret

    def get_finterp(self, pset = None):
        xp = self.times
        yp = self.frequencies
        zp = self.interpolate(xp)
        if pset in ('abs',):
            zp = np.abs(zp)
        if zp.dtype == complex:
            return interp2d_complex(xp, yp, zp)
        else:
            return interp2d(xp, yp, zp)

    def plot_spectrum(self, times, freqs,
                      fsave,
                      figsize = None, 
                      cmaptype = 'jet', pcolorbins = 100,
                      xlabel = None, ylabel = None,
                      xlim = None, ylim = None,
                      yticks = None,
                      title = None):
        # plot setting
        if figsize is None:
            figsize = (12, 7)
        cmap = plt.get_cmap(cmaptype)
        if ylim is None:
            ylim = [self.frequencies[0], self.frequencies[-1]]
        if yticks is None:
            yticksval = np.logspace(np.log10(ylim[0]), np.log10(ylim[1]), 5)
            yticks = (yticksval, ['%.1f'%_freq for _freq in yticksval])
        if title is None:
            title = self._info

        x = times
        y = freqs
        z = self.get_finterp(pset = 'abs')(x,y)
        if xlabel is None:
            idx_tpeak_0, idx_fpeak_0 = get_2D_argpeak(z)
            tpeak = '%.2f'%x[idx_tpeak_0]
            fpeak = '%.1f'%y[idx_fpeak_0]
            snrpeak = '%.3f'%z[idx_fpeak_0, idx_tpeak_0]
            xlabel = f'loudest snr = {snrpeak}, at geocent gps = {tpeak}, f = {fpeak}'

        levels = MaxNLocator(nbins=pcolorbins).tick_values(z.min(), z.max())
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        fig = plt.figure(figsize = figsize)
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(x, y, z, cmap = cmap, norm = norm)
        fig.colorbar(im, ax=ax)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.yscale('log')
        if isinstance(yticks, tuple):
            plt.yticks(*yticks)
        plt.savefig(fsave ,dpi = 200)
        plt.close()

    def plot_spectrum_with_track(self,
                                 tmpl, gps_trigger, fsave,
                                 figsize = None, 
                                 cmaptype = 'jet', pcolorbins = 100,
                                 xlabel = None, ylabel = None,
                                 yticks = None,
                                 title = None):
        # Track
        track_x, track_y = tmpl.get_track(gps_trigger)
        ntrack = len(track_x)
        if ntrack > 1e4:
            fs_plot = int(self.fs * (1e4 / ntrack))
            track_x = resample(track_x, tmpl.fs, fs_plot)
            track_y = resample(track_y, tmpl.fs, fs_plot)
        else:
            fs_plot = self.fs
        # plot setting
        if figsize is None:
            figsize = (12, 14)
        else:
            figsize = (figsize[0], 2*figsize[1])
        cmap = plt.get_cmap(cmaptype)
        ylim = (self.frequencies[0], self.frequencies[-1])
        xlim = (max(track_x[0] - 0.5, self.trange[0]), min(track_x[1] + 0.5, self.trange[1]))
        
        if yticks is None:
            yticksval = np.logspace(np.log10(ylim[0]), np.log10(ylim[1]), 5)
            yticks = (yticksval, ['%.1f'%_freq for _freq in yticksval])
        
        if title is None:
            title = self._info

        x = np.arange(xlim[0], xlim[1], 1./fs_plot)
        y = np.logspace(np.log10(ylim[0]), np.log10(ylim[1]), 500)
        func = self.get_finterp(pset = 'abs')
        z = func(x,y)
        track_spec = func(track_x, track_y)
        track_trace = track_spec[np.arange(ntrack), np.arange(ntrack)]
        
        if self.size > 1e4:
            time_back = scipy_resample(self.times, 1e4)
        else:
            time_back = self.times
        track_back_spec = func(time_back, track_y)
        track_back = np.median(track_back_spec, axis = 1)
        significance = np.sum(track_trace) / np.median(z) / ntrack
        
        if xlabel is None:
            xlabel = f'track significance = {significance}'

        levels = MaxNLocator(nbins=pcolorbins).tick_values(z.min(), z.max())
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        fig = plt.figure(figsize = figsize)
        gs = gridspec.GridSpec(10, 1)
        gs.update(left=0.05, right=0.95, wspace=0.02)
        plt.title(title)
        ax1 = plt.subplot(gs[:6,0])
        ax2 = plt.subplot(gs[7:,0])
        im = ax1.pcolormesh(x, y, z, cmap = cmap, norm = norm)
        fig.colorbar(im, ax=ax1)
        ax1.plot(track_x, track_y, '-', color='#ba7b00', zorder=3, lw=1.5)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        ax1.set_ylim(ylim)
        ax1.set_xlim(xlim)
        ax1.set_yscale('log')
        if isinstance(yticks, tuple):
            ax1.set_yticks(*yticks)
        # Plot track evolution
        ax2.plot(track_x, track_trace, label = 'track SNR')
        ax2.plot(track_x, track_back, label = 'Background median SNR')
        ax2.set_xlabel(f'gps [s]')
        ax2.set_ylabel(f'SNR')
        ax2.legend()
        plt.savefig(fsave ,dpi = 200)
        plt.close()
        return significance


    def calc_integrate_track(self, track_x, track_y,
                    dt_idx = 10):
        dt_idx = int(dt_idx)
        integ = np.zeros(len(self.frequencies))
        freqgrad = np.gradient(self.frequencies)
        time_split = np.zeros(len(self.frequencies))
        for i, freq in enumerate(self.frequencies):
            deltafreq = track_y - freq
            idx = np.where( np.abs(deltafreq) == np.min(np.abs(deltafreq)) )[0][0]
            gps = track_x[idx]
            deltatime = self.x + self.epoch[i] - gps
            dtrange = dt_idx / freq
            idx_times = np.where( np.abs(deltatime) < dtrange )[0]
            integ[i] = np.max(np.abs(self._array[i, idx_times]))
            time_split[i] = gps
        ret = np.sum(integ * np.gradient(time_split) * freqgrad)
        return ret / (time_split[-1] - time_split[0]) / (self.frequencies[-1] - self.frequencies[0])
        

def get_2D_argpeak(matrix):
    arg = np.where(matrix == np.max(matrix))
    return arg[1][0], arg[0][0]

