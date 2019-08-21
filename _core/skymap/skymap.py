#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:24:24 2019

@author: drizl
"""

import numpy as np
import matplotlib
import warnings
import six
from .pix import vec2pix

max_nside = 1 << 29
UNSEEN = 999 # CPP
pi = np.pi
dtor = pi / 180.0


# def pix2vec(nside, ipix, nest=False):
#     # CPP
#     pass



def isnpixok(npix):
    nside = np.sqrt(np.asarray(npix) / 12.0)
    return nside == np.floor(nside)

def isnsideok(nside, nest=False):
    if hasattr(nside, "__len__"):
        if not isinstance(nside, np.ndarray):
            nside = np.asarray(nside)
        is_nside_ok = (
            (nside == nside.astype(np.int)) & (nside > 0) & (nside <= max_nside)
        )
        if nest:
            is_nside_ok &= (nside.astype(np.int) & (nside.astype(np.int) - 1)) == 0
    else:
        is_nside_ok = nside == int(nside) and 0 < nside <= max_nside
        if nest:
            is_nside_ok = is_nside_ok and (int(nside) & (int(nside) - 1)) == 0
    return is_nside_ok

def nside2npix(nside):
    return 12 * nside * nside

def get_min_valid_nside(npix):
    order = 0.5 * np.log2(npix / 12.0)
    return 1 << int(np.ceil(order))


def ma_to_array(m):
    try:
        return m.filled()
    except AttributeError:
        try:
            return np.array([mm.filled() for mm in m])
        except AttributeError:
            pass
    return m

def get_map_size(m):
    if isinstance(m, dict):
        if "nside" in m:
            return nside2npix(m["nside"])
        elif hasattr(m, "nside"):
            return nside2npix(m.nside)
        else:
            nside = get_min_valid_nside(max(m.keys()) + 1)
            return nside2npix(nside)
    else:
        if isnpixok(len(m)):
            return len(m)
        else:
            raise ValueError("Wrong pixel number (it is not 12*nside**2)")

def mask_bad(m, badval=UNSEEN, rtol=1.0e-5, atol=1.0e-8):
    m = np.asarray(m)
    atol = np.absolute(atol)
    rtol = np.absolute(rtol)
    return np.absolute(m - badval) <= atol + rtol * np.absolute(badval)


def maptype(m):
    if not hasattr(m, "__len__"):
        raise TypeError("input map is a scalar")
    if len(m) == 0:
        raise TypeError("input map has length zero")

    try:
        npix = len(m[0])
    except TypeError:
        npix = None

    if npix is not None:
        for mm in m[1:]:
            if len(mm) != npix:
                raise TypeError("input maps have different npix")
        if isnpixok(len(m[0])):
            return len(m)
        else:
            raise TypeError("bad number of pixels")
    else:
        if isnpixok(len(m)):
            return 0
        else:
            raise TypeError("bad number of pixels")


def npix2nside(npix):
    nside = np.sqrt(np.asarray(npix) / 12.0)
    if nside != np.floor(nside):
        raise ValueError('Incorrect npix')
    return int(np.sqrt(npix / 12.0))

def get_nside(m):
    typ = maptype(m)
    if typ == 0:
        return npix2nside(len(m))
    else:
        return npix2nside(len(m[0]))

def check_nside(nside, nest=False):
    """Raises exception is nside is not valid"""
    if not np.all(isnsideok(nside, nest=nest)):
        raise ValueError(
            "%s is not a valid nside parameter (must be a power of 2, less than 2**30)"
            % str(nside)
        )

# def is_ma(m):
#     return hasattr(m, "filled") or hasattr(m[0], "filled")

# def fit_dipole(m, nest=False, bad=UNSEEN, gal_cut=0):
#     m = ma_to_array(m)
#     m = np.asarray(m)
#     npix = m.size
#     nside = npix2nside(npix)
#     if nside > 128:
#         bunchsize = npix // 24
#     else:
#         bunchsize = npix
#     aa = np.zeros((4, 4), dtype=np.float64)
#     v = np.zeros(4, dtype=np.float64)
#     for ibunch in range(npix // bunchsize):
#         ipix = np.arange(ibunch * bunchsize, (ibunch + 1) * bunchsize)
#         ipix = ipix[(m.flat[ipix] != bad) & (np.isfinite(m.flat[ipix]))]
#         x, y, z = pix2vec(nside, ipix, nest)
#         if gal_cut > 0:
#             w = np.abs(z) >= np.sin(gal_cut * np.pi / 180)
#             ipix = ipix[w]
#             x = x[w]
#             y = y[w]
#             z = z[w]
#             del w
#         aa[0, 0] += ipix.size
#         aa[1, 0] += x.sum()
#         aa[2, 0] += y.sum()
#         aa[3, 0] += z.sum()
#         aa[1, 1] += (x ** 2).sum()
#         aa[2, 1] += (x * y).sum()
#         aa[3, 1] += (x * z).sum()
#         aa[2, 2] += (y ** 2).sum()
#         aa[3, 2] += (y * z).sum()
#         aa[3, 3] += (z ** 2).sum()
#         v[0] += m.flat[ipix].sum()
#         v[1] += (m.flat[ipix] * x).sum()
#         v[2] += (m.flat[ipix] * y).sum()
#         v[3] += (m.flat[ipix] * z).sum()
#     aa[0, 1] = aa[1, 0]
#     aa[0, 2] = aa[2, 0]
#     aa[0, 3] = aa[3, 0]
#     aa[1, 2] = aa[2, 1]
#     aa[1, 3] = aa[3, 1]
#     aa[2, 3] = aa[3, 2]
#     res = np.dot(np.linalg.inv(aa), v)
#     mono = res[0]
#     dipole = res[1:4]
#     return mono, dipole

# def ma(m, badval=UNSEEN, rtol=1e-5, atol=1e-8, copy=True):
#     return np.ma.masked_values(np.array(m), badval, rtol=rtol, atol=atol, copy=copy)


# def remove_dipole(
#     m, nest=False, bad=UNSEEN, gal_cut=0, fitval=False, copy=True, verbose=True
# ):
#     input_ma = is_ma(m)
#     m = ma_to_array(m)
#     m = np.array(m, copy=copy)
#     npix = m.size
#     nside = npix2nside(npix)
#     if nside > 128:
#         bunchsize = npix // 24
#     else:
#         bunchsize = npix
#     mono, dipole = fit_dipole(m, nest=nest, bad=bad, gal_cut=gal_cut)
#     for ibunch in range(npix // bunchsize):
#         ipix = np.arange(ibunch * bunchsize, (ibunch + 1) * bunchsize)
#         ipix = ipix[(m.flat[ipix] != bad) & (np.isfinite(m.flat[ipix]))]
#         x, y, z = pix2vec(nside, ipix, nest)
#         m.flat[ipix] -= dipole[0] * x
#         m.flat[ipix] -= dipole[1] * y
#         m.flat[ipix] -= dipole[2] * z
#         m.flat[ipix] -= mono
#     if verbose:
#         from . import rotator as R

#         lon, lat = R.vec2dir(dipole, lonlat=True)
#         amp = np.sqrt((dipole * dipole).sum())
#         print(
#             "monopole: {0:g}  dipole: lon: {1:g}, lat: {2:g}, amp: {3:g}".format(
#                 mono, lon, lat, amp
#             )
#         )
#     if is_ma:
#         m = ma(m)
#     if fitval:
#         return m, mono, dipole
#     else:
#         return m

# def fit_monopole(m, nest=False, bad=UNSEEN, gal_cut=0):
#     m = ma_to_array(m)
#     m = np.asarray(m)
#     npix = m.size
#     nside = npix2nside(npix)
#     if nside > 128:
#         bunchsize = npix // 24
#     else:
#         bunchsize = npix
#     aa = v = 0.0
#     for ibunch in range(npix // bunchsize):
#         ipix = np.arange(ibunch * bunchsize, (ibunch + 1) * bunchsize)
#         ipix = ipix[(m.flat[ipix] != bad) & (np.isfinite(m.flat[ipix]))]
#         x, y, z = pix2vec(nside, ipix, nest)
#         if gal_cut > 0:
#             w = np.abs(z) >= np.sin(gal_cut * np.pi / 180)
#             ipix = ipix[w]
#             x = x[w]
#             y = y[w]
#             z = z[w]
#             del w
#         aa += ipix.size
#         v += m.flat[ipix].sum()
#     mono = v / aa
#     return mono


# def remove_monopole(
#     m, nest=False, bad=UNSEEN, gal_cut=0, fitval=False, copy=True, verbose=True
# ):
#     input_ma = is_ma(m)
#     m = ma_to_array(m)
#     m = np.array(m, copy=copy)
#     npix = m.size
#     nside = npix2nside(npix)
#     if nside > 128:
#         bunchsize = npix // 24
#     else:
#         bunchsize = npix
#     mono = fit_monopole(m, nest=nest, bad=bad, gal_cut=gal_cut)
#     for ibunch in range(npix // bunchsize):
#         ipix = np.arange(ibunch * bunchsize, (ibunch + 1) * bunchsize)
#         ipix = ipix[(m.flat[ipix] != bad) & (np.isfinite(m.flat[ipix]))]
#         x, y, z = pix2vec(nside, ipix, nest)
#         m.flat[ipix] -= mono
#     if verbose:
#         print("monopole: {0:g}".format(mono))
#     if input_ma:
#         m = ma(m)
#     if fitval:
#         return m, mono
#     else:
#         return m

def lonlat2thetaphi(lon, lat):
    return np.pi / 2.0 - np.radians(lat), np.radians(lon)


# def get_interp_val(m, theta, phi, nest=False, lonlat=False):
#     m2 = m.ravel()
#     nside = npix2nside(m2.size)
#     if lonlat:
#         theta, phi = lonlat2thetaphi(theta, phi)
#     if nest:
#         r = _get_interpol_nest(nside, theta, phi)
#     else:
#         r = _get_interpol_ring(nside, theta, phi)
#     p = np.array(r[0:4])
#     w = np.array(r[4:8])
#     del r
#     return np.sum(m2[p] * w, 0)


def get_color_table(vmin, vmax, val, cmap=None, norm=None):
    # Create color table
    newcmap = create_colormap(cmap)
    if type(norm) is str:
        if norm.lower().startswith("log"):
            norm = LogNorm2(clip=False)
        elif norm.lower().startswith("hist"):
            norm = HistEqNorm(clip=False)
        else:
            norm = None
    if norm is None:
        norm = LinNorm2(clip=False)

    norm.vmin = vmin
    norm.vmax = vmax
    norm.autoscale_None(val)

    return newcmap, norm

def create_colormap(cmap):
    if type(cmap) == str:
        cmap0 = matplotlib.cm.get_cmap(cmap)
    elif type(cmap) in [
        matplotlib.colors.LinearSegmentedColormap,
        matplotlib.colors.ListedColormap,
    ]:
        cmap0 = cmap
    else:
        cmap0 = matplotlib.cm.get_cmap(matplotlib.rcParams["image.cmap"])
    if hasattr(cmap0, "_segmentdata"):
        newcm = matplotlib.colors.LinearSegmentedColormap(
            "newcm", cmap0._segmentdata, cmap0.N
        )
    else:
        newcm = cmap0
    newcm.set_over(newcm(1.0))
    newcm.set_under("w")
    newcm.set_bad("gray")
    return newcm

class LogNorm2(matplotlib.colors.Normalize):
    """
    Normalize a given value to the 0-1 range on a log scale
    """

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = "array"
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = "scalar"
            val = np.ma.array([value]).astype(np.float)

        val = np.ma.masked_where(np.isinf(val.data), val)

        self.autoscale_None(val)
        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin <= 0:
            raise ValueError("values must all be positive")
        elif vmin == vmax:
            return type(value)(0.0 * np.asarray(value))
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax), mask=mask)
            result = (np.ma.log(val) - np.log(vmin)) / (np.log(vmax) - np.log(vmin))
            result.data[result.data < 0] = 0.0
            result.data[result.data > 1] = 1.0
            result[np.isinf(val.data)] = -np.inf
            if result.mask is not np.ma.nomask:
                result.mask[np.isinf(val.data)] = False
        if vtype == "scalar":
            result = result[0]
        return result

    def autoscale_None(self, A):
        " autoscale only None-valued vmin or vmax"
        if self.vmin is None or self.vmax is None:
            val = np.ma.masked_where(np.isinf(A.data), A)
            matplotlib.colors.Normalize.autoscale_None(self, val)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = float(self.vmin), float(self.vmax)

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin * np.ma.power((vmax / vmin), val)
        else:
            return vmin * np.pow((vmax / vmin), value)

class HistEqNorm(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, clip=False):
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)
        self.xval = None
        self.yval = None

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = "array"
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = "scalar"
            val = np.ma.array([value]).astype(np.float)

        self.autoscale_None(val)

        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            return 0.0 * val
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax), mask=mask)
            result = np.ma.array(
                np.interp(val, self.xval, self.yval), mask=np.ma.getmask(val)
            )
            result[np.isinf(val.data)] = -np.inf
        if vtype == "scalar":
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        if matplotlib.cbook.iterable(value):
            vtype = "array"
            val = np.ma.array(value)
        else:
            vtype = "scalar"
            val = np.ma.array([value])
        result = np.ma.array(
            self._lininterp(val, self.yval, self.xval), mask=np.ma.getmask(val)
        )
        result[np.isinf(val.data)] = -np.inf
        if vtype == "scalar":
            result = result[0]
        return result

    def autoscale_None(self, val):
        changed = False
        if self.vmin is None:
            self.vmin = val.min()
            changed = True
        if self.vmax is None:
            self.vmax = val.max()
            changed = True
        if changed or self.xval is None or self.yval is None:
            self._set_xyvals(val)

    def autoscale(self, val):
        self.vmin = val.min()
        self.vmax = val.max()
        self._set_xyvals(val)

    def _set_xyvals(self, val):
        data = np.ma.asarray(val).ravel()
        w = np.isinf(data.data)
        if data.mask is not np.ma.nomask:
            w = w | data.mask
        data2 = data.data[~w]
        if data2.size < 3:
            self.yval = np.array([0, 1], dtype=np.float)
            self.xval = np.array([self.vmin, self.vmax], dtype=np.float)
            return
        bins = min(data2.size // 20, 5000)
        if bins < 3:
            bins = data2.size
        try:
            # for numpy 1.1, use new bins format (left and right edges)
            hist, bins = np.histogram(
                data2, bins=bins, range=(self.vmin, self.vmax), new=True
            )
        except TypeError:
            # for numpy <= 1.0 or numpy >= 1.2, no new keyword
            hist, bins = np.histogram(data2, bins=bins, range=(self.vmin, self.vmax))
        if bins.size == hist.size + 1:
            # new bins format, remove last point
            bins = bins[:-1]
        hist = hist.astype(np.float) / np.float(hist.sum())
        self.yval = np.concatenate([[0.0], hist.cumsum(), [1.0]])
        self.xval = np.concatenate(
            [[self.vmin], bins + 0.5 * (bins[1] - bins[0]), [self.vmax]]
        )

    def _lininterp(self, x, X, Y):
        if hasattr(x, "__len__"):
            xtype = "array"
            xx = np.asarray(x).astype(np.float)
        else:
            xtype = "scalar"
            xx = np.asarray([x]).astype(np.float)
        idx = X.searchsorted(xx)
        yy = xx * 0
        yy[idx > len(X) - 1] = Y[-1]  # over
        yy[idx <= 0] = Y[0]  # under
        wok = np.where((idx > 0) & (idx < len(X)))  # the good ones
        iok = idx[wok]
        yywok = Y[iok - 1] + (
            (Y[iok] - Y[iok - 1]) / (X[iok] - X[iok - 1]) * (xx[wok] - X[iok - 1])
        )
        w = np.where(((X[iok] - X[iok - 1]) == 0))  # where are the nan ?
        yywok[w] = Y[iok[w] - 1]  # replace by previous value
        wl = np.where(xx[wok] == X[0])
        yywok[wl] = Y[0]
        wh = np.where(xx[wok] == X[-1])
        yywok[wh] = Y[-1]
        yy[wok] = yywok
        if xtype == "scalar":
            yy = yy[0]
        return yy
    
class LinNorm2(matplotlib.colors.Normalize):
    """
    Normalize a given value to the 0-1 range on a lin scale
    """

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = "array"
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = "scalar"
            val = np.ma.array([value]).astype(np.float)

        winf = np.isinf(val.data)
        val = np.ma.masked_where(winf, val)

        self.autoscale_None(val)
        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            return type(value)(0.0 * np.asarray(value))
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax), mask=mask)
            result = (val - vmin) * (1.0 / (vmax - vmin))
            result.data[result.data < 0] = 0.0
            result.data[result.data > 1] = 1.0
            result[winf] = -np.inf
            if result.mask is not np.ma.nomask:
                result.mask[winf] = False
        if vtype == "scalar":
            result = result[0]
        return result

    def autoscale_None(self, A):
        " autoscale only None-valued vmin or vmax"
        if self.vmin is None or self.vmax is None:
            val = np.ma.masked_where(np.isinf(A.data), A)
            matplotlib.colors.Normalize.autoscale_None(self, val)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = float(self.vmin), float(self.vmax)

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin + (vmax - vmin) * val
        else:
            return vmin + (vmax - vmin) * value

#--------------------------------------------------------#
#                                                        #
#                                                        #
#                                                        #
#--------------------------------------------------------#
        
        
def mollview(
    map=None,
    fig=None,
    rot=None,
    coord=None,
    unit="",
    xsize=800,
    title="Mollweide view",
    nest=False,
    min=None,
    max=None,
    flip="astro",
    gal_cut=0,
    format="%g",
    format2="%g",
    cbar=True,
    cmap=None,
    notext=False,
    norm=None,
    hold=False,
    margins=None,
    sub=None,
    nlocs=2,
    return_projected_map=False,
):
    # Create the figure
    import pylab

    # Ensure that the nside is valid
    nside = get_nside(map)
    check_nside(nside, nest=nest)

    if not (hold or sub):
        f = pylab.figure(fig, figsize=(8.5, 5.4))
        extent = (0.02, 0.05, 0.96, 0.9)
    elif hold:
        f = pylab.gcf()
        left, bottom, right, top = np.array(f.gca().get_position()).ravel()
        extent = (left, bottom, right - left, top - bottom)
        f.delaxes(f.gca())
    else:  # using subplot syntax
        f = pylab.gcf()
        if hasattr(sub, "__len__"):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub // 100, (sub % 100) // 10, (sub % 10)
        if idx < 1 or idx > ncols * nrows:
            raise ValueError("Wrong values for sub: %d, %d, %d" % (nrows, ncols, idx))
        c, r = (idx - 1) % ncols, (idx - 1) // ncols
        if not margins:
            margins = (0.01, 0.0, 0.0, 0.02)
        extent = (
            c * 1.0 / ncols + margins[0],
            1.0 - (r + 1) * 1.0 / nrows + margins[1],
            1.0 / ncols - margins[2] - margins[0],
            1.0 / nrows - margins[3] - margins[1],
        )
        extent = (
            extent[0] + margins[0],
            extent[1] + margins[1],
            extent[2] - margins[2] - margins[0],
            extent[3] - margins[3] - margins[1],
        )
        # extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    # f=pylab.figure(fig,figsize=(8.5,5.4))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12) + np.inf
            cbar = False
        map = ma_to_array(map)
        ax = HpxMollweideAxes(
            f, extent, coord=coord, rot=rot, format=format2, flipconv=flip
        )
        f.add_axes(ax)
        img = ax.projmap(
            map,
            nest=nest,
            xsize=xsize,
            coord=coord,
            vmin=min,
            vmax=max,
            cmap=cmap,
            norm=norm,
        )
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0, 1, im.cmap.N + 1))
            v = np.linspace(im.norm.vmin, im.norm.vmax, im.cmap.N)
            if matplotlib.__version__ >= "0.91.0":
                cb = f.colorbar(
                    im,
                    ax=ax,
                    orientation="horizontal",
                    shrink=0.5,
                    aspect=25,
                    ticks=BoundaryLocator(nlocs, norm),
                    pad=0.05,
                    fraction=0.1,
                    boundaries=b,
                    values=v,
                    format=format,
                )
            else:
                # for older matplotlib versions, no ax kwarg
                cb = f.colorbar(
                    im,
                    orientation="horizontal",
                    shrink=0.5,
                    aspect=25,
                    ticks=BoundaryLocator(nlocs, norm),
                    pad=0.05,
                    fraction=0.1,
                    boundaries=b,
                    values=v,
                    format=format,
                )
            cb.solids.set_rasterized(True)
        ax.set_title(title)
        if not notext:
            ax.text(
                0.86,
                0.05,
                ax.proj.coordsysstr,
                fontsize=14,
                fontweight="bold",
                transform=ax.transAxes,
            )
        if cbar:
            cb.ax.text(
                0.5,
                -1.0,
                unit,
                fontsize=14,
                transform=cb.ax.transAxes,
                ha="center",
                va="center",
            )
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            # pylab.show()
    if return_projected_map:
        return img


def graticule(dpar=None, dmer=None, coord=None, local=None, **kwds):
    import pylab

    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if len(f.get_axes()) == 0:
            ax = HpxMollweideAxes(f, (0.02, 0.05, 0.96, 0.9), coord=coord)
            f.add_axes(ax)
            ax.text(
                0.86,
                0.05,
                ax.proj.coordsysstr,
                fontsize=14,
                fontweight="bold",
                transform=ax.transAxes,
            )
        for ax in f.get_axes():
            if isinstance(ax, SphericalProjAxes):
                ax.graticule(dpar=dpar, dmer=dmer, coord=coord, local=local, **kwds)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()


class BoundaryLocator(matplotlib.ticker.Locator):
    def __init__(self, N=2, norm=None):
        if N < 2:
            raise ValueError("Number of locs must be greater than 1")
        self.Nlocs = N
        self.norm = norm

    def __call__(self):
        if matplotlib.__version__ < "0.98":
            vmin, vmax = self.viewInterval.get_bounds()
        else:
            vmin, vmax = self.axis.get_view_interval()
        if self.norm == "log":
            locs = np.log10(vmin) + np.arange(self.Nlocs) * (
                np.log10(vmax) - np.log10(vmin)
            ) / (self.Nlocs - 1.0)
            locs = 10 ** (locs)
        else:
            locs = vmin + np.arange(self.Nlocs) * (vmax - vmin) / (self.Nlocs - 1.0)
        return locs

    def autoscale(self):
        self.verify_intervals()
        vmin, vmax = self.dataInterval.get_bounds()
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        return vmin, vmax



    
class SphericalProjAxes(matplotlib.axes.Axes):

    def __init__(self, ProjClass, *args, **kwds):
        if not issubclass(ProjClass, SphericalProj):
            raise TypeError(
                "First argument must be a SphericalProj class " "(or derived from)"
            )
        self.proj = ProjClass(
            rot=kwds.pop("rot", None),
            coord=kwds.pop("coord", None),
            flipconv=kwds.pop("flipconv", None),
            **kwds.pop("arrayinfo", {})
        )
        kwds.setdefault("format", "%g")
        kwds.setdefault("coordprec", 2)
        kwds["aspect"] = "equal"
        super(SphericalProjAxes, self).__init__(*args, **kwds)
        self.axis("off")
        self.set_autoscale_on(False)
        xmin, xmax, ymin, ymax = self.proj.get_extent()
        self.set_xlim(xmin, xmax)
        self.set_ylim(ymin, ymax)
        dx, dy = self.proj.ang2xy(pi / 2.0, 1.0 * dtor, direct=True)
        self._segment_threshold = 16.0 * np.sqrt(dx ** 2 + dy ** 2)
        self._segment_step_rad = 0.1 * pi / 180
        self._do_border = True
        self._gratdef = {}
        self._gratdef["local"] = False
        self._gratdef["dpar"] = 30.0
        
    def set_format(self, f):
        """Set the format string for value display
        """
        self._format = f
        return f

    def set_coordprec(self, n):
        """Set the number of digits after floating point for coord display.
        """
        self._coordprec = n

    def projmap(
        self,
        map,
        vec2pix_func,
        vmin=None,
        vmax=None,
        badval=UNSEEN,
        cmap=None,
        norm=None,
        rot=None,
        coord=None,
        **kwds
    ):
        img = self.proj.projmap(map, vec2pix_func, rot=rot, coord=coord)
        w = ~(np.isnan(img) | np.isinf(img) | mask_bad(img, badval=badval))
        try:
            if vmin is None:
                vmin = img[w].min()
        except ValueError:
            vmin = 0.0
        try:
            if vmax is None:
                vmax = img[w].max()
        except ValueError:
            vmax = 0.0
        if vmin > vmax:
            vmin = vmax
        if vmin == vmax:
            vmin -= 1.0
            vmax += 1.0
        cm, nn = get_color_table(vmin, vmax, img[w], cmap=cmap, norm=norm)
        ext = self.proj.get_extent()
        img = np.ma.masked_values(img, badval)
        aximg = self.imshow(
            img,
            extent=ext,
            cmap=cm,
            norm=nn,
            interpolation="nearest",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            **kwds
        )
        xmin, xmax, ymin, ymax = self.proj.get_extent()
        self.set_xlim(xmin, xmax)
        self.set_ylim(ymin, ymax)
        return img
    
    def graticule(
        self, dpar=None, dmer=None, coord=None, local=None, verbose=True, **kwds
    ):
        """Draw a graticule.

        Input:
         - dpar: angular separation between parallels in degree
         - dmer: angular separation between meridians in degree
         - coord: coordinate system of the graticule ('G', 'E' or 'C')
         - local: if True, no rotation performed at all
        """
        gratargs = (dpar, dmer, coord, local)
        gratkwds = kwds
        if dpar is None:
            dpar = self._gratdef["dpar"]
        if local is None:
            local = self._gratdef["local"]
        if dmer is None:
            dmer = dpar
        dpar = abs(dpar) * dtor
        dmer = abs(dmer) * dtor
        if not local:
            vec = dir2vec(self.proj.get_center())
            vec0 = Rotator(coord=self.proj.mkcoord(coord=coord)).I(vec)
        else:
            vec = (1, 0, 0)
            vec0 = (1, 0, 0)
        u_pmin, u_pmax = kwds.pop("pmax", None), kwds.pop("pmin", None)
        u_mmin, u_mmax = kwds.pop("mmin", None), kwds.pop("mmax", None)
        if u_pmin:
            u_pmin = (pi / 2.0 - u_pmin * dtor) % pi
        if u_pmax:
            u_pmax = (pi / 2.0 - u_pmax * dtor) % pi
        if u_mmin:
            u_mmin = (((u_mmin + 180.0) % 360) - 180) * dtor
        if u_mmax:
            u_mmax = (((u_mmax + 180.0) % 360) - 180) * dtor
        pmin, pmax = self.get_parallel_interval(vec0)
        mmin, mmax = self.get_meridian_interval(vec0)
        if u_pmin:
            pmin = u_pmin
        if u_pmax:
            pmax = u_pmax
        if u_mmin:
            mmin = u_mmin
        if u_mmax:
            mmax = u_pmax
        if verbose:
            print(
                "{0} {1} {2} {3}".format(
                    pmin / dtor, pmax / dtor, mmin / dtor, mmax / dtor
                )
            )
        if not kwds.pop("force", False):
            dpar, dmer = self._get_interv_graticule(
                pmin, pmax, dpar, mmin, mmax, dmer, verbose=verbose
            )
        theta_list = np.around(np.arange(pmin, pmax + 0.5 * dpar, dpar) / dpar) * dpar
        phi_list = np.around(np.arange(mmin, mmax + 0.5 * dmer, dmer) / dmer) * dmer
        theta = np.arange(
            pmin, pmax, min((pmax - pmin) / 100.0, self._segment_step_rad)
        )
        phi = np.arange(mmin, mmax, min((mmax - mmin) / 100.0, self._segment_step_rad))
        equator = False
        gratlines = []
        kwds.setdefault("lw", 1)
        kwds.setdefault("color", "k")
        for t in theta_list:
            if abs(t - pi / 2.0) < 1.0e-10:
                fmt = "-"
                equator = True
            elif abs(t) < 1.0e-10:  # special case: north pole
                t = 1.0e-10
                fmt = "-"
            elif abs(t - pi) < 1.0e-10:  # special case: south pole
                t = pi - 1.0e-10
                fmt = "-"
            else:
                fmt = ":"
            gratlines.append(
                self.projplot(
                    phi * 0.0 + t, phi, fmt, coord=coord, direct=local, **kwds
                )
            )
        if not equator and pmin <= pi / 2.0 and pi / 2 <= pmax:
            gratlines.append(
                self.projplot(
                    phi * 0.0 + pi / 2.0, phi, "-", coord=coord, direct=local, **kwds
                )
            )
        for p in phi_list:
            if abs(p) < 1.0e-10:
                fmt = "-"
            else:
                fmt = ":"
            gratlines.append(
                self.projplot(
                    theta, theta * 0.0 + p, fmt, coord=coord, direct=local, **kwds
                )
            )
        # Now the borders (only useful for full sky projection)
        if hasattr(self, "_do_border") and self._do_border:
            theta = np.arange(0, 181) * dtor
            gratlines.append(
                self.projplot(theta, theta * 0 - pi, "-k", lw=1, direct=True)
            )
            gratlines.append(
                self.projplot(theta, theta * 0 + 0.9999 * pi, "-k", lw=1, direct=True)
            )
            phi = np.arange(-180, 180) * dtor
            gratlines.append(
                self.projplot(phi * 0 + 1.0e-10, phi, "-k", lw=1, direct=True)
            )
            gratlines.append(
                self.projplot(phi * 0 + pi - 1.0e-10, phi, "-k", lw=1, direct=True)
            )
        if hasattr(self, "_graticules"):
            self._graticules.append((gratargs, gratkwds, gratlines))
        else:
            self._graticules = [(gratargs, gratkwds, gratlines)]
        return dpar, dmer

    def get_parallel_interval(self, vx, vy=None, vz=None):
        """Get the min and max value of theta of the parallel to cover the
        field of view.

        Input:
          - the normalized vector of the direction of the center of the
            projection, in the reference frame of the graticule.
        Return:
          - vmin,vmax : between 0 and pi, vmin<vmax, the interval of theta
                        for the parallels crossing the field of view
        """
        if vy is None and vz is None:
            vx, vy, vz = vx
        elif vy is None or vz is None:
            raise ValueError("Both vy and vz must be given or both not given")
        a = np.arccos(vz)
        fov = self.proj.get_fov()
        vmin = max(0.0, a - fov / 2.0)
        vmax = min(pi, a + fov / 2.0)
        return vmin, vmax

    def get_meridian_interval(self, vx, vy=None, vz=None):
        """Get the min and max value of phi of the meridians to cover the field
        of view.

        Input:
          - the normalized vector of the direction of the center of the
            projection, in the reference frame of the graticule.
        Return:
          - vmin,vmax : the interval of phi for the
                        meridians crossing the field of view.
        """
        if vy is None and vz is None:
            vx, vy, vz = vx
        elif vy is None or vz is None:
            raise ValueError("Both vy and vz must be given or both not given")
        fov = self.proj.get_fov()
        th = np.arccos(vz)
        if th <= fov / 2.0:  # test whether north pole is visible
            return -np.pi, np.pi
        if abs(th - pi) <= fov / 2.0:  # test whether south pole is visible
            return -np.pi, np.pi
        sth = np.sin(th)
        phi0 = np.arctan2(vy, vx)
        return phi0 - fov / sth / 2.0, phi0 + fov / sth / 2.0

    def _get_interv_graticule(self, pmin, pmax, dpar, mmin, mmax, dmer, verbose=True):
        def set_prec(d, n, nn=2):
            arcmin = False
            if d / n < 1.0:
                d *= 60
                arcmin = True
                nn = 1
            x = d / n
            y = nn * x
            ex = np.floor(np.log10(y))
            z = np.around(y / 10 ** ex) * 10 ** ex / nn
            if arcmin:
                z = 1.0 / np.around(60.0 / z)
            return z

        max_n_par = 18
        max_n_mer = 36
        n_par = (pmax - pmin) / dpar
        n_mer = (mmax - mmin) / dmer
        if n_par > max_n_par:
            dpar = set_prec((pmax - pmin) / dtor, max_n_par / 2) * dtor
        if n_mer > max_n_mer:
            dmer = set_prec((mmax - mmin) / dtor, max_n_mer / 2, nn=1) * dtor
        if dmer / dpar < 0.2 or dmer / dpar > 5.0:
            dmer = dpar = max(dmer, dpar)
        vdeg = int(np.floor(np.around(dpar / dtor, 10)))
        varcmin = (dpar / dtor - vdeg) * 60.0
        if verbose:
            print(
                "The interval between parallels is {0:d} deg {1:.2f}'.".format(
                    vdeg, varcmin
                )
            )
        vdeg = int(np.floor(np.around(dmer / dtor, 10)))
        varcmin = (dmer / dtor - vdeg) * 60.0
        if verbose:
            print(
                "The interval between meridians is {0:d} deg {1:.2f}'.".format(
                    vdeg, varcmin
                )
            )
        return dpar, dmer

    def projplot(self, *args, **kwds):
        fmt = None
        if len(args) < 1:
            raise ValueError("No argument given")
        if len(args) == 1:
            theta, phi = np.asarray(args[0])
        elif len(args) == 2:
            if type(args[1]) is str:
                fmt = args[1]
                theta, phi = np.asarray(args[0])
            else:
                theta, phi = np.asarray(args[0]), np.asarray(args[1])
        elif len(args) == 3:
            if type(args[2]) is not str:
                raise TypeError("Third argument must be a string")
            else:
                theta, phi = np.asarray(args[0]), np.asarray(args[1])
                fmt = args[2]
        else:
            raise TypeError("Three args maximum")
        rot = kwds.pop("rot", None)
        if rot is not None:
            rot = np.array(np.atleast_1d(rot), copy=1)
            rot.resize(3)
            rot[1] = rot[1] - 90.0
        coord = self.proj.mkcoord(kwds.pop("coord", None))[::-1]
        lonlat = kwds.pop("lonlat", False)
        vec = dir2vec(theta, phi, lonlat=lonlat)
        vec = (Rotator(rot=rot, coord=coord, eulertype="Y")).I(vec)
        x, y = self.proj.vec2xy(vec, direct=kwds.pop("direct", False))
        x, y = self._make_segment(
            x, y, threshold=kwds.pop("threshold", self._segment_threshold)
        )
        thelines = []
        for xx, yy in zip(x, y):
            if fmt is not None:
                try:  # works in matplotlib 1.3 and earlier
                    linestyle, marker, color = matplotlib.axes._process_plot_format(fmt)
                except:  # matplotlib 1.4 and later
                    linestyle, marker, color = matplotlib.axes._axes._process_plot_format(
                        fmt
                    )
                kwds.setdefault("linestyle", linestyle)
                kwds.setdefault("marker", marker)
                if color is not None:
                    kwds.setdefault("color", color)
            l = matplotlib.lines.Line2D(xx, yy, **kwds)
            self.add_line(l)
            thelines.append(l)
        return thelines

    def _make_segment(self, x, y, threshold=None):
        if threshold is None:
            threshold = self._segment_threshold
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        d2 = np.sqrt((np.roll(x, 1) - x) ** 2 + (np.roll(y, 1) - y) ** 2)
        w = np.where(d2 > threshold)[0]
        # w=w[w!=0]
        xx = []
        yy = []
        if len(w) == 1:
            x = np.roll(x, -w[0])
            y = np.roll(y, -w[0])
            xx.append(x)
            yy.append(y)
        elif len(w) >= 2:
            xx.append(x[0 : w[0]])
            yy.append(y[0 : w[0]])
            for i in six.moves.xrange(len(w) - 1):
                xx.append(x[w[i] : w[i + 1]])
                yy.append(y[w[i] : w[i + 1]])
            xx.append(x[w[-1] :])
            yy.append(y[w[-1] :])
        else:
            xx.append(x)
            yy.append(y)
        return xx, yy



class MollweideAxes(SphericalProjAxes):

    def __init__(self, *args, **kwds):
        kwds.setdefault("coordprec", 2)
        super(MollweideAxes, self).__init__(MollweideProj, *args, **kwds)
        self.set_xlim(-2.01, 2.01)
        self.set_ylim(-1.01, 1.01)

    def projmap(self, map, vec2pix_func, xsize=800, **kwds):
        self.proj.set_proj_plane_info(xsize=xsize)
        img = super(MollweideAxes, self).projmap(map, vec2pix_func, **kwds)
        self.set_xlim(-2.01, 2.01)
        self.set_ylim(-1.01, 1.01)
        return img


class HpxMollweideAxes(MollweideAxes):
    def projmap(self, map, nest=False, **kwds):
        nside = npix2nside(get_map_size(map))
        f = lambda x, y, z: vec2pix(nside, x, y, z, nest=nest)
        return super(HpxMollweideAxes, self).projmap(map, f, **kwds)



#--------------------------------------------------------#
#                                                        #
#                                                        #
#                                                        #
#--------------------------------------------------------#




class SphericalProj(object):
    """
    This class defines functions for spherical projection.
    
    This class contains class method for spherical projection computation. It 
    should not be instantiated. It should be inherited from and methods should
    be overloaded for desired projection.
    """

    name = "None"

    def __init__(self, rot=None, coord=None, flipconv=None, **kwds):
        self.rotator = Rotator(rot=rot, coord=None, eulertype="ZYX")
        self.coordsys = Rotator(coord=coord).coordout
        self.coordsysstr = Rotator(coord=coord).coordoutstr
        self.set_flip(flipconv)
        self.set_proj_plane_info(**kwds)

    def set_proj_plane_info(self, **kwds):
        allNone = True
        for v in kwds.values():
            if v is not None:
                allNone = False
        if not allNone:
            self._arrayinfo = dict(kwds)
        else:
            self._arrayinfo = None

    def get_proj_plane_info(self):
        return self._arrayinfo

    arrayinfo = property(
        get_proj_plane_info, doc="Dictionary with information on the projection array"
    )


    def __eq__(self, a):
        if type(a) is not type(self):
            return False
        return (self.rotator == a.rotator) and (self.coordsys == a.coordsys)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        pass

    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        pass

    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        pass

    def xy2vec(self, x, y=None, direct=False):
        pass

    def xy2ij(self, x, y=None):
        pass

    def ij2xy(self, i=None, j=None):
        pass

    def projmap(self, map, vec2pix_func, rot=None, coord=None):
        x, y = self.ij2xy()
        if np.__version__ >= "1.1":
            matype = np.ma.core.MaskedArray
        else:
            matype = np.ma.array
        if type(x) is matype and x.mask is not np.ma.nomask:
            w = x.mask == False
        else:
            w = slice(None)
        img = np.zeros(x.shape, np.float64) - np.inf
        vec = self.xy2vec(np.asarray(x[w]), np.asarray(y[w]))
        vec = (Rotator(rot=rot, coord=self.mkcoord(coord))).I(vec)
        pix = vec2pix_func(vec[0], vec[1], vec[2])
        # support masked array for map, or a dictionnary (for explicit pixelisation)
        if isinstance(map, matype) and map.mask is not np.ma.nomask:
            mpix = map[pix]
            mpix[map.mask[pix]] = UNSEEN
        elif isinstance(map, dict):
            is_pix_seen = np.in1d(pix, map.keys()).reshape(pix.shape)
            is_pix_unseen = ~is_pix_seen
            mpix = np.zeros_like(img[w])
            mpix[is_pix_unseen] = UNSEEN
            pix_seen = pix[is_pix_seen]
            iterable = (map[p] for p in pix_seen)
            mpix[is_pix_seen] = np.fromiter(iterable, mpix.dtype, count=pix_seen.size)
        else:
            mpix = map[pix]
        img[w] = mpix
        return img

    def set_flip(self, flipconv):
        """flipconv is either 'astro' or 'geo'. None will be default.
        
        With 'astro', east is toward left and west toward right. 
        It is the opposite for 'geo'
        """
        if flipconv is None:
            flipconv = "astro"  # default
        if flipconv == "astro":
            self._flip = -1
        elif flipconv == "geo":
            self._flip = 1
        else:
            raise ValueError("flipconv must be 'astro', 'geo' or None for default.")

    def get_extent(self):
        pass

    def get_fov(self):
        return 2.0 * pi

    def get_center(self, lonlat=False):
        lon, lat = np.asarray(self.rotator.rots[0][0:2]) * 180 / pi
        if lonlat:
            return lon, lat
        else:
            return pi / 2.0 - lat * dtor, lon * dtor

    def mkcoord(self, coord):
        if self.coordsys is None:
            return (coord, coord)
        elif coord is None:
            return (self.coordsys, self.coordsys)
        elif type(coord) is str:
            return (coord, self.coordsys)
        else:
            return (tuple(coord)[0], self.coordsys)
        




class MollweideProj(SphericalProj):
    """This class provides class methods for Mollweide projection.
    """

    name = "Mollweide"
    __molldata = []

    def __init__(self, rot=None, coord=None, xsize=800, **kwds):
        self.__initialise_data()
        super(MollweideProj, self).__init__(rot=rot, coord=coord, xsize=xsize, **kwds)

    def set_proj_plane_info(self, xsize):
        super(MollweideProj, self).set_proj_plane_info(xsize=xsize)

    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct:
            theta, phi = vec2dir(self.rotator(vx, vy, vz))
        else:
            theta, phi = vec2dir(vx, vy, vz)
        flip = self._flip
        X, Y = MollweideProj.__molldata
        # set phi in [-pi,pi]
        phi = (phi + pi) % (2 * pi) - pi
        lat = pi / 2.0 - theta
        A = MollweideProj.__lininterp(X, Y, lat)
        x = flip * 2.0 / pi * phi * np.cos(A)
        y = np.sin(A)
        return x, y

    # vec2xy.__doc__ = SphericalProj.vec2xy.__doc__ % (name, name)

    def xy2vec(self, x, y=None, direct=False):
        flip = self._flip
        if y is None:
            x, y = x
        mask = np.asarray(x) ** 2 / 4.0 + np.asarray(y) ** 2 > 1.0
        w = np.where(mask == False)
        if not mask.any():
            mask = np.ma.nomask
        if not hasattr(x, "__len__"):
            if mask is not np.ma.nomask:
                return np.nan, np.nan, np.nan
            else:
                s = np.sqrt((1 - y) * (1 + y))
                a = np.arcsin(y)
                z = 2.0 / pi * (a + y * s)
                phi = flip * pi / 2.0 * x / np.maximum(s, 1.0e-6)
                sz = np.sqrt((1 - z) * (1 + z))
                vec = sz * np.cos(phi), sz * np.sin(phi), z
                if not direct:
                    return self.rotator.I(vec)
                else:
                    return vec
        else:
            vec = (
                np.zeros(x.shape) + np.nan,
                np.zeros(x.shape) + np.nan,
                np.zeros(x.shape) + np.nan,
            )
            s = np.sqrt((1 - y[w]) * (1 + y[w]))
            a = np.arcsin(y[w])
            vec[2][w] = 2.0 / pi * (a + y[w] * s)
            phi = flip * pi / 2.0 * x[w] / np.maximum(s, 1.0e-6)
            sz = np.sqrt((1 - vec[2][w]) * (1 + vec[2][w]))
            vec[0][w] = sz * np.cos(phi)
            vec[1][w] = sz * np.sin(phi)
            if not direct:
                return self.rotator.I(vec)
            else:
                return vec

    # xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name, name)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        return self.vec2xy(dir2vec(theta, phi, lonlat=lonlat), direct=direct)

    # ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name, name)

    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        vec = self.xy2vec(x, y, direct=direct)
        return vec2dir(vec, lonlat=lonlat)

    # xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name, name)

    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError(
                "No projection plane array information defined for " "this projector"
            )
        xsize = int(self.arrayinfo["xsize"])
        ysize = xsize // 2
        if y is None:
            x, y = x
        xc, yc = (xsize - 1.0) / 2.0, (ysize - 1.0) / 2.0
        if hasattr(x, "__len__"):
            j = np.around(x * xc / 2.0 + xc).astype(np.long)
            i = np.around(yc + y * yc).astype(np.long)
            mask = x ** 2 / 4.0 + y ** 2 > 1.0
            if not mask.any():
                mask = np.ma.nomask
            j = np.ma.array(j, mask=mask)
            i = np.ma.array(i, mask=mask)
        else:
            if x ** 2 / 4.0 + y ** 2 > 1.0:
                i, j = np.nan, np.nan
            else:
                j = np.around(x * xc / 2.0 + xc).astype(np.long)
                i = np.around(yc + y * yc).astype(np.long)
        return i, j

    # xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name, name)

    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError(
                "No projection plane array information defined for " "this projector"
            )
        xsize = int(self.arrayinfo["xsize"])
        ysize = xsize // 2
        xc, yc = (xsize - 1.0) / 2.0, (ysize - 1.0) / 2.0
        if i is None and j is None:
            idx = np.outer(np.arange(ysize), np.ones(xsize))
            y = (idx - yc) / yc
            idx = np.outer(np.ones(ysize), np.arange(xsize))
            x = 2.0 * (idx - xc) / xc
            mask = x ** 2 / 4.0 + y ** 2 > 1.0
            if not mask.any():
                mask = np.ma.nomask
            x = np.ma.array(x, mask=mask)
            y = np.ma.array(y, mask=mask)
        elif i is not None and j is not None:
            y = (np.asarray(i) - yc) / yc
            x = 2.0 * (np.asarray(j) - xc) / xc
            if x ** 2 / 4.0 + y ** 2 > 1.0:
                x, y = np.nan, np.nan
        elif i is not None and j is None:
            i, j = i
            y = (np.asarray(i) - yc) / yc
            x = 2.0 * (np.asarray(j) - xc) / xc
            if x ** 2 / 4.0 + y ** 2 > 1.0:
                x, y = np.nan, np.nan
        else:
            raise TypeError("i and j must be both given or both not given")
        return x, y

    # ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name, name)

    def get_extent(self):
        return (-2.0, 2.0, -1.0, 1.0)

    @staticmethod
    def __initialise_data():
        if len(MollweideProj.__molldata) == 0:
            X = (np.arange(1.0, 180.0, 1.0) - 90.0) * dtor
            Y = MollweideProj.__findRoot(
                MollweideProj.__fmoll, MollweideProj.__dfmoll, X.copy(), X, niter=10
            )
            X = np.concatenate([[-pi / 2], X, [pi / 2]])
            Y = np.concatenate([[-pi / 2], Y, [pi / 2]])
            MollweideProj.__molldata.append(X)
            MollweideProj.__molldata.append(Y)
        return

    @staticmethod
    def __findRoot(f, df, x0, argsf=None, argsdf=None, niter=100):
        x = x0
        niter = min(abs(niter), 1000)
        i = 0
        while i < niter:
            dx = -f(x, argsf) / df(x, argsdf)
            x += dx
            i += 1
        return x

    @staticmethod
    def __fmoll(x, args):
        return 2.0 * x + np.sin(2.0 * x) - pi * np.sin(args)

    @staticmethod
    def __dfmoll(x, args):
        return 2.0 * (1.0 + np.cos(2.0 * x))

    @staticmethod
    def __lininterp(X, Y, x):
        idx = X.searchsorted(x)
        y = Y[idx - 1] + (Y[idx] - Y[idx - 1]) / (X[idx] - X[idx - 1]) * (
            x - X[idx - 1]
        )
        return y


#--------------------------------------------------------#
#                                                        #
#                                                        #
#                                                        #
#--------------------------------------------------------#
coordname = {"G": "Galactic", "E": "Ecliptic", "C": "Equatorial"}


def dir2vec(theta, phi=None, lonlat=False):
    if phi is None:
        theta, phi = theta
    if lonlat:
        lon, lat = theta, phi
        theta, phi = np.pi / 2.0 - np.radians(lat), np.radians(lon)
    ct, st, cp, sp = np.cos(theta), np.sin(theta), np.cos(phi), np.sin(phi)
    vec = np.empty((3, ct.size), np.float64)
    vec[0, :] = st * cp
    vec[1, :] = st * sp
    vec[2, :] = ct
    return vec.squeeze()

def vec2dir(vec, vy=None, vz=None, lonlat=False):
    if np.any(np.isnan(vec)):
        return np.nan, np.nan
    if vy is None and vz is None:
        vx, vy, vz = vec
    elif vy is not None and vz is not None:
        vx = vec
    else:
        raise TypeError("You must either give both vy and vz or none of them")
    r = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    ang = np.empty((2, r.size))
    ang[0, :] = np.arccos(vz / r)
    ang[1, :] = np.arctan2(vy, vx)
    if lonlat:
        ang = np.degrees(ang)
        np.negative(ang[0, :], ang[0, :])
        ang[0, :] += 90.0
        return ang[::-1, :].squeeze()
    else:
        return ang.squeeze()

def check_coord(c):
    if c is None:
        return c
    if not isinstance(c, six.string_types):
        raise TypeError(
            "Coordinate must be a string (G[alactic],"
            " E[cliptic], C[elestial]"
            " or Equatorial=Celestial)"
        )
    if c[0].upper() == "G":
        x = "G"
    elif c[0].upper() == "E" and c != "Equatorial":
        x = "E"
    elif c[0].upper() == "C" or c == "Equatorial":
        x = "C"
    else:
        raise ValueError(
            "Wrong coordinate (either G[alactic],"
            " E[cliptic], C[elestial]"
            " or Equatorial=Celestial)"
        )
    return x
    

def normalise_rot(rot, deg=False):
    if deg:
        convert = np.pi / 180.0
    else:
        convert = 1.0
    if rot is None:
        rot = np.zeros(3)
    else:
        rot = np.array(rot, np.float64).flatten() * convert
        rot.resize(3, refcheck=False)
    return rot

def normalise_coord(coord):
    coord_norm = []
    if coord is None:
        coord = (None, None)
    coord = tuple(coord)
    if len(coord) > 2:
        raise TypeError(
            "Coordinate must be a string (G[alactic],"
            " E[cliptic] or C[elestial])"
            " or a sequence of 2 strings"
        )
    for x in coord:
        coord_norm.append(check_coord(x))
    if len(coord_norm) < 2:
        coord_norm.append(coord_norm[0])
    return tuple(coord_norm)

class ConsistencyWarning(Warning):
    """Warns for a problem in the consistency of data
    """

    pass

def euler_matrix_new(a1, a2, a3, X=True, Y=False, ZYX=False, deg=False):
    t_k = 0
    if ZYX:
        t_k = t_k + 1
    # if X:   t_k = t_k + 1
    if Y:
        t_k = t_k + 1
    if t_k > 1:
        raise ValueError("Choose either X, Y or ZYX convention")

    convert = 1.0
    if deg:
        convert = np.pi / 180.0

    c1 = np.cos(a1 * convert)
    s1 = np.sin(a1 * convert)
    c2 = np.cos(a2 * convert)
    s2 = np.sin(a2 * convert)
    c3 = np.cos(a3 * convert)
    s3 = np.sin(a3 * convert)

    if ZYX:
        m1 = np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])  # around   z

        m2 = np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]])  # around   y

        m3 = np.array([[1, 0, 0], [0, c3, -s3], [0, s3, c3]])  # around   x

    elif Y:
        m1 = np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])  # around   z

        m2 = np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]])  # around   y

        m3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])  # around   z

    else:
        m1 = np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])  # around   z

        m2 = np.array([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])  # around   x

        m3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])  # around   z

    M = np.dot(m3.T, np.dot(m2.T, m1.T))

    return M


def get_rotation_matrix(rot, deg=False, eulertype="ZYX"):
    rot = normalise_rot(rot, deg=deg)
    if not np.allclose(rot, np.zeros(3), rtol=0.0, atol=1.0e-15):
        do_rot = True
    else:
        do_rot = False
    if eulertype == "X":
        matrot = euler_matrix_new(rot[0], -rot[1], rot[2], X=True)
    elif eulertype == "Y":
        matrot = euler_matrix_new(rot[0], -rot[1], rot[2], Y=True)
    else:
        matrot = euler_matrix_new(rot[0], -rot[1], rot[2], ZYX=True)

    return matrot, do_rot, rot

def get_coordconv_matrix(coord):
    coord_norm = normalise_coord(coord)

    if coord_norm[0] == coord_norm[1]:
        matconv = np.identity(3)
        do_conv = False
    else:
        eps = 23.452294 - 0.0130125 - 1.63889e-6 + 5.02778e-7
        eps = eps * np.pi / 180.0

        # ecliptic to galactic
        e2g = np.array(
            [
                [-0.054882486, -0.993821033, -0.096476249],
                [0.494116468, -0.110993846, 0.862281440],
                [-0.867661702, -0.000346354, 0.497154957],
            ]
        )

        # ecliptic to equatorial
        e2q = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, np.cos(eps), -1.0 * np.sin(eps)],
                [0.0, np.sin(eps), np.cos(eps)],
            ]
        )

        # galactic to ecliptic
        g2e = np.linalg.inv(e2g)

        # galactic to equatorial
        g2q = np.dot(e2q, g2e)

        # equatorial to ecliptic
        q2e = np.linalg.inv(e2q)

        # equatorial to galactic
        q2g = np.dot(e2g, q2e)

        if coord_norm == ("E", "G"):
            matconv = e2g
        elif coord_norm == ("G", "E"):
            matconv = g2e
        elif coord_norm == ("E", "C"):
            matconv = e2q
        elif coord_norm == ("C", "E"):
            matconv = q2e
        elif coord_norm == ("C", "G"):
            matconv = q2g
        elif coord_norm == ("G", "C"):
            matconv = g2q
        else:
            raise ValueError("Wrong coord transform :", coord_norm)
        do_conv = True

    return matconv, do_conv, coord_norm

def rotateVector(rotmat, vec, vy=None, vz=None, do_rot=True):
    if vy is None and vz is None:
        if do_rot:
            return np.tensordot(rotmat, vec, axes=(1, 0))
        else:
            return vec
    elif vy is not None and vz is not None:
        if do_rot:
            return np.tensordot(rotmat, np.array([vec, vy, vz]), axes=(1, 0))
        else:
            return vec, vy, vz
    else:
        raise TypeError("You must give either vec only or vec, vy " "and vz parameters")


def rotateDirection(rotmat, theta, phi=None, do_rot=True, lonlat=False):
    vx, vy, vz = rotateVector(rotmat, dir2vec(theta, phi, lonlat=lonlat), do_rot=do_rot)
    return vec2dir(vx, vy, vz, lonlat=lonlat)



#----------------CLASS--------------#
class Rotator(object):

    ErrMessWrongPar = (
        "rot and coord must be single elements or " "sequence of same size."
    )

    def __init__(self, rot=None, coord=None, inv=None, deg=True, eulertype="ZYX"):
        rot_is_seq = hasattr(rot, "__len__") and hasattr(rot[0], "__len__")
        coord_is_seq = (
            hasattr(coord, "__len__")
            and hasattr(coord[0], "__len__")
            and type(coord[0]) is not str
        )
        if rot_is_seq and coord_is_seq:
            if len(rot) != len(coord):
                raise ValueError(Rotator.ErrMessWrongPar)
            else:
                rots = rot
                coords = coord
        elif (rot_is_seq or coord_is_seq) and (rot is not None and coord is not None):
            raise ValueError(Rotator.ErrMessWrongPar)
        else:
            rots = [rot]
            coords = [coord]
        inv_is_seq = hasattr(inv, "__len__")
        if inv_is_seq:
            if len(inv) != len(rots):
                raise ValueError("inv must have same length as rot and/or coord")
            invs = inv
        else:
            invs = [inv] * len(rots)
        # check the argument and normalize them
        if eulertype in ["ZYX", "X", "Y"]:
            self._eultype = eulertype
        else:
            self._eultype = "ZYX"
        self._rots = []
        self._coords = []
        self._invs = []
        for r, c, i in zip(rots, coords, invs):
            rn = normalise_rot(r, deg=deg)
            #            if self._eultype in ['X','Y']:
            #                rn[1] = -rn[1]
            cn = normalise_coord(c)
            self._rots.append(rn)  # append(rn) or insert(0, rn) ?
            self._coords.append(cn)  # append(cn) or insert(0, cn) ?
            self._invs.append(bool(i))
        if not self.consistent:
            warnings.warn(
                "The chain of coord system rotations is not consistent",
                category=ConsistencyWarning,
            )
        self._update_matrix()

    def _update_matrix(self):
        self._matrix = np.identity(3)
        self._do_rotation = False
        for r, c, i in zip(self._rots, self._coords, self._invs):
            rotmat, do_rot, rotnorm = get_rotation_matrix(r, eulertype=self._eultype)
            convmat, do_conv, coordnorm = get_coordconv_matrix(c)
            r = np.dot(rotmat, convmat)
            if i:
                r = r.T
            self._matrix = np.dot(self._matrix, r)
            self._do_rotation = self._do_rotation or (do_rot or do_conv)

    def _is_coords_consistent(self):
        for c, i in zip(self._coords, self._invs):
            break
        for cnext, inext in zip(self._coords[1:], self._invs[1:]):
            if c[i] != cnext[not inext]:
                return False
            c, i = cnext, inext
        return True

    consistent = property(
        _is_coords_consistent, doc="consistency of the coords transform chain"
    )

    def __eq__(self, a):
        if type(a) is not type(self):
            return False
        # compare the _rots
        v = [np.allclose(x, y, rtol=0, atol=1e-15) for x, y in zip(self._rots, a._rots)]
        return (
            np.array(v).all()
            and (self._coords == a._coords)
            and (self._invs == a._invs)
        )

    def __call__(self, *args, **kwds):
        if kwds.pop("inv", False):
            m = self._matrix.T
        else:
            m = self._matrix
        lonlat = kwds.pop("lonlat", False)
        if len(args) == 1:
            arg = args[0]
            if not hasattr(arg, "__len__") or len(arg) < 2 or len(arg) > 3:
                raise TypeError("Argument must be a sequence of 2 or 3 " "elements")
            if len(arg) == 2:
                return rotateDirection(
                    m, arg[0], arg[1], self._do_rotation, lonlat=lonlat
                )
            else:
                return rotateVector(m, arg[0], arg[1], arg[2], self._do_rotation)
        elif len(args) == 2:
            return rotateDirection(
                m, args[0], args[1], self._do_rotation, lonlat=lonlat
            )
        elif len(args) == 3:
            return rotateVector(m, args[0], args[1], args[2], self._do_rotation)
        else:
            raise TypeError("Either 1, 2 or 3 arguments accepted")

    def __mul__(self, a):
        """Composition of rotation.
        """
        if not isinstance(a, Rotator):
            raise TypeError(
                "A Rotator can only multiply another Rotator "
                "(composition of rotations)"
            )
        rots = self._rots + a._rots
        coords = self._coords + a._coords
        invs = self._invs + a._invs
        return Rotator(rot=rots, coord=coords, inv=invs, deg=False)

    def __rmul__(self, b):
        if not isinstance(b, Rotator):
            raise TypeError(
                "A Rotator can only be multiplied by another Rotator "
                "(composition of rotations)"
            )
        rots = b._rots + self._rots
        coords = b._coords + self._coords
        invs = self._invs + b._invs
        return Rotator(rot=rots, coord=coords, inv=invs, deg=False)

    def __nonzero__(self):
        return self._do_rotation

    def get_inverse(self):
        rots = self._rots[::-1]
        coords = self._coords[::-1]
        invs = [not i for i in self._invs[::-1]]
        return Rotator(rot=rots, coord=coords, inv=invs, deg=False)

    # I = property(get_inverse,doc='Return a new rotator representing the '
    #             'inverse rotation')

    def I(self, *args, **kwds):
        """Rotate the given vector or direction using the inverse matrix.
        rot.I(vec) <==> rot(vec,inv=True)
        """
        kwds["inv"] = True
        return self.__call__(*args, **kwds)

    @property
    def mat(self):
        """The matrix representing the rotation.
        """
        return np.matrix(self._matrix)

    @property
    def coordin(self):
        """The input coordinate system.
        """
        if not self.consistent:
            return None
        for c, i in zip(self._coords, self._invs):
            pass
        return c[i]

    @property
    def coordout(self):
        """The output coordinate system.
        """
        if not self.consistent:
            return None
        for c, i in zip(self._coords, self._invs):
            pass
        return c[not i]

    @property
    def coordinstr(self):
        """The input coordinate system in str.
        """
        return coordname.get(self.coordin, "")

    @property
    def coordoutstr(self):
        """The output coordinate system in str.
        """
        return coordname.get(self.coordout, "")

    @property
    def rots(self):
        """The sequence of rots defining the rotation.
        """
        return self._rots

    @property
    def coords(self):
        """The sequence of coords defining the rotation.
        """
        return self._coords

    def do_rot(self, i):
        """Returns True if rotation is not (close to) identity.
        """
        return not np.allclose(self.rots[i], np.zeros(3), rtol=0.0, atol=1.0e-15)

    def angle_ref(self, *args, **kwds):
        R = self
        lonlat = kwds.get("lonlat", False)
        inv = kwds.get("inv", False)
        if len(args) == 1:
            arg = args[0]
            if not hasattr(arg, "__len__") or len(arg) < 2 or len(arg) > 3:
                raise TypeError("Argument must be a sequence of 2 or 3 " "elements")
            if len(arg) == 2:
                v = dir2vec(arg[0], arg[1], lonlat=lonlat)
            else:
                v = arg
        elif len(args) == 2:
            v = dir2vec(args[0], args[1], lonlat=lonlat)
        elif len(args) == 3:
            v = args
        else:
            raise TypeError("Either 1, 2 or 3 arguments accepted")
        vp = R(v, inv=inv)
        north_pole = R([0.0, 0.0, 1.0], inv=inv)
        sinalpha = north_pole[0] * vp[1] - north_pole[1] * vp[0]
        cosalpha = north_pole[2] - vp[2] * np.dot(north_pole, vp)
        return np.arctan2(sinalpha, cosalpha)

    # def rotate_alm(self, alm, lmax=None, mmax=None):

    #     rotated_alm = alm.copy()  # rotate_alm works inplace
    #     rotate_alm(rotated_alm, matrix=self.mat, lmax=lmax, mmax=mmax)
    #     return rotated_alm

    # def rotate_map_alms(self, m, use_pixel_weights=True, lmax=None, mmax=None):
    #     alm = map2alm(
    #         m, use_pixel_weights=use_pixel_weights, lmax=lmax, mmax=mmax
    #     )
    #     rotated_alm = self.rotate_alm(alm, lmax=lmax, mmax=mmax)
    #     return sphtfunc.alm2map(
    #         rotated_alm, lmax=lmax, mmax=mmax, nside=pixelfunc.get_nside(m)
    #     )

    # def rotate_map_pixel(self, m):
    #     if maptype(m) == 0:  # a single map is converted to a list
    #         m = [m]
    #     npix = len(m[0])
    #     nside = npix2nside(npix)
    #     theta_pix_center, phi_pix_center = pix2ang(
    #         nside=nside, ipix=np.arange(npix)
    #     )

    #     # Rotate the pixels center of the new reference frame to the original frame
    #     theta_pix_center_rot, phi_pix_center_rot = self.I(
    #         theta_pix_center, phi_pix_center
    #     )

    #     # Interpolate the original map to the pixels centers in the new ref frame
    #     m_rotated = [
    #         get_interp_val(each, theta_pix_center_rot, phi_pix_center_rot)
    #         for each in m
    #     ]

    #     # Rotate polarization
    #     if len(m_rotated) > 1:
    #         # Create a complex map from QU  and apply the rotation in psi due to the rotation
    #         # Slice from the end of the array so that it works both for QU and IQU
    #         L_map = (m_rotated[-2] + m_rotated[-1] * 1j) * np.exp(
    #             1j * 2 * self.angle_ref(theta_pix_center_rot, phi_pix_center_rot)
    #         )

    #         # Overwrite the Q and U maps with the correct values
    #         m_rotated[-2] = np.real(L_map)
    #         m_rotated[-1] = np.imag(L_map)
    #     else:
    #         m_rotated = m_rotated[0]

    #     return m_rotated

    def __repr__(self):
        return (
            "[ "
            + ", ".join([str(self._coords), str(self._rots), str(self._invs)])
            + " ]"
        )

    __str__ = __repr__


#--------------------------------------------------------#
#                                                        #
#                                                        #
#                                                        #
#--------------------------------------------------------#
# MAX_NSIDE = (
#     8192
# )  # The maximum nside up to which most operations (e.g. map2alm) will work


# def check_max_nside(nside):
#     if nside > MAX_NSIDE:
#         raise ValueError(
#             "nside {nside} of map cannot be larger than "
#             "MAX_NSIDE {max_nside}".format(nside=nside, max_nside=MAX_NSIDE)
#         )

#     return 0

    
# def map2alm(
#     maps,
#     lmax=None,
#     mmax=None,
#     iter=3,
#     pol=True,
#     use_weights=False,
#     datapath=None,
#     gal_cut=0,
#     use_pixel_weights=False,
# ):
#     maps = ma_to_array(maps)
#     info = maptype(maps)
#     nside = pixelfunc.get_nside(maps)
#     check_max_nside(nside)

#     if use_pixel_weights:
#         if use_weights:
#             raise RuntimeError("Either use pixel or ring weights")
#         with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
#             "remote_timeout", 30
#         ):
#             pixel_weights_filename = data.get_pkg_data_filename(
#                 "full_weights/healpix_full_weights_nside_%04d.fits" % nside,
#                 package="healpy",
#             )
#     else:
#         pixel_weights_filename = None

#     if pol or info in (0, 1):
#         alms = _sphtools.map2alm(
#             maps,
#             niter=iter,
#             datapath=datapath,
#             use_weights=use_weights,
#             lmax=lmax,
#             mmax=mmax,
#             gal_cut=gal_cut,
#             pixel_weights_filename=pixel_weights_filename,
#         )
#     else:
#         # info >= 2 and pol is False : spin 0 spht for each map
#         alms = [
#             _sphtools.map2alm(
#                 mm,
#                 niter=iter,
#                 datapath=datapath,
#                 use_weights=use_weights,
#                 lmax=lmax,
#                 mmax=mmax,
#                 gal_cut=gal_cut,
#                 pixel_weights_filename=pixel_weights_filename,
#             )
#             for mm in maps
#         ]
#     return np.array(alms)

