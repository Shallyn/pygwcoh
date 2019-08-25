#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 20:00:24 2019

@author: drizl
"""

import numpy as np

#[ 1.52911759,  0.78550497,  1.57079633,  0.05103658,  3.09055608]
#[ 0.        ,  0.78539816,  1.61988371,  0.78539816,  0.78539816]

def alert(cond, msg):
    if cond:
        raise ValueError(msg)

def isqrt(x):
    return np.sqrt(x + 0.5)

def nside2order(nside):
    if not isinstance(nside, int):
        return False
    x = np.log2(nside)
    if x - int(x) > 1e-2:
        return False
    else:
        return int(x)

def fmodulo(v1, v2):
    if (v1>=0):
        if v1 < v2:
            return v1
        else:
            return np.mod(v1, v2)
    tmp=np.mod(v1, v2)
    if tmp == v2:
        return 0
    else:
        return tmp

def nside2npix(nside):
    return 12 * nside * nside

def npix2nside(npix):
    return int(np.sqrt(npix / 12))

def loc2pix(z, phi, sth, have_sth, Hpix):
    za = abs(z)
    tt = fmodulo(phi*2 / np.pi, 4.0) # in [0,4)
    alert((tt < 0) or (tt >= 4), "Error: tt:{}".format(tt))
    nside = Hpix.nside
    order = Hpix.order
    ncap = Hpix.ncap
    npix = Hpix.npix
    if Hpix.scheme == 'ring':
        if za<=2/3:
            nl4 = 4*nside
            temp1 = nside*(0.5+tt)
            temp2 = nside*z*0.75
            jp = int(temp1-temp2)
            jm = int(temp1+temp2)
            ir = nside + 1 + jp - jm
            kshift = 1 - (ir & 1)
            t1 = jp+jm-nside+kshift+1+nl4+nl4
            if order > 0:
                ip = (t1>>1)&(nl4-1)
            else:
                ip = ((t1>>1)%nl4)
            return ncap + (ir-1)*nl4 + ip
        else:
            tp = tt-int(tt)
            if (za<0.99) and ( not have_sth):
                tmp = nside*np.sqrt(3*(1-za))
            else:
                tmp = nside*sth/np.sqrt((1.+za)/3.)
            jp = int(tp*tmp)
            jm = int((1.0-tp)*tmp)
            ir = jp+jm+1
            ip = int(tt*ir)
            alert((ip<0) or (ip>=4*ir),"must not happen: ip = {}".format(ip))
            if z > 0:
                return 2*ir*(ir-1) + ip
            else:
                return npix - 2*ir*(ir+1) + ip
    else:
        raise TypeError('Unsupported scheme: nest')
        # if za <= 2/3:
        #     temp1 = nside*(0.5+tt)
        #     temp2 = nside*(z*0.75)
        #     jp = int(temp1-temp2)
        #     jm = int(temp1+temp2)
        #     ifp = jp >> order
        #     ifm = jm >> order
        #     if ifp==ifm:
        #         face_num = (ifp|4)
        #     elif ifp<ifm:
        #         face_num = ifp
        #     else:
        #         face_num = ifm+8
        #     ix = jm & (nside-1)
        #     iy = nside - (jp & (nside-1)) - 1
        #     return xyf2nest(ix,iy,face_num)
        # else:
        #     ntt = min(3,int(tt))
        #     tp = tt-ntt
        #     if (za<0.99) and ( not have_sth):
        #         tmp = nside*np.sqrt(3*(1-za))
        #     else:
        #         tmp = nside*sth/np.sqrt((1.+za)/3.)
        #     jp = int(tp*tmp)
        #     jm = int((1.0-tp)*tmp)
        #     jp=min(jp,nside-1)
        #     jm=min(jm,nside-1)
        #     if z>=0:
        #         xyf2nest(nside-jm -1,nside-jp-1,ntt, order)
        #     else:
        #         xyf2nest(jp,jm,ntt+8, order)


def vec2pix(nside, x, y, z, nest = False):
    if nest:
        scheme = 'nest'
    else:
        scheme = 'ring'
    Hpix = HealPix(nside, scheme)
    if not isinstance(x, list) and not isinstance(x, np.ndarray):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        
    N = len(x)
    pix = []

    for i in range(N):
        xi = x[i]
        yi = y[i]
        zi = z[i]
        xl = 1./np.sqrt(xi**2 + yi**2 + zi**2)
        phi = np.arctan2(yi,xi);
        nz = zi*xl
        if (abs(nz)>0.99):
            ret = loc2pix (nz,phi,np.sqrt(xi*xi+yi*yi)*xl,True, Hpix)
        else:
            ret = loc2pix (nz,phi,0,False, Hpix)
        pix.append(ret)
    return np.array(pix)
    

def pix2ang(nside, ipix):
    order = nside2order(nside)
    if not order and order != 0:
        raise ValueError('nside should be power of 2')
    nside = 2 ** order
    npface = nside**2
    ncap = (npface - nside) << 1
    npix = 12 * npface
    fact2 = 4 / npix
    fact1 = (nside << 1) * fact2
    
    list_theta = []
    list_phi = []
    if not isinstance(ipix, list) and not isinstance(ipix, np.ndarray):
        ipix = [ipix]
    for pix in ipix:
        have_sth = False
        sth = 0
        if pix < ncap: #counted from North pole
            iring = (1 + int(isqrt(1 + 2*pix))) >> 1 
            iphi = (pix + 1) - 2 * iring * (iring - 1)
            tmp = iring ** 2 * fact2
            z = 1 - tmp
            if z > 0.99:
                have_sth = True
                sth = np.sqrt(tmp*(2 - tmp))
            phi = (iphi - 0.5) * np.pi / iring / 2
        elif pix < (npix - ncap): # Equatorial region
            nl4 = 4 * nside
            ip = pix - ncap
            if order >=0:
                tmp = ip >> int(order + 2)
            else:
                tmp = ip / nl4
            iring = tmp + nside
            iphi = ip - nl4 * tmp + 1
            if (iring + nside) & 1:
                fodd = 1
            else:
                fodd = 0.5
            z = (2 * nside - iring) * fact1
            phi = (iphi - fodd) * np.pi * 0.75 * fact1
        else: #  South Polar cap
            ip = npix - pix
            iring = (1 + int(isqrt(2 * ip - 1))) >> 1
            iphi = 4 * iring + 1 - (ip - 2 * iring*(iring - 1))
            tmp = (iring * iring) * fact2
            z = tmp - 1
            if z < -0.99:
                sth = np.sqrt(tmp * (2 - tmp))
                have_sth = True
            phi = (iphi - 0.5) * np.pi / iring / 2
        if have_sth:
            ret_theta = np.arctan2(sth , z)
            ret_phi = phi
        else:
            ret_theta = np.arccos(z)
            ret_phi = phi
        if ret_theta < 0:
            ret_theta = np.pi + ret_theta
        list_theta.append(ret_theta)
        list_phi.append(ret_phi)
    return np.array(list_theta), np.array(list_phi)

class HealPix(object):
    def __init__(self, nside, scheme):
        self.order = nside2order(nside)
        if not self.order and self.order != 0:
            raise ValueError('nside should be power of 2')
        self.nside = 2 ** self.order
        self.npface = nside**2
        self.ncap = (self.npface - nside) << 1
        self.npix = 12 * self.npface
        self.fact2 = 4 / self.npix
        self.fact1 = (nside << 1) * self.fact2
        self.scheme = scheme
    
    def pix2ang(self):
        list_theta = []
        list_phi = []
        for pix in ipix:
            have_sth = False
            sth = 0
            if pix < self.ncap:
                iring = (1 + int(isqrt(1 + 2*pix))) >> 1
                iphi = (pix + 1) - 2 * iring * (iring - 1)
                tmp = iring ** 2 * self.fact2
                z = 1 - tmp
                if z > 0.99:
                    have_sth = True
                    sth = np.sqrt(tmp*(2 - tmp))
                phi = (iphi - 0.5) * np.pi / iring / 2
            elif pix < (self.npix - self.ncap):
                nl4 = 4 * nside
                ip = pix - self.ncap
                if self.order >=0:
                    tmp = ip >> int(self.order + 2)
                else:
                    tmp = ip / nl4
                iring = tmp + nside
                iphi = ip - nl4 * tmp + 1
                if (iring + nside) & 1:
                    fodd = 1
                else:
                    fodd = 0.5
                z = (2 * nside - iring) * self.fact1
                phi = (iphi - fodd) * np.pi * 0.75 * self.fact1
            else:
                ip = self.npix - pix
                iring = (1 + int(isqrt(2 * ip - 1))) >> 1
                iphi = 4 * iring + 1 - (ip - 2 * iring*(iring - 1))
                tmp = (iring * iring) * self.fact2
                z = tmp - 1
                if z < -0.99:
                    sth = np.sqrt(tmp * (2 - tmp))
                    have_sth = True
                phi = (iphi - 0.5) * np.pi / iring / 2
            if have_sth:
                ret_theta = np.arctan(sth / z)
                ret_phi = phi
            else:
                ret_theta = np.arccos(z)
                ret_phi = phi
            if ret_theta < 0:
                ret_theta = np.pi + ret_theta
            list_theta.append(ret_theta)
            list_phi.append(ret_phi)
        return list_theta, list_phi

        

        