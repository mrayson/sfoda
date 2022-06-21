"""
Fast implementation of xyz-enu transformation.

This follows http://wiki.cython.org/tutorials/numpy#Addingtypes

Reference: 'ADCP Coordinate Transformation', a booklet
produced by Teledyne RD Instruments and available from their
web site.  Version P/N 951-6079-00 (January 2008).

"""

cdef extern from "math.h":
    double cos(double)
    double sin(double)
    double tan(double)
    double atan(double)

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t


@cython.boundscheck(False)
def _heading_rotate(np.ndarray[DTYPE_t, ndim=3] xyz,
                    np.ndarray[DTYPE_t, ndim=1] heading):
    """
    Inner calculation for the case of heading-only;
    various combinations of dimensions are handled by a wrapper.
    Masks are handled by _heading_rotate_m.

    We are now assuming that there is a heading value
    for each element in the first dimension of xyz.
    """
    cdef np.ndarray[DTYPE_t, ndim=3] enu
    cdef unsigned int i, j
    cdef unsigned int nt, nd
    cdef double to_rad = np.pi / 180.0
    cdef double h, ch, sh
    enu = xyz.copy()
    nt, nd = xyz.shape[0], xyz.shape[1]  # xyz.shape[:2] did not work

    for i in range(nt):
        h = heading[i] * to_rad
        ch = cos(h)
        sh = sin(h)
        for j in range(nd):
            enu[i,j,0] = xyz[i,j,0] * ch + xyz[i,j,1] * sh
            enu[i,j,1] = xyz[i,j,0] * (-sh) + xyz[i,j,1] * ch

    return enu

@cython.boundscheck(False)
def _heading_rotate_m(np.ndarray[DTYPE_t, ndim=3] xyz,
                      np.ndarray[DTYPE_t, ndim=1] heading,
                      np.ndarray[np.int8_t, ndim=3] xmask,
                      np.ndarray[np.int8_t, ndim=1] hmask):
    """
    Inner calculation for the case of heading-only, with full
    masking support, including output of a resulting mask.

    """
    cdef np.ndarray[DTYPE_t, ndim=3] enu
    cdef np.ndarray[np.int8_t, ndim=3] outmask
    cdef unsigned int i, j
    cdef unsigned int nt, nd
    cdef double to_rad = np.pi / 180.0
    cdef double h, ch, sh
    enu = xyz.copy()
    outmask = xmask.copy()
    nt, nd = xyz.shape[0], xyz.shape[1]  # xyz.shape[:2] did not work

    for i in range(nt):
        if hmask[i]:
            for j in range(nd):
                outmask[i, j, 0] = outmask[i, j, 1] = 1
            continue
        h = heading[i] * to_rad
        ch = cos(h)
        sh = sin(h)
        for j in range(nd):
            if xmask[i, j, 0] or xmask[i, j, 1]:
                outmask[i, j, 0] = outmask[i, j, 1] = 1
                continue
            enu[i,j,0] = xyz[i,j,0] * ch + xyz[i,j,1] * sh
            enu[i,j,1] = xyz[i,j,0] * (-sh) + xyz[i,j,1] * ch

    return enu, outmask

@cython.boundscheck(False)
def _hpr_rotate(np.ndarray[DTYPE_t, ndim=3] xyz,
                    np.ndarray[DTYPE_t, ndim=1] heading,
                    np.ndarray[DTYPE_t, ndim=1] pitch,
                    np.ndarray[DTYPE_t, ndim=1] roll,
                    orientation,
                    gimbal):
    """
    Instrument coordinates xyz transformed by heading, pitch, roll

    Version for the case with no mask.
    """
    cdef np.ndarray[DTYPE_t, ndim=3] enu
    cdef unsigned int i, j
    cdef unsigned int nt, nd
    cdef double to_rad = np.pi / 180.0
    cdef double h, ch, sh
    cdef double p, cp, sp
    cdef double r, cr, sr
    cdef double m00, m01, m02, m10, m11, m12, m20, m21, m22
    cdef int adj_pitch = (not gimbal)
    cdef int orient_up = (not (orientation == "down"))
    enu = xyz.copy()
    nt, nd = xyz.shape[0], xyz.shape[1]  # xyz.shape[:2] did not work

    for i in range(nt):
        h = heading[i] * to_rad
        ch = cos(h)
        sh = sin(h)
        r = roll[i] * to_rad
        if orient_up:    # roll by additional 180 degrees to flip it up
            cr = -cos(r)
            sr = -sin(r)
        else:
            cr = cos(r)
            sr = sin(r)
        p = pitch[i] * to_rad
        if adj_pitch:
            p = atan(tan(p) * cos(r))
        cp = cos(p)
        sp = sin(p)

        m00 = ch*cr + sh*sp*sr
        m01 = sh*cp
        m02 = ch*sr - sh*sp*cr
        m10 = -sh*cr + ch*sp*sr
        m11 = ch*cp
        m12 = -sh*sr - ch*sp*cr
        m20 = -cp*sr
        m21 = sp
        m22 = cp*cr

        for j in range(nd):
            enu[i,j,0] = xyz[i,j,0] * m00 + xyz[i,j,1] * m01 + xyz[i,j,2] * m02
            enu[i,j,1] = xyz[i,j,0] * m10 + xyz[i,j,1] * m11 + xyz[i,j,2] * m12
            enu[i,j,2] = xyz[i,j,0] * m20 + xyz[i,j,1] * m21 + xyz[i,j,2] * m22

    return enu

@cython.boundscheck(False)
def _hpr_rotate_m(np.ndarray[DTYPE_t, ndim=3] xyz,
                    np.ndarray[DTYPE_t, ndim=1] heading,
                    np.ndarray[DTYPE_t, ndim=1] pitch,
                    np.ndarray[DTYPE_t, ndim=1] roll,
                    np.ndarray[np.int8_t, ndim=3] xmask,
                    np.ndarray[np.int8_t, ndim=1] hmask,
                    orientation,
                    gimbal):
    """
    Instrument coordinates xyz transformed by heading, pitch, roll

    Version for the case with no mask.
    """
    cdef np.ndarray[DTYPE_t, ndim=3] enu
    cdef np.ndarray[np.int8_t, ndim=3] outmask
    cdef unsigned int i, j
    cdef unsigned int nt, nd
    cdef double to_rad = np.pi / 180.0
    cdef double h, ch, sh
    cdef double p, cp, sp
    cdef double r, cr, sr
    cdef double m00, m01, m02, m10, m11, m12, m20, m21, m22
    cdef int adj_pitch = (not gimbal)
    cdef int orient_up = (not (orientation == "down"))
    enu = xyz.copy()
    outmask = xmask.copy()
    nt, nd = xyz.shape[0], xyz.shape[1]  # xyz.shape[:2] did not work

    for i in range(nt):
        if hmask[i]:
            for j in range(nd):
                outmask[i, j, 0] = outmask[i, j, 1] = 1
            continue
        h = heading[i] * to_rad
        ch = cos(h)
        sh = sin(h)
        r = roll[i] * to_rad
        if orient_up:    # roll by additional 180 degrees to flip it up
            cr = -cos(r)
            sr = -sin(r)
        else:
            cr = cos(r)
            sr = sin(r)
        p = pitch[i] * to_rad
        if adj_pitch:
            p = atan(tan(p) * cos(r))
        cp = cos(p)
        sp = sin(p)

        m00 = ch*cr + sh*sp*sr
        m01 = sh*cp
        m02 = ch*sr - sh*sp*cr
        m10 = -sh*cr + ch*sp*sr
        m11 = ch*cp
        m12 = -sh*sr - ch*sp*cr
        m20 = -cp*sr
        m21 = sp
        m22 = cp*cr

        for j in range(nd):
            if xmask[i,j,0] or xmask[i,j,1] or xmask[i,j,2]:
                outmask[i,j,0] = outmask[i,j,1] = outmask[i,j,2] = 1
                continue
            enu[i,j,0] = xyz[i,j,0] * m00 + xyz[i,j,1] * m01 + xyz[i,j,2] * m02
            enu[i,j,1] = xyz[i,j,0] * m10 + xyz[i,j,1] * m11 + xyz[i,j,2] * m12
            enu[i,j,2] = xyz[i,j,0] * m20 + xyz[i,j,1] * m21 + xyz[i,j,2] * m22

    return enu, outmask



