"""
Fast 1D linear interpolation using numba

Note that this code will linearly extrapolate out of bounds points. 

Code "borrowed" from here: https://gist.github.com/shoyer/c8bb9e4d689b297a170296426a08aef7

Usage:
```
yout = interp1d_numbar(xout, xin, yin)
```

"""

from numba import guvectorize

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 


def _interp1d(xnew, xvals, yvals, ynew):
    i = 0
    N = len(xvals)
    if xnew[0] < xvals[0]:
        x_a = 0.0
        y_a = 0.0
        x_b = xvals[0]
        y_b = yvals[0]
    else:
        while xnew[0] >= xvals[i] and i < N:
            i += 1
        if xnew[0] == xvals[i]:
            ynew[0] = yvals[i]
            return
        if i == N:
            i = N-1
        x_a = xvals[i-1]
        y_a = yvals[i-1]
        x_b = xvals[i]
        y_b = yvals[i]
    slope = (xnew[0] - x_a)/(x_b - x_a)
    ynew[0] = slope * (y_b-y_a) + y_a
    return

  
interp1d_numba = guvectorize(
    ['float64[:], float64[:], float64[:], float64[:]'],
    "(),(n),(n) -> ()", nopython=True)(_interp1d)
