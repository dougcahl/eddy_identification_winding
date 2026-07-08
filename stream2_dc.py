"""Port of stream2_dc.m — custom streamline integrator (D. Cahl, USC 2023).

Unlike matplotlib's streamplot (used by the v1.0 Python code), this integrator
has NO visited-cell mask, so circling streamlines keep going instead of being
terminated when they re-enter a cell. This is what allows the stricter
winding/closure thresholds of the newer MATLAB repo to work.

Velocity at the current position is an inverse-distance weighted average of
the 4 surrounding grid nodes (NaN nodes are skipped); the step is a forward
Euler step of dr * cell-size along the local flow direction.
"""

import numpy as np

try:  # numba-compiled kernel (LLVM); bit-identical, ~100x faster
    from kernels_numba import stream2_kernel as _kernel
except Exception:
    _kernel = None


def stream2_dc(x, y, u, v, startx, starty, dr=0.1, N=10000):
    if _kernel is not None:
        X, Y, status = _kernel(
            np.ascontiguousarray(x, dtype=np.float64),
            np.ascontiguousarray(y, dtype=np.float64),
            np.ascontiguousarray(u, dtype=np.float64),
            np.ascontiguousarray(v, dtype=np.float64),
            float(startx), float(starty), float(dr), int(N))
        if status == 1:
            raise IndexError('streamline left the grid')
        return X, Y
    return stream2_dc_py(x, y, u, v, startx, starty, dr, N)


def stream2_dc_py(x, y, u, v, startx, starty, dr=0.1, N=10000):
    """Integrate a single streamline from (startx, starty).

    x, y, u, v : 2-D arrays, x/y a (possibly slightly irregular) plaid grid.
                 Axis 0 ("row") varies in x, axis 1 ("col") varies in y,
                 matching the MATLAB [lats,lons]=meshgrid(lat,lon) layout.
    dr         : step size as a fraction of a cell (MATLAB default 0.1)
    N          : maximum number of vertices (MATLAB default 10000)

    Returns (X, Y) length-N arrays padded with NaN after the last point.
    Raises IndexError when the trajectory leaves the grid (the MATLAB version
    errors there too and the caller catches it).
    """
    X = np.full(N, np.nan)
    Y = np.full(N, np.nan)
    X[0] = startx
    Y[0] = starty

    nrow, ncol = x.shape
    d = (x - startx) ** 2 + (y - starty) ** 2
    row, col = np.unravel_index(np.argmin(d), d.shape)
    row = int(row)
    col = int(col)

    # plain nested lists: much faster scalar access than numpy indexing
    xL = x.tolist(); yL = y.tolist()
    uL = u.tolist(); vL = v.tolist()

    sx = float(startx)
    sy = float(starty)

    for i in range(1, N):
        x1 = xL[row][col]; y1 = yL[row][col]
        u1 = uL[row][col]; v1 = vL[row][col]
        row2 = row + 1 if sx > x1 else row - 1
        col2 = col + 1 if sy > y1 else col - 1
        if row2 < 0 or row2 >= nrow or col2 < 0 or col2 >= ncol:
            raise IndexError('streamline left the grid')

        x2 = xL[row2][col]; y2 = yL[row2][col]
        x3 = xL[row][col2]; y3 = yL[row][col2]
        x4 = xL[row2][col2]; y4 = yL[row2][col2]
        u2 = uL[row2][col]; v2 = vL[row2][col]
        u3 = uL[row][col2]; v3 = vL[row][col2]
        u4 = uL[row2][col2]; v4 = vL[row2][col2]

        r1 = ((sx - x1) ** 2 + (sy - y1) ** 2) ** 0.5
        r2 = ((sx - x2) ** 2 + (sy - y2) ** 2) ** 0.5
        r3 = ((sx - x3) ** 2 + (sy - y3) ** 2) ** 0.5
        r4 = ((sx - x4) ** 2 + (sy - y4) ** 2) ** 0.5

        if r1 == 0:
            u0, v0 = u1, v1
        elif r2 == 0:
            u0, v0 = u2, v2
        elif r3 == 0:
            u0, v0 = u3, v3
        elif r4 == 0:
            u0, v0 = u4, v4
        else:
            n1 = u1 == u1  # not NaN
            n2 = u2 == u2
            n3 = u3 == u3
            n4 = u4 == u4
            if not (n1 or n2 or n3 or n4):
                break
            norm = (n1 / r1) + (n2 / r2) + (n3 / r3) + (n4 / r4)
            u0 = ((u1 / r1 if n1 else 0.0) + (u2 / r2 if n2 else 0.0)
                  + (u3 / r3 if n3 else 0.0) + (u4 / r4 if n4 else 0.0)) / norm
            v0 = ((v1 / r1 if n1 else 0.0) + (v2 / r2 if n2 else 0.0)
                  + (v3 / r3 if n3 else 0.0) + (v4 / r4 if n4 else 0.0)) / norm

        U = (u0 ** 2 + v0 ** 2) ** 0.5
        if U == 0 or U != U:  # stalled flow: MATLAB degenerates to NaNs here
            break
        ddx = abs(x2 - x1)
        dt = ddx / U * dr

        sx = sx + dt * u0
        sy = sy + dt * v0
        X[i] = sx
        Y[i] = sy

        r1 = ((sx - x1) ** 2 + (sy - y1) ** 2) ** 0.5
        r2 = ((sx - x2) ** 2 + (sy - y2) ** 2) ** 0.5
        r3 = ((sx - x3) ** 2 + (sy - y3) ** 2) ** 0.5
        r4 = ((sx - x4) ** 2 + (sy - y4) ** 2) ** 0.5
        if r2 < r1 and r2 < r3 and r2 < r4:
            row = row2
        elif r3 < r1 and r3 < r2 and r3 < r4:
            col = col2
        elif r4 < r1 and r4 < r2 and r4 < r3:
            row = row2
            col = col2

    return X, Y
