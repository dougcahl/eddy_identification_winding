"""Numba-compiled kernels for the hot loops (stream2_dc integration and the
winding-angle scan). Compiled once per machine (cache=True) via LLVM —
the same code generation path a Rust/C implementation would use, with
strict IEEE float64 semantics (fastmath is off), so results are
bit-identical to the pure-Python fallbacks.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def stream2_kernel(x, y, u, v, startx, starty, dr, N):
    """Port of stream2_dc.m. Returns (X, Y, status); status 1 means the
    trajectory left the grid (MATLAB errors there; caller raises)."""
    X = np.full(N, np.nan)
    Y = np.full(N, np.nan)
    X[0] = startx
    Y[0] = starty

    nrow, ncol = x.shape

    # nearest grid node to the start
    best = 1e300
    row = 0
    col = 0
    for r in range(nrow):
        for c in range(ncol):
            d = (x[r, c] - startx) ** 2 + (y[r, c] - starty) ** 2
            if d < best:
                best = d
                row = r
                col = c

    sx = startx
    sy = starty

    for i in range(1, N):
        x1 = x[row, col]; y1 = y[row, col]
        u1 = u[row, col]; v1 = v[row, col]
        row2 = row + 1 if sx > x1 else row - 1
        col2 = col + 1 if sy > y1 else col - 1
        if row2 < 0 or row2 >= nrow or col2 < 0 or col2 >= ncol:
            return X, Y, 1  # left the grid

        x2 = x[row2, col]; y2 = y[row2, col]
        x3 = x[row, col2]; y3 = y[row, col2]
        x4 = x[row2, col2]; y4 = y[row2, col2]
        u2 = u[row2, col]; v2 = v[row2, col]
        u3 = u[row, col2]; v3 = v[row, col2]
        u4 = u[row2, col2]; v4 = v[row2, col2]

        r1 = ((sx - x1) ** 2 + (sy - y1) ** 2) ** 0.5
        r2 = ((sx - x2) ** 2 + (sy - y2) ** 2) ** 0.5
        r3 = ((sx - x3) ** 2 + (sy - y3) ** 2) ** 0.5
        r4 = ((sx - x4) ** 2 + (sy - y4) ** 2) ** 0.5

        if r1 == 0:
            u0 = u1; v0 = v1
        elif r2 == 0:
            u0 = u2; v0 = v2
        elif r3 == 0:
            u0 = u3; v0 = v3
        elif r4 == 0:
            u0 = u4; v0 = v4
        else:
            n1 = u1 == u1
            n2 = u2 == u2
            n3 = u3 == u3
            n4 = u4 == u4
            if not (n1 or n2 or n3 or n4):
                break
            norm = 0.0
            u0 = 0.0
            v0 = 0.0
            if n1:
                norm += 1.0 / r1; u0 += u1 / r1; v0 += v1 / r1
            if n2:
                norm += 1.0 / r2; u0 += u2 / r2; v0 += v2 / r2
            if n3:
                norm += 1.0 / r3; u0 += u3 / r3; v0 += v3 / r3
            if n4:
                norm += 1.0 / r4; u0 += u4 / r4; v0 += v4 / r4
            u0 = u0 / norm
            v0 = v0 / norm

        U = (u0 * u0 + v0 * v0) ** 0.5
        if U == 0 or U != U:
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

    return X, Y, 0


@njit(cache=True)
def winding_kernel(sx, sy, winding_thres, d_thres, baddir_thres):
    """Winding-angle scan of one streamline (same update order as
    eddy_subroutine.m). Returns (eddy, ii, winding, dir1)."""
    eddy = 0
    winding = 0.0
    dir1 = 0.0
    baddir = 0.0
    d1 = np.inf
    sx0 = sx[0]
    sy0 = sy[0]
    ang0 = np.degrees(np.arctan2(sy[1] - sy0, sx[1] - sx0))
    n = len(sx)
    ii = 1
    for ii in range(1, n - 1):
        dx = sx[ii + 1] - sx[ii]
        dy = sy[ii + 1] - sy[ii]
        ang1 = np.degrees(np.arctan2(dy, dx))
        ang_diff = ang1 - ang0
        if ang_diff > 180:
            ang_diff -= 360
        if ang_diff < -180:
            ang_diff += 360
        if dir1 == 0:
            dir1 = float((ang_diff > 0) - (ang_diff < 0))

        if ii > 1:
            newdir = float((ang_diff > 0) - (ang_diff < 0))
            if newdir != dir1:
                baddir += ang_diff
                if abs(baddir) > baddir_thres:
                    break
            else:
                baddir = 0.0

        ang0 = ang1
        winding += ang_diff

        if abs(winding) > winding_thres:
            dx = sx[ii] - sx0
            dy = sy[ii] - sy0
            d = (dx * dx + dy * dy) ** 0.5
            if eddy == 1 and d > d1:
                break
            if d < d_thres:
                eddy = 1
                d1 = d

    return eddy, ii, winding, dir1
