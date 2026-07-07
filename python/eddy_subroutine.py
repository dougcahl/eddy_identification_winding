"""Python port of the current MATLAB repo's eddy_subroutine.m
(dougcahl/eddy_identification_winding, D. Cahl, USC 2023).

Features ported from the MATLAB upgrade (missing from the v1.0 Python code):
  - stream2_dc integrator (no visited-cell termination), so the stricter
    thresholds (winding > 300 deg, closure < 10 km) work
  - Cahl et al. 2023 clustering: streamlines belong to the same eddy when a
    center lies inside the other's streamline polygon (inpolygon test) or the
    polygons overlap; falls back to a center-distance test
  - consistent UTM zone handling (geog2utm_nodisp / utm2ll ports)
  - ellipse fit with the sqrt(2*eig) axis scaling
  - per-streamline lat/lon saved to the .mat results (save_streams)
  - map figure with speed field, unit-vector quiver, ETOPO1 bathymetry
    contours, and land/coastline (from the ETOPO1 land mask)
  - results written as .mat with the same variable names as the MATLAB repo
"""

import math

import numpy as np
from matplotlib.path import Path
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
from scipy.io import savemat

from stream2_dc import stream2_dc
from utm_dc import geog2utm_nodisp, utm2ll, auto_zone
from etopo1_reader import read_etopo1

try:  # numba-compiled winding scan; bit-identical to the loop below
    from kernels_numba import winding_kernel as _winding_kernel
except Exception:
    _winding_kernel = None


def _inpolygon(px, py, polyx, polyy):
    """MATLAB inpolygon equivalent. px/py scalar or arrays -> bool or bool array."""
    poly = Path(np.column_stack([np.asarray(polyx).ravel(), np.asarray(polyy).ravel()]))
    pts = np.column_stack([np.atleast_1d(px).ravel(), np.atleast_1d(py).ravel()])
    res = poly.contains_points(pts)
    return bool(res[0]) if res.size == 1 else res


def _winding_py(stline_x, stline_y, winding_thres, d_thres, baddir_thres):
    """Pure-Python winding-angle scan (fallback when numba is unavailable);
    same update order as the new MATLAB code."""
    sxL = stline_x.tolist()
    syL = stline_y.tolist()
    eddy = 0
    winding = 0.0
    dir1 = 0.0
    baddir = 0.0
    d1 = math.inf
    sx0 = sxL[0]; sy0 = syL[0]
    ang0 = math.degrees(math.atan2(syL[1] - sy0, sxL[1] - sx0))
    ii = 1
    for ii in range(1, len(sxL) - 1):
        dx = sxL[ii + 1] - sxL[ii]
        dy = syL[ii + 1] - syL[ii]
        ang1 = math.degrees(math.atan2(dy, dx))
        ang_diff = ang1 - ang0
        if ang_diff > 180:
            ang_diff -= 360
        if ang_diff < -180:
            ang_diff += 360
        if dir1 == 0:
            dir1 = (ang_diff > 0) - (ang_diff < 0)

        if ii > 1:
            newdir = (ang_diff > 0) - (ang_diff < 0)
            if newdir != dir1:
                baddir += ang_diff
                if abs(baddir) > baddir_thres:
                    break
            else:
                baddir = 0.0

        ang0 = ang1
        winding += ang_diff

        if abs(winding) > winding_thres:
            dx = sxL[ii] - sx0
            dy = syL[ii] - sy0
            d = math.hypot(dx, dy)
            if eddy and d > d1:
                break
            if d < d_thres:
                eddy = 1
                d1 = d

    return eddy, ii, winding, dir1


def _base_map(ax, plt, params, lons, lats, u, v, lon_min, lon_max, lat_min, lat_max):
    """Speed field + quiver + bathymetry + land patch, mirroring the m_map figure."""
    if params['plt_quiver']:
        U = np.sqrt(u ** 2 + v ** 2)
        pc = ax.pcolormesh(lons, lats, U, shading='auto', zorder=1)
        with np.errstate(invalid='ignore'):
            ax.quiver(lons, lats, u / U, v / U, color='k', scale=55,
                      width=0.0012, zorder=3)
        cb = plt.colorbar(pc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label('Speed (m/s)')

    try:
        elev, elon, elat = read_etopo1(lon_min - 0.05, lon_max + 0.05,
                                       lat_min - 0.05, lat_max + 0.05)
        if params['bathy']:
            cs = ax.contour(elon, elat, elev, levels=sorted(params['levels']),
                            colors='k', linewidths=0.5, linestyles='-', zorder=2)
            ax.clabel(cs, fontsize=7, fmt='%g')
        # land patch + coastline from the ETOPO1 land mask (m_gshhs equivalent)
        land = np.ma.masked_less(elev, 0)
        ax.contourf(elon, elat, land, levels=[0, elev.max() + 1],
                    colors=[(0.8, 0.8, 0.8)], zorder=4)
        ax.contour(elon, elat, elev, levels=[0], colors='k',
                   linewidths=0.7, zorder=5)
    except FileNotFoundError as e:
        print('warning: %s (skipping bathymetry/coastline)' % e)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect(1 / np.cos(np.deg2rad((lat_min + lat_max) / 2)))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')


def eddy_subroutine(params, u, v, time, lons, lats, lon=None, lat=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    save_data_name = params['save_data_name']
    save_fig_name = params['save_fig_name']
    save_streams = params['save_streams']
    dij = params['dij']

    if params['lonlatlimits']:
        lon_min = params['lon_min']; lon_max = params['lon_max']
        lat_min = params['lat_min']; lat_max = params['lat_max']
        jlon = (lon >= lon_min) & (lon <= lon_max)
        jlat = (lat >= lat_min) & (lat <= lat_max)
        lon1 = lon[jlon]; lat1 = lat[jlat]
        # arrays are (nlon, nlat), matching MATLAB [lats,lons]=meshgrid(lat,lon)
        lons = lons[np.ix_(jlon, jlat)]
        lats = lats[np.ix_(jlon, jlat)]
        u = u[np.ix_(jlon, jlat)]
        v = v[np.ix_(jlon, jlat)]
    else:
        lon_min = float(np.floor(10 * np.nanmin(lons)) / 10)
        lon_max = float(np.ceil(10 * np.nanmax(lons)) / 10)
        lat_min = float(np.floor(10 * np.nanmin(lats)) / 10)
        lat_max = float(np.ceil(10 * np.nanmax(lats)) / 10)

    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    winding_thres = params['winding_thres']
    d_thres = params['d_thres']
    baddir_thres = params['baddir_thres']
    param_center_dist = params['param_center_dist']
    new_dist_thres = params['new_dist_thres']
    dlon = params['radar_km_resolution'] / 6
    dlat = dlon

    # ---------------- identification loop ----------------
    i_streams = 0
    streamsx = []; streamsy = []
    streamsx_center = []; streamsy_center = []
    streamslat = []; streamslon = []
    streamslat_center = []; streamslon_center = []
    winding_dir = []; angular_vel = []

    # whole-grid UTM conversion, reused for every window whose auto-zone matches
    # (geog2utm output depends only on the zone, so slicing is exact)
    gzone = auto_zone(lons)
    gxf, gyf, _ = geog2utm_nodisp(lons.ravel(), lats.ravel(), gzone)
    GX = gxf.reshape(lons.shape)
    GY = gyf.reshape(lons.shape)

    # fast_omega=1 (default): angular velocity interpolated on ONE whole-grid
    # triangulation per timestep instead of a fresh one per eddy streamline
    # (~200 ms each). Affects only the omega diagnostic, at the level of the
    # existing scatteredInterpolant-vs-griddata difference; set 0 for the
    # per-window triangulation that matches the pure-Python baseline bit-exactly.
    fast_omega = params.get('fast_omega', 1)
    g_interp = None

    # column-major order, matching MATLAB find(~isnan(u)) so eddy numbering agrees
    cols_t, rows_t = np.where(~np.isnan(u.T))
    rows, cols = rows_t, cols_t
    for k in range(len(rows)):
        ik0 = rows[k]; jk0 = cols[k]
        start_lon = lons[ik0, jk0]
        start_lat = lats[ik0, jk0]

        ik = slice(max(0, ik0 - dij), min(lons.shape[0], ik0 + dij + 1))
        jk = slice(max(0, jk0 - dij), min(lons.shape[1], jk0 + dij + 1))
        lons1 = lons[ik, jk]; lats1 = lats[ik, jk]
        u1 = u[ik, jk]; v1 = v[ik, jk]

        UTMzone = auto_zone(lons1)
        if UTMzone == gzone:
            x1 = GX[ik, jk]
            y1 = GY[ik, jk]
        else:  # window straddles a zone boundary differently: compute exactly
            x1f, y1f, UTMzone = geog2utm_nodisp(lons1.ravel(), lats1.ravel(), UTMzone)
            x1 = x1f.reshape(lons1.shape)
            y1 = y1f.reshape(lons1.shape)
        x1f = x1.ravel(); y1f = y1.ravel()
        startx, starty, UTMzone = geog2utm_nodisp(start_lon, start_lat, UTMzone)
        startx = float(startx[0]); starty = float(starty[0])

        try:
            X, Y = stream2_dc(x1, y1, u1, v1, startx, starty, 0.2, 2000)
        except IndexError:
            continue

        if len(X) < 20 or np.isnan(X[9]):
            continue

        nn = np.where(np.isnan(X))[0]
        if nn.size:
            stline_x = X[:nn[0]]
            stline_y = Y[:nn[0]]
        else:
            stline_x = X
            stline_y = Y
        if len(stline_x) < 3:
            continue

        if _winding_kernel is not None:
            eddy, ii, winding, dir1 = _winding_kernel(
                np.ascontiguousarray(stline_x), np.ascontiguousarray(stline_y),
                float(winding_thres), float(d_thres), float(baddir_thres))
        else:
            eddy, ii, winding, dir1 = _winding_py(
                stline_x, stline_y, winding_thres, d_thres, baddir_thres)

        if eddy:
            i_streams += 1
            sx = stline_x[:ii + 1]
            sy = stline_y[:ii + 1]
            streamsx.append(sx)
            streamsy.append(sy)
            streamsx_center.append(np.mean(sx))
            streamsy_center.append(np.mean(sy))
            slat, slon = utm2ll(sx * 1000, sy * 1000, UTMzone)
            streamslat.append(slat)
            streamslon.append(slon)
            clat, clon = utm2ll(streamsx_center[-1] * 1000,
                                streamsy_center[-1] * 1000, UTMzone)
            streamslat_center.append(float(clat[0]))
            streamslon_center.append(float(clon[0]))
            winding_dir.append(dir1)

            # angular velocity in deg/sec (griddata == LinearNDInterpolator on
            # a Delaunay triangulation)
            if fast_omega:
                if g_interp is None:
                    gtri = Delaunay(np.column_stack([gxf, gyf]))
                    g_interp = (LinearNDInterpolator(gtri, u.ravel()),
                                LinearNDInterpolator(gtri, v.ravel()))
                u2 = g_interp[0](sx, sy)
                v2 = g_interp[1](sx, sy)
            else:
                tri = Delaunay(np.column_stack([x1f, y1f]))
                u2 = LinearNDInterpolator(tri, u1.ravel())(sx, sy)
                v2 = LinearNDInterpolator(tri, v1.ravel())(sx, sy)
            vv = np.sqrt(u2 ** 2 + v2 ** 2)[:-1]
            dr = np.hypot(np.diff(sx), np.diff(sy))
            with np.errstate(divide='ignore', invalid='ignore'):
                dt = np.sum(dr / vv * 1000)
            angular_vel.append(winding / dt)

    # ---------------- streamline clustering ----------------
    eddy_nums = np.arange(i_streams)
    for i1 in range(i_streams - 1):
        if eddy_nums[i1] != i1:
            continue
        for j1 in range(i1 + 1, i_streams):
            if new_dist_thres:
                # Cahl et al., 2023: same eddy when one center is inside the
                # other's streamline polygon, or the polygons overlap
                in1 = _inpolygon(streamslon_center[j1], streamslat_center[j1],
                                 streamslon[i1], streamslat[i1])
                in2 = _inpolygon(streamslon_center[i1], streamslat_center[i1],
                                 streamslon[j1], streamslat[j1])
                in11 = _inpolygon(streamslon[j1], streamslat[j1],
                                  streamslon[i1], streamslat[i1])
                in22 = _inpolygon(streamslon[i1], streamslat[i1],
                                  streamslon[j1], streamslat[j1])
                if in1 or in2 or np.sum(in11) > 1 or np.sum(in22) > 1:
                    eddy_nums[j1] = i1
                else:
                    dxa, dya, _ = geog2utm_nodisp(
                        streamslon_center[i1], streamslat_center[i1],
                        streamslon_center[j1], streamslat_center[j1])
                    if np.hypot(dxa[0], dya[0]) < param_center_dist:
                        eddy_nums[j1] = i1
            else:  # Sadarjoen and Post 2000
                dxa, dya, _ = geog2utm_nodisp(
                    streamslon_center[i1], streamslat_center[i1],
                    streamslon_center[j1], streamslat_center[j1])
                d = np.hypot(dxa[0], dya[0])
                if d < param_center_dist and winding_dir[i1] == winding_dir[j1]:
                    eddy_nums[j1] = i1

    # ---------------- figure ----------------
    h = None
    ax = None
    if params['plt_debug'] and not params['plt_each_eddy']:
        h = plt.figure(figsize=(13, 11))
        ax = h.add_subplot(111)
        _base_map(ax, plt, params, lons, lats, u, v,
                  lon_min, lon_max, lat_min, lat_max)

    # ---------------- eddy parameters ----------------
    uniq = np.unique(eddy_nums) if i_streams else np.array([], dtype=int)
    eddies = 0
    eddy_center_lon = []; eddy_center_lat = []
    eddy_dir = []; eddy_angular_vel = []; eddy_streamlines = []
    eddy_length_x = []; eddy_length_y = []; eddy_ellipse_theta = []
    eddy_streamlines_lat = []; eddy_streamlines_lon = []

    for j1 in range(len(uniq)):
        k = np.where(eddy_nums == uniq[j1])[0]
        eddies += 1
        ptsx = np.concatenate([streamslon[kk] for kk in k])
        ptsy = np.concatenate([streamslat[kk] for kk in k])
        if save_streams:
            eddy_streamlines_lon.append([streamslon[kk] for kk in k])
            eddy_streamlines_lat.append([streamslat[kk] for kk in k])

        px, py, UTMzone = geog2utm_nodisp(ptsx, ptsy)
        cx = np.mean(px); cy = np.mean(py)
        clat, clon = utm2ll(cx * 1000, cy * 1000, UTMzone)
        eddy_center_lon.append(float(clon[0]))
        eddy_center_lat.append(float(clat[0]))
        eddy_dir.append(winding_dir[k[0]])
        eddy_angular_vel.append(np.mean([angular_vel[kk] for kk in k]))
        eddy_streamlines.append(len(k))

        M = np.cov(px, py)
        l, dvec = np.linalg.eigh(M)
        dx1 = dvec[0, 0]; dy1 = dvec[1, 0]
        l1 = np.sqrt(2 * abs(l[0]))
        l2 = np.sqrt(2 * abs(l[1]))
        ellipse_theta = np.degrees(np.arctan2(dy1, dx1))
        eddy_ellipse_theta.append(ellipse_theta)
        eddy_length_x.append(l1)
        eddy_length_y.append(l2)

        # streamline and ellipse plots
        if params['plt_debug'] and len(k) > params['plt_min_streams']:
            if params['plt_each_eddy']:
                h = plt.figure(figsize=(13, 11))
                ax = h.add_subplot(111)
                _base_map(ax, plt, params, lons, lats, u, v,
                          eddy_center_lon[-1] - dlon, eddy_center_lon[-1] + dlon,
                          eddy_center_lat[-1] - dlat, eddy_center_lat[-1] + dlat)

            for kk in k:
                ax.plot(streamslon[kk], streamslat[kk], 'm', linewidth=0.4, zorder=6)
                ax.plot(streamslon_center[kk], streamslat_center[kk], 'rx',
                        markersize=4, zorder=7)

            it = np.arange(1, 361)
            ex = l1 * np.cos(np.deg2rad(it))
            ey = l2 * np.sin(np.deg2rad(it))
            th = np.deg2rad(ellipse_theta)
            exr = ex * np.cos(th) - ey * np.sin(th)
            eyr = ex * np.sin(th) + ey * np.cos(th)
            ax.plot(eddy_center_lon[-1], eddy_center_lat[-1], 'kx',
                    markersize=12, markeredgewidth=2, zorder=8)
            ccx, ccy, UTMzone = geog2utm_nodisp(eddy_center_lon[-1], eddy_center_lat[-1])
            elat, elon2 = utm2ll((ccx[0] + exr) * 1000, (ccy[0] + eyr) * 1000, UTMzone)
            ax.plot(np.append(elon2, elon2[0]), np.append(elat, elat[0]),
                    color=(0.85, 0.85, 0.85), linewidth=1.5, zorder=8)

            if params['plt_debug_sv'] and params['plt_each_eddy']:
                h.savefig('%s_eddy%d.png' % (save_fig_name, j1 + 1),
                          dpi=params['print_res'], bbox_inches='tight')
                if params['plt_close']:
                    plt.close(h)

    # ---------------- save data (same variable names as MATLAB) ----------------
    n_ed = eddies
    maxk = max((len(s) for s in eddy_streamlines_lon), default=0)
    cell_lon = np.empty((n_ed, maxk), dtype=object)
    cell_lat = np.empty((n_ed, maxk), dtype=object)
    for a in range(n_ed):
        for b in range(maxk):
            cell_lon[a, b] = np.zeros((0, 0))
            cell_lat[a, b] = np.zeros((0, 0))
    for a in range(n_ed):
        if save_streams and a < len(eddy_streamlines_lon):
            for b, s in enumerate(eddy_streamlines_lon[a]):
                cell_lon[a, b] = np.asarray(s).reshape(-1, 1)
            for b, s in enumerate(eddy_streamlines_lat[a]):
                cell_lat[a, b] = np.asarray(s).reshape(-1, 1)

    savemat(save_data_name + '.mat', {
        'eddy_center_lat': np.asarray(eddy_center_lat),
        'eddy_center_lon': np.asarray(eddy_center_lon),
        'eddy_dir': np.asarray(eddy_dir, dtype=float),
        'eddy_angular_vel': np.asarray(eddy_angular_vel),
        'eddy_streamlines': np.asarray(eddy_streamlines, dtype=float),
        'eddy_length_x': np.asarray(eddy_length_x),
        'eddy_length_y': np.asarray(eddy_length_y),
        'eddy_ellipse_theta': np.asarray(eddy_ellipse_theta),
        'eddy_streamlines_lat': cell_lat,
        'eddy_streamlines_lon': cell_lon,
        'UTMzone': UTMzone if i_streams else 0,
        'time': np.asarray(time, dtype=float),
    })

    # ---------------- save plot ----------------
    if params['plt_debug_sv'] and params['plt_debug'] and not params['plt_each_eddy'] and h is not None:
        h.savefig(save_fig_name + '.png', dpi=params['print_res'], bbox_inches='tight')
        if params['plt_close']:
            plt.close(h)

    print('%d eddies found' % eddies)
    for i_eddy in range(eddies):
        print('eddy #%d has %d streamlines' % (i_eddy + 1, eddy_streamlines[i_eddy]))
    print(' ')
    return eddies, eddy_streamlines
