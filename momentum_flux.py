#!/usr/bin/env python3
"""Optional analysis module: eddy momentum flux and EKE from HF radar
surface currents, following Cahl (2023), "HF Radar: Shining a Light on Ocean
Currents", Chapter 4 (University of South Carolina).

Method (thesis section 4.3.2):
  1. Rotate (u, v) into a shelf-break reference frame set by the local
     isobath bearing B (compass degrees): the along-shelf component points
     downstream along the isobath and the cross-shelf component points
     offshore (positive flux = momentum leaving the shelf, e.g. toward the
     Gulf Stream in Long Bay, SC; Lee et al. 1991 convention).
  2. Reynolds-decompose against a time mean: u = ubar + u', v = vbar + v'.
  3. Momentum flux = < u'_cross * v'_along > (reported in cm^2 s^-2) and
     EKE = sqrt(< u'^2 + v'^2 >) (cm/s), per grid point.
  4. Sample the flux along a shelf-break segment for an along-shelf profile.

VALIDITY — read before using (all discussed in the thesis):
  * These are second-order statistics. The thesis used a FULL YEAR of data
    with yearly/monthly means. Records shorter than ~1 month give unstable,
    eddy-dominated "means"; this module warns loudly below 30 days.
  * Currents must be SUBTIDAL (detided / low-pass filtered): otherwise
    <u'v'> reflects tidal ellipse orientation, not eddy momentum transport.
    A spectral check warns when semidiurnal energy looks significant; an
    optional boxcar low-pass (lowpass_hours) is provided.
  * The rotation bearing is only meaningful along a roughly straight isobath
    segment (thesis: 56 deg N at the 100 m isobath in Long Bay, valid over
    ~40 km). Use isobath_bearing() to estimate it from ETOPO1 elsewhere.
  * HF radar u,v errors are GDOP-correlated, which biases the covariance;
    the thesis limited the calculation to where both radars' look directions
    are within 30 deg of boresight. Supply a quality mask if you have one.
  * This is a SURFACE momentum flux signature, not a depth-integrated
    transport, and may include wind-driven contributions.
"""

import os
import warnings

import numpy as np
from netCDF4 import Dataset

from etopo1_reader import read_etopo1

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root


# ---------------------------------------------------------------- helpers
def isobath_bearing(depth, near_lonlat, lon_range, lat_range, fit_km=40):
    """Estimate the local isobath bearing (compass degrees of the along-
    isobath direction) from ETOPO1 near a point, by PCA of the depth-contour
    vertices within fit_km of near_lonlat. Returns bearing in [0, 180)."""
    import matplotlib.pyplot as plt
    elev, elon, elat = read_etopo1(lon_range[0], lon_range[1],
                                   lat_range[0], lat_range[1])
    fig = plt.figure()
    cs = plt.contour(elon, elat, elev, levels=[-abs(depth)])
    segs = [p.vertices for p in cs.get_paths()]
    plt.close(fig)
    if not segs:
        raise ValueError('no %g m isobath found in the given range' % depth)
    pts = np.vstack(segs)
    km_per_deg = 111.32
    dx = (pts[:, 0] - near_lonlat[0]) * km_per_deg * np.cos(np.deg2rad(near_lonlat[1]))
    dy = (pts[:, 1] - near_lonlat[1]) * km_per_deg
    keep = np.hypot(dx, dy) < fit_km
    if keep.sum() < 10:
        raise ValueError('isobath does not pass within %g km of the point' % fit_km)
    x = dx[keep] - dx[keep].mean()
    y = dy[keep] - dy[keep].mean()
    _, vecs = np.linalg.eigh(np.cov(x, y))
    ex, ey = vecs[:, -1]  # principal direction (largest eigenvalue)
    bearing = np.degrees(np.arctan2(ex, ey)) % 180
    return bearing


def _tidal_energy_fraction(u, v, dt_hours):
    """Fraction of (spatial-mean speed) variance in the 10-14 h band."""
    spd = np.sqrt(np.nanmean(u, axis=(0, 1)) ** 2
                  + np.nanmean(v, axis=(0, 1)) ** 2)
    spd = spd - np.nanmean(spd)
    spd = np.nan_to_num(spd)
    F = np.fft.rfft(spd * np.hanning(len(spd)))
    freq = np.fft.rfftfreq(len(spd), d=dt_hours)  # cycles per hour
    P = np.abs(F) ** 2
    band = (freq > 1 / 14.0) & (freq < 1 / 10.0)
    return float(P[band].sum() / P[1:].sum()) if P[1:].sum() > 0 else 0.0


def momentum_flux(lons, lats, u, v, time_datenum, bearing_deg,
                  lowpass_hours=None, min_valid_frac=0.5, quality_mask=None):
    """Compute momentum flux and EKE fields.

    lons, lats : 2-D grids;  u, v : (ni, nj, nt) in m/s;
    time_datenum : (nt,) MATLAB datenum;  bearing_deg : compass bearing of
    the along-isobath (downstream) direction.

    Returns dict with ubar, vbar (m/s), flux (cm^2/s^2), eke (cm/s),
    record_days, tidal_fraction, warnings (list of strings).
    """
    u = np.asarray(u, dtype=float).copy()
    v = np.asarray(v, dtype=float).copy()
    t = np.asarray(time_datenum, dtype=float)
    notes = []

    record_days = float(t[-1] - t[0])
    dt_hours = float(np.median(np.diff(t)) * 24)
    if record_days < 30:
        msg = ('SHORT RECORD: %.1f days. Momentum flux/EKE are second-order '
               'statistics; the thesis used a full year. Treat these results '
               'as a demonstration, not stable statistics.' % record_days)
        warnings.warn(msg)
        notes.append(msg)

    tf = _tidal_energy_fraction(u, v, dt_hours)
    if tf > 0.2 and lowpass_hours is None:
        msg = ('POSSIBLE TIDAL ENERGY: %.0f%% of the domain-mean speed '
               'variance sits in the 10-14 h band. Detide or pass '
               'lowpass_hours (e.g. 33) or <u\'v\'> will reflect tidal '
               'ellipses, not eddies.' % (100 * tf))
        warnings.warn(msg)
        notes.append(msg)

    if lowpass_hours:
        w = max(1, int(round(lowpass_hours / dt_hours)))
        kernel = np.ones(w) / w
        shape = u.shape
        u2 = u.reshape(-1, shape[2])
        v2 = v.reshape(-1, shape[2])
        for arr in (u2, v2):
            for r in range(arr.shape[0]):
                row = arr[r]
                if np.all(np.isnan(row)):
                    continue
                filled = np.where(np.isnan(row), np.nanmean(row), row)
                sm = np.convolve(filled, kernel, mode='same')
                sm[np.isnan(row)] = np.nan
                arr[r] = sm
        u = u2.reshape(shape)
        v = v2.reshape(shape)
        notes.append('boxcar low-pass applied: %g h' % lowpass_hours)

    if quality_mask is not None:
        u[~quality_mask] = np.nan
        v[~quality_mask] = np.nan

    # rotation: e_along points downstream along the isobath (compass bearing),
    # e_cross points offshore (90 deg clockwise from e_along)
    Br = np.deg2rad(bearing_deg)
    e_along = (np.sin(Br), np.cos(Br))
    e_cross = (e_along[1], -e_along[0])
    uc = u * e_cross[0] + v * e_cross[1]
    va = u * e_along[0] + v * e_along[1]

    valid = np.isfinite(u).sum(axis=2)
    enough = valid >= min_valid_frac * u.shape[2]

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        ubar = np.nanmean(u, axis=2)
        vbar = np.nanmean(v, axis=2)
        ucp = uc - np.nanmean(uc, axis=2, keepdims=True)
        vap = va - np.nanmean(va, axis=2, keepdims=True)
        flux = np.nanmean(ucp * vap, axis=2) * 1e4          # cm^2 s^-2
        eke = np.sqrt(np.nanmean(u ** 2 + v ** 2, axis=2)
                      - np.nanmean(u, axis=2) ** 2
                      - np.nanmean(v, axis=2) ** 2) * 100   # cm/s

    for arr in (ubar, vbar, flux, eke):
        arr[~enough] = np.nan

    return {'ubar': ubar, 'vbar': vbar, 'flux': flux, 'eke': eke,
            'record_days': record_days, 'tidal_fraction': tf,
            'bearing_deg': float(bearing_deg), 'warnings': notes}


def segment_profile(lons, lats, field, p1, p2, n=40):
    """Sample a field along the segment p1->p2 (lon, lat) at n points by
    nearest grid node. Returns (seg_lon, seg_lat, dist_km, values)."""
    seg_lon = np.linspace(p1[0], p2[0], n)
    seg_lat = np.linspace(p1[1], p2[1], n)
    km_per_deg = 111.32
    dist = np.hypot((seg_lon - seg_lon[0]) * km_per_deg
                    * np.cos(np.deg2rad(seg_lat.mean())),
                    (seg_lat - seg_lat[0]) * km_per_deg)
    vals = np.full(n, np.nan)
    for k in range(n):
        d2 = (lons - seg_lon[k]) ** 2 + (lats - seg_lat[k]) ** 2
        i, j = np.unravel_index(np.argmin(d2), d2.shape)
        vals[k] = field[i, j]
    return seg_lon, seg_lat, dist, vals


def save_nc(path, lons, lats, res, seg=None):
    with Dataset(path, 'w', format='NETCDF4') as nc:
        nc.title = 'eddy momentum flux and EKE (Cahl 2023, thesis ch. 4 method)'
        nc.bearing_deg = res['bearing_deg']
        nc.record_days = res['record_days']
        nc.tidal_fraction = res['tidal_fraction']
        nc.warnings = ' | '.join(res['warnings']) if res['warnings'] else 'none'
        nc.createDimension('i', lons.shape[0])
        nc.createDimension('j', lons.shape[1])
        for name, units, data in [
                ('lon', 'degrees_east', lons), ('lat', 'degrees_north', lats),
                ('ubar', 'm s-1', res['ubar']), ('vbar', 'm s-1', res['vbar']),
                ('momentum_flux', 'cm2 s-2', res['flux']),
                ('eke', 'cm s-1', res['eke'])]:
            vv = nc.createVariable(name, 'f8', ('i', 'j'), zlib=True,
                                   fill_value=np.nan)
            vv.units = units
            vv[:] = data
        if seg is not None:
            nc.createDimension('seg_pt', len(seg[0]))
            for name, units, data in [
                    ('seg_lon', 'degrees_east', seg[0]),
                    ('seg_lat', 'degrees_north', seg[1]),
                    ('seg_dist', 'km', seg[2]),
                    ('seg_flux', 'cm2 s-2', seg[3])]:
                vv = nc.createVariable(name, 'f8', ('seg_pt',))
                vv.units = units
                vv[:] = data
    print('saved %s' % path)


# ---------------------------------------------------------------- demo
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from eddy_uvdata_loop import load_data
    from results_io import load_tracks_nc

    lons, lats, u1, v1, time1 = load_data()

    # thesis values for Long Bay, SC: 100 m isobath bearing 56 deg N and the
    # ~40 km shelf-break segment where both radars look near boresight
    P1 = (-78.28, 32.78)
    P2 = (-77.90, 32.99)
    THETA = 56.0
    try:
        est = isobath_bearing(100, ((P1[0] + P2[0]) / 2, (P1[1] + P2[1]) / 2),
                              (lons.min(), lons.max()), (lats.min(), lats.max()))
        print('isobath_bearing(100 m) estimate: %.1f degN (thesis uses %.0f)'
              % (est, THETA))
    except Exception as e:
        print('isobath bearing estimate failed: %s' % e)

    res = momentum_flux(lons, lats, u1, v1, time1, THETA)
    for w in res['warnings']:
        print('WARNING:', w)
    print('record: %.1f days | tidal fraction: %.0f%%'
          % (res['record_days'], 100 * res['tidal_fraction']))

    seg = segment_profile(lons, lats, res['flux'], P1, P2)
    print('segment mean flux: %.1f cm2/s2 (thesis 2013 yearly mean ~100)'
          % np.nanmean(seg[3]))
    os.makedirs('data/results', exist_ok=True)
    save_nc('data/results/momentum_flux.nc', lons, lats, res, seg)

    # ------------------------------------------------ figures
    lon_min = float(np.floor(10 * np.nanmin(lons)) / 10)
    lon_max = float(np.ceil(10 * np.nanmax(lons)) / 10)
    lat_min = float(np.floor(10 * np.nanmin(lats)) / 10)
    lat_max = float(np.ceil(10 * np.nanmax(lats)) / 10)
    try:
        elev, elon, elat = read_etopo1(lon_min - .05, lon_max + .05,
                                       lat_min - .05, lat_max + .05)
    except FileNotFoundError:
        elev = None

    def base(ax):
        if elev is not None:
            ax.contour(elon, elat, elev,
                       levels=[-1000, -500, -250, -100, -50, -25, -10],
                       colors=[(.55, .55, .55)], linewidths=0.5, zorder=2)
            land = np.ma.masked_less(elev, 0)
            ax.contourf(elon, elat, land, levels=[0, elev.max() + 1],
                        colors=[(.8, .8, .8)], zorder=3)
            ax.contour(elon, elat, elev, levels=[0], colors='k',
                       linewidths=0.7, zorder=4)
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_aspect(1 / np.cos(np.deg2rad((lat_min + lat_max) / 2)))

    short = res['record_days'] < 30
    tag = ('DEMONSTRATION — %.1f-day record (thesis used a full year)'
           % res['record_days']) if short else ''

    fig, axs = plt.subplots(2, 2, figsize=(17, 14))
    ax = axs[0, 0]
    spd = np.hypot(res['ubar'], res['vbar'])
    pc = ax.pcolormesh(lons, lats, spd, shading='auto', zorder=1)
    with np.errstate(invalid='ignore'):
        ax.quiver(lons[::3, ::3], lats[::3, ::3],
                  (res['ubar'] / spd)[::3, ::3], (res['vbar'] / spd)[::3, ::3],
                  color='k', scale=40, width=0.0018, zorder=2)
    plt.colorbar(pc, ax=ax, fraction=0.045).set_label('Mean speed (m/s)')
    base(ax)
    ax.set_title('(a) Mean surface currents')

    ax = axs[0, 1]
    vmax = np.nanpercentile(np.abs(res['flux']), 98)
    pc = ax.pcolormesh(lons, lats, res['flux'], cmap='RdBu_r',
                       vmin=-vmax, vmax=vmax, shading='auto', zorder=1)
    plt.colorbar(pc, ax=ax, fraction=0.045).set_label(
        "Momentum flux <u'v'> (cm$^2$ s$^{-2}$)")
    ax.plot([P1[0], P2[0]], [P1[1], P2[1]], 'k-', lw=2.5, zorder=5)
    base(ax)
    ax.set_title("(b) Momentum flux, rotated %g$^\\circ$N (positive = offshore)"
                 % res['bearing_deg'])

    ax = axs[1, 0]
    pc = ax.pcolormesh(lons, lats, res['eke'], shading='auto', zorder=1)
    plt.colorbar(pc, ax=ax, fraction=0.045).set_label('EKE (cm/s)')
    base(ax)
    ax.set_title("(c) Eddy kinetic energy $\\sqrt{<u'^2+v'^2>}$")

    ax = axs[1, 1]
    ax.plot(seg[2], seg[3], 'b.-')
    ax.axhline(0, color='k', lw=0.8)
    ax.axhline(np.nanmean(seg[3]), color='k', ls='--', lw=1,
               label='segment mean %.0f cm$^2$s$^{-2}$' % np.nanmean(seg[3]))
    ax.set_xlabel('Distance along shelf-break segment (km, SW to NE)')
    ax.set_ylabel("Momentum flux (cm$^2$ s$^{-2}$)")
    ax.grid(alpha=0.4)
    ax.legend()
    ax.set_title('(d) Flux along the 100 m isobath segment')

    fig.suptitle('Momentum flux and EKE — Long Bay, SC (data2)\n%s' % tag,
                 fontsize=13, color='crimson' if short else 'black')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig('data/results_figs/momentum_flux_maps.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # flux + eddy tracks overlay (thesis section 4.4.4)
    fig, ax = plt.subplots(figsize=(11, 9))
    pc = ax.pcolormesh(lons, lats, res['flux'], cmap='RdBu_r',
                       vmin=-vmax, vmax=vmax, shading='auto', zorder=1)
    plt.colorbar(pc, ax=ax, fraction=0.045).set_label(
        "Momentum flux <u'v'> (cm$^2$ s$^{-2}$)")
    base(ax)
    ax.plot([P1[0], P2[0]], [P1[1], P2[1]], 'k-', lw=2.5, zorder=5)
    try:
        tr = load_tracks_nc('data/results/data2_tracks.nc')
        for latc, lonc, om in zip(tr['lat_center'], tr['lon_center'],
                                  tr['omega']):
            if len(latc) < 12:
                continue
            c = 'b' if np.nanmean(om) > 0 else 'r'
            ax.plot(lonc, latc, '-' + c, lw=1.6, zorder=6)
            ax.plot(lonc[0], latc[0], 'o', mfc='none', mec=(0, .8, 0),
                    mew=2, zorder=7)
            ax.plot(lonc[-1], latc[-1], 'kx', ms=9, mew=2, zorder=7)
    except Exception as e:
        print('no tracks overlaid: %s' % e)
    ax.set_title('Momentum flux with eddy tracks (blue CCW, red CW)\n%s' % tag,
                 color='crimson' if short else 'black')
    fig.savefig('data/results_figs/momentum_flux_eddies.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    print('figures saved to data/results_figs/')
    print('MOMENTUM_FLUX_DONE')
