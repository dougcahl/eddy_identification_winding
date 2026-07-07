# eddy_identification_winding — Python v2

Python port of the current MATLAB repo
([dougcahl/eddy_identification_winding](https://github.com/dougcahl/eddy_identification_winding)),
bringing the MATLAB upgrades to the v1.0 Python code.

## What's new vs the v1.0 Python code

| Feature | v1.0 Python | this v2 port |
|---|---|---|
| Streamline integrator | modified matplotlib `streamplot` (terminates when a trajectory re-enters a visited cell — circling streamlines get cut short) | `stream2_dc.py`, port of `stream2_dc.m`: inverse-distance-weighted 4-node velocity, Euler steps, no cell mask |
| Thresholds | needs lenient 270° winding / 25 km closure to compensate | works with the repo defaults: 300° winding / 10 km closure |
| Clustering | Sadarjoen & Post 2000 (center distance + same rotation) | **Cahl et al. 2023**: `inpolygon` tests — same eddy when one streamline's center lies inside the other's polygon or the polygons overlap; center-distance fallback (`new_dist_thres = 1`); set `0` for the old method |
| Ellipse fit | `sqrt(eig)` axes | `sqrt(2*eig)` axes (as in the MATLAB repo) |
| Map figure | quiver only | speed field, unit-vector quiver, **ETOPO1 bathymetry contours** with labels, **land patch + coastline** |
| Streamline output | not saved | per-streamline lat/lon saved in the `.mat` (`save_streams`) |
| Results format | NetCDF | `.mat` with the same variable names as the MATLAB repo (`eddy_center_lat`, `eddy_dir`, `eddy_streamlines_lat`, ...) — directly comparable/loadable in MATLAB |
| UTM | `utm` pip package (WGS84 auto-zones) | ports of `geog2utm_nodisp.m` / `utm2ll.m` with the exact MATLAB constants |

## Performance

The hot loops (`stream2_dc` integration and the winding-angle scan) have
numba-compiled kernels in `kernels_numba.py` (LLVM, strict IEEE — results are
bit-identical to the pure-Python fallbacks, which are used automatically when
numba is not installed). Angular velocity is interpolated on one whole-grid
Delaunay triangulation per timestep (`fast_omega=1`, default); set
`params['fast_omega'] = 0` for per-window triangulation (matches the
pure-Python baseline bit-exactly, ~50 s/timestep instead of ~1.5 s).
A/B tested on 6 timesteps (`test_fast_omega.py`): detection geometry —
eddy counts, centers, directions, streamline counts and coordinates, ellipse
fits — is **bit-identical** either way; only the omega diagnostic shifts
(relative differences up to ~10% on individual eddies, absolute ≲ 1e-4 deg/s,
occasionally ±1 NaN), smaller than the existing MATLAB↔Python interpolant
difference in absolute terms.

data2 example (90×95 grid, ~2,380 valid points, one timestep, single core):

| Implementation | per timestep |
|---|---|
| MATLAB repo (R2024b JIT) | ~10–27 s |
| Python v2, pure Python | ~395 s |
| Python v2 + numba kernels | ~52 s |
| Python v2 + numba + fast_omega | **~1.5 s** |

The full 101-timestep pipeline runs in ~32 s on 14 cores (identification
27.9 s + tracking 2.4 s + analysis 1.3 s; was ~24 min pure Python, MATLAB
~36 min single-core). Verified against MATLAB with `verify_vs_matlab.py`:
all 191 eddies match with centers within 10⁻¹³ degrees and identical tracks.

## Requirements

Python 3.10+ with `numpy`, `scipy`, `matplotlib`, `netCDF4`; `numba` optional
but strongly recommended (pure-Python fallback is ~40x slower per point).

Bathymetry/coastline need the ETOPO1 binary `etopo1_ice_g_i2.bin`
(same file m_map uses). Put it in `./etopo1/` or edit `SEARCH_PATHS` in
`etopo1_reader.py`. Without it the code still runs and simply skips
bathymetry and coastline. The coastline is drawn from the ETOPO1 land mask
(1 arc-minute), standing in for the GSHHS coastline that m_map uses.

## Run the example

```
python eddy_uvdata.py
```

Runs the bundled Delaware Bay sample (`data/202004030900_hfr_usegc_2km_rtv_uwls_NDBC.nc`,
US East Coast 2 km HF radar, 2020-04-03 09:00 UTC) with the same subregion and
parameters as `eddy_uvdata.m`. Outputs:

- `data/results/<fn>.mat` — eddy centers, rotation, angular velocity, ellipse
  fit, per-streamline lat/lon
- `data/results_figs/<fn>.png` — map figure

## Tracking pipeline (multi-timestep example)

```
python eddy_uvdata_loop.py       # identification for all 101 timesteps (parallel)
python eddy_tracking.py          # link identifications into tracks
python analyze_eddy_tracks.py    # track map + statistics figures
```

## Data formats — NetCDF only

The pipeline is fully NetCDF-native; no `.mat` files are read or written.

- **Input**: `data/data2.nc` (or any lon/lat/u/v/time source). Convert a
  MATLAB `data2.mat` with `python convert_data.py data`.
- **Results**: per-timestep `data2_<i>.nc` and `data2_tracks.nc`, written by
  `results_io.py` — per-eddy variables plus ragged streamline/track series as
  concatenated arrays with `*_start`/`*_len` index vectors.
- Time is CF-style `days since 1970-01-01`, with the exact original MATLAB
  datenum kept in `matlab_datenum`.
- `convert_data.py results` imports MATLAB-produced result `.mat` files to
  the same `.nc` schema, and `verify_vs_matlab.py` compares a Python run
  against a MATLAB run file-by-file.

## Files

- `eddy_uvdata.py` — single-file example driver (port of `eddy_uvdata.m`)
- `eddy_uvdata_loop.py` — multi-timestep identification, parallel over cores (port of `eddy_uvdata_loop.m`)
- `eddy_tracking.py` — eddy tracking with the Cahl 2023 ellipse/inpolygon continuation test (port of `eddy_tracking.m`)
- `analyze_eddy_tracks.py` — track map + histogram statistics (port of `analyze_eddy_tracks.m`)
- `eddy_subroutine.py` — identification + clustering + plotting + save (port of `eddy_subroutine.m`)
- `stream2_dc.py` — streamline integrator (port of `stream2_dc.m`)
- `utm_dc.py` — `geog2utm_nodisp` / `utm2ll` ports
- `etopo1_reader.py` — ETOPO1 subregion reader (replaces `m_etopo2`)
- `results_io.py` — native NetCDF read/write for results and tracks
- `convert_data.py` — import MATLAB `.mat` data/results to the `.nc` schema
- `verify_vs_matlab.py` — file-by-file comparison of a Python run vs a MATLAB run
- `kernels_numba.py` — numba (LLVM) kernels for the integrator and winding scan
