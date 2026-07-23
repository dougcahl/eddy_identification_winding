#!/usr/bin/env python3
"""Generate a detection figure for every timestep of data2 and build a
browsable HTML gallery grouped by number of eddies.

Output: data/gallery/index.html (open in a browser) + data/gallery/t###.png
"""

import os
import glob
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta

os.chdir(os.path.dirname(os.path.abspath(__file__)))

OUT = 'data/gallery'


def run_one(i):
    from eddy_uvdata_loop import load_data, make_params
    from eddy_subroutine import eddy_subroutine
    lons, lats, u1, v1, time1 = load_data()
    params = make_params()
    params['plt_debug'] = 1
    params['plt_close'] = 1
    params['print_res'] = 120
    params['save_data_name'] = 'data/testing/gal_%d' % i
    params['save_fig_name'] = '%s/t%03d' % (OUT, i)
    eddy_subroutine(params, u1[:, :, i - 1], v1[:, :, i - 1], time1[i - 1],
                    lons, lats)
    return i


def build_index():
    from results_io import load_ident_nc
    entries = []
    for i in range(1, 102):
        r = load_ident_nc('data/results/data2_%d.nc' % i)
        streams = sorted((int(s) for s in r['eddy_streamlines']), reverse=True)
        dt = datetime(1970, 1, 1) + timedelta(days=r['time'] - 719529.0)
        entries.append((i, len(streams), streams, dt))

    groups = {}
    for e in entries:
        groups.setdefault(e[1], []).append(e)

    html = ['''<!doctype html><html><head><meta charset="utf-8">
<title>Long Bay eddy detections — all 101 timesteps</title>
<style>
body { font-family: "Segoe UI", system-ui, sans-serif; background: #fafbfc;
       color: #14232e; margin: 0; padding: 2rem 1.5rem; }
h1 { font-size: 1.5rem; }
h2 { font-size: 1.15rem; margin: 2.2rem 0 0.8rem; color: #1b6e8c;
     border-bottom: 1px solid #dce4ea; padding-bottom: 0.3rem; }
.nav { position: sticky; top: 0; background: #fafbfcee; padding: 0.6rem 0;
       z-index: 10; border-bottom: 1px solid #dce4ea; }
.nav a { margin-right: 1.2rem; color: #1b6e8c; font-weight: 600;
         text-decoration: none; }
.nav a:hover { text-decoration: underline; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(430px, 1fr));
        gap: 1rem; }
.card { background: #fff; border: 1px solid #dce4ea; border-radius: 6px;
        padding: 0.6rem; }
.card img { width: 100%; height: auto; display: block; border-radius: 3px; }
.cap { font-size: 0.82rem; color: #55687a; margin-top: 0.45rem; }
.cap b { color: #14232e; }
</style></head><body>
<h1>Long Bay eddy detections — all 101 timesteps</h1>
<p style="color:#55687a">data2 (3 km WERA HF radar, Jan 3&ndash;5 2013, half-hour steps),
Python winding-angle pipeline. Grouped by number of detected eddies; click any
map to open the full-size image.</p>
<div class="nav">''']
    counts = sorted(groups, reverse=True)
    for n in counts:
        html.append('<a href="#g%d">%d eddies (%d)</a>' % (n, n, len(groups[n])))
    html.append('</div>')

    for n in counts:
        html.append('<h2 id="g%d">%d eddies — %d timesteps</h2>' % (n, n, len(groups[n])))
        html.append('<div class="grid">')
        for i, _, streams, dt in groups[n]:
            fn = 't%03d.png' % i
            html.append(
                '<div class="card"><a href="%s" target="_blank">'
                '<img loading="lazy" src="%s" alt="timestep %d"></a>'
                '<div class="cap"><b>t%d</b> &middot; %s UTC &middot; '
                'streamlines: %s</div></div>'
                % (fn, fn, i, i, dt.strftime('%Y-%m-%d %H:%M'),
                   ', '.join(map(str, streams))))
        html.append('</div>')
    html.append('</body></html>')

    with open(os.path.join(OUT, 'index.html'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(html))
    print('gallery index written: %s/index.html' % OUT)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    os.makedirs('data/testing', exist_ok=True)
    workers = max(1, (os.cpu_count() or 4) - 2)
    print('rendering 101 figures on %d workers' % workers)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        done = 0
        for i in ex.map(run_one, range(1, 102)):
            done += 1
            if done % 10 == 0:
                print('%d/101 rendered' % done, flush=True)
    for f in glob.glob('data/testing/gal_*.nc'):
        os.remove(f)
    build_index()
    print('GALLERY_DONE')
