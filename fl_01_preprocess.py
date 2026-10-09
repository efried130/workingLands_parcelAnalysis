# ---
# jupytext:
#   text_representation:
#     extension: .py
#     format_name: percent
# ---

# %% [markdown]
# # `fl_01_preprocess.py` — build the statewide cache. Run once.
#
# The source downloads are national (CSB) or statewide (NHD, NWI, HUC12). This
# script cuts each one down to the state and stores the result as one layer in
# **`clipped_layers/<st>_state_layers.gpkg`** (`ky_` for Kentucky, `tn_` for
# Tennessee - the state is `STATE_FIPS` in §0 of the config).
#
# Once that file exists, **this script never has to run again.**
# `fl_02_statewide.py` reads the cache and never opens a source file.
#
# ```
# raw_data/  ──▶  fl_01_preprocess.py  ──▶  <st>_state_layers.gpkg  ──▶  fl_02_statewide.py
#  (GB)            (slow, once)              (one file, ~900 MB)      (fast, repeatable)
# ```
#
# **Every switch this uses is in `fl_config.py`.** The one that matters here is
# `RUN_PREPROCESS`: `'auto'` builds only the layers that are missing.
#
# **What this is NOT.** It does not clip to a county, and it does not build a
# parcel table. The statewide analysis measures the whole state and finds
# parcels per county itself, so both of those were cost with no benefit.

# %% PP0 - setup
from fl_common import *            # noqa: F401,F403 - and fl_config through it

describe()
ensure_dirs()
check_raw_dir()
setup_matplotlib()

# %% [markdown]
# ## Step 1 — Find every input, and say what is missing
#
# Resolves every pattern in the registry to a real path **before** anything
# slow runs. A missing download fails here, in seconds, with the name of the
# file — not forty minutes later inside a spatial overlay.
#
# A missing *optional* file is reported with the exact consequence of its
# absence, so a degraded run is never mistaken for a complete one.
#
# **Validation:** `_preprocess/validation/input_inventory.csv`

# %% PP1 - the input inventory
rule('INPUT INVENTORY')

_inv = []
for _name, _spec in DATASETS.items():
    _path = resolve_dataset(_spec, _name)
    _spec['path'] = _path
    _ok = _path is not None and os.path.exists(_path)
    _layers = ''
    if _ok:
        _open = open_path(_path, verbose=True)
        if str(_open).lower().endswith(('.gpkg', '.gdb')):
            try:
                _all = sorted(n for n, _ in list_layers(_open))
                # does the file actually contain the LAYER this row asks for?
                # a present file with the wrong layer name fails later, not here
                _want = resolve_layer(_path, _spec)
                _mark = ('' if not _want else
                         ' [OK]' if _want in _all else f' [NO {_want!r}]')
                _layers = (', '.join(_all))[:52] + _mark
            except Exception:
                _layers = '(unreadable)'
    _inv.append({
        'dataset': _name,
        'found': 'yes' if _ok else ('MISSING' if _spec['required']
                                    else 'absent'),
        'MB': (round(os.path.getsize(_path) / 1e6, 1)
               if _ok and os.path.isfile(_path) else None),
        'file': os.path.basename(_path) if _ok else '(no match)',
        'layers': _layers})

inventory = pd.DataFrame(_inv)
print(inventory.to_string(index=False))
inventory.to_csv(os.path.join(P_VAL, 'input_inventory.csv'), index=False)

_missing_req = [r['dataset'] for r in _inv if r['found'] == 'MISSING']
_missing_opt = [r['dataset'] for r in _inv if r['found'] == 'absent']

if _missing_opt:
    print('\nOPTIONAL INPUTS ABSENT - here is exactly what that costs:')
    for _n in _missing_opt:
        print(f'  {_n:<12} {DATASETS[_n]["degrades"]}')

if _missing_req or _missing_opt:
    print('\nPATTERNS THAT MATCHED NOTHING:')
    for _n in _missing_req + _missing_opt:
        _p = DATASETS[_n]['pattern']
        print(f'  {_n:<12} {_p if isinstance(_p, str) else " | ".join(_p)}')
    show_raw_dir()
    print('\nIf you can SEE the file above, the pattern is wrong, not the '
          'data.\nEdit that row in DATASETS (fl_config.py §8) - the match is '
          'recursive and\ncase-insensitive, so the file name alone is usually '
          'enough.')

assert not _missing_req, (
    'required inputs missing: ' + ', '.join(_missing_req)
    + f'\nPut them in {RAW_DIR}, or fix the pattern in fl_config.py §8 using '
      f'the listing above.')
print(f'\nall required inputs present '
      f'({len(DATASETS) - len(_missing_opt)} of {len(DATASETS)} resolved)')

# %% [markdown]
# ## Step 2 — The study area
#
# Built first, because everything else is cut to it. Three pieces:
#
# | | what | how |
# |---|---|---|
# | `counties` | one row per county | dissolved up from the county-subdivision file |
# | `state_geom` | the state | union of the counties |
# | `pa_geo` | the priority area | **trimmed to the state**, and cached for reporting only |
#
# **Trimming the priority area to the state is the single most important line
# here.** The file ships covering many states; used untrimmed, every "precise"
# operation below runs against a multi-state geometry and does nothing useful.

# %% PP2 - the study area
rule('STUDY AREA')

_have = cache_layers()
_from_cache = RUN_PREPROCESS is False or (
    RUN_PREPROCESS == 'auto' and {'counties', 'cousub'} <= _have)

if _from_cache:
    counties = load_cached('counties')
    cousub = load_cached('cousub', required=False)
else:
    cousub = read_source('cousub')
    counties = cousub.dissolve(by='COUNTYFP', as_index=False)[
        ['COUNTYFP', 'geometry']]
    print(f'  {"counties":<12} {len(counties):>9,} dissolved from cousub')

# cousub carries the town names Pass C uses to decide whether an owner mails
# to an address INSIDE the county. Without it that split falls back to nothing
# and the absentee share is quietly wrong, so it is read from source whenever
# the cache does not have it.
if cousub is None:
    cousub = read_source('cousub')

state_geom = counties.geometry.union_all()
STATE_AC = state_geom.area / AC
print(f'  state area   {STATE_AC:,.0f} ac ({STATE_AC / 640:,.0f} sq mi) '
      f'across {len(counties)} counties')

# The priority area, cut to the state. REPORTING ONLY: it is cached so the
# statewide run can fill `pa_share` and `in_pa`, decide where results are
# reported (REPORT_WITHIN_PA) and - if you ever set COUNTY_EXTENT_MODE to a
# pa_* value - decide county extents. It never decides what ground
# preprocessing measures.
#
# THE FILE'S PRESENCE IS THE SWITCH. Like every other required=False dataset
# in the registry, the PA is used if raw_data/ holds a file matching the 'pa'
# pattern and ignored if it does not. There is no USE_PA_LAYER any more: it
# was a second on/off switch that no other dataset had, and fl_02 never
# honoured it anyway. To turn the priority area off, take the file out of
# raw_data/ and rebuild pa_geo.
# The priority area is the layer most likely to be REPLACED rather than
# re-downloaded - a new operational boundary arrives and the old file is
# overwritten. So it is rebuilt whenever the source file has changed, not only
# when it is missing from the cache. Nothing else in preprocessing depends on
# it: every resource layer is cut to the STATE, so a new priority area costs
# one small re-read and nothing else.
_PA_WHY = stale_reason('pa_geo', DATASETS['pa'])
_PA_STALE = ('pa_geo' not in _have
             or RUN_PREPROCESS is True
             or _PA_WHY is not None)
pa_geo = (load_cached('pa_geo', required=False)
          if 'pa_geo' in _have and not _PA_STALE else None)
if 'pa_geo' in _have and _PA_STALE:
    print(f'  {"pa_geo":<12} '
          f'{_PA_WHY or "RUN_PREPROCESS = True"} - rebuilding')
if pa_geo is None and DATASETS['pa'].get('path'):
    _pa = read_source('pa')
    if _pa is not None and len(_pa):
        print(f'  {"pa":<12} {len(_pa):>9,} features '
              f'({_pa.geometry.area.sum() / AC:,.0f} ac, multi-state)')
        pa_geo = gpd.clip(_pa, gpd.GeoSeries([state_geom], crs=WORKING_CRS),
                          keep_geom_type=True)
        pa_geo = pa_geo[~pa_geo.geometry.is_empty].reset_index(drop=True)
        _pac = pa_geo.geometry.area.sum() / AC
        print(f'  {"pa_geo":<12} {len(pa_geo):>9,} features, {_pac:,.0f} ac '
              f'inside the state ({100 * _pac / STATE_AC:.0f}%)')
elif pa_geo is None:
    # Not an error, and the reason matters: no FILE matched, as opposed to a
    # file that was read and turned out to hold nothing in this state.
    print(f'  {"pa":<12} no priority-area file matched in {RAW_DIR} - '
          f'pa_share and in_pa\n               will be 0 for every '
          f'{COUNTY_WORD.lower()}. Drop one in and re-run to enable them.')

# THE GROUND THE CACHE COVERS is the whole state. Simplified ONCE, up front: a
# spatial test costs vertices on BOTH sides, and the state outline is the
# complicated one.
MASK = state_geom.simplify(MASK_SIMPLIFY_M)
print(f'  clip mask    {shapely.get_num_coordinates(state_geom):,} vertices '
      f'-> {shapely.get_num_coordinates(MASK):,} after '
      f'simplify({MASK_SIMPLIFY_M} m)')
_b = state_geom.bounds
AOI_BBOX = (_b[0] - MI_M, _b[1] - MI_M, _b[2] + MI_M, _b[3] + MI_M)

# %% [markdown]
# ## Step 3 — Build the cached layers
#
# The slow cell, and the reason `RUN_PREPROCESS` exists. Each source file is
# read with its row filter and bounding box **pushed down to the driver**, cut
# to the state, given `acres` and `miles` columns, and written as one layer.
#
# On a first full run expect tens of minutes, most of it reading CSB. After
# that, `'auto'` finds every layer already cached and this takes seconds.

# %% PP3 - build the layers
rule('BUILDING THE CACHE')


MANIFEST = read_cache_manifest()


def build_layer(name, force=False):
    """Read one registered dataset, cut it to the state, cache it."""
    spec = DATASETS[name]
    if not force and name in cache_layers():
        _why = stale_reason(name, spec, MANIFEST)
        if _why:
            print(f'  {name:<12} {_why} - rebuilding')
        else:
            g = load_cached(name, required=spec['required'])
            # A REGISTRY THAT NOW ASKS FOR MORE COLUMNS IS A STALE LAYER
            # TOO, and layer_is_stale cannot see it: the SOURCE FILE has
            # not changed, only which of its fields we want. Raising
            # CROP_SEQ_YEARS from 1 to 3 would otherwise leave a
            # one-year cache in place and silently give no sequence.
            _short = []
            _ns = newest_spec(spec)
            if g is not None and _ns is not None:
                _al, _pat, _n = _ns
                if len(cdl_year_cols(g, _pat)) < _n:
                    _short = [f'{_al} years']
            if g is not None and len(g) and not _short:
                return g
            print(f'  {name:<12} '
                  + (f'cached copy is short of {", ".join(_short)} '
                     f'the registry now asks for - rebuilding'
                     if _short else 'cached copy is EMPTY - rebuilding'))

    g = read_source(name, bbox=AOI_BBOX)
    if g is None or not len(g):
        return g

    if spec['stage'] == 'layer':
        # simplify BEFORE the spatial op, not after: otherwise the test pays
        # the full vertex cost and the saving arrives too late to help
        g['geometry'] = g.geometry.simplify(SIMPLIFY_M)
        g = g[g.geometry.is_valid & ~g.geometry.is_empty]
        n0 = len(g)
        if CLIP_MODE == 'select':
            g = g[g.geometry.intersects(MASK)].copy()
            print(f'  {"":<12} in the state: {n0:,} -> {len(g):,} '
                  f'(whole features)')
        else:
            _pre = g
            g = gpd.clip(g, gpd.GeoSeries([MASK], crs=WORKING_CRS))
            g = g[~g.geometry.is_empty & g.geometry.notna()]
            g = keep_dimension(g, like=_pre, label=name)
            print(f'  {"":<12} clipped to the state: {n0:,} -> {len(g):,}')
        g = repair(g).reset_index(drop=True)

    # GeoPackage field names are case-insensitive, so a source 'ACRES' column
    # collides with the 'acres' computed below. Rename the source rather than
    # drop it - NWI's own acreage is worth keeping for comparison.
    for c in list(g.columns):
        if c != 'geometry' and c.lower() in ('acres', 'miles'):
            g = g.rename(columns={c: f'src_{c.lower()}'})
    g['acres'] = g.geometry.area / AC
    g['miles'] = g.geometry.length / MI_M

    if len(g):
        gpkg_safe(g).to_file(CACHE, layer=name, driver='GPKG')
        MANIFEST[name] = {'source': source_fingerprint(spec.get('path')),
                          'features': int(len(g)),
                          'built_at': datetime.datetime.now().isoformat(
                              timespec='seconds')}
        write_cache_manifest(MANIFEST)
        print(f'  {"":<12} cached as layer {name!r}')
    return g


# `pa` and `cousub` are handled above; `watersheds` is a boundary layer that
# is cached whole because it is context for the maps, not an analysis input.
_todo = [n for n, s in DATASETS.items() if n not in ('pa', 'cousub')]

if RUN_PREPROCESS is False:
    print('RUN_PREPROCESS = False - loading the cache, building nothing.')
    _build = []
elif RUN_PREPROCESS is True:
    print('RUN_PREPROCESS = True - rebuilding every layer.')
    _build = _todo
else:
    _build = [n for n in _todo if n not in cache_layers()]
    print(f'RUN_PREPROCESS = auto - '
          f'{len(cache_layers() & set(_todo))} layer(s) cached, '
          f'building {len(_build)}: {", ".join(_build) or "(none)"}')

_layers = {n: build_layer(n, force=n in _build) for n in _todo}

csb        = _layers.get('csb')
nhd_flow   = _layers.get('nhd_flow')
nhd_wb     = _layers.get('nhd_wb')
nhd_area   = _layers.get('nhd_area')
wetlands   = _layers.get('wetlands')
huc12      = _layers.get('huc12')
watersheds = _layers.get('watersheds')

# The three boundary layers the statewide script needs, written into the SAME
# file. One cache means no second path to get wrong, and no "which file has
# pa_geo in it" question six weeks from now.
for _nm, _g, _src in [('counties', counties, 'cousub'),
                      ('cousub', cousub, 'cousub'),
                      ('pa_geo', pa_geo, 'pa')]:
    if _g is None or not len(_g):
        continue
    if _nm in cache_layers() and not _build and not (
            _nm == 'pa_geo' and _PA_STALE):
        continue
    gpkg_safe(_g).to_file(CACHE, layer=_nm, driver='GPKG')
    MANIFEST[_nm] = {'source': source_fingerprint(DATASETS[_src].get('path')),
                     'features': int(len(_g)),
                     'built_at': datetime.datetime.now().isoformat(
                         timespec='seconds')}
    write_cache_manifest(MANIFEST)
    print(f'  {_nm:<12} cached as layer {_nm!r}')

print('\nin memory:')
for _n in ['counties', 'cousub', 'pa_geo', 'watersheds', 'csb', 'nhd_flow',
           'nhd_wb', 'nhd_area', 'wetlands', 'huc12']:
    _g = globals().get(_n)
    print(f'  {_n:<12} {"-" if _g is None else f"{len(_g):,} features"}')

if csb is not None:
    print(f'\ncropland: {csb["acres"].sum():,.0f} ac across {len(csb):,} '
          f'fields' + ('' if 'CDL' in csb.columns
                       else '   (NO CDL COLUMN - crop tables will be skipped)'))
if nhd_flow is not None:
    print(f'streams : {nhd_flow["miles"].sum():,.0f} mi across '
          f'{len(nhd_flow):,} reaches')

# %% [markdown]
# ## Step 4 — Validation: did preprocessing lose anything?
#
# Clipping and filtering are destructive by design, so the only honest thing to
# do is measure what they removed. Three checks, all written to
# `_preprocess/validation/`.

# %% PP4 - validation
rule('PREPROCESSING VALIDATION')

# CHECK 1. Every layer the statewide analysis needs is present and non-empty.
# This is the check that catches a GeoPackage which was written but did not
# survive - the failure that used to surface eleven lines into the next script
# as a confusing KeyError.
_req = {'counties': counties, 'csb': csb, 'nhd_flow': nhd_flow,
        'huc12': huc12}
_rows = []
for _n in list(_req) + ['cousub', 'pa_geo', 'nhd_wb', 'nhd_area', 'wetlands',
                        'watersheds']:
    _g = _req.get(_n, globals().get(_n))
    _rows.append({
        'layer': _n, 'required': _n in _req, 'present': _g is not None,
        'features': 0 if _g is None else len(_g),
        'empty_geom': 0 if _g is None else int(_g.geometry.is_empty.sum()),
        'invalid': 0 if _g is None else int((~_g.geometry.is_valid).sum())})
layer_check = pd.DataFrame(_rows)
print(layer_check.to_string(index=False))
layer_check.to_csv(os.path.join(P_VAL, 'layer_check.csv'), index=False)

for _r in _rows:
    if _r['required'] and not _r['present']:
        raise AssertionError(f'{_r["layer"]} is required and absent')
    if _r['required'] and _r['features'] == 0:
        raise AssertionError(
            f'{_r["layer"]} came back EMPTY. Check that the source overlaps '
            f'the state,\n  and that its `where` filter in fl_config.py §8 is '
            f'right for this vintage.')

# CHECK 2. Does each layer REACH the state? A layer that survives the cut but
# sits in one corner means the wrong tile was downloaded - which no row count
# will tell you.
print('\ncoverage of the state:')
_cov = []
for _n in ['csb', 'nhd_flow', 'nhd_wb', 'nhd_area', 'wetlands', 'huc12']:
    _g = globals().get(_n)
    if _g is None or not len(_g):
        continue
    _bb = shp_box(*_g.total_bounds)
    _cov.append({'layer': _n,
                 'bbox_covers_state_pct': round(
                     100 * _bb.intersection(state_geom).area
                     / state_geom.area, 1)})
coverage = pd.DataFrame(_cov)
print(coverage.to_string(index=False))
coverage.to_csv(os.path.join(P_VAL, 'state_coverage.csv'), index=False)
if len(coverage) and (coverage['bbox_covers_state_pct'] < 95).any():
    print('  WARNING: a layer does not span the state. Fine if the resource '
          'genuinely\n  is not there statewide; a wrong-tile bug if it is.')

# CHECK 3. Do the counties the cache reaches look like the whole state? This
# is the number the statewide run turns into `source_cover`, and a county the
# cache does not reach is reported as no_source_data rather than as empty.
if csb is not None and len(csb):
    shapely.prepare(state_geom)
    _cg = np.asarray(counties.geometry.values, dtype=object)
    _hit = np.zeros(len(counties), bool)
    _tree = shapely.STRtree(np.asarray(csb.geometry.values, dtype=object))
    _qi, _si = _tree.query(_cg, predicate='intersects')
    _hit[np.unique(_qi)] = True
    print(f'\ncounties with cropland in the cache: {int(_hit.sum())} of '
          f'{len(counties)}')
    if _hit.sum() < len(counties):
        print('  The rest have no CSB fields. That can be real (a fully '
              'forested county)\n  or a clipped source - the statewide run '
              'reports them as no_source_data\n  and keeps them out of every '
              'average either way.')

# %% [markdown]
# ## Step 5 — One picture of what got cached
#
# So a wrong clip is **visible** rather than inferred from a table.
#
# **Figure:** `_preprocess/figures/cached_inputs.png`

# %% PP5 - the overview figure
fig, ax = plt.subplots(figsize=(13, 9))
ax.set_xticks([]); ax.set_yticks([])
for _s in ax.spines.values():
    _s.set_visible(False)
ax.set_aspect('equal')

_leg = []
if csb is not None and len(csb):
    csb.plot(ax=ax, facecolor=CROP_C, edgecolor='none', zorder=2)
    _leg.append(Patch(facecolor=CROP_C,
                      label=f'Cropland — {csb["acres"].sum():,.0f} ac'))
if wetlands is not None and len(wetlands):
    wetlands.plot(ax=ax, facecolor=WET_C, edgecolor='none', alpha=0.6,
                  zorder=3)
    _leg.append(Patch(facecolor=WET_C, alpha=0.6,
                      label=f'NWI wetlands — {len(wetlands):,}'))
if nhd_wb is not None and len(nhd_wb):
    nhd_wb.plot(ax=ax, facecolor=WB_C, edgecolor='none', zorder=4)
    _leg.append(Patch(facecolor=WB_C,
                      label=f'NHD waterbodies — {len(nhd_wb):,}'))
if nhd_flow is not None and len(nhd_flow):
    nhd_flow.plot(ax=ax, color=FLOW_C, linewidth=0.12, alpha=0.7, zorder=5)
    _leg.append(Line2D([], [], color=FLOW_C, lw=1.4,
                       label=f'NHD flowlines — {nhd_flow["miles"].sum():,.0f} mi'))
if huc12 is not None and len(huc12):
    huc12.boundary.plot(ax=ax, color=MUTED, linewidth=0.25, zorder=6)
    _leg.append(Line2D([], [], color=MUTED, lw=1.0,
                       label=f'HUC12 subwatersheds — {len(huc12):,}'))
if pa_geo is not None and len(pa_geo):
    pa_geo.boundary.plot(ax=ax, color=OVER_C, linewidth=1.4, zorder=7)
    _leg.append(Line2D([], [], color=OVER_C, lw=1.6,
                       label='Priority area (reporting only)'))
counties.boundary.plot(ax=ax, color=INK2, linewidth=0.4, zorder=8)
_leg.append(Line2D([], [], color=INK2, lw=1.0,
                   label=f'Counties — {len(counties)}'))

# blank sky for the legend, so it covers no data
_y0, _y1 = ax.get_ylim()
ax.set_ylim(_y0, _y1 + 0.24 * (_y1 - _y0))
_lg = ax.legend(handles=_leg, loc='upper left', frameon=True,
                fontsize=FS['legend'], borderpad=0.8, labelspacing=0.5,
                title='Cached layers', title_fontsize=FS['legend_title'])
_lg.get_frame().set_facecolor(SURFACE)
_lg.get_frame().set_edgecolor(MUTED)
_lg.get_title().set_color(INK2)
_lg.get_title().set_fontweight('bold')
ax.set_title(f'Preprocessed inputs — {STATE_NAME}',
             fontsize=FS['fig_title'], loc='left', pad=12)
ax.set_xlabel(f'{os.path.basename(CACHE)} · '
              f'{os.path.getsize(CACHE) / 1e6:,.0f} MB · '
              f'{STATE_AC / 640:,.0f} sq mi',
              fontsize=FS['caption'], color=INK2)
_p = os.path.join(P_FIG, 'cached_inputs.png')
fig.savefig(_p, dpi=FIG_DPI, bbox_inches='tight')
show_close(fig)
print(f'wrote {_p}')

# %% PP6 - done
rule('PREPROCESSING COMPLETE')
print(f'  cache   {CACHE}')
print(f'          {os.path.getsize(CACHE) / 1e6:,.1f} MB')
print(f'          layers: {", ".join(sorted(cache_layers()))}')
print(f'\n  Set RUN_PREPROCESS = False in fl_config.py to skip all of the '
      f'above next time.')
print(f'  Next: python fl_02_statewide.py')

# If the cache was built on scratch disk, get it somewhere durable NOW - not
# at the end of the NEXT script, six hours later, on a machine that might have
# rebooted in between.
if SCRATCH_DIR and os.path.abspath(CACHE).startswith(
        os.path.abspath(SCRATCH_DIR)):
    os.makedirs(FINAL_CLIPPED_DIR, exist_ok=True)
    _final_cache = os.path.join(FINAL_CLIPPED_DIR, CACHE_NAME)
    _jobs = [(CACHE, _final_cache, 'cache')]
    # THE SIDECAR GOES WITH IT. Copying the GeoPackage alone is what
    # made a restored cache rebuild all nine layers: the durable copy
    # had no record of what built it. It is a few KB against ~900 MB.
    if os.path.exists(CACHE_MANIFEST):
        _jobs.append((CACHE_MANIFEST,
                      os.path.splitext(_final_cache)[0] + '.cache.json',
                      'cache manifest'))
    copy_to_final(_jobs)
