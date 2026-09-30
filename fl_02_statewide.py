# ---
# jupytext:
#   text_representation:
#     extension: .py
#     format_name: percent
# ---

# %% [markdown]
# # `fl_02_statewide.py` — the statewide riparian-hotspot analysis
#
# Reads the cache `fl_01_preprocess.py` built and produces the whole
# deliverable. **It never opens a source file**, which is why it is repeatable
# in seconds once Pass A is done.
#
# **Three passes.**
#
# | pass | runs | what it does | cost |
# |---|---|---|---|
# | **A** | once | one statewide density surface: corridor → riparian cropland → 1-mile grid → 3-mile band density | hours |
# | **B** | 120×, cheap | each county's slice of that surface, its zones, its rank, its row in the manifest | seconds |
# | **C** | N×, moderate | the parcel work, for the counties that have a parcel file | minutes each |
#
# **Why Pass A can exist at all.** The hotspot method never looks at a parcel —
# zones come from cropland and water alone — and it never looks at a county
# except to decide the extent. So running it on one geometry or on a state is
# the same code with a different polygon. That is also what makes the
# free-vs-proprietary parcel comparison a real experiment: **the zone map costs
# nothing in parcel licensing.**
#
# **What running it statewide fixes.** In a county-by-county run the 3-mile
# neighbourhood is cut off at the county line, so every border cell is measured
# over less land than an interior one and reads low for an administrative
# reason. One statewide surface removes that for every interior county.
#
# ---
#
# **Every switch is in `fl_config.py`.** The two that decide what this run
# costs:
#
# * `RUN_SURFACE` — leave on `'auto'`. Once `state_surface.gpkg` exists,
#   changing `BASIS`, the floor, a colour or a figure re-runs in **seconds**.
#   Set it to `True` only when the cache, `RIP_BUFFER_FT`, `CELL_MI`,
#   `BAND_MI` or `TILE_KM` changes.
# * `COUNTIES` — a short list of FIPS runs the whole pipeline in minutes, which
#   is how you test a change before committing to 120 counties.

# %% S0 - setup
from fl_common import *            # noqa: F401,F403 - and fl_config through it

describe()
ensure_dirs()
setup_matplotlib()

if not os.path.exists(CACHE):
    raise FileNotFoundError(
        f'no cache at {CACHE}\n'
        f'  Run fl_01_preprocess.py first. This script reads only what\n'
        f'  preprocessing leaves in clipped_layers/.')


# %% [markdown]
# ## Pass A — one statewide surface
#
# Everything in this part runs **once**, whatever the threshold is set to. It
# produces `_state/output/state_surface.gpkg`, which holds the 1-mile grid with
# its 3-mile neighbourhood density already computed. Thresholds, zones,
# counties and figures all read that file, so trying a different basis or floor
# costs seconds, not hours.

# %% B1 - read the cache and work out what it actually covers
rule('PASS A - INPUTS')

# The cache has to hold the ANALYSIS layers, not merely exist. A file with
# only `parcels_enriched` in it is a Section 2 that did not finish, or a
# GeoPackage that lost its layers on the way to disk - and the old check, which
# only asked whether the file had ANY layers, let that through to fail eleven
# lines later on a confusing KeyError.
REQUIRED_LAYERS = {'counties', 'csb', 'nhd_flow', 'huc12'}
OPTIONAL_LAYERS = {'cousub', 'nhd_wb', 'nhd_area', 'wetlands', 'watersheds',
                   'pa_geo'}
_have = cache_layers()
_missing = REQUIRED_LAYERS - _have
if _missing:
    _p = os.path.abspath(CACHE)
    while _p != os.path.dirname(_p) and not os.path.isdir(_p):
        _p = os.path.dirname(_p)
    print(f'\n{"=" * 74}\nTHE CACHE IS NOT USABLE\n{"=" * 74}')
    print(f'  wanted     : {CACHE}')
    print(f'  exists?    : {os.path.exists(CACHE)}')
    if os.path.exists(CACHE):
        _sz = os.path.getsize(CACHE) / 1e6
        _mt = datetime.datetime.fromtimestamp(
            os.path.getmtime(CACHE)).isoformat(timespec="seconds")
        print(f'  size       : {_sz:,.1f} MB    last written: {_mt}')
    print(f'  layers in it   : {", ".join(sorted(_have)) or "(none)"}')
    print(f'  layers REQUIRED: {", ".join(sorted(REQUIRED_LAYERS))}')
    print(f'  MISSING        : {", ".join(sorted(_missing))}')
    print(f'  also expected  : '
          f'{", ".join(sorted(OPTIONAL_LAYERS - _have)) or "(all present)"}')
    if not os.path.exists(CACHE):
        print(f'\n  deepest folder that does exist: {_p}')
        try:
            _k = sorted(os.listdir(_p))[:20]
            print(f'  it contains: {", ".join(_k) if _k else "(nothing)"}')
        except Exception:
            pass
    if _have and not _have & REQUIRED_LAYERS:
        # the damaging case: the file is there and has SOMETHING in it, so it
        # looks like a cache, but nothing the analysis needs survived
        print('\n  The file exists and has layers, but none of the ones the\n'
              '  analysis needs. That is a preprocessing run that did not\n'
              '  finish, or one whose writes did not survive. Re-run\n'
              '  fl_01_preprocess.py with RUN_PREPROCESS = True and watch for\n'
              '  a "cached as layer" line for every REQUIRED layer above.')
    raise FileNotFoundError(
        f'{CACHE} is missing {", ".join(sorted(_missing))}.\n'
        f'  If the path looks wrong, fix PROJECT_DIR / CACHE_NAME in '
        f'fl_config.py section 1.\n'
        f'  If it looks right, run fl_01_preprocess.py - this script never '
        f'opens a\n  source file, it reads what preprocessing leaves behind.')
print(f'cache layers: {", ".join(sorted(_have))}')

counties = load_cached('counties')
cousub   = load_cached('cousub', required=False)
pa_geo    = load_cached('pa_geo', required=False)
csb      = load_cached('csb')
nhd_flow = load_cached('nhd_flow')
nhd_wb   = load_cached('nhd_wb', required=False)
nhd_area = load_cached('nhd_area', required=False)
wetlands = load_cached('wetlands', required=False)
huc12    = load_cached('huc12', required=False)

counties['FIPS'] = counties['COUNTYFP'].astype(str).str.zfill(3)
counties = counties.dissolve(by='FIPS', as_index=False)[['FIPS', 'geometry']]
counties['county_ac'] = counties.geometry.area / AC

# Step 1. County names, for whatever state this is. Nothing MEASURED depends on
# them - they are labels, folder names and the string Pass C matches parcel
# files against - so the cascade is cheapest-first and never fatal:
#
#   cousub, if this state's file carries a county-name column (most do not)
#   CSB's CNTY, already in the cache and free
#   a TIGER county file in raw_data/, if one was downloaded
#   the Census national county list, fetched once and cached (COUNTY_NAME_FETCH)
#   COUNTY_NAMES from the config, for an offline machine or a correction
#   'County 047', so a nameless county still gets a row and a folder
_names, _name_src = {}, []
_cf = county_name_field(cousub)
if _cf is not None:
    _cn = (cousub[['COUNTYFP', _cf]].dropna().astype(str)
           .assign(COUNTYFP=lambda d: d['COUNTYFP'].str.zfill(3))
           .drop_duplicates('COUNTYFP'))
    _names = dict(zip(_cn['COUNTYFP'], _cn[_cf].str.title()))
    if _names:
        _name_src.append(f'cousub.{_cf} ({len(_names)})')
if {'CNTYFIPS', 'CNTY'} <= set(csb.columns):
    _n = (csb[['CNTYFIPS', 'CNTY']].dropna().astype(str)
          .assign(CNTYFIPS=lambda d: d['CNTYFIPS'].str.zfill(3))
          .drop_duplicates('CNTYFIPS'))
    _n0 = len(_names)
    for _k, _v in zip(_n['CNTYFIPS'], _n['CNTY'].str.title()):
        _names.setdefault(_k, _v)
    if len(_names) > _n0:
        _name_src.append(f'CSB.CNTY (+{len(_names) - _n0})')
# CSB is clipped to the cache scope, so it names only the counties the cache
# reaches. A TIGER county file fills the rest if one was downloaded.
_hit = sorted(glob.glob(os.path.join(RAW_DIR, '**', 'tl_*county*.shp'),
                        recursive=True))
for _h in _hit:
    try:
        _c = gpd.read_file(_h)
    except Exception:
        continue
    if 'COUNTYFP' not in _c.columns or 'NAME' not in _c.columns:
        continue
    if 'STATEFP' in _c.columns:
        _c = _c[_c['STATEFP'].astype(str).str.zfill(2) == STATE_FIPS]
    _n0 = len(_names)
    for _k, _v in zip(_c['COUNTYFP'].astype(str).str.zfill(3), _c['NAME']):
        _names.setdefault(_k, str(_v))
    if len(_names) > _n0:
        _name_src.append(f'{os.path.basename(_h)} (+{len(_names) - _n0})')
    break

# Anything still unnamed: the Census list (fetched once, then cached), then
# whatever you put in COUNTY_NAMES.
_need = set(counties['FIPS']) - set(_names)
if _need and COUNTY_NAME_FETCH:
    _n0 = len(_names)
    for _k, _v in census_county_names().items():
        _names.setdefault(_k, _v)
    if len(_names) > _n0:
        _name_src.append(f'Census list (+{len(_names) - _n0})')
for _k, _v in (COUNTY_NAMES or {}).items():
    _names[str(_k).zfill(3)] = _v          # your table WINS - it is the override

counties['NAME'] = counties['FIPS'].map(_names).fillna(
    COUNTY_WORD + ' ' + counties['FIPS'])
_unnamed = int(counties['NAME'].str.startswith(COUNTY_WORD + ' ').sum())
print(f'counties    : {len(counties)}  '
      f'[names: {", ".join(_name_src) if _name_src else "none found"}]'
      + (f' ({_unnamed} without a name)' if _unnamed else ''))

STATE_GEOM = counties.geometry.union_all()
STATE_AC = STATE_GEOM.area / AC
print(f'state area  : {STATE_AC:,.0f} ac ({STATE_AC / 640:,.0f} sq mi)')

# Step 2. What ground do the cached layers actually cover? The notebook's cache
# is scoped to the priority area, so 68 counties have no cropland or water in
# it. Reporting those as "no hotspots" would be a lie; they are reported as
# uncovered, and are kept out of the state and region averages.
# ONE CACHE, AND IT IS STATEWIDE BY CONSTRUCTION. fl_01_preprocess.py cuts
# every layer to the state and writes them all - including `pa_geo` - into the
# same file, so there is no second path to get wrong and no "which file has the
# priority area in it" question six weeks from now.
#
# `pa_geo` here is REPORTING ONLY: it fills in_pa / pa_share, and it decides
# county extents only if COUNTY_EXTENT_MODE is set to a pa_* value. Reporting
# INSIDE the priority area is a different switch - REPORT_WITHIN_PA, applied
# in B6, below - and that one costs nothing.
# It never decides what ground was measured.
SOURCE_GEOM = STATE_GEOM
SOURCE_SCOPE = 'state'          # carried into state_surface.json as provenance
SOURCE_AC = SOURCE_GEOM.area / AC
print(f'cache scope : statewide - {SOURCE_AC:,.0f} ac (100% of the state)')
if pa_geo is not None and len(pa_geo):
    print(f'priority area: {pa_geo.geometry.area.sum() / AC:,.0f} ac '
          f'({100 * pa_geo.geometry.area.sum() / AC / STATE_AC:.0f}% of the '
          f'state) - reporting only')
else:
    print('priority area: not in the cache - in_pa and pa_share will be 0.\n'
          '               The PA is enabled by being there: put a file '
          'matching the \'pa\'\n               pattern in raw_data/ and '
          're-run fl_01_preprocess.py.')

shapely.prepare(SOURCE_GEOM)
_cg = np.asarray(counties.geometry.values, dtype=object)
counties['cache_ac'] = shapely.area(
    shapely.intersection(_cg, SOURCE_GEOM)) / AC
counties['cache_share'] = (counties['cache_ac'] /
                           counties['county_ac']).round(4)
print(f'counties with any ground in the cache: '
      f'{int((counties["cache_share"] > 0.001).sum())} of {len(counties)}; '
      f'{int((counties["cache_share"] >= 0.999).sum())} wholly inside it')
# Whether a county can be ANALYSED is not this number - it is whether the
# cache covers the county's own analysed EXTENT, which for a county clipped to
# the priority area is a different question. That is settled in B2, once the
# extents exist.


# %% B2 - the per-county extents and the analysed ground
# Step 3. Each county's extent, by the rule COUNTY_EXTENT_MODE sets. The basis
# is stored per county because it changes what county_avg divides by, and two
# counties on different bases are not directly comparable on that number.
_PA_GEOM = (pa_geo.geometry.union_all() if pa_geo is not None and len(pa_geo)
            else None)
if _PA_GEOM is not None:
    shapely.prepare(_PA_GEOM)
    _pa_ac = shapely.area(shapely.intersection(_cg, _PA_GEOM)) / AC
else:
    _pa_ac = np.zeros(len(counties))
counties['pa_ac'] = _pa_ac
counties['pa_share'] = (counties['pa_ac'] / counties['county_ac']).round(4)
# `intersects` would say yes to a county that only TOUCHES the boundary, and
# its "extent" would be a zero-area line. The share is the honest test.
counties['in_pa'] = counties['pa_share'] >= MIN_PA_SHARE
if _PA_GEOM is not None:
    _edge = int(((counties['pa_share'] > 0) & ~counties['in_pa']).sum())
    if _edge:
        print(f'  {_edge} county/counties overlap the priority area by less '
              f'than {MIN_PA_SHARE:.0%} - treated as outside it')

_ext, _basis = [], []
for _i, _row in counties.iterrows():
    _g = _row.geometry
    if COUNTY_EXTENT_MODE == 'whole_county' or _PA_GEOM is None:
        _ext.append(_g)
        _basis.append('whole county')
    elif _row['in_pa']:
        _ext.append(_g.intersection(_PA_GEOM))
        _basis.append('county n priority area')
    elif COUNTY_EXTENT_MODE == 'pa_where_available':
        # LEGACY - see fl_config.py section 3. Kept so a run made under
        # it still reproduces; the supported way to report on a priority
        # area is REPORT_WITHIN_PA, in B6, which needs no rebuild.
        _ext.append(_g)
        _basis.append('whole county')
    else:
        _ext.append(shp_box(0, 0, 0, 0))
        _basis.append('excluded - outside the priority area')
counties['extent'] = _ext
counties['extent_basis'] = _basis
counties['extent_ac'] = area_ac(counties['extent'].values)
counties = counties[counties['extent_ac'] > 0].reset_index(drop=True)

# Step 3a. Coverage is measured against the county's own EXTENT, not against
# the whole county. A county clipped to the priority area has an extent that
# lies inside the cache by construction and is fully analysable, even though
# most of the county is not in the cache at all. Measuring against the county
# would blank a real result; measuring against the extent asks the question
# that matters - was every acre this county is judged on actually looked at?
_eg0 = np.asarray(counties['extent'].values, dtype=object)
counties['source_ac'] = shapely.area(
    shapely.intersection(_eg0, SOURCE_GEOM)) / AC
counties['source_cover'] = (counties['source_ac'] /
                            counties['extent_ac']).round(4)
counties['covered'] = counties['source_cover'] >= MIN_SOURCE_COVER
print(f'\nanalysable  : {int(counties["covered"].sum())} of {len(counties)} '
      f'counties have at least {MIN_SOURCE_COVER:.0%} of their extent inside '
      f'the cache')
_part = int(((counties['source_cover'] > 0.001) & ~counties['covered']).sum())
if _part:
    print(f'              {_part} more are PARTLY covered - reported as '
          f'no_source_data, because a\n              county average over ground '
          f'that was never measured is not an average')
if not counties['covered'].all():
    print('  A county the cache does not reach is reported as no_source_data '
          'and kept\n  OUT of every average. If that surprises you, check the '
          'coverage table in\n  _preprocess/validation/ - a source layer that '
          'was clipped before you got it\n  looks exactly like a county with '
          'no resource in it.')

if COUNTIES != 'all':
    _want = {str(c).zfill(3) for c in COUNTIES}
    counties = counties[counties['FIPS'].isin(_want)].reset_index(drop=True)
    print(f'COUNTIES filter: {len(counties)} of the requested {len(_want)} found')

print(f'\nextent basis:')
print(counties.groupby('extent_basis')['FIPS'].size().to_string())

# Step 4. The ground Pass A actually measures: the union of the county extents,
# cut to what the cache covers. Measuring outside the cache would mix real
# zeros with missing data, and nothing downstream could tell them apart.
ANALYSIS_GEOM = union_all(
    np.asarray(counties['extent'].values, dtype=object)).intersection(SOURCE_GEOM)
ANALYSIS_AC = ANALYSIS_GEOM.area / AC
# Each county's share of that ground. This - not the county extent - is what
# zones are attributed against, so a county the cache never reached cannot be
# credited with a zone that merely laps over its boundary.
counties['analysed_geom'] = shapely.intersection(
    np.asarray(counties['extent'].values, dtype=object), ANALYSIS_GEOM)
print(f'\nanalysed ground: {ANALYSIS_AC:,.0f} ac '
      f'({ANALYSIS_AC / 640:,.0f} sq mi), '
      f'{100 * ANALYSIS_AC / STATE_AC:.0f}% of {STATE_NAME}')


# %% B3 - the water corridor and the riparian cropland, one tile at a time
# Step 5. The corridor is every stream, river, waterbody and wetland buffered
# by RIP_BUFFER_FT and unioned, so a stream running through a wetland is
# counted once. Statewide that is one geometry with millions of vertices, and
# building it whole is what would put this on a bigger machine.
#
# Instead the state is cut into TILE_KM squares. Each tile buffers only the
# water within a HALO_M margin of itself, overlays the cropland, and keeps the
# result clipped to the tile alone. The halo is fifty times the buffer width,
# so a corridor near a tile edge is complete before the clip removes the part
# belonging to the neighbour - the tiling changes how many PIECES a field is
# split into, never how many acres it has.
def _tiles(geom, km):
    x0, y0, x1, y1 = geom.bounds
    step = km * 1000.0
    out = []
    for i in range(int(np.ceil((x1 - x0) / step))):
        for j in range(int(np.ceil((y1 - y0) / step))):
            out.append(shp_box(x0 + i * step, y0 + j * step,
                               x0 + (i + 1) * step, y0 + (j + 1) * step))
    keep = [t for t in out if geom.intersects(t)]
    return keep


def corridor_for(bbox, verbose=False):
    """Per-source corridors plus their union, for the water inside one bbox."""
    parts = {}
    for name in USE_SOURCES:
        spec = WATER_SOURCES[name]
        src = globals().get(spec['var'])
        if src is None:
            continue
        g = src.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        g = g[(~g.geometry.is_empty) & (g.geometry.area + g.geometry.length > 0)]
        if name == 'wetland' and WETLAND_EXCLUDE and len(g):
            g = drop_nwi_types(g, verbose=verbose)
        if not len(g):
            continue
        parts[name] = corridor_one(g, spec['kind'], spec['interior'], RIP_BUFFER_FT)
    if not parts:
        return {}, None
    return parts, union_all(np.asarray(list(parts.values()), dtype=object),
                            grid_size=0.01)


# %% B3a - what the surface on disk was built from
# Pass A is the only expensive stage, so 'auto' reuses state_surface.gpkg when
# it is there. That is only safe if the surface was built under the SAME
# assumptions as this run, and several settings that change it are nowhere near
# Pass A in the file. The worst of them is the county filter: it is applied in
# B2, ABOVE, so `FL_COUNTIES=047,225` - the two-county test the README
# recommends - writes a two-county surface. Reusing that for a statewide run
# would report the whole state out of 2/120ths of it and never say so.
#
# So the surface carries a signature of everything it depends on, and a
# mismatch is either rebuilt ('auto') or refused (False). Note what is NOT in
# here: BASIS, the floor, the palette, TOP_N_CONTACTS, anything in Pass B or C.
# Those are downstream of the surface, and changing them must stay free.
# WHICH CACHE LAYERS PASS A ACTUALLY READS. `pa_geo` is deliberately not one
# of them: under COUNTY_EXTENT_MODE = 'whole_county' it fills in_pa and
# pa_share and decides nothing the surface measures, and under the B6
# reporting mask it decides only what is REPORTED, which is downstream of the
# surface either way.
_SURFACE_LAYERS = ('counties', 'cousub', 'csb', 'nhd_flow', 'nhd_wb',
                   'nhd_area', 'wetlands', 'huc12')


def _cache_fingerprint():
    """The cache by what Pass A reads out of it, rather than by its size.

    THE BUG THIS REPLACES. The signature used to carry the GeoPackage's byte
    count, which moves whenever ANY layer inside it is rewritten - including
    pa_geo, which Pass A never opens. So dropping in a new priority area read
    as a rebuilt cache and cost a five-hour Pass A for a layer that changes no
    measurement.

    The preprocessing manifest already records, per layer, the source file it
    was built from (path, mtime, size) and how many features it produced.
    Hashing those for the layers above is both a STRONGER identity than a byte
    count - a same-size edit to CSB moves it, where the size would not - and a
    NARROWER one, because pa_geo is not in it. `built_at` is deliberately
    excluded: re-running preprocessing over unchanged sources rewrites the
    timestamp and rebuilds nothing, and that must not cost a Pass A either.

    With no manifest there is nothing to hash, so it falls back to the size,
    which is what earlier runs recorded.
    """
    man = read_cache_manifest()
    rows = {k: {'source': man[k].get('source'),
                'features': man[k].get('features')}
            for k in _SURFACE_LAYERS if k in man}
    if not rows:
        return {'by': 'size',
                'bytes': (os.path.getsize(CACHE) if os.path.exists(CACHE)
                          else None)}
    return {'by': 'layers', 'n': len(rows),
            'sha': hashlib.sha1(json.dumps(rows, sort_keys=True,
                                           default=str).encode()
                                ).hexdigest()[:16]}


def surface_signature():
    """Every setting the density surface is a function of."""
    _f = sorted(counties['FIPS'].astype(str))
    sig = {
        # THE CACHE BY WHAT PASS A READS OUT OF IT - not by where it is,
        # and not by how big it is. Copying it to a fast local disk is the
        # SAME cache, so the path is out; storing it once threw away a
        # finished Pass A. Its SIZE is out for the same kind of reason one
        # level down - see _cache_fingerprint, above.
        'cache': os.path.basename(CACHE),
        'cache_layers': _cache_fingerprint(),
        'source_scope': SOURCE_SCOPE,
        'county_extent_mode': COUNTY_EXTENT_MODE,
        'min_source_cover': float(MIN_SOURCE_COVER),
        'n_counties': len(_f),
        'counties': hashlib.sha1(','.join(_f).encode()).hexdigest()[:12],
        'analysis_ac': round(float(ANALYSIS_AC), 1),
        'cell_mi': float(CELL_MI),
        'band_mi': float(BAND_MI),
        'tile_km': float(TILE_KM),
        'buffer_ft': float(RIP_BUFFER_FT),
        'sources': ','.join(sorted(USE_SOURCES)),
        'wetland_exclude': ','.join(sorted(WETLAND_EXCLUDE or [])),
    }
    # The priority area only enters the surface when it decides an extent.
    # Under 'whole_county' a new PA file changes the reporting and nothing
    # else, and a rebuild would be five hours spent on no difference.
    if COUNTY_EXTENT_MODE != 'whole_county' and _PA_GEOM is not None:
        sig['min_pa_share'] = float(MIN_PA_SHARE)
        # The priority area by its SHAPE, not by its file. Re-downloading the
        # same boundary, or re-saving it, changes the file's mtime and size and
        # nothing that matters; area, perimeter, bounds and part count together
        # do not move unless the extent really moved. A five-hour rebuild
        # should cost a real edit, not a touch.
        sig['pa_extent'] = [round(_PA_GEOM.area / AC, 1),
                            round(_PA_GEOM.length / MI_M, 3),
                            int(shapely.get_num_geometries(_PA_GEOM))] + [
                                round(_v, 1) for _v in _PA_GEOM.bounds]
    return sig


_SIG_TOL = {'analysis_ac': 1e-3}        # relative; float noise is not a change

# Keys stored as a file path by an older run, and compared by basename now.
# A surface built before the cache was identified by name rather than by its
# full path would otherwise read as a mismatch on every machine but the one
# that built it.
_SIG_BASENAME = {'cache'}


def surface_diff(stored, now):
    """(differs, not_recorded). A key the older run never wrote cannot be
    checked - that is a caveat, not a mismatch."""
    differs, missing = [], []
    for k, v in now.items():
        if k not in stored:
            missing.append(k)
            continue
        w = stored[k]
        if k in _SIG_BASENAME and isinstance(w, str):
            w = os.path.basename(w)
        if (k in _SIG_TOL and isinstance(w, (int, float))
                and isinstance(v, (int, float))
                and not isinstance(w, bool)):
            if abs(float(v) - float(w)) <= _SIG_TOL[k] * max(abs(float(v)), 1):
                continue
        elif w == v:
            continue
        differs.append((k, w, v))
    return differs, missing


_SURF_META, _SURF_DIFF, _SURF_MISSING = None, [], []
if os.path.exists(SURFACE_GPKG) and os.path.exists(SURFACE_JSON):
    with open(SURFACE_JSON) as _f:
        _SURF_META = json.load(_f)
    # Older runs wrote the fields flat; newer ones nest them under
    # 'signature'. Merge, so a surface built by either is comparable.
    _stored = dict(_SURF_META)
    _stored.update(_SURF_META.get('signature') or {})
    _SURF_DIFF, _SURF_MISSING = surface_diff(_stored, surface_signature())
elif os.path.exists(SURFACE_GPKG):
    print(f'\nNOTE: {os.path.basename(SURFACE_GPKG)} is on disk but its '
          f'state_surface.json is not.\n      Nothing records what it was '
          f'built from, so it is rebuilt rather than trusted.')

_want_surface = RUN_SURFACE is True or (
    RUN_SURFACE == 'auto' and _SURF_META is None)
_sig_forced = False                      # did a mismatch force the rebuild?

if _SURF_DIFF and not _want_surface:
    _w = max(len(k) for k, _o, _n in _SURF_DIFF)
    print(f'\n{"!" * 74}\nTHE SURFACE ON DISK WAS BUILT UNDER DIFFERENT '
          f'SETTINGS.\n  built {_SURF_META.get("built_at", "?")}\n')
    for _k, _old, _new in _SURF_DIFF:
        print(f'  {_k:<{_w}}  surface: {_old!r}\n  {"":<{_w}}  this run: {_new!r}')
    print('!' * 74)
    if not SURFACE_GUARD:
        print('SURFACE_GUARD = False - reusing it anyway. Every number below '
              'is\nmeasured on the OLD ground.')
    elif RUN_SURFACE == 'auto':
        print("RUN_SURFACE = 'auto' - rebuilding Pass A on the settings above.\n"
              'Set RUN_SURFACE = False to stop instead, or SURFACE_GUARD = '
              'False to\nreuse the old surface regardless.')
        _want_surface, _sig_forced = True, True
    else:
        raise RuntimeError(
            'RUN_SURFACE = False asks to reuse state_surface.gpkg, but it was '
            'built under the settings listed above. Set RUN_SURFACE = True to '
            'rebuild it, or SURFACE_GUARD = False if you know why the '
            'signature differs and want it reused anyway.')
elif _SURF_MISSING and not _want_surface:
    print(f'\nNOTE: this surface predates the signature check; '
          f'{len(_SURF_MISSING)} of its\n      inputs '
          f'({", ".join(_SURF_MISSING[:4])}'
          f'{", ..." if len(_SURF_MISSING) > 4 else ""}) are not recorded and '
          f'cannot be verified.\n      Rebuild it once (RUN_SURFACE = True) '
          f'and every later run is checked.')

if not _want_surface:
    print(f'\nRUN_SURFACE = {RUN_SURFACE!r} and {os.path.basename(SURFACE_GPKG)} '
          f'exists - loading it.')
    grid = gpd.read_file(SURFACE_GPKG, layer='grid').to_crs(WORKING_CRS)
    rip_crop = gpd.read_file(SURFACE_GPKG, layer='rip_crop').to_crs(WORKING_CRS)
    _meta = _SURF_META
    RIP_TOTAL = _meta['rip_total_ac']
    CORRIDOR_AC = _meta['corridor_ac']
    print(f'  grid {len(grid):,} cells | rip_crop {len(rip_crop):,} pieces | '
          f'{RIP_TOTAL:,.0f} ac')
else:
    rule(f'PASS A - WATER CORRIDOR AT {RIP_BUFFER_FT:.0f} FT, STATEWIDE')
    import time as _t
    _t0 = _t.time()

    _CSB_COLS = (['CDL', 'geometry'] if 'CDL' in csb.columns else ['geometry'])
    if 'CDL' not in csb.columns:
        print('  NOTE: no CDL column on csb - crop identity is unavailable')
    _csb_all = csb[_CSB_COLS].copy()
    _csb_all['_fid'] = np.arange(len(_csb_all))

    TILES = _tiles(ANALYSIS_GEOM, TILE_KM)
    print(f'  {len(TILES)} tiles of {TILE_KM:g} km, halo {HALO_M:,.0f} m')

    # Step 5a. Each tile's result is written as it is finished. A Colab session
    # that disconnects an hour in otherwise costs the whole of Pass A; with this
    # it costs one tile. Re-running picks up where it stopped, and the folder is
    # safe to delete once state_surface.gpkg exists.
    # RESUME_TILES and TILE_DIR come from fl_config.py; the checkpoints go
    # to the durable folder, so a reboot four hours in costs one tile.
    if not RESUME_TILES:
        for _f in glob.glob(os.path.join(TILE_DIR, 'tile_*')):
            os.remove(_f)

    # A checkpoint is a piece of ONE extent cut at ONE buffer. If the settings
    # moved since these were written - either because the signature check just
    # forced this rebuild, or because a run was interrupted and something was
    # changed before it resumed - resuming would stitch two different runs
    # together into a surface that is internally inconsistent and looks fine.
    # The tiles carry the same signature the surface does.
    _now_sig = surface_signature()
    _old_sig = None
    if os.path.exists(TILE_SIG):
        try:
            with open(TILE_SIG) as _f:
                _old_sig = json.load(_f)
        except (ValueError, OSError):
            _old_sig = None
    _tile_diff = (surface_diff(_old_sig, _now_sig)[0] if _old_sig is not None
                  else ([('(no signature recorded)', '?', '?')]
                        if glob.glob(os.path.join(TILE_DIR, 'tile_*.json'))
                        else []))
    if _tile_diff and SURFACE_GUARD:
        _n_rm = 0
        for _f in glob.glob(os.path.join(TILE_DIR, 'tile_*')):
            os.remove(_f)
            _n_rm += 1
        if _n_rm:
            print(f'  cleared {_n_rm} tile checkpoints - they were cut under '
                  f'different settings\n  ('
                  f'{", ".join(str(_k) for _k, _o, _n in _tile_diff[:4])})')
    elif _tile_diff:
        print('  SURFACE_GUARD = False - resuming tiles that were cut under '
              'different settings')
    os.makedirs(TILE_DIR, exist_ok=True)
    with open(TILE_SIG, 'w') as _f:
        json.dump(_now_sig, _f, indent=2)

    _done = len(glob.glob(os.path.join(TILE_DIR, 'tile_*.json')))
    if _done:
        print(f'  resuming: {_done} of {len(TILES)} tiles already on disk in '
              f'{os.path.basename(TILE_DIR)}/\n'
              f'  (set RESUME_TILES = False in B3 to start over)')

    _pieces, _corr_ac, _src_ac = [], 0.0, {}
    for _k, _tile in enumerate(TILES, 1):
        _tj = os.path.join(TILE_DIR, f'tile_{_k:04d}.json')
        _tg = os.path.join(TILE_DIR, f'tile_{_k:04d}.gpkg')
        if os.path.exists(_tj):
            with open(_tj) as _f:
                _meta_t = json.load(_f)
            _corr_ac += _meta_t['corridor_ac']
            for _n, _v in _meta_t['src_ac'].items():
                _src_ac[_n] = _src_ac.get(_n, 0.0) + _v
            if _meta_t['n_pieces']:
                _pieces.append(gpd.read_file(_tg).to_crs(WORKING_CRS))
            continue

        _cell = _tile.intersection(ANALYSIS_GEOM)
        _here_corr, _here_src, _ov = 0.0, {}, None
        if not _cell.is_empty:
            _hb = _tile.buffer(HALO_M).bounds
            _parts, _corr = corridor_for(_hb)
            _corr_here = (None if _corr is None
                          else _corr.intersection(_cell))
            if _corr_here is not None and not _corr_here.is_empty:
                _here_corr = _corr_here.area / AC
                for _n, _g in _parts.items():
                    _here_src[_n] = _g.intersection(_cell).area / AC
                _cs = _csb_all.cx[_hb[0]:_hb[2], _hb[1]:_hb[3]]
                if len(_cs):
                    _ov = fast_overlay(_cs, _corr_here, f'tile {_k}')
                    _ov = clip_to(_ov, _cell, f'tile {_k}') if len(_ov) else _ov
                    if len(_ov):
                        # which source each piece sits on, decided in this tile
                        _rp = np.asarray(
                            _ov.geometry.representative_point().values,
                            dtype=object)
                        for _n, _g in _parts.items():
                            shapely.prepare(_g)
                            _ov[f'on_{_n}'] = shapely.within(_rp, _g)
                    else:
                        _ov = None

        _corr_ac += _here_corr
        for _n, _v in _here_src.items():
            _src_ac[_n] = _src_ac.get(_n, 0.0) + _v
        if _ov is not None and len(_ov):
            gpkg_safe(_ov).to_file(_tg, layer='rip_crop', driver='GPKG')
            _pieces.append(_ov)
        with open(_tj, 'w') as _f:
            json.dump({'corridor_ac': _here_corr, 'src_ac': _here_src,
                       'n_pieces': int(0 if _ov is None else len(_ov))}, _f)

        _el = _t.time() - _t0
        print(f'    tile {_k:>4}/{len(TILES)}  '
              f'{sum(len(p) for p in _pieces):>9,} pieces  '
              f'{_el:>6.0f}s elapsed  '
              f'~{_el / max(_k - _done, 1) * (len(TILES) - _k) / 60:>5.1f} min '
              f'left', flush=True)

    assert _pieces, 'no cropland on water anywhere in the analysed ground'
    rip_crop = pd.concat(_pieces, ignore_index=True)
    rip_crop = gpd.GeoDataFrame(rip_crop, geometry='geometry', crs=WORKING_CRS)
    for _n in USE_SOURCES:
        if f'on_{_n}' not in rip_crop.columns:
            rip_crop[f'on_{_n}'] = False
        rip_crop[f'on_{_n}'] = rip_crop[f'on_{_n}'].fillna(False).astype(bool)
    rip_crop['ac'] = rip_crop.geometry.area / AC
    _flags = rip_crop[[f'on_{n}' for n in USE_SOURCES]]
    rip_crop['water_source'] = [
        '+'.join([n for n, v in zip(USE_SOURCES, row) if v]) or 'none'
        for row in _flags.to_numpy()]
    RIP_TOTAL = rip_crop['ac'].sum()
    CORRIDOR_AC = _corr_ac

    print(f'\n  corridor      {CORRIDOR_AC:>12,.0f} ac')
    for _n in USE_SOURCES:
        if _n in _src_ac:
            print(f'    {_n:<11} {_src_ac[_n]:>12,.0f} ac (before the union)')
    print(f'  riparian crop {RIP_TOTAL:>12,.0f} ac across '
          f'{rip_crop["_fid"].nunique():,} fields, '
          f'{len(rip_crop):,} pieces')
    print(f'  {_t.time() - _t0:.0f}s')
    print(rip_crop.groupby('water_source')['ac'].agg(['size', 'sum'])
          .rename(columns={'size': 'polygons', 'sum': 'acres'}).round(0)
          .sort_values('acres', ascending=False).head(10).to_string())


# %% B4 - the 1-mile grid
# Step 6. One lattice over the whole analysed ground, built once. `land_ac` is
# the part of each cell inside that ground - edge cells are partial, and using
# the full cell area would understate their density.
if _want_surface:
    rule('PASS A - THE 1-MILE GRID')
    _x0, _y0, _x1, _y1 = ANALYSIS_GEOM.bounds
    _nx = int(np.ceil((_x1 - _x0) / CELL_M))
    _ny = int(np.ceil((_y1 - _y0) / CELL_M))
    print(f'  lattice {_nx:,} x {_ny:,} = {_nx * _ny:,} candidate cells')

    _boxes = np.array([shp_box(_x0 + i * CELL_M, _y0 + j * CELL_M,
                               _x0 + (i + 1) * CELL_M, _y0 + (j + 1) * CELL_M)
                       for i in range(_nx) for j in range(_ny)], dtype=object)
    shapely.prepare(ANALYSIS_GEOM)
    _keep = shapely.intersects(ANALYSIS_GEOM, _boxes)
    _boxes = _boxes[_keep]
    grid = gpd.GeoDataFrame(geometry=_boxes, crs=WORKING_CRS)
    print(f'  {len(grid):,} cells intersect the analysed ground')

    # interior cells are already their own answer; only edge cells are cut
    _full = shapely.contains_properly(ANALYSIS_GEOM, _boxes)
    _land = np.full(len(_boxes), CELL_M ** 2 / AC)
    _land[~_full] = shapely.area(
        shapely.intersection(_boxes[~_full], ANALYSIS_GEOM)) / AC
    grid['land_ac'] = _land
    print(f'  {int(_full.sum()):,} whole cells, {int((~_full).sum()):,} cut by '
          f'the boundary')
    grid = grid[grid['land_ac'] > 0].reset_index(drop=True)

    # Step 7. Allocate the riparian cropland to cells by exact intersection.
    # The pair list comes from an STRtree rather than gpd.overlay: identical
    # arithmetic, but it does not build a GeoDataFrame of every intersection
    # piece, which at this scale is the difference between minutes and hours.
    # B5b checks it against gpd.overlay on a sample.
    _gg = np.asarray(grid.geometry.values, dtype=object)
    _rg = np.asarray(rip_crop.geometry.values, dtype=object)
    _tree = shapely.STRtree(_gg)
    _ri, _ci = _tree.query(_rg, predicate='intersects')
    _inter_ac = shapely.area(shapely.intersection(_rg[_ri], _gg[_ci])) / AC
    grid['rip_ac'] = (pd.Series(_inter_ac).groupby(_ci).sum()
                      .reindex(range(len(grid))).fillna(0.0).to_numpy())
    print(f'  {len(_ri):,} cropland-piece x cell intersections -> '
          f'{grid["rip_ac"].sum():,.0f} ac allocated '
          f'(of {RIP_TOTAL:,.0f} ac total)')


# %% B5 - the 3-mile neighbourhood density
# Step 8. THE LINE THIS WHOLE SCRIPT EXISTS FOR. The neighbourhood is taken
# over the statewide lattice, so a cell near a county line sees its real
# 3-mile surroundings instead of the truncated piece that happens to fall
# inside its own county. Density is the ratio of sums - neighbourhood acres of
# riparian cropland over neighbourhood acres of land - never the mean of the
# cells' own ratios.
if _want_surface:
    _ctr = np.c_[grid.geometry.centroid.x, grid.geometry.centroid.y]
    _nb = cKDTree(_ctr).query_ball_point(_ctr, r=BAND_M)
    _lens = np.fromiter(map(len, _nb), int, len(_nb))
    _flat = np.concatenate(_nb)
    _starts = np.r_[0, np.cumsum(_lens)[:-1]]

    _rip = grid['rip_ac'].to_numpy(float)
    _lnd = grid['land_ac'].to_numpy(float)
    grid['n_nbhd']       = _lens
    grid['nbhd_rip_ac']  = band_sums(_rip, _flat, _starts)
    grid['nbhd_land_ac'] = band_sums(_lnd, _flat, _starts)
    grid['rip_per_sqmi'] = (grid['nbhd_rip_ac'] / grid['nbhd_land_ac'] * 640).round(2)
    grid['own_per_sqmi'] = (grid['rip_ac'] / grid['land_ac'] * 640).round(2)
    print(f'  neighbourhood: {_lens.mean():.1f} cells on average, '
          f'{_lens.min()} to {_lens.max()}')

    # the reduceat sums must equal the plain per-cell sums
    _rs = np.random.default_rng(0).choice(len(grid), size=min(200, len(grid)),
                                          replace=False)
    _chk = max(abs(_rip[_nb[i]].sum() - grid['nbhd_rip_ac'].iloc[i])
               for i in _rs)
    assert _chk < 1e-6, f'band_sums disagrees with the per-cell sum by {_chk}'
    print(f'  CHECK band sums vs per-cell sums on {len(_rs)} cells: '
          f'max difference {_chk:.2e} ac')

    # how much of each cell's band is real ground rather than off the edge
    _full_band = np.pi * (BAND_MI ** 2) * 640
    grid['band_cover'] = (grid['nbhd_land_ac'] / _full_band).clip(0, 1).round(3)

    grid.index.name = 'cell'
    gpkg_safe(grid.reset_index()).to_file(SURFACE_GPKG, layer='grid',
                                          driver='GPKG')
    gpkg_safe(rip_crop).to_file(SURFACE_GPKG, layer='rip_crop', driver='GPKG')
    with open(SURFACE_JSON, 'w') as _f:
        json.dump({'built_at': datetime.datetime.now().isoformat(timespec='seconds'),
                   'cache': CACHE, 'source_scope': SOURCE_SCOPE,
                   'county_extent_mode': COUNTY_EXTENT_MODE,
                   'analysis_ac': ANALYSIS_AC, 'corridor_ac': CORRIDOR_AC,
                   'rip_total_ac': float(RIP_TOTAL), 'n_cells': int(len(grid)),
                   'cell_mi': CELL_MI, 'band_mi': BAND_MI,
                   'buffer_ft': RIP_BUFFER_FT,
                   # what a later run checks itself against before reusing this
                   'signature': surface_signature()}, _f, indent=2)
    print(f'\n  surface -> {SURFACE_GPKG}')

    # The per-tile checkpoints have done their job now that the surface
    # exists. KEEP_TILES = True (the default) is insurance against a crash in
    # a LATER stage: the surface is there, but if you had to rebuild it you
    # would otherwise re-tile from zero.
    if not KEEP_TILES:
        _n_rm = 0
        for _f in glob.glob(os.path.join(TILE_DIR, 'tile_*')):
            os.remove(_f)
            _n_rm += 1
        if _n_rm:
            print(f'  removed {_n_rm} tile checkpoints (KEEP_TILES = True '
                  f'keeps them)')

if 'cell' in grid.columns:
    grid = grid.set_index('cell')
grid.index.name = 'cell'


# %% B5b - does the fast allocation agree with gpd.overlay?
# Step 9. The STRtree allocation in B4 replaces gpd.overlay for speed. That is
# only legitimate if it gives the same acres, so it is checked against
# gpd.overlay on one county's worth of cells rather than asserted.
_sample_fips = counties.loc[counties['covered'], 'FIPS']
if len(_sample_fips) and _want_surface:
    _sg = counties.loc[counties['FIPS'] == _sample_fips.iloc[0],
                       'extent'].iloc[0]
    _sub = grid[grid.geometry.intersects(_sg)]
    if len(_sub):
        # every piece that touches THESE CELLS, not every piece in the county:
        # a cell on the county line also holds cropland from the other side,
        # and leaving that out would make the control wrong rather than the
        # thing it is controlling.
        _sb = _sub.total_bounds
        _sr = rip_crop.cx[_sb[0]:_sb[2], _sb[1]:_sb[3]]
        _ov = gpd.overlay(_sr[['geometry']],
                          _sub[['geometry']].reset_index()[['cell', 'geometry']],
                          how='intersection', keep_geom_type=True)
        _ov['ac'] = _ov.geometry.area / AC
        _ref = _ov.groupby('cell')['ac'].sum().reindex(_sub.index).fillna(0.0)
        _d = float((_ref - _sub['rip_ac']).abs().max())
        print(f'\nCHECK STRtree allocation vs gpd.overlay on county '
              f'{_sample_fips.iloc[0]} ({len(_sub):,} cells):')
        print(f'  largest per-cell difference {_d:.6f} ac '
              f'({"PASS" if _d < 1e-6 else "FAIL"})')
        assert _d < 1e-6, 'the fast allocation does not match gpd.overlay'
        with open(os.path.join(S_VAL, 'equivalence_checks.txt'), 'w') as _f:
            _f.write(f'STRtree cell allocation vs gpd.overlay, county '
                     f'{_sample_fips.iloc[0]}, {len(_sub)} cells: '
                     f'max per-cell difference {_d:.3e} ac\n')


# %% B6 - the reporting mask: report inside the priority area, keep the surface
# WHAT WAS MEASURED IS NOT WHERE IT IS REPORTED, and conflating those two is
# what makes a priority-area run expensive.
#
# COUNTY_EXTENT_MODE decides what Pass A MEASURES. Setting it to 'pa_only'
# moves analysis_ac, n_counties and the county hash, so the signature no
# longer matches and Pass A rebuilds - hours - and the surface that comes back
# is WORSE: the 1-mile lattice is then laid only over the priority area, so
# cKDTree finds fewer neighbours near its edge and every edge cell's 3-mile
# band is measured over a truncated piece of ground. That is exactly the
# defect C7 quantifies for county lines, reintroduced against a boundary that
# is programmatic rather than ecological.
#
# REPORT_WITHIN_PA decides where the results are REPORTED. It runs HERE, after
# the surface has loaded and after every surface_signature() call, so it is
# structurally invisible to the signature and costs seconds.
#
# WHAT BECOMES PRIORITY-AREA GROUND: land_ac, rip_ac, own_per_sqmi, RIP_TOTAL,
# the county extents and ANALYSIS_GEOM - and therefore all three references,
# the floor sweep, the hot rule, the zones, the manifest and every figure,
# none of which needed an edit because every one of them reads those.
#
# WHAT DOES NOT MOVE: rip_per_sqmi and the nbhd_* columns. The 3-mile
# neighbourhood density stays exactly as Pass A measured it statewide, so a
# cell just inside the boundary still sees its real surroundings. That is the
# entire point of doing this here instead of rebuilding, and `band_in_report`
# records how much of each cell's band lies inside the mask, so the cells that
# lean on outside ground are a column rather than a footnote.
#
# NOT COMPARABLE WITH A STATEWIDE RUN. The averages and the derived floor are
# computed over different ground - which is what was asked for - so
# `report_mask` is stamped on every zone row and every manifest row.
REPORT_MASK = None
if globals().get('REPORT_WITHIN_PA'):
    if _PA_GEOM is None:
        raise ValueError(
            'REPORT_WITHIN_PA = True but there is no priority area to mask '
            'with:\n  pa_geo is not in the cache. Check the "priority '
            'area:" line in PASS A -\n  INPUTS above. The PA is enabled by '
            'being there, so put a file matching\n  the \'pa\' pattern in '
            'raw_data/ and re-run fl_01_preprocess.py; that costs\n  one '
            'small layer rebuild and leaves every resource layer alone.')
    if COUNTY_EXTENT_MODE != 'whole_county':
        print(f'\nNOTE: COUNTY_EXTENT_MODE = {COUNTY_EXTENT_MODE!r} has '
              f'already cut the extents to the\n      priority area, so the '
              f'mask has little left to do - and that run had to\n      '
              f'rebuild Pass A to get there. The cheap combination is '
              f"'whole_county'\n      plus REPORT_WITHIN_PA = True: one "
              f'statewide surface, reported inside the PA.')
    REPORT_MASK = 'pa'

if REPORT_MASK == 'pa':
    rule('REPORTING MASK - THE PRIORITY AREA')
    _m_cells0, _m_ac0, _m_rip0, _m_cty0 = (len(grid), ANALYSIS_AC, RIP_TOTAL,
                                           len(counties))

    # Step 1. The counties FIRST, so the reported ground is the union of what
    # the counties actually hold - the same order B2 builds it in. A county
    # under MIN_PA_SHARE is dropped rather than clipped to a sliver: that is
    # what MIN_PA_SHARE is for, and a county average over one square mile is
    # not an average.
    counties = counties[counties['in_pa']].reset_index(drop=True)
    if not len(counties):
        raise ValueError(
            f'no {COUNTY_WORD.lower()} overlaps the priority area by at least '
            f'{MIN_PA_SHARE:.0%} (MIN_PA_SHARE),\n  so there is nothing to '
            f'report on. Check that the PA layer is really this state\'s.')
    counties['extent'] = shapely.intersection(
        np.asarray(counties['extent'].values, dtype=object), _PA_GEOM)
    counties['extent_ac'] = area_ac(counties['extent'].values)
    counties = counties[counties['extent_ac'] > 0].reset_index(drop=True)
    counties['extent_basis'] = 'county n priority area (reporting mask)'
    _m_eg = np.asarray(counties['extent'].values, dtype=object)
    counties['source_ac'] = shapely.area(
        shapely.intersection(_m_eg, SOURCE_GEOM)) / AC
    counties['source_cover'] = (counties['source_ac'] /
                                counties['extent_ac']).round(4)
    counties['covered'] = counties['source_cover'] >= MIN_SOURCE_COVER

    # Step 2. The reported ground, by the same expression B2 uses for the
    # analysed ground - union of the extents, cut to what the cache covers.
    ANALYSIS_GEOM = union_all(_m_eg).intersection(SOURCE_GEOM)
    ANALYSIS_AC = ANALYSIS_GEOM.area / AC
    counties['analysed_geom'] = shapely.intersection(_m_eg, ANALYSIS_GEOM)
    shapely.prepare(ANALYSIS_GEOM)

    # Step 3. How much of each cell's BAND falls inside the mask - computed
    # BEFORE any cell is dropped, because a cell outside the mask still
    # contributes ground to a neighbour's band. That is precisely why the
    # density is not being recomputed, and this is the column that says so.
    _m_gg = np.asarray(grid.geometry.values, dtype=object)
    if 'land_ac_measured' not in grid.columns:
        grid['land_ac_measured'] = grid['land_ac']
    _m_full = shapely.contains_properly(ANALYSIS_GEOM, _m_gg)
    _m_land = np.full(len(grid), CELL_M ** 2 / AC)
    _m_land[~_m_full] = shapely.area(
        shapely.intersection(_m_gg[~_m_full], ANALYSIS_GEOM)) / AC
    _m_ctr = np.c_[grid.geometry.centroid.x, grid.geometry.centroid.y]
    _m_nb = cKDTree(_m_ctr).query_ball_point(_m_ctr, r=BAND_M)
    _m_lens = np.fromiter(map(len, _m_nb), int, len(_m_nb))
    _m_flat = np.concatenate(_m_nb)
    _m_starts = np.r_[0, np.cumsum(_m_lens)[:-1]]
    _m_all = band_sums(grid['land_ac_measured'].to_numpy(float), _m_flat,
                       _m_starts)
    grid['band_in_report'] = np.round(
        np.where(_m_all > 0, band_sums(_m_land, _m_flat, _m_starts) / _m_all,
                 0.0), 3)

    # Step 4. The cells. `land_ac` becomes the part inside the mask; what Pass
    # A measured is kept beside it as `land_ac_measured`, so the two can be
    # compared rather than confused.
    grid['land_ac'] = np.round(_m_land, 4)
    grid = grid[grid['land_ac'] > 0].copy()

    # Step 5. The cropland, and every per-cell number that is a function of
    # it. `rip_per_sqmi` and the nbhd_* columns are deliberately NOT here.
    rip_crop = clip_to(rip_crop, ANALYSIS_GEOM, 'priority area')
    rip_crop['ac'] = rip_crop.geometry.area / AC
    RIP_TOTAL = float(rip_crop['ac'].sum())
    _m_gg2 = np.asarray(grid.geometry.values, dtype=object)
    _m_rg2 = np.asarray(rip_crop.geometry.values, dtype=object)
    _m_ri, _m_ci = shapely.STRtree(_m_gg2).query(_m_rg2, predicate='intersects')
    grid['rip_ac'] = (
        pd.Series(shapely.area(shapely.intersection(_m_rg2[_m_ri],
                                                    _m_gg2[_m_ci])) / AC)
        .groupby(_m_ci).sum()
        .reindex(range(len(grid))).fillna(0.0).to_numpy())
    grid['own_per_sqmi'] = (grid['rip_ac'] /
                            grid['land_ac'].replace(0, np.nan) * 640).round(2)

    print(f'  priority area : {_PA_GEOM.area / AC / 640:,.0f} sq mi')
    print(f'  reported on   : {ANALYSIS_AC / 640:,.0f} sq mi, '
          f'{100 * ANALYSIS_AC / max(_m_ac0, 1e-9):.0f}% of the '
          f'{_m_ac0 / 640:,.0f} sq mi Pass A measured')
    print(f'  {COUNTY_WORDS:<14}: {len(counties)} of {_m_cty0}')
    print(f'  cells         : {len(grid):,} of {_m_cells0:,}')
    print(f'  riparian crop : {RIP_TOTAL:,.0f} ac of {_m_rip0:,.0f} '
          f'({100 * RIP_TOTAL / max(_m_rip0, 1e-9):.0f}%)')
    _m_lean = int((grid['band_in_report'] < 0.95).sum())
    print(f'\n  THE DENSITY SURFACE IS UNCHANGED - only what is reported on '
          f'it has moved.\n  {_m_lean:,} cells '
          f'({_m_lean / max(len(grid), 1):.0%}) draw part of their 3-mile '
          f'band from ground outside\n  the mask, which is deliberate: a cell '
          f'just inside the boundary still sees its\n  real surroundings '
          f'instead of a band cut at a programmatic line. '
          f'band_in_report\n  carries the share for every cell.')
    print(f'\n  The averages and the floor below are derived over '
          f'priority-area ground and are\n  NOT comparable with a statewide '
          f"run. report_mask = 'pa' is on every zone row\n  and every "
          f'manifest row so the two cannot be read as one table.')


# %% [markdown]
# ## The threshold
#
# One absolute number per cell — `rip_per_sqmi`, riparian-crop acres per square
# mile over its 3-mile neighbourhood — and three references to measure it
# against. All three are carried as columns, so which one decides "hot" is a
# reporting choice made here rather than a parameter baked into the surface.
#
# | column | reference | the question it answers |
# |---|---|---|
# | `ratio_county` | the county's own average | where in this county should this biologist look? |
# | `ratio_state`  | the state's average | which counties deserve a biologist at all? |
# | `ratio_region` | the HUC8's average | where, allowing for the fact that the Purchase and the Bluegrass are different places? |
# | `p_state` | — | the cell's statewide percentile |
#
# A ratio on its own always finds a top tail, so a county with almost no
# riparian cropland manufactures "hotspots" out of three fields on a creek.
# The absolute FLOOR is what stops that, and it is derived from the statewide
# distribution and swept rather than picked.

# %% C1 - which county and which region each cell belongs to
rule('ATTRIBUTION - COUNTY AND REGION')

_ctr_pts = gpd.GeoSeries(grid.geometry.centroid.values, crs=WORKING_CRS)
_cty = gpd.GeoDataFrame(counties[['FIPS', 'NAME']].copy(),
                        geometry=counties['extent'].values, crs=WORKING_CRS)
_j = gpd.sjoin(gpd.GeoDataFrame(geometry=_ctr_pts.values, crs=WORKING_CRS),
               _cty, predicate='within', how='left')
_j = _j[~_j.index.duplicated(keep='first')]
grid['FIPS'] = _j['FIPS'].to_numpy()
grid['county'] = _j['NAME'].to_numpy()
_orphan = grid['FIPS'].isna()
if _orphan.any():
    # a centroid on a boundary can miss every polygon; give it the nearest
    _nn = gpd.sjoin_nearest(
        gpd.GeoDataFrame(geometry=_ctr_pts[_orphan.to_numpy()].values,
                         crs=WORKING_CRS), _cty, how='left')
    _nn = _nn[~_nn.index.duplicated(keep='first')]
    grid.loc[_orphan, 'FIPS'] = _nn['FIPS'].to_numpy()
    grid.loc[_orphan, 'county'] = _nn['NAME'].to_numpy()
    print(f'  {int(_orphan.sum())} cell(s) on a boundary assigned to the '
          f'nearest county')
print(f'  {grid["FIPS"].nunique()} counties hold at least one cell')

# Step 1. Regions are HUC8s, the first eight digits of the cached HUC12 codes.
# No extra download, and it is the unit the water program already thinks in.
HAVE_REGIONS = huc12 is not None and len(huc12)
if HAVE_REGIONS:
    _hc = next((c for c in huc12.columns if c.lower() == 'huc12'), None)
    HAVE_REGIONS = _hc is not None
if HAVE_REGIONS:
    _h = huc12[[_hc, 'geometry']].copy()
    _h['HUC8'] = _h[_hc].astype(str).str.replace(r'\D', '', regex=True).str[:8]
    _h = _h[_h['HUC8'].str.len() == 8]
    regions = _h.dissolve(by='HUC8', as_index=False)[['HUC8', 'geometry']]
    _rg = np.asarray(regions.geometry.values, dtype=object)
    regions['region_ac_full'] = shapely.area(_rg) / AC
    regions['region_ac'] = shapely.area(
        shapely.intersection(_rg, ANALYSIS_GEOM)) / AC
    regions['region_cover'] = (regions['region_ac'] /
                               regions['region_ac_full']).round(3)
    _jr = gpd.sjoin(gpd.GeoDataFrame(geometry=_ctr_pts.values, crs=WORKING_CRS),
                    regions[['HUC8', 'geometry']], predicate='within', how='left')
    _jr = _jr[~_jr.index.duplicated(keep='first')]
    grid['HUC8'] = _jr['HUC8'].to_numpy()
    print(f'  {len(regions)} HUC8 regions, '
          f'{grid["HUC8"].notna().sum():,} of {len(grid):,} cells inside one')
else:
    regions = gpd.GeoDataFrame({'HUC8': [], 'region_ac': [], 'region_cover': []},
                               geometry=[], crs=WORKING_CRS)
    grid['HUC8'] = np.nan
    print('  no HUC12 layer in the cache - the region basis falls back to '
          'the state')


# %% C2 - the three references
# Step 2. Each reference is a ratio of sums: riparian-crop acres over land
# acres, times 640. Never a mean of the cells' own ratios, which would weight a
# tiny edge cell the same as a whole one.
rule('REFERENCES - COUNTY, STATE, REGION')

# Step 2a. The references are AREAL ratios, taken against the ground itself
# rather than against the cells that happen to cover it. That is the
# notebook's own definition - RIP_TOTAL / EXTENT_AC * 640 - and keeping it
# matters for two reasons: the number stays comparable with a single-county
# run, and the threshold stops depending on where the lattice lines fall. A
# cell straddling a county line is owned by one county but holds cropland from
# both, and on a boundary stream that was enough to move a county average by
# 14%.
def areal_avg(geoms):
    """Riparian crop acres per square mile inside each geometry, exactly."""
    g = np.asarray(geoms, dtype=object)
    rg = np.asarray(rip_crop.geometry.values, dtype=object)
    gi, ci = shapely.STRtree(g).query(rg, predicate='intersects')
    ia = shapely.area(shapely.intersection(rg[gi], g[ci])) / AC
    rip = (pd.Series(ia).groupby(ci).sum()
           .reindex(range(len(g))).fillna(0.0).to_numpy())
    land = shapely.area(g) / AC
    return rip, land, np.where(land > 0, rip / np.where(land > 0, land, 1) * 640,
                               np.nan)


AVG_STATE = RIP_TOTAL / ANALYSIS_AC * 640
print(f'  state average : {AVG_STATE:.2f} riparian crop ac per sq mi '
      f'over {ANALYSIS_AC / 640:,.0f} sq mi')

counties['rip_ac_exact'], counties['analysed_ac'], _ca = areal_avg(
    counties['analysed_geom'].values)
counties['avg_county'] = _ca
_cty_stat = grid.groupby('FIPS').agg(rip_ac=('rip_ac', 'sum'),
                                     land_ac=('land_ac', 'sum'),
                                     n_cells=('rip_ac', 'size'))
_cty_stat['avg_county'] = counties.set_index('FIPS')['avg_county']
grid['avg_county'] = grid['FIPS'].map(_cty_stat['avg_county']).astype(float)

if HAVE_REGIONS:
    regions['analysed_geom'] = shapely.intersection(
        np.asarray(regions.geometry.values, dtype=object), ANALYSIS_GEOM)
    regions['rip_ac_exact'], _rl, regions['avg_region'] = areal_avg(
        regions['analysed_geom'].values)
    _reg_stat = grid.dropna(subset=['HUC8']).groupby('HUC8').agg(
        rip_ac=('rip_ac', 'sum'), land_ac=('land_ac', 'sum'),
        n_cells=('rip_ac', 'size'))
    _reg_stat['avg_region'] = regions.set_index('HUC8')['avg_region']
    _reg_stat = _reg_stat.join(regions.set_index('HUC8')['region_cover'])
    # a region that is mostly outside the analysed ground, or too small to be
    # its own sample, is not a reference - those cells use the state
    _reg_stat['usable'] = ((_reg_stat['n_cells'] >= REGION_MIN_CELLS) &
                           (_reg_stat['region_cover'].fillna(1)
                            >= REGION_MIN_COVER))
    grid['avg_region'] = grid['HUC8'].map(
        _reg_stat.loc[_reg_stat['usable'], 'avg_region']).astype(float)
    grid['region_ref'] = np.where(grid['avg_region'].notna(), 'HUC8', 'state')
    grid['avg_region'] = grid['avg_region'].fillna(AVG_STATE)
    print(f'  regions       : {int(_reg_stat["usable"].sum())} of '
          f'{len(_reg_stat)} HUC8s are usable as their own reference '
          f'(>= {REGION_MIN_CELLS} cells, >= {REGION_MIN_COVER:.0%} covered);')
    print(f'                  the rest fall back to the state average')
else:
    grid['avg_region'] = AVG_STATE
    grid['region_ref'] = 'state'

grid['avg_state'] = AVG_STATE
for _b in ('county', 'state', 'region'):
    grid[f'ratio_{_b}'] = (grid['rip_per_sqmi'] /
                           grid[f'avg_{_b}'].replace(0, np.nan)).round(3)
grid['p_state'] = grid['rip_per_sqmi'].rank(pct=True).round(4)

print(f'\n  how much the reference moves between counties:')
_sp = (_cty_stat['avg_county'].describe()[['min', '25%', '50%', '75%', 'max']]
       .round(2))
print(_sp.to_frame('ac per sq mi').to_string())
print(f'  A cell at {_sp["50%"] * HOT_MULTIPLE:.1f} ac/sq mi is HOT in the '
      f'median county and\n  not hot in the top quartile. That spread is the '
      f'reason for the floor.')


# %% C3 - the floor sweep
# Step 3. Sweep the absolute cut across the statewide distribution and report,
# at each value: how many counties keep a zone, how much ground is hot, how
# much of the resource that ground holds, and the concentration - the share of
# the resource divided by the share of the area. Concentration is the number
# worth preserving; the acre figure is just where it happens to land.
rule('FLOOR SWEEP - DERIVING THE ABSOLUTE CUT')

_ref = grid[f'avg_{BASIS}'].to_numpy(float)
_dens = grid['rip_per_sqmi'].to_numpy(float)
_own = grid['own_per_sqmi'].to_numpy(float)
_land = grid['land_ac'].to_numpy(float)
_ripv = grid['rip_ac'].to_numpy(float)
_ratio_ok = (_dens >= HOT_MULTIPLE * _ref) & (_own >= _ref)

_FIPS_ARR = grid['FIPS'].to_numpy()
_LAND_TOT = _land.sum()
_RIP_TOT = max(_ripv.sum(), 1e-9)

# The cells are a lattice, so zones are connected components of a boolean
# image and cost nothing to count. Doing it this way means the sweep can
# report the number of ZONES and the number of counties that keep a real zone
# - MIN_ZONE_CELLS applied - rather than counting counties that merely hold a
# hot cell, which is a weaker claim wearing the same name.
# 4-connectivity, because two cells touching only at a CORNER union into a
# MultiPolygon in C4 and are therefore two zones there too.
_gb = np.array([g.bounds for g in grid.geometry])
_CI = np.round((_gb[:, 0] - _gb[:, 0].min()) / CELL_M).astype(int)
_CJ = np.round((_gb[:, 1] - _gb[:, 1].min()) / CELL_M).astype(int)
_CSHAPE = (_CJ.max() + 1, _CI.max() + 1)


def zone_stats(h):
    """(number of zones, mask of cells inside one) for a hot-cell mask."""
    img = np.zeros(_CSHAPE, bool)
    img[_CJ[h], _CI[h]] = True
    lab, n = ndimage.label(img)
    if not n:
        return 0, np.zeros(len(h), bool)
    sizes = np.bincount(lab.ravel())
    keep = np.zeros(sizes.size, bool)
    keep[MIN_ZONE_CELLS and np.flatnonzero(sizes >= MIN_ZONE_CELLS)] = True
    keep[0] = False                                   # label 0 is background
    return int(keep.sum()), keep[lab[_CJ, _CI]] & h


def sweep_at(floor_value, pctl=np.nan):
    """Apply the whole rule at one floor and report what it keeps.

    The same function feeds the REPORTED table and every knee rule, so they
    can never disagree about what a floor does - which they would the moment
    one of them was reimplemented."""
    h = _ratio_ok & (_dens >= floor_value)
    ha, hr = _land[h].sum(), _ripv[h].sum()
    pa = 100 * ha / _LAND_TOT
    pr = 100 * hr / _RIP_TOT
    nz, inz = zone_stats(h)
    return {'pctl': pctl, 'floor_ac_per_sqmi': round(float(floor_value), 2),
            'hot_cells': int(h.sum()), 'hot_sqmi': round(ha / 640, 1),
            'hot_land_ac': round(float(ha), 1), 'hot_rip_ac': round(float(hr), 1),
            'pct_area_hot': round(pa, 2), 'pct_resource_hot': round(pr, 1),
            'concentration': round(pr / pa, 2) if pa > 0 else np.nan,
            'n_zones': nz,
            'counties_with_hot_cells': int(pd.Series(_FIPS_ARR[h]).nunique()),
            'counties_with_zones': int(pd.Series(_FIPS_ARR[inz]).nunique())}


def _floor_at(p):
    return 0.0 if p <= 0 else float(np.percentile(_dens, p))


floor_sweep = pd.DataFrame([sweep_at(_floor_at(_p), _p)
                            for _p in FLOOR_SWEEP_PCTLS])

# %% [markdown]
# ### Step 3a — where the floor should sit, found rather than asserted
#
# Raising the floor always buys **concentration** (the designated ground holds
# a larger share of the resource than it does of the area) and always costs
# **capture** (some riparian cropland ends up outside the map). Every rule
# below is a different answer to "when does that trade stop paying", and the
# run computes **all of them, every time**, so the choice is made against
# numbers rather than in the abstract. `FLOOR_KNEE_RULE` picks which one
# `FLOOR_MODE = 'knee'` adopts.
#
# | rule | cost of raising the floor | parameters | what it protects |
# |---|---|---|---|
# | `'resource'` | riparian cropland left outside the map | `FLOOR_KNEE_WEIGHT` | the resource |
# | `'counties'` | counties that keep no zone | `FLOOR_KNEE_WEIGHT` | statewide delivery |
# | `'curve'` | — (geometry of the capture curve) | none | nothing; it is the shape |
# | `'marginal'` | the density of the ground being thrown away | `FLOOR_KNEE_MARGINAL_X` | the exchange rate |
#
# **Zone count is deliberately not a cost.** Raising the floor can knock out a
# single bridging cell and SPLIT one zone into two, so `n_zones` moves both
# ways as the floor rises. An argmax over a non-monotone curve finds the place
# a bridge happened to break, not a knee. It is reported in the table because
# it is worth seeing; it is not optimised against.

# %% C3b - the knee rules
_kp = np.round(np.arange(0.0, FLOOR_KNEE_MAX_PCTL + 1e-9, FLOOR_KNEE_GRID), 4)
knee_curve = pd.DataFrame([sweep_at(_floor_at(_p), _p) for _p in _kp])

_C = knee_curve['concentration'].to_numpy(float)     # gain, every rule
_R = knee_curve['pct_resource_hot'].to_numpy(float)  # capture
_A = knee_curve['pct_area_hot'].to_numpy(float)      # area designated
_NC = knee_curve['counties_with_zones'].to_numpy(float)
_HR = knee_curve['hot_rip_ac'].to_numpy(float)
_HA = knee_curve['hot_land_ac'].to_numpy(float)
_OKC = np.isfinite(_C) & (knee_curve['hot_cells'].to_numpy() > 0)
_SPAN = (float(np.nanmax(_C[_OKC])) - float(np.nanmin(_C[_OKC]))
         if _OKC.any() else 0.0)


def _norm(v, lo=None, hi=None):
    lo = float(np.nanmin(v[_OKC])) if lo is None else lo
    hi = float(np.nanmax(v[_OKC])) if hi is None else hi
    return (v - lo) / max(hi - lo, 1e-9)


def _pick(net, note):
    """argmax takes the FIRST maximum, i.e. the smallest floor that ties -
    the floor should be the least cut that earns its keep."""
    net = np.where(_OKC, net, -np.inf)
    i = int(np.argmax(net))
    return {'pctl': float(_kp[i]), 'floor_ac_per_sqmi': float(_floor_at(_kp[i])),
            'note': note, 'net': net}


def knee_by(rule):
    """Each rule returns the percentile it would adopt and why, plus the
    per-step curve it argued from (for floor_knee.csv)."""
    gain = _norm(_C)

    if rule == 'resource':
        # cost = the share of the state's riparian cropland that raising the
        # floor pushes outside the map. Monotone, smooth, and the number a
        # reviewer actually asks about.
        cost = _norm(-_R, -float(np.nanmax(_R[_OKC])), -float(np.nanmin(_R[_OKC])))
        return _pick(gain - FLOOR_KNEE_WEIGHT * cost,
                     f'concentration gained vs resource left out, weight '
                     f'{FLOOR_KNEE_WEIGHT:g}'), gain, cost

    if rule == 'counties':
        # cost = counties that end up with no zone at all. Administrative
        # rather than ecological, and coarse: with few counties one of them is
        # a large fraction of the whole range.
        n0, n1 = float(_NC[_OKC][0]), float(np.nanmin(_NC[_OKC]))
        cost = ((n0 - _NC) / max(n0 - n1, 1e-9) if n0 > n1
                else np.zeros_like(_NC))
        return _pick(gain - FLOOR_KNEE_WEIGHT * cost,
                     f'concentration gained vs counties losing their last '
                     f'zone, weight {FLOOR_KNEE_WEIGHT:g}'), gain, cost

    if rule == 'curve':
        # The capture curve, with no weight to argue about. Normalise resource
        # captured and area designated to 0-1 and take the point furthest
        # above the diagonal - the classic knee, and identical to maximising
        # (captured - designated) in normalised units.
        rn, an = _norm(_R), _norm(_A)
        return _pick(rn - an,
                     'furthest point above the chord of the capture curve '
                     '(no parameters)'), rn, an

    if rule == 'marginal':
        # The exchange rate, in the units everything else is already in.
        # Between two adjacent floors, the ground that LEAVES the map has its
        # own density: dropped riparian acres over dropped land acres. Keep
        # raising the floor while what leaves is worse than the state average;
        # stop when it would start costing better-than-average ground.
        d_rip = np.r_[np.diff(-_HR), 0.0]
        d_land = np.r_[np.diff(-_HA), 0.0]
        raw = np.where(d_land > 1e-9, d_rip / np.maximum(d_land, 1e-9) * 640,
                       np.nan)
        # It is a ratio of two small differences, so it spikes wherever a step
        # drops almost no ground. A centred rolling median keeps the level and
        # drops the spikes; without it one noisy step picks the floor.
        marg = (pd.Series(raw).rolling(FLOOR_KNEE_MARGINAL_SMOOTH,
                                       center=True, min_periods=1)
                .median().to_numpy()
                if FLOOR_KNEE_MARGINAL_SMOOTH > 1 else raw)
        knee_curve['marginal_raw_ac_per_sqmi'] = np.round(raw, 2)
        thr = FLOOR_KNEE_MARGINAL_X * AVG_STATE
        bad = np.isfinite(marg) & (marg >= thr) & _OKC
        i = int(np.flatnonzero(bad)[0]) if bad.any() else int(
            np.flatnonzero(_OKC)[-1])
        return ({'pctl': float(_kp[i]),
                 'floor_ac_per_sqmi': float(_floor_at(_kp[i])),
                 'note': f'the ground the next step would drop reaches '
                         f'{FLOOR_KNEE_MARGINAL_X:g}x the state average '
                         f'({thr:.2f} ac/sq mi)',
                 'net': marg}, marg, np.full_like(marg, thr))

    raise ValueError(f'unknown FLOOR_KNEE_RULE: {rule!r}')


def knee_plateau(rule, pick):
    """(lo, hi) - the floors this rule cannot tell apart from its pick.

    For the three argmax rules it is the CONTIGUOUS run of swept floors around
    the pick whose net is within FLOOR_KNEE_PLATEAU_TOL of the best, measured
    against that rule's own net range. Contiguous on purpose: a scatter of
    equally good points elsewhere on the axis is a different statement from a
    flat optimum, and only the second one means "the argmax is arbitrary here".

    'marginal' is a crossing rather than an argmax. The true crossing lies
    somewhere between the last swept floor under the threshold and the first
    one over it, so that bracket is its range - which is also why it is
    usually the tightest of the four.
    """
    if not np.isfinite(pick.get('pctl', np.nan)):
        return (np.nan, np.nan)
    i = int(np.argmin(np.abs(_kp - pick['pctl'])))
    if rule == 'marginal':
        j = min(i + 1, len(_kp) - 1)
        return (float(_floor_at(_kp[i])), float(_floor_at(_kp[j])))
    net = np.asarray(pick['net'], dtype=float)
    ok = np.isfinite(net) & _OKC
    if not ok.any() or not ok[i]:
        return (np.nan, np.nan)
    span = float(net[i]) - float(np.nanmin(net[ok]))
    near = ok & (net >= float(net[i]) - FLOOR_KNEE_PLATEAU_TOL * max(span, 1e-12))
    lo = hi = i
    while lo - 1 >= 0 and near[lo - 1]:
        lo -= 1
    while hi + 1 < len(net) and near[hi + 1]:
        hi += 1
    return (float(_floor_at(_kp[lo])), float(_floor_at(_kp[hi])))


KNEE_RULES = ('resource', 'counties', 'curve', 'marginal')
KNEE_OK = bool(_OKC.sum() >= 3 and _SPAN >= FLOOR_KNEE_MIN_RANGE)
knee_picks, _rule_rows = {}, []
if KNEE_OK:
    for _r in KNEE_RULES:
        _pk, _g, _c = knee_by(_r)
        knee_picks[_r] = _pk
        if _r == 'marginal':
            # this rule has no normalised gain/cost pair - what it argues
            # from is a density in ac/sq mi against a threshold in ac/sq mi,
            # so the columns say that rather than wearing borrowed names
            knee_curve['marginal_ac_per_sqmi'] = np.round(_g, 2)
            knee_curve['marginal_threshold'] = np.round(_c, 2)
        else:
            knee_curve[f'{_r}_gain'] = np.round(_g, 4)
            knee_curve[f'{_r}_cost'] = np.round(_c, 4)
            knee_curve[f'{_r}_net'] = np.round(
                np.where(np.isfinite(_pk['net']), _pk['net'], np.nan), 4)
        _row = sweep_at(_pk['floor_ac_per_sqmi'], _pk['pctl'])
        _plo, _phi = knee_plateau(_r, _pk)
        _row.update({'rule': _r, 'adopted': _r == FLOOR_KNEE_RULE,
                     'plateau_lo': round(_plo, 2), 'plateau_hi': round(_phi, 2),
                     'why': _pk['note']})
        _rule_rows.append(_row)

KNEE_PCTL = knee_picks.get(FLOOR_KNEE_RULE, {}).get('pctl', np.nan)
KNEE_FLOOR = knee_picks.get(FLOOR_KNEE_RULE, {}).get('floor_ac_per_sqmi',
                                                     np.nan)
KNEE_OK = KNEE_OK and np.isfinite(KNEE_FLOOR)

knee_rules = pd.DataFrame(_rule_rows)
if len(knee_rules):
    _RCOLS = ['rule', 'adopted', 'pctl', 'floor_ac_per_sqmi',
              'plateau_lo', 'plateau_hi',
              'pct_resource_hot', 'pct_area_hot', 'concentration', 'n_zones',
              'counties_with_zones', 'hot_sqmi', 'why']
    knee_rules = knee_rules[_RCOLS]
    knee_rules.to_csv(os.path.join(S_VAL, 'floor_knee_rules.csv'), index=False)
KNEE_PLATEAU = (np.nan, np.nan)
if KNEE_OK and len(knee_rules) and (knee_rules['rule'] == FLOOR_KNEE_RULE).any():
    _kpl = knee_rules.loc[knee_rules['rule'] == FLOOR_KNEE_RULE].iloc[0]
    KNEE_PLATEAU = (float(_kpl['plateau_lo']), float(_kpl['plateau_hi']))

knee_curve['is_knee'] = knee_curve['pctl'].eq(KNEE_PCTL)
knee_curve.to_csv(os.path.join(S_VAL, 'floor_knee.csv'), index=False)

# the adopted knee's own row joins the reported table, flagged, so one CSV
# answers "what was adopted, what was recommended, and what does each keep?"
floor_sweep['adopted'] = False
floor_sweep['knee'] = False
if KNEE_OK and not np.isclose(KNEE_PCTL, FLOOR_SWEEP_PCTLS, atol=1e-9).any():
    _kr = sweep_at(KNEE_FLOOR, KNEE_PCTL)
    _kr.update({'adopted': False, 'knee': True})
    floor_sweep = pd.concat([floor_sweep, pd.DataFrame([_kr])],
                            ignore_index=True).sort_values('pctl')
elif KNEE_OK:
    floor_sweep.loc[np.isclose(floor_sweep['pctl'], KNEE_PCTL), 'knee'] = True

# Step 3b. Choose the floor.
if FLOOR_MODE == 'off':
    FLOOR = 0.0
elif FLOOR_MODE == 'absolute':
    FLOOR = float(FLOOR_AC_PER_SQMI)
elif FLOOR_MODE == 'knee':
    FLOOR = KNEE_FLOOR if KNEE_OK else float(np.percentile(_dens, FLOOR_PCTL))
else:
    FLOOR = float(np.percentile(_dens, FLOOR_PCTL))

_adopted_pctl = (100.0 * float((_dens < FLOOR).mean()) if FLOOR > 0 else 0.0)
# one flag, on the LOWEST percentile that produces this floor - several
# percentiles give floor 0.00 when a quarter of the state holds nothing
floor_sweep = floor_sweep.reset_index(drop=True)
_near = (floor_sweep['floor_ac_per_sqmi'] - round(FLOOR, 2)).abs()
if len(_near) and float(_near.min()) < 0.005:
    floor_sweep.loc[int(_near.idxmin()), 'adopted'] = True
floor_sweep.round(3).to_csv(os.path.join(S_VAL, 'floor_sweep.csv'),
                            index=False)
print(floor_sweep.to_string(index=False))
print('\n  concentration = share of the riparian cropland divided by share of '
      'the area.\n  Raising the floor always buys concentration and always '
      'costs capture; every rule\n  below is a different answer to when that '
      'trade stops paying.')

if KNEE_OK:
    print(f'\n  WHERE EACH RULE WOULD PUT THE FLOOR  (searched on a '
          f'{FLOOR_KNEE_GRID:g}-percentile grid)')
    print(knee_rules.drop(columns=['why']).to_string(index=False))
    for _, _rr in knee_rules.iterrows():
        print(f'    {_rr["rule"]:<9} {_rr["why"]}')
    print(f'    -> {os.path.join(S_VAL, "floor_knee_rules.csv")} and '
          f'{os.path.join(S_VAL, "floor_knee.csv")}')
    # Agreement between rules is the real signal. Four different arguments
    # landing in the same place is worth more than any one of them; a wide
    # spread means the sweep has no knee worth the name and the floor is a
    # judgement call that should be made and recorded as one.
    _fl = knee_rules['floor_ac_per_sqmi'].to_numpy(float)
    _lo, _hi = float(np.nanmin(_fl)), float(np.nanmax(_fl))
    if _hi <= 1.25 * max(_lo, 1e-9):
        print(f'    The four rules agree to within '
              f'{_hi / max(_lo, 1e-9):.2f}x - a real knee.')
    else:
        print(f'    The rules span {_lo:.2f} to {_hi:.2f} ac/sq mi '
              f'({_hi / max(_lo, 1e-9):.1f}x). They are answering different '
              f'questions and\n    this sweep has no single knee - pick the '
              f'question, and say which in the brief.')
    # HOW SHARP each pick is. A rule whose plateau spans a 1.5x range of
    # floors is not choosing between them on any evidence the rule contains;
    # saying so is the difference between a derived number and a derived
    # number that is also defensible to two decimal places.
    print(f'\n  HOW SHARP IS EACH PICK  (floors each rule scores within '
          f'{FLOOR_KNEE_PLATEAU_TOL:.0%} of its own best)')
    for _, _rr in knee_rules.iterrows():
        _lo, _hi = float(_rr['plateau_lo']), float(_rr['plateau_hi'])
        if not (np.isfinite(_lo) and np.isfinite(_hi)):
            continue
        _wide = _hi > 1.5 * max(_lo, 1e-9)
        _tail = ('the bracket the crossing falls in' if _rr['rule'] == 'marginal'
                 else '  <- FLAT: the argmax is choosing between near-equals'
                 if _wide else '  <- sharp')
        print(f'    {_rr["rule"]:<9} pick {_rr["floor_ac_per_sqmi"]:6.2f}   '
              f'plateau {_lo:6.2f} - {_hi:6.2f}   {_tail}')
    # WHERE THE PLATEAUS OVERLAP. Two rules rating the same floor near-optimal
    # is a stronger argument than either rule's argmax, because it does not
    # depend on which cost term you believe.
    _grid_fl = np.array([_floor_at(_p) for _p in _kp], dtype=float)
    _cov = np.zeros(len(_kp), dtype=int)
    for _, _rr in knee_rules.iterrows():
        _lo, _hi = float(_rr['plateau_lo']), float(_rr['plateau_hi'])
        if np.isfinite(_lo) and np.isfinite(_hi):
            _cov += ((_grid_fl >= _lo - 1e-9) &
                     (_grid_fl <= _hi + 1e-9)).astype(int)
    if _cov.max() >= 2:
        _sel = _grid_fl[_cov == _cov.max()]
        _names = [_rr['rule'] for _, _rr in knee_rules.iterrows()
                  if np.isfinite(float(_rr['plateau_lo']))
                  and float(_rr['plateau_lo']) - 1e-9 <= _sel.min()
                  and float(_rr['plateau_hi']) + 1e-9 >= _sel.max()]
        print(f'    {_cov.max()} of {len(knee_rules)} plateaus overlap at '
              f'{_sel.min():.2f} - {_sel.max():.2f} ac/sq mi'
              + (f'  ({", ".join(_names)})' if _names else ''))
    else:
        print('    No two plateaus overlap - the rules are not merely '
              'imprecise, they disagree.')

    _nspan = int(knee_curve['counties_with_zones'].max()
                 - knee_curve['counties_with_zones'].min())
    if FLOOR_KNEE_RULE == 'counties' and _nspan and _nspan < 10:
        print(f'    NOTE: counties_with_zones spans only {_nspan} value(s), '
              f'so one county is {1 / _nspan:.0%} of\n    the cost term - '
              f"the 'counties' rule is coarse on this run.")
else:
    print(f'\n  NO KNEE: concentration moves by only {_SPAN:.3f} across the '
          f'whole sweep, so there is\n  no trade-off to find. Any floor in '
          f'this range does much the same thing.')

_FLOOR_WHY = (
    f'{FLOOR_PCTL:g}th percentile of the statewide cell distribution'
    if FLOOR_MODE == 'percentile' else
    f'the knee, {FLOOR_KNEE_RULE!r} rule'
    if FLOOR_MODE == 'knee' and KNEE_OK else
    f'FLOOR_MODE = {FLOOR_MODE!r} found no knee, so FLOOR_PCTL = '
    f'{FLOOR_PCTL:g} was used instead' if FLOOR_MODE == 'knee' else
    'set by hand' if FLOOR_MODE == 'absolute' else 'no floor')
print(f'\n  FLOOR = {FLOOR:.2f} riparian crop ac per sq mi  ({_FLOOR_WHY})')
print(f'  state average is {AVG_STATE:.2f}, so the floor is '
      f'{FLOOR / max(AVG_STATE, 1e-9):.2f}x it'
      + (f' and sits at the {_adopted_pctl:.0f}th percentile of the cells.'
         if FLOOR > 0 else '.'))
if KNEE_OK and FLOOR_MODE != 'knee' and FLOOR > 0:
    _gap = FLOOR / max(KNEE_FLOOR, 1e-9)
    print(f"  The {FLOOR_KNEE_RULE!r} knee is at {KNEE_FLOOR:.2f}; the "
          f'adopted floor is {_gap:.2f}x that. '
          + ('Close enough that the choice barely matters.'
             if 0.8 <= _gap <= 1.25 else
             "Set FLOOR_MODE = 'knee' to adopt it instead."))

# %% C4 - hot cells and zones
# Step 4. THE RULE. Both conditions, plus the own-density condition the
# notebook already had: the 3-mile window covers ~28 sq mi, so a cell can clear
# the neighbourhood test purely by sitting NEAR a concentration while holding
# almost none itself. Requiring its own density to clear the reference too
# means a hot cell IS a concentration rather than a neighbour of one.
#
#   BASIS = 'county' and FLOOR_MODE = 'off' reduce this to the notebook's rule
#   exactly:  (vs_county >= HOT_MULTIPLE) & (own_per_sqmi >= county_avg)
rule('HOT CELLS AND ZONES')

grid['ratio'] = grid[f'ratio_{BASIS}']
grid['reference'] = grid[f'avg_{BASIS}']
grid['hot'] = ((grid['ratio'] >= HOT_MULTIPLE) &
               (grid['own_per_sqmi'] >= grid['reference']) &
               (grid['rip_per_sqmi'] >= FLOOR))
grid['threshold_applied'] = np.maximum(
    HOT_MULTIPLE * grid['reference'].to_numpy(float), FLOOR).round(2)

_nh = int(grid['hot'].sum())
print(f'  basis                : {BASIS}')
print(f'  ratio condition      : neighbourhood density >= {HOT_MULTIPLE:g}x its '
      f'{BASIS} average')
print(f'  floor condition      : neighbourhood density >= {FLOOR:.2f} ac/sq mi')
print(f'  hot cells            : {_nh:,} of {len(grid):,} '
      f'({_nh / len(grid):.1%} of the analysed ground)')
print(f'  riparian cropland in them: '
      f'{grid.loc[grid["hot"], "rip_ac"].sum():,.0f} ac of '
      f'{grid["rip_ac"].sum():,.0f} '
      f'({100 * grid.loc[grid["hot"], "rip_ac"].sum() / max(grid["rip_ac"].sum(), 1e-9):.0f}%)')
if _nh:
    _cr = ((100 * grid.loc[grid['hot'], 'rip_ac'].sum() / max(grid['rip_ac'].sum(), 1e-9)) /
           (100 * grid.loc[grid['hot'], 'land_ac'].sum() / grid['land_ac'].sum()))
    print(f'  concentration        : {_cr:.1f}x')

# Step 5. Contiguous hot cells become zones.
hot = grid[grid['hot']]
if len(hot):
    _merged = union_all(np.asarray(hot.geometry.values, dtype=object))
    _pieces = (list(_merged.geoms) if _merged.geom_type.startswith('Multi')
               else [_merged])
    zones = gpd.GeoDataFrame(geometry=_pieces, crs=WORKING_CRS)
    zones['n_cells'] = (zones.geometry.area / (CELL_M ** 2)).round().astype(int)
    zones = zones[zones['n_cells'] >= MIN_ZONE_CELLS].reset_index(drop=True)
else:
    zones = gpd.GeoDataFrame(geometry=[], crs=WORKING_CRS)

if len(zones):
    _zg = np.asarray(zones.geometry.values, dtype=object)
    zones['zone_ac'] = shapely.area(
        shapely.intersection(_zg, ANALYSIS_GEOM)) / AC
    _rg = np.asarray(rip_crop.geometry.values, dtype=object)
    _zt = shapely.STRtree(_zg)
    _ri, _zi = _zt.query(_rg, predicate='intersects')
    _ia = shapely.area(shapely.intersection(_rg[_ri], _zg[_zi])) / AC
    zones['rip_crop_ac'] = (pd.Series(_ia).groupby(_zi).sum()
                            .reindex(range(len(zones))).fillna(0.0).to_numpy())
    zones['rip_per_sqmi'] = (zones['rip_crop_ac'] /
                             (zones['zone_ac'] / 640)).round(1)
    zones['pct_inside'] = (zones['zone_ac'] /
                           (zones['n_cells'] * CELL_M ** 2 / AC)).round(2)
    zones['vs_state'] = (zones['rip_per_sqmi'] / AVG_STATE).round(2)
    zones = zones.sort_values('rip_crop_ac', ascending=False).reset_index(drop=True)
    zones.index.name = 'zone_id'
    print(f'\n  {len(zones):,} zones of at least {MIN_ZONE_CELLS} cells, '
          f'{zones["zone_ac"].sum() / 640:,.0f} sq mi, '
          f'{zones["rip_crop_ac"].sum():,.0f} ac of riparian cropland')
else:
    print('\n  no zones met the threshold anywhere - lower HOT_MULTIPLE or the '
          'floor')


# %% C5 - zones cross county lines, so say so
# Step 6. With one statewide surface a zone is free to straddle a boundary, and
# many do. Each zone is kept WHOLE - cutting it would report two twenty-square-
# mile zones where there is one of forty - and split into per-county rows in
# zone_county_split.csv. A shared zone is useful information for a biologist:
# it names who to call next door.
if len(zones):
    _zg = np.asarray(zones.geometry.values, dtype=object)
    _eg = np.asarray(counties['analysed_geom'].values, dtype=object)
    _et = shapely.STRtree(_eg)
    _zi, _ci = _et.query(_zg, predicate='intersects')
    _sa = shapely.area(shapely.intersection(_zg[_zi], _eg[_ci])) / AC
    split = pd.DataFrame({'zone_id': _zi, 'FIPS': counties['FIPS'].to_numpy()[_ci],
                          'NAME': counties['NAME'].to_numpy()[_ci],
                          'zone_ac_in_county': _sa})
    # A zone is a union of whole 1-mile cells, so it laps a little way over
    # the boundary of the county it belongs to. MIN_ZONE_COUNTY_AC (config) is
    # the cut that calls that overhang rather than a share of the zone.
    split = split[split['zone_ac_in_county'] >= MIN_ZONE_COUNTY_AC].copy()

    # the riparian cropland inside each zone-county piece
    _rga = np.asarray(rip_crop.geometry.values, dtype=object)
    _rt = shapely.STRtree(_rga)
    _row_of = {f: i for i, f in enumerate(counties['FIPS'].to_numpy())}
    _pieces_ac = []
    for _zid, _fips in zip(split['zone_id'].to_numpy(), split['FIPS'].to_numpy()):
        _g = _zg[int(_zid)].intersection(_eg[_row_of[_fips]])
        _hits = _rt.query(_g, predicate='intersects')
        _pieces_ac.append(
            float(shapely.area(shapely.intersection(_rga[_hits], _g)).sum() / AC)
            if len(_hits) else 0.0)
    split['rip_crop_ac_in_county'] = _pieces_ac
    split['pct_of_zone'] = (split['zone_ac_in_county'] /
                            split['zone_id'].map(zones['zone_ac'])).round(3)
    split = split.sort_values(['zone_id', 'zone_ac_in_county'],
                              ascending=[True, False]).reset_index(drop=True)
    split.round(2).to_csv(os.path.join(S_OUT, 'zone_county_split.csv'),
                          index=False)

    _nz = split.groupby('zone_id')['FIPS'].nunique()
    zones['n_counties'] = zones.index.map(_nz).fillna(1).astype(int)
    _primary = split.sort_values('zone_ac_in_county', ascending=False) \
                    .drop_duplicates('zone_id').set_index('zone_id')
    zones['primary_FIPS'] = zones.index.map(_primary['FIPS'])
    zones['primary_county'] = zones.index.map(_primary['NAME'])
    zones['counties'] = zones.index.map(
        split.groupby('zone_id')['NAME'].apply(lambda s: ', '.join(sorted(set(s)))))
    _shared = int((zones['n_counties'] > 1).sum())
    print(f'  {_shared:,} of {len(zones):,} zones cross a county line '
          f'({_shared / max(len(zones), 1):.0%})')
    print(f'  -> {os.path.join(S_OUT, "zone_county_split.csv")}')
    gpkg_safe(zones.reset_index()).to_file(SURFACE_GPKG, layer='zones',
                                           driver='GPKG')
else:
    split = pd.DataFrame(columns=['zone_id', 'FIPS', 'NAME',
                                  'zone_ac_in_county',
                                  'rip_crop_ac_in_county', 'pct_of_zone'])


# %% C5b - the water statistics a biologist asks for, on every zone
# Step 6b. The zone table so far says how big a zone is and how much riparian
# cropland is in it. That is not what gets asked in the room. The questions are
# "which creek is this", "how much stream is in there", "is it wetland or is it
# river", and "is there another one next door". Those are the same fields the
# single-county run put on a PARCEL (n_stream_segments, rip_ac_<source>,
# adjacency) - this puts them on a ZONE, which is the unit that exists whether
# or not anyone has bought parcel data.
#
# Everything here is measured against the SAME geometry the zone acres come
# from, so a stream mile and an acre in the same row describe the same ground.
_NHD_STATS = [
    ('stream',    'nhd_flow', ['permanent_identifier', 'reachcode'], 'line'),
    ('river',     'nhd_area', ['permanent_identifier'], 'polygon'),
    ('waterbody', 'nhd_wb',   ['permanent_identifier'], 'polygon'),
    ('wetland',   'wetlands', ['GLOBALID', 'OBJECTID', 'WETLAND_ID'], 'polygon'),
]
_TREES = {}


def _src_tree(var):
    """One STRtree per source layer, built once and reused. Building the
    statewide flowline index is seconds; building it per zone would be hours."""
    if var not in _TREES:
        g = globals().get(var)
        if g is None or not len(g):
            _TREES[var] = (None, None, None)
        else:
            arr = np.asarray(g.geometry.values, dtype=object)
            _TREES[var] = (arr, shapely.STRtree(arr), g)
    return _TREES[var]


def water_stats(geoms, names=True):
    """Per-geometry water statistics: riparian cropland split by source, NHD
    feature counts, stream miles inside, and the named creeks.

    Counts are of features that INTERSECT the geometry; miles and acres are
    CLIPPED to it. Those are different questions and the column names say
    which - n_stream_segments answers "how many reaches touch this", stream_mi
    answers "how much channel is inside it"."""
    g = np.asarray(geoms, dtype=object)
    out = pd.DataFrame(index=pd.RangeIndex(len(g)))
    if not len(g):
        return out

    # riparian cropland, split by the water it sits on
    rg = np.asarray(rip_crop.geometry.values, dtype=object)
    _ri, _gi = shapely.STRtree(g).query(rg, predicate='intersects')
    _ia = shapely.area(shapely.intersection(rg[_ri], g[_gi])) / AC
    for _n in USE_SOURCES:
        _c = f'on_{_n}'
        _m = (rip_crop[_c].to_numpy(bool)[_ri] if _c in rip_crop.columns
              else np.zeros(len(_ri), bool))
        out[f'rip_ac_{_n}'] = (pd.Series(_ia[_m]).groupby(_gi[_m]).sum()
                               .reindex(range(len(g))).fillna(0.0)
                               .to_numpy().round(1))

    # NHD / NWI features
    for _n, _var, _pref, _kind in _NHD_STATS:
        _arr, _tree, _gdf = _src_tree(_var)
        _cnt = {'stream': 'n_stream_segments', 'river': 'n_rivers',
                'waterbody': 'n_waterbodies',
                'wetland': 'n_wetlands'}.get(_n, f'n_{_n}s')
        if _tree is None:
            out[_cnt] = 0
            if _kind == 'line':
                out[f'{_n}_mi'] = 0.0
            continue
        # STRtree.query(array) returns (index into `g`, index into the tree)
        _qi, _si = _tree.query(g, predicate='intersects')
        _idc = next((c for c in _pref if c in _gdf.columns), None)
        if _idc is not None:
            _key = _gdf[_idc].astype(str).to_numpy()[_si]
            out[_cnt] = (pd.DataFrame({'q': _qi, 'k': _key})
                         .groupby('q')['k'].nunique()
                         .reindex(range(len(g))).fillna(0).astype(int).to_numpy())
        else:
            out[_cnt] = (pd.Series(np.ones(len(_qi))).groupby(_qi).sum()
                         .reindex(range(len(g))).fillna(0).astype(int).to_numpy())
        if _kind == 'line':
            _ln = shapely.length(
                shapely.intersection(_arr[_si], g[_qi])) / MI_M
            out[f'{_n}_mi'] = (pd.Series(_ln).groupby(_qi).sum()
                               .reindex(range(len(g))).fillna(0.0)
                               .to_numpy().round(2))
            if names and 'gnis_name' in _gdf.columns:
                _nm = _gdf['gnis_name'].astype(str).to_numpy()[_si]
                _ok = ~pd.isna(_nm) & (_nm != 'nan') & (_nm != '')
                _t = (pd.DataFrame({'q': _qi[_ok], 'nm': _nm[_ok],
                                    'mi': _ln[_ok]})
                      .groupby(['q', 'nm'])['mi'].sum().reset_index()
                      .sort_values(['q', 'mi'], ascending=[True, False]))
                out['named_waters'] = (
                    _t.groupby('q')['nm'].apply(lambda s: ', '.join(s.head(6)))
                    .reindex(range(len(g))).fillna('').to_numpy())
    return out


# ---- what is growing in a zone, and which subwatersheds it sits in ---
# Step 6b-ii. PASS A/B ONLY, AND PARCEL-INVARIANT BY CONSTRUCTION. Everything
# below reads `rip_crop` (Pass A), `huc12` (the cache) and three tables in
# fl_config - CDL_GROUPS, GRP_ORDER, CDL_NAMES. It never opens a parcel file,
# never imports anything Pass C defines, and writes only into the _state
# folder. A county with no parcel licence gets these columns in full, which is
# the whole point: the zone map is the deliverable that costs nothing per
# county, and "what is growing there" is part of the zone, not part of a
# parcel.
#
# The FIELD NAMES match the ones Pass C puts on a parcel - rip_<group>_ac,
# rip_cover, rip_cdl_dominant, huc12 - so a zone row and a parcel row can be
# read with one legend. That is a shared DEFINITION, not shared code: both
# sides derive their grouping from CDL_GROUPS in the config, independently,
# so neither can run without the other and neither can drift from it.
#
# TWO DIFFERENCES FROM THE PARCEL FORM, both forced by scale:
#
#   A parcel sits in ONE subwatershed. A zone is tens of square miles and
#   routinely spans several, so membership is a LIST (`n_huc12`, `huc12_ids`)
#   and the largest share is named on its own as `primary_huc12`.
#
#   A parcel's crop string is short enough to carry whole. A zone can hold
#   thirty classes, so `rip_cdl_acres` carries the largest ZONE_CDL_TOP and
#   the full split goes to zone_cover_by_cdl.csv.
#
# THE CDL ACRES ADD UP AND THE WATER ACRES DO NOT, and that is not an
# inconsistency to reconcile. `rip_ac_<source>` counts a piece once per
# corridor it sits in, so a field on a stream running through a wetland is in
# both columns and the union is `rip_crop_ac`. A piece carries exactly one CDL
# code, so the `rip_<group>_ac` columns sum to `rip_crop_ac` instead.
#
# An unlisted code is a row crop. That is what CDL's long tail of minor crops
# actually is, and it is the same default the config's grouping implies.
_ZONE_CDL_GRP = {c: g for g, codes in CDL_GROUPS.items() for c in codes}

# A zone laps into the subwatershed next door as well as into the county next
# door, and for the same reason - it is a union of whole 1-mile cells. Smaller
# cut than MIN_ZONE_COUNTY_AC because a HUC12 is 20-40 sq mi where a county is
# ~300: one cell of a HUC12 is a real share of it, not an overhang.
MIN_ZONE_HUC12_AC = float(globals().get('MIN_ZONE_HUC12_AC', 64.0))
ZONE_CDL_TOP = int(globals().get('ZONE_CDL_TOP', 8))
ZONE_HUC12_TOP = int(globals().get('ZONE_HUC12_TOP', 8))


def cover_huc_stats(geoms, names=True):
    """Cover mix and subwatershed membership, per geometry.

    Returns (wide, cdl_long, huc_long) - one row per geometry, plus the two
    unrolled tables the digest columns summarise. Acres are CLIPPED to the
    geometry, so a zone's cover acres describe that zone's ground and nothing
    beyond it; HUC12 counts are of subwatersheds the geometry INTERSECTS by
    more than MIN_ZONE_HUC12_AC."""
    g = np.asarray(geoms, dtype=object)
    n = len(g)
    wide = pd.DataFrame(index=pd.RangeIndex(n))
    cdl_long = pd.DataFrame(columns=['row', 'CDL', 'cdl_name', 'grp',
                                     'rip_crop_ac'])
    huc_long = pd.DataFrame(columns=['row', 'huc12', 'huc12_name',
                                     'ac_in_huc12', 'pct_of_zone'])
    if not n:
        return wide, cdl_long, huc_long

    # --- what is growing there, by CDL class ------------------------------
    for _g in GRP_ORDER:
        wide[f'rip_{_g}_ac'] = 0.0
    wide['rip_cover'] = 'none'
    wide['rip_n_cdl'] = 0
    wide['rip_cdl_dominant'] = pd.NA
    wide['rip_cdl_dominant_name'] = ''
    wide['rip_cdl_acres'] = ''

    if 'CDL' in rip_crop.columns:
        rg = np.asarray(rip_crop.geometry.values, dtype=object)
        _ri, _gi = shapely.STRtree(g).query(rg, predicate='intersects')
        if len(_ri):
            _ia = shapely.area(shapely.intersection(rg[_ri], g[_gi])) / AC
            _d = pd.DataFrame({
                'row': _gi,
                'CDL': pd.to_numeric(pd.Series(rip_crop['CDL'].to_numpy()[_ri]),
                                     errors='coerce'),
                'rip_crop_ac': _ia})
            _d = _d[_d['rip_crop_ac'] > 0]
        else:
            _d = pd.DataFrame(columns=['row', 'CDL', 'rip_crop_ac'])
        if len(_d):
            cdl_long = (_d.groupby(['row', 'CDL'], dropna=False)['rip_crop_ac']
                        .sum().reset_index())
            cdl_long['CDL'] = cdl_long['CDL'].astype('Int64')
            cdl_long['cdl_name'] = [
                CDL_NAMES.get(int(c), f'CDL {int(c)}') if pd.notna(c)
                else 'unclassified' for c in cdl_long['CDL']]
            cdl_long['grp'] = [
                _ZONE_CDL_GRP.get(int(c), 'row_crop') if pd.notna(c)
                else 'row_crop' for c in cdl_long['CDL']]
            cdl_long = cdl_long.sort_values(['row', 'rip_crop_ac'],
                                            ascending=[True, False])

            _gp = cdl_long.pivot_table(index='row', columns='grp',
                                       values='rip_crop_ac', aggfunc='sum',
                                       fill_value=0.0)
            for _g in GRP_ORDER:
                if _g in _gp.columns:
                    wide[f'rip_{_g}_ac'] = (_gp[_g].reindex(range(n))
                                            .fillna(0.0).round(1).to_numpy())
            _gc = [f'rip_{x}_ac' for x in GRP_ORDER]
            wide['rip_cover'] = (wide[_gc].idxmax(axis=1)
                                 .str.replace('rip_', '', regex=False)
                                 .str.replace('_ac', '', regex=False)
                                 .where(wide[_gc].max(axis=1) > 0, 'none'))
            wide['rip_n_cdl'] = (cdl_long.groupby('row').size()
                                 .reindex(range(n)).fillna(0).astype(int)
                                 .to_numpy())
            _dom = cdl_long.drop_duplicates('row').set_index('row')
            wide['rip_cdl_dominant'] = _dom['CDL'].reindex(range(n)).to_numpy()
            wide['rip_cdl_dominant_name'] = (
                _dom['cdl_name'].reindex(range(n)).fillna('').to_numpy())
            # built as a column and grouped as a Series, so this does not
            # depend on DataFrameGroupBy.apply and its include_groups churn
            _lab = (cdl_long['cdl_name'] + ':'
                    + cdl_long['rip_crop_ac'].round(0).astype(int).astype(str))
            wide['rip_cdl_acres'] = (
                _lab.groupby(cdl_long['row'])
                .apply(lambda s: ', '.join(s.head(ZONE_CDL_TOP)))
                .reindex(range(n)).fillna('').to_numpy())

    # --- which subwatersheds it sits in -----------------------------------
    wide['n_huc12'] = 0
    wide['primary_huc12'] = ''
    if names:
        wide['primary_huc12_name'] = ''
    wide['pct_in_primary_huc12'] = np.nan
    wide['huc12_ids'] = ''

    _hc = (next((c for c in huc12.columns if c.lower() == 'huc12'), None)
           if huc12 is not None and len(huc12) else None)
    if _hc:
        _h = huc12[huc12[_hc].astype(str).str.len() == 12]
        _hn = next((c for c in _h.columns
                    if c.lower() in ('name', 'huc12_name')), None)
        if len(_h):
            hg = np.asarray(_h.geometry.values, dtype=object)
            _hi, _gi2 = shapely.STRtree(g).query(hg, predicate='intersects')
            if len(_hi):
                _ha = shapely.area(shapely.intersection(hg[_hi],
                                                        g[_gi2])) / AC
                huc_long = pd.DataFrame({
                    'row': _gi2,
                    'huc12': _h[_hc].astype(str).to_numpy()[_hi],
                    'huc12_name': (_h[_hn].astype(str).to_numpy()[_hi]
                                   if _hn else ''),
                    'ac_in_huc12': _ha})
                huc_long = huc_long[
                    huc_long['ac_in_huc12'] >= MIN_ZONE_HUC12_AC].copy()
        if len(huc_long):
            _tot = pd.Series(shapely.area(g) / AC)
            huc_long['pct_of_zone'] = (
                huc_long['ac_in_huc12']
                / huc_long['row'].map(_tot).replace(0, np.nan)).round(3)
            huc_long = huc_long.sort_values(['row', 'ac_in_huc12'],
                                            ascending=[True, False])
            wide['n_huc12'] = (huc_long.groupby('row').size()
                               .reindex(range(n)).fillna(0).astype(int)
                               .to_numpy())
            _p = huc_long.drop_duplicates('row').set_index('row')
            wide['primary_huc12'] = (_p['huc12'].reindex(range(n))
                                     .fillna('').to_numpy())
            if names:
                wide['primary_huc12_name'] = (_p['huc12_name']
                                              .reindex(range(n))
                                              .fillna('').to_numpy())
            wide['pct_in_primary_huc12'] = (_p['pct_of_zone']
                                            .reindex(range(n)).to_numpy())
            wide['huc12_ids'] = (
                huc_long.groupby('row')['huc12']
                .apply(lambda s: ', '.join(s.head(ZONE_HUC12_TOP)))
                .reindex(range(n)).fillna('').to_numpy())
    return wide, cdl_long, huc_long


if len(zones):
    rule('ZONE WATER STATISTICS')
    _zg = np.asarray(zones.geometry.values, dtype=object)
    _ws = water_stats(_zg)
    for _c in _ws.columns:
        zones[_c] = _ws[_c].to_numpy()

    # the cover mix and the subwatersheds, in the same fields Pass C
    # puts on a parcel - but measured here, off the surface, so a
    # county with no parcel licence gets them too. The two long tables
    # are the unrolled version of the digest columns, so the summary
    # stays narrow and the per-class acreage is still checkable.
    _cs, _cdl_long, _huc_long = cover_huc_stats(_zg)
    for _c in _cs.columns:
        zones[_c] = _cs[_c].to_numpy()
    ZONE_CDL_CSV = os.path.join(S_OUT, 'zone_cover_by_cdl.csv')
    ZONE_HUC_CSV = os.path.join(S_OUT, 'zone_huc12.csv')
    if len(_cdl_long):
        (_cdl_long.rename(columns={'row': 'zone_id'})
         [['zone_id', 'CDL', 'cdl_name', 'grp', 'rip_crop_ac']]
         .round(2).to_csv(ZONE_CDL_CSV, index=False))
    if len(_huc_long):
        (_huc_long.rename(columns={'row': 'zone_id'})
         .round(2).to_csv(ZONE_HUC_CSV, index=False))

    # Step 6c. ZONE ADJACENCY - the same idea as the parcel-level
    # `adj_parcel_on_list`, one level up. Zones are unions of whole cells and
    # are disjoint by construction, so "adjacent" has to mean "within working
    # distance", and the honest distance is the neighbourhood radius the
    # analysis already uses: two zones a 3-mile band apart are one trip.
    _zt2 = shapely.STRtree(_zg)
    _qi2, _si2 = _zt2.query(shapely.buffer(_zg, ZONE_NEIGHBOR_MI * MI_M),
                            predicate='intersects')
    _self = _qi2 == _si2
    zones['n_zones_within_band'] = (
        pd.Series(np.ones((~_self).sum())).groupby(_qi2[~_self]).sum()
        .reindex(range(len(zones))).fillna(0).astype(int).to_numpy())
    _d = shapely.distance(_zg[_qi2[~_self]], _zg[_si2[~_self]]) / MI_M
    zones['nearest_zone_mi'] = (
        pd.Series(_d).groupby(_qi2[~_self]).min()
        .reindex(range(len(zones))).round(2).to_numpy())

    # the reference each zone was actually judged against, and the cut it had
    # to clear - so a zone row can be read without the manifest beside it
    _zc2 = grid[grid['hot']].copy()
    _zj = gpd.sjoin(
        gpd.GeoDataFrame(geometry=gpd.GeoSeries(_zc2.geometry.centroid.values,
                                                crs=WORKING_CRS)),
        gpd.GeoDataFrame(zones[[]].reset_index()[['zone_id']],
                         geometry=_zg, crs=WORKING_CRS),
        predicate='within', how='left')
    _zc2['zone_id'] = _zj['zone_id'].to_numpy()
    _zref = _zc2.dropna(subset=['zone_id']).groupby('zone_id').agg(
        reference_used=('reference', 'median'),
        threshold_applied=('threshold_applied', 'median'),
        max_ratio=('ratio', 'max'), median_ratio=('ratio', 'median'))
    for _c in _zref.columns:
        zones[_c] = zones.index.map(_zref[_c]).astype(float).round(2)
    zones['basis'] = BASIS
    zones['report_mask'] = REPORT_MASK or 'none'
    zones['floor_ac_per_sqmi'] = round(FLOOR, 2)

    ZONE_CSV = os.path.join(S_OUT, 'zone_summary.csv')
    _ZCOLS = (['primary_FIPS', 'primary_county', 'counties', 'n_counties',
               'n_cells', 'zone_ac', 'pct_inside', 'rip_crop_ac',
               'rip_per_sqmi', 'vs_state', 'basis', 'report_mask',
               'reference_used',
               'floor_ac_per_sqmi', 'threshold_applied', 'median_ratio',
               'max_ratio']
              + [f'rip_ac_{n}' for n in USE_SOURCES]
              + ['n_stream_segments', 'stream_mi', 'n_rivers',
                 'n_waterbodies', 'n_wetlands', 'named_waters']
              + ['n_huc12', 'primary_huc12', 'primary_huc12_name',
                 'pct_in_primary_huc12', 'huc12_ids']
              + [f'rip_{g}_ac' for g in GRP_ORDER]
              + ['rip_cover', 'rip_cdl_dominant', 'rip_cdl_dominant_name',
                 'rip_n_cdl', 'rip_cdl_acres']
              + ['n_zones_within_band', 'nearest_zone_mi'])
    _ZCOLS = [c for c in _ZCOLS if c in zones.columns]
    zones.reset_index()[['zone_id'] + _ZCOLS].round(2).to_csv(ZONE_CSV,
                                                              index=False)
    gpkg_safe(zones.reset_index()).to_file(SURFACE_GPKG, layer='zones',
                                           driver='GPKG')

    print(f'  riparian cropland by the water it sits on, across all '
          f'{len(zones)} zones:')
    for _n in USE_SOURCES:
        _c = f'rip_ac_{_n}'
        if _c in zones.columns:
            print(f'    {_n:<11} {zones[_c].sum():>10,.0f} ac   '
                  f'{int((zones[_c] > 0).sum()):>3} zones have some')
    print(f'  stream inside the zones : {zones.get("stream_mi", pd.Series(dtype=float)).sum():,.0f} mi '
          f'across {int(zones.get("n_stream_segments", pd.Series(dtype=int)).sum()):,} reach intersections')
    print(f'  zones with another zone within {ZONE_NEIGHBOR_MI:g} mi: '
          f'{int((zones["n_zones_within_band"] > 0).sum())} of {len(zones)}')

    # WHAT IS ACTUALLY GROWING THERE. These acres DO sum to rip_crop_ac,
    # unlike the per-source acres above - a piece has one CDL code and
    # can sit in two corridors.
    _gc0 = [f'rip_{g}_ac' for g in GRP_ORDER if f'rip_{g}_ac' in zones]
    if _gc0:
        _tot0 = float(zones[_gc0].to_numpy().sum())
        print(f'  cover of the {zones["rip_crop_ac"].sum():,.0f} ac of '
              f'riparian cropland inside the zones:')
        for _g in GRP_ORDER:
            _c = f'rip_{_g}_ac'
            if _c in zones.columns and zones[_c].sum() > 0:
                print(f'    {_g:<11} {zones[_c].sum():>10,.0f} ac '
                      f'{100 * zones[_c].sum() / max(_tot0, 1e-9):>5.1f}%   '
                      f'{int((zones[_c] > 0).sum()):>3} zones have some')
        _dm = zones['rip_cdl_dominant_name'].replace('', np.nan).dropna()
        if len(_dm):
            _top = _dm.value_counts().head(5)
            _lbl = ', '.join(f'{k} ({v})' for k, v in _top.items())
            print(f'  largest single class in a zone: {_lbl}')
    if 'n_huc12' in zones.columns and zones['n_huc12'].sum():
        _multi = int((zones['n_huc12'] > 1).sum())
        print(f'  subwatersheds  : {_multi} of {len(zones)} zones '
              f'({_multi / max(len(zones), 1):.0%}) span more than one HUC12; '
              f'the median\n                   zone touches '
              f'{zones["n_huc12"].median():.0f} and the largest touches '
              f'{int(zones["n_huc12"].max())}, over '
              f'{int(zones["n_huc12"].sum()):,} zone x HUC12\n'
              f'                   memberships in all '
              f'(overhang under {fmt(MIN_ZONE_HUC12_AC)} ac dropped)')
    print(f'  -> {ZONE_CSV}')
    if len(_cdl_long):
        print(f'  -> {ZONE_CDL_CSV}   one row per zone x CDL class')
    if len(_huc_long):
        print(f'  -> {ZONE_HUC_CSV}        one row per zone x HUC12')

    # Step 6d. The same fields on each zone x county PIECE, so the county packet
    # and the state table cannot disagree about who has what.
    if len(split):
        _eg2 = np.asarray(counties['analysed_geom'].values, dtype=object)
        _row_of2 = {f: i for i, f in enumerate(counties['FIPS'].to_numpy())}
        _pieces = [_zg[int(z)].intersection(_eg2[_row_of2[f]])
                   for z, f in zip(split['zone_id'], split['FIPS'])]
        _ps = water_stats(_pieces, names=False)
        for _c in _ps.columns:
            split[_c] = _ps[_c].to_numpy()
        # cover and subwatershed for the PIECE, not for the whole zone -
        # a zone half in Union and half in Henderson does not have one
        # crop mix, and the county packet should carry its own.
        _cs2, _, _ = cover_huc_stats(_pieces, names=False)
        for _c in _cs2.columns:
            split[_c] = _cs2[_c].to_numpy()
        split['rip_per_sqmi_in_county'] = (
            split['rip_crop_ac_in_county'] /
            (split['zone_ac_in_county'] / 640).replace(0, np.nan)).round(1)
        split.round(2).to_csv(os.path.join(S_OUT, 'zone_county_split.csv'),
                              index=False)
        print(f'  zone_county_split.csv now carries the same water columns '
              f'per county piece')


# %% C6 - what the basis choice actually changes
# Step 7. The same surface, thresholded three ways. This is the table that
# makes the choice arguable in a meeting instead of asserted in a footnote.
rule('BASIS COMPARISON - THE SAME SURFACE, THREE REFERENCES')

_cmp = []
for _b in ('county', 'state', 'region'):
    _r = grid[f'avg_{_b}'].to_numpy(float)
    for _fl, _lab in [(0.0, 'no floor'), (FLOOR, f'floor {FLOOR:.1f}')]:
        _h = ((grid['rip_per_sqmi'].to_numpy() >= HOT_MULTIPLE * _r) &
              (_own >= _r) & (grid['rip_per_sqmi'].to_numpy() >= _fl))
        _pa = 100 * _land[_h].sum() / _land.sum()
        _pr = 100 * _ripv[_h].sum() / max(_ripv.sum(), 1e-9)
        # DOES THE FLOOR ACTUALLY DO ANYTHING ON THIS BASIS? The ratio test
        # already demands density >= HOT_MULTIPLE x reference. Wherever that
        # product is above the floor, the floor is implied and removes nothing
        # - the two rows then come out IDENTICAL, which reads like a bug and is
        # not one. This column says how many cells the floor removed, so an
        # identical pair is explained on its own row.
        _ratio_only = ((grid['rip_per_sqmi'].to_numpy() >= HOT_MULTIPLE * _r) &
                       (_own >= _r))
        _cmp.append({'basis': _b, 'floor': _lab, 'hot_cells': int(_h.sum()),
                     'hot_sqmi': round(_land[_h].sum() / 640, 1),
                     'pct_area_hot': round(_pa, 2),
                     'pct_resource_hot': round(_pr, 1),
                     'concentration': round(_pr / _pa, 2) if _pa else np.nan,
                     'counties_with_hot': int(
                         pd.Series(grid['FIPS'].to_numpy()[_h]).nunique()),
                     'cells_cut_by_floor': int(_ratio_only.sum() - _h.sum()),
                     'min_ratio_threshold': round(
                         float(HOT_MULTIPLE * np.nanmin(_r)), 2),
                     'floor_binds': bool(_fl > 0 and
                                         (_ratio_only.sum() - _h.sum()) > 0)})
basis_cmp = pd.DataFrame(_cmp)
basis_cmp.to_csv(os.path.join(S_VAL, 'basis_comparison.csv'), index=False)
print(basis_cmp.to_string(index=False))
print(f'\n  The "county / no floor" row is the notebook\'s rule applied to all '
      f'{grid["FIPS"].nunique()}\n  counties at once. Compare its '
      f'counties_with_hot with the chosen row: the gap is\n  the set of '
      f'counties whose hotspots exist only because they were measured\n'
      f'  against themselves.')

# Step 7a. WHY TWO ROWS CAN BE IDENTICAL. The floor is an absolute number; the
# ratio test is a number relative to a reference. On a basis whose reference is
# high enough, HOT_MULTIPLE x reference already exceeds the floor everywhere,
# so the floor never removes a cell and the two rows match to the digit. That
# is a real finding about the basis, not a failed computation, and it is said
# here rather than left for someone to notice in the CSV.
_dead = basis_cmp[(basis_cmp['floor'] != 'no floor') &
                  (~basis_cmp['floor_binds'])]
if FLOOR > 0 and len(_dead):
    print(f'\n  THE FLOOR IS INERT ON: {", ".join(_dead["basis"])}')
    for _, _dr in _dead.iterrows():
        print(f'    {_dr["basis"]:<7} lowest ratio cut anywhere is '
              f'{_dr["min_ratio_threshold"]:.2f} ac/sq mi, which is already '
              f'above the\n            {FLOOR:.2f} floor - so the floor removes '
              f'0 cells and the two rows are identical.')
    if 'state' in set(_dead['basis']):
        print(f'    For the STATE basis that is arithmetic, not coincidence: '
              f'one reference\n    ({AVG_STATE:.2f}) applies everywhere, so the '
              f'ratio cut is a single number\n    {HOT_MULTIPLE:g} x '
              f'{AVG_STATE:.2f} = {HOT_MULTIPLE * AVG_STATE:.2f} ac/sq mi. The '
              f'floor only ever binds when it\n    exceeds that - i.e. when '
              f'FLOOR_PCTL puts the floor above {HOT_MULTIPLE:g}x the state\n'
              f'    average. Here the floor is {FLOOR / max(AVG_STATE, 1e-9):.2f}x '
              f'it, and {FLOOR / max(AVG_STATE, 1e-9):.2f} < {HOT_MULTIPLE:g}.')
    print(f'    The floor exists for the county and region bases, where a low '
          f'local\n    reference makes {HOT_MULTIPLE:g}x cheap to clear. It is '
          f'doing its job there:')
    _live = basis_cmp[(basis_cmp['floor'] != 'no floor') &
                      (basis_cmp['floor_binds'])]
    for _, _lr in _live.iterrows():
        print(f'      {_lr["basis"]:<7} floor removes '
              f'{_lr["cells_cut_by_floor"]:,} cells')


# %% C7 - what the statewide neighbourhood fixed
# Step 8. In the county-by-county version the 3-mile band stops at the county
# line, so a border cell is averaged over less land than an interior one. This
# measures how many cells that was, and how far their density moved.
rule('BAND TRUNCATION - WHAT THE STATEWIDE SURFACE FIXED')

_own_cty = grid['FIPS'].to_numpy()
_nb_cty_rip = np.zeros(len(grid))
_nb_cty_land = np.zeros(len(grid))
# Under a reporting mask the grid was re-cut in B6, AFTER B5 built these
# neighbour lists, so _flat indexes cells that are no longer in it. The
# check compares a statewide band against a county-only band; with the
# mask on, neither side is the thing it was built to compare.
_same = (_own_cty[_flat] == np.repeat(_own_cty, _lens)
         if _want_surface and REPORT_MASK is None else None)
if _same is not None:
    _nb_cty_rip = np.add.reduceat(np.where(_same, _ripv[_flat], 0.0), _starts)
    _nb_cty_land = np.add.reduceat(np.where(_same, _land[_flat], 0.0), _starts)
    _cty_only = np.where(_nb_cty_land > 0, _nb_cty_rip / _nb_cty_land * 640, 0.0)
    grid['rip_per_sqmi_county_band'] = np.round(_cty_only, 2)
    _cut = _nb_cty_land < 0.999 * grid['nbhd_land_ac'].to_numpy()
    _shift = np.abs(_cty_only - _dens)
    _trunc = pd.DataFrame({
        'cells_total': [len(grid)],
        'cells_with_truncated_band': [int(_cut.sum())],
        'pct_truncated': [round(100 * _cut.mean(), 1)],
        'median_density_shift': [round(float(np.median(_shift[_cut])), 2) if _cut.any() else 0.0],
        'p90_density_shift': [round(float(np.percentile(_shift[_cut], 90)), 2) if _cut.any() else 0.0],
        'max_density_shift': [round(float(_shift.max()), 2)]})
    _trunc.to_csv(os.path.join(S_VAL, 'band_truncation.csv'), index=False)
    print(_trunc.to_string(index=False))
    print(f'\n  {int(_cut.sum()):,} cells ({100 * _cut.mean():.0f}%) sit close '
          f'enough to a county line that a\n  county-only band would have '
          f'measured them over less land. Those are the cells\n  the statewide '
          f'surface changes; everything in the interior is identical.')
elif REPORT_MASK is not None:
    print('  skipped: the reporting mask re-cut the grid after Pass A, so the '
          'neighbour\n  lists no longer index it. Run once with '
          'REPORT_WITHIN_PA = False to measure this.')
else:
    print('  surface was loaded from disk, so the neighbour lists are not in '
          'memory -\n  set RUN_SURFACE = True to recompute this check.')


# %% [markdown]
# ## Pass B — the county attribution and the manifest
#
# Cheap, because the surface already exists. Each county gets its slice of the
# grid, the zones that touch it, a folder, and a row in
# `county_manifest.csv`.
#
# **The manifest is the run ledger.** A 120-county run dies partway at some
# point — a bad geometry, a disconnected Drive, a session timeout — and the
# manifest is what makes that resumable instead of restartable. It also carries
# the two things that stop a statewide table being misread: `extent_basis`
# (whether this county was measured against county ∩ priority area or against
# the whole county) and `source_cover` (whether the cached layers reach it at
# all). A county with `status = no_source_data` is not a county without
# hotspots.

# %% D1 - assemble the manifest
rule('PASS B - COUNTY ATTRIBUTION')

RUN_AT = datetime.datetime.now().isoformat(timespec='seconds')

_g = grid.reset_index()
_cell_cty = _g.groupby('FIPS').agg(
    n_cells=('cell', 'size'), land_ac=('land_ac', 'sum'),
    rip_ac=('rip_ac', 'sum'), hot_cells=('hot', 'sum'),
    median_density=('rip_per_sqmi', 'median'),
    p_state_median=('p_state', 'median'))
_hot_cty = _g[_g['hot']].groupby('FIPS').agg(
    hot_land_ac=('land_ac', 'sum'), hot_rip_ac=('rip_ac', 'sum'))
_cell_cty = _cell_cty.join(_hot_cty).fillna({'hot_land_ac': 0.0,
                                             'hot_rip_ac': 0.0})

if len(zones):
    _zc = split.groupby('FIPS').agg(
        n_zones=('zone_id', 'nunique'),
        zone_ac=('zone_ac_in_county', 'sum'),
        zone_rip_ac=('rip_crop_ac_in_county', 'sum'))
    _shared_cty = (split.merge(zones[['n_counties']], left_on='zone_id',
                               right_index=True)
                   .query('n_counties > 1').groupby('FIPS')['zone_id'].nunique())
else:
    _zc = pd.DataFrame(columns=['n_zones', 'zone_ac', 'zone_rip_ac'])
    _shared_cty = pd.Series(dtype=int)

man = counties[['FIPS', 'NAME', 'county_ac', 'in_pa', 'pa_share',
                'extent_basis', 'extent_ac', 'source_cover', 'covered',
                'analysed_ac', 'rip_ac_exact', 'avg_county']].copy()
# A spreadsheet turns '003' into 3, so carry a text id that survives the trip.
# It is also the county folder name, which makes the two easy to line up.
man['county_id'] = (STATE_FIPS + man['FIPS'] + '_' +
                    man['NAME'].map(safe_name))
man = man.merge(_cell_cty, left_on='FIPS', right_index=True, how='left')
man = man.merge(_zc, left_on='FIPS', right_index=True, how='left')
man['n_zones_shared'] = man['FIPS'].map(_shared_cty).fillna(0).astype(int)
man[['n_cells', 'hot_cells', 'n_zones']] = man[
    ['n_cells', 'hot_cells', 'n_zones']].fillna(0).astype(int)
for _c in ('land_ac', 'rip_ac', 'hot_land_ac', 'hot_rip_ac', 'zone_ac',
           'zone_rip_ac'):
    man[_c] = man[_c].fillna(0.0)

# A county the cache never reached has no measurements, and a zero would read
# as one. Blank is the honest value, and it keeps those counties out of every
# ranking rather than tied at the bottom of it.
man.loc[~man['covered'], ['land_ac', 'rip_ac', 'rip_ac_exact', 'analysed_ac',
                          'avg_county', 'hot_land_ac', 'hot_rip_ac',
                          'zone_ac', 'zone_rip_ac']] = np.nan
man.loc[~man['covered'], ['n_cells', 'hot_cells', 'n_zones']] = 0

# Two totals, deliberately both present. `rip_crop_ac` is the exact acreage
# inside the county's analysed extent - the number to quote. `cell_rip_ac` is
# what the cells attributed to this county hold, which is what the hot-cell
# percentages below are a share of. They differ by the lattice overhang at the
# county line, and printing only one of them is how that difference turns into
# an argument later.
man['analysed_sqmi'] = (man['analysed_ac'] / 640).round(1)
man['cell_rip_ac']   = man['rip_ac'].round(0)
man['cell_sqmi']     = (man['land_ac'] / 640).round(1)
man['avg_county']    = man['avg_county'].astype(float).round(2)
man['avg_state']     = round(AVG_STATE, 2)
man['avg_region']    = man['FIPS'].map(
    _g.groupby('FIPS')['avg_region'].median()).astype(float).round(2)
man['reference_used'] = man['FIPS'].map(
    _g.groupby('FIPS')['reference'].median()).astype(float).round(2)
man['threshold_applied'] = np.where(
    man['reference_used'].notna(),
    np.maximum(HOT_MULTIPLE * man['reference_used'].fillna(0), FLOOR),
    np.nan).round(2)
man['floor_ac_per_sqmi'] = round(FLOOR, 2)
man['basis'] = BASIS
man['report_mask'] = REPORT_MASK or 'none'

man['hot_sqmi']          = (man['hot_land_ac'] / 640).round(1)
man['pct_area_hot']      = (100 * man['hot_land_ac'] /
                            man['land_ac'].replace(0, np.nan)).round(1)
man['pct_resource_hot']  = (100 * man['hot_rip_ac'] /
                            man['rip_ac'].replace(0, np.nan)).round(1)
man['concentration']     = (man['pct_resource_hot'] /
                            man['pct_area_hot'].replace(0, np.nan)).round(2)
man['zone_sqmi']         = (man['zone_ac'] / 640).round(1)
man['rip_crop_ac']       = man['rip_ac_exact'].round(0)
man['pct_state_resource'] = (100 * man['rip_ac_exact'] /
                             max(RIP_TOTAL, 1e-9)).round(2)

# Step 1. Status. A county with no zones is a RESULT - it tells the programme
# where not to spend time - so the reason is recorded rather than left as a
# blank row.
def _status(r):
    if not r['covered']:
        return 'no_source_data'
    if r['n_cells'] == 0:
        return 'no_cells'
    if not (r['rip_ac'] > 0):
        return 'no_riparian_cropland'
    if r['n_zones'] > 0:
        return 'zones'
    if r['hot_cells'] > 0:
        return 'hot_cells_but_no_zone'
    if r['median_density'] < FLOOR:
        return 'no_zones_below_floor'
    return 'no_zones_below_ratio'


man['status'] = man.apply(_status, axis=1)


def _rank(s):
    """Rank high to low, with the unrankable - a county the cache never
    reached - sorted last rather than crashing the cast."""
    r = s.rank(ascending=False, method='min')
    return r.fillna(len(s) + 1).astype(int)


man['rank_by_resource'] = _rank(man['rip_ac_exact'])
man['rank_by_density'] = _rank(man['avg_county'])
man['rank_by_zone_resource'] = _rank(man['zone_rip_ac'])
man['has_parcels'] = False
man['parcel_path'] = ''
man['run_at'] = RUN_AT
man['error'] = ''

# Step 2. Which counties already have a parcel file. Parcels are bought per
# county, so this column is the shopping list: it is what Pass C - the
# notebook's own parcel work - can run on today.
_cand = [p for p in glob.glob(os.path.join(RAW_DIR, '**', '*'),
                              recursive=True)
         if p.lower().endswith(PARCEL_EXTENSIONS)]
_cand_base = [os.path.basename(p).lower() for p in _cand]
for _i, _f in man['FIPS'].items():
    _nm = safe_name(man.at[_i, 'NAME']).lower()
    _hits = [p for p, b in zip(_cand, _cand_base) if _nm and _nm in b]
    if _hits:
        man.at[_i, 'has_parcels'] = True
        man.at[_i, 'parcel_path'] = _hits[0]

for _c in ('rip_crop_ac', 'zone_rip_ac', 'hot_rip_ac', 'county_ac',
           'extent_ac', 'analysed_ac'):
    man[_c] = man[_c].round(0)

MAN_COLS = ['county_id', 'FIPS', 'NAME', 'status', 'extent_basis', 'pa_share',
            'source_cover', 'county_ac', 'extent_ac', 'analysed_sqmi',
            'cell_sqmi', 'rip_crop_ac', 'cell_rip_ac', 'pct_state_resource',
            'avg_county', 'avg_region',
            'avg_state', 'basis', 'report_mask', 'reference_used',
            'floor_ac_per_sqmi',
            'threshold_applied', 'n_cells', 'hot_cells', 'hot_sqmi',
            'pct_area_hot', 'pct_resource_hot', 'concentration', 'n_zones',
            'n_zones_shared', 'zone_sqmi', 'zone_rip_ac', 'rank_by_resource',
            'rank_by_density', 'rank_by_zone_resource', 'has_parcels',
            'parcel_path', 'run_at', 'error']
man = man[MAN_COLS].sort_values('rank_by_zone_resource').reset_index(drop=True)
man.to_csv(MANIFEST_CSV, index=False)

print(f'  {len(man)} counties -> {MANIFEST_CSV}')
print('\n  status:')
print(man.groupby('status').size().sort_values(ascending=False)
      .to_frame('counties').to_string())
print(f'\n  top 15 by riparian cropland inside a zone:')
print(man.loc[man['n_zones'] > 0,
              ['FIPS', 'NAME', 'analysed_sqmi', 'rip_crop_ac', 'avg_county',
               'threshold_applied', 'n_zones', 'zone_sqmi', 'zone_rip_ac',
               'pct_resource_hot']]
      .head(15).to_string(index=False))


# %% D2 - one folder per county
# Step 3. County first, because the county is the delivery unit: a biologist
# gets one folder. `cells.gpkg` carries every column the statewide surface has,
# so a county packet can be re-thresholded on its own without going back to the
# state file.
if WRITE_COUNTY_FILES:
    rule('PASS B - COUNTY FOLDERS')
    _made = 0
    for _, _r in man.iterrows():
        _dir = os.path.join(C_DIR, _r['county_id'])
        os.makedirs(os.path.join(_dir, 'output'), exist_ok=True)
        os.makedirs(os.path.join(_dir, 'figures'), exist_ok=True)
        _cells = grid[grid['FIPS'] == _r['FIPS']]
        if len(_cells):
            gpkg_safe(_cells.reset_index()).to_file(
                os.path.join(_dir, 'output', 'cells.gpkg'), layer='cells',
                driver='GPKG')
        if len(zones) and _r['n_zones']:
            _ids = split.loc[split['FIPS'] == _r['FIPS'], 'zone_id'].unique()
            _zz = zones.loc[zones.index.isin(_ids)].reset_index()
            gpkg_safe(_zz).to_file(os.path.join(_dir, 'output', 'zones.gpkg'),
                                   layer='zones', driver='GPKG')
        with open(os.path.join(_dir, 'output', 'summary.json'), 'w') as _f:
            json.dump({k: (None if pd.isna(v) else
                           (v.item() if hasattr(v, 'item') else v))
                       for k, v in _r.items()}, _f, indent=2, default=str)
        _made += 1
    print(f'  {_made} county folders under {C_DIR}')

# Step 4. The statewide grid and the county summary, in one place.
gpkg_safe(grid.reset_index()).to_file(SURFACE_GPKG, layer='grid_scored',
                                      driver='GPKG')
man.to_csv(os.path.join(S_OUT, 'county_summary.csv'), index=False)
print(f'  scored grid -> {SURFACE_GPKG} (layer "grid_scored")')


# %% [markdown]
# ## Figures
#
# Two state maps, two justification charts, one comparison, and one map per
# county. The colour is doing a job in each: **one blue ramp, light to dark,
# for density** (a magnitude), **one orange for hot** (a state, not a series),
# and grey for everything that is context. Nothing is coloured by rank, so a
# county does not change colour when another county is added.
#
# The density classes are the same everywhere — the state map, the comparison
# panels and all 120 county maps — which is the whole point of having a floor:
# the legend means the same thing in Ballard as in Christian.

# %% E0 - the density classes and the drawing helpers
# Every colour, type size and DPI is in fl_config.py section 7. What is derived
# HERE is everything that depends on THIS RUN's data or threshold, which is why
# it cannot live in a config file.
DENS_CMAP = ListedColormap(SEQ)
DENS_CMAP.set_bad(NEUTRAL)

# Step 1. Density classes, fixed once from the statewide distribution so every
# map in the run shares a legend. Breaks are percentiles of the cells that hold
# any riparian cropland at all - classing on the full distribution would spend
# four of six classes on zero.
_pos = grid.loc[grid['rip_per_sqmi'] > 0, 'rip_per_sqmi']
if len(_pos) > 10:
    DENS_BREAKS = [0.0] + [float(np.percentile(_pos, p))
                           for p in (20, 40, 60, 80, 95)] + [float(_pos.max())]
else:
    DENS_BREAKS = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, max(float(grid['rip_per_sqmi'].max()), 21.0)]
DENS_BREAKS = sorted(set(round(b, 2) for b in DENS_BREAKS))
while len(DENS_BREAKS) < len(SEQ) + 1:
    DENS_BREAKS.append(DENS_BREAKS[-1] + 1.0)
DENS_NORM = BoundaryNorm(DENS_BREAKS, DENS_CMAP.N)
# enough decimals that two classes never print the same label
_dp = 0
while _dp < 3 and any(round(DENS_BREAKS[i], _dp) == round(DENS_BREAKS[i + 1], _dp)
                      for i in range(len(SEQ))):
    _dp += 1
DENS_LABELS = [f'{DENS_BREAKS[i]:,.{_dp}f} – {DENS_BREAKS[i + 1]:,.{_dp}f}'
               for i in range(len(SEQ))]
# Step 1b. The RATIO scale used by the county sheets. It is derived from this
# run's own distribution and from the threshold this run actually applied - not
# a fixed 0-4 - so the colour bar says something true about these data. It is
# computed ONCE, statewide, so every county sheet still shares one scale.
# It is ANCHORED ON THE THRESHOLD THIS RUN ACTUALLY APPLIED, not on a number
# carried over from somewhere else: vmax is a fixed multiple of HOT_MULTIPLE,
# so the black threshold line always sits at the same height on the bar and
# raising HOT_MULTIPLE re-scales every map with it. Anchoring on a percentile
# instead was tried and is worse - one heavy-tailed county drags vmax up and
# squashes the range everyone else lives in.
RATIO_VMAX = float(np.round(HOT_MULTIPLE * RATIO_HEADROOM, 1))
_RATIO_SAT = 100 * float((grid['ratio'].to_numpy(float) > RATIO_VMAX).mean())

# Step 1c. THE COUNTY SHEETS ARE DRAWN AGAINST THE CUT, NOT AGAINST THE
# REFERENCE. A ratio against the county or region average is a number whose
# meaning changes from sheet to sheet: where the local reference is very low,
# a cell can sit at 4x it and still hold a third of what the floor demands.
# Lawrence County is the case in point - region average 0.6, floor 7.03 - and
# its sheet came out solid red with zero zones, which reads as a bug and is
# arithmetic.
#
# Dividing by `threshold_applied` - max(HOT_MULTIPLE x reference, FLOOR), the
# cut that cell actually had to clear - fixes the meaning: 1.0 is the line,
# everywhere, on every sheet. The within-county pattern is untouched, because
# for most counties this is a division by a constant.
grid['vs_cut'] = (grid['rip_per_sqmi'] /
                  grid['threshold_applied'].replace(0, np.nan))
CUT_VMAX = 2.0            # 0 to 2x the cut, so the 1.0 line sits mid-bar
_CUT_SAT = 100 * float((grid['vs_cut'].to_numpy(float) > CUT_VMAX).mean())
# how many counties are FLOOR-limited: the floor, not the ratio, is what they
# had to clear. Their sheets carry a line saying so.
_FLOOR_LIMITED = set(
    grid.loc[HOT_MULTIPLE * grid['reference'].to_numpy(float) < FLOOR,
             'FIPS'].dropna().unique())
print(f'county sheets: coloured against the cut each county had to clear '
      f'(0 to {CUT_VMAX:g}x,\n                1.0 = the cut); '
      f'{_CUT_SAT:.1f}% of cells sit above the top of that bar.')
print(f'               {len(_FLOOR_LIMITED)} counties are FLOOR-limited '
      f'({HOT_MULTIPLE:g}x their reference is below the\n               '
      f'{fmt(FLOOR)} floor), and their sheets say so.')
print(f'ratio scale   : 0 to {RATIO_VMAX:g}x the {BASIS} average '
      f'({RATIO_HEADROOM:g}x the {HOT_MULTIPLE:g}x threshold, which is marked '
      f'on the bar);\n                {_RATIO_SAT:.1f}% of cells sit above it '
      f'and are drawn at the top colour')
print('density classes (riparian crop ac per sq mi, 3-mile neighbourhood):')
for _l, _c in zip(DENS_LABELS, SEQ):
    print(f'  {_l:>18}   {_c}')


def cell_raster(sub, col):
    """Cells are equal squares on a lattice, so draw them as an image rather
    than as tens of thousands of polygons. Same geometry, a fraction of the
    time, and no hairline seams between neighbours."""
    gb = np.array([g.bounds for g in sub.geometry])
    gx0, gy0 = gb[:, 0].min(), gb[:, 1].min()
    ii = np.round((gb[:, 0] - gx0) / CELL_M).astype(int)
    jj = np.round((gb[:, 1] - gy0) / CELL_M).astype(int)
    arr = np.full((jj.max() + 1, ii.max() + 1), np.nan)
    arr[jj, ii] = np.asarray(sub[col], dtype=float)
    ext = [gx0, gx0 + (ii.max() + 1) * CELL_M,
           gy0, gy0 + (jj.max() + 1) * CELL_M]
    return arr, ext


def fmt(v, width=1):
    """Enough decimals to be true. A floor of 1.12 must not print as '1'."""
    return f'{v:,.0f}' if abs(v) >= 10 else f'{v:,.{width + 1}f}'.rstrip('0').rstrip('.')


def map_figsize(width_in, bounds, pad_in=1.6, n_panels=1):
    """Height that matches the ground, so an equal-aspect map does not ship
    with two inches of white above and below it."""
    x0, y0, x1, y1 = bounds
    panel_w = width_in / n_panels
    return (width_in, max(3.0, panel_w * (y1 - y0) / max(x1 - x0, 1e-9) + pad_in))


def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect('equal')
    return ax


def legend_headroom(ax, frac=0.22):
    """Make blank sky above the map for the legend to sit in.

    A legend placed inside a map always covers something; on a wide state the
    upper left is the Ohio River counties, which is exactly the ground with
    the most riparian cropland in the state. Extending the y limit instead
    creates space that is genuinely empty - the axes keep an equal aspect, so
    the map is not stretched, the box simply gets shorter inside a taller
    frame. That is the only way to put a legend top-left and be certain it
    covers no data."""
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + frac * (y1 - y0))
    return ax


def frame(lg):
    lg.get_frame().set_facecolor(SURFACE)
    lg.get_frame().set_edgecolor(MUTED)
    lg.get_frame().set_linewidth(0.6)
    lg.get_title().set_color(INK2)
    lg.get_title().set_fontweight('bold')
    return lg


def density_legend(ax, loc='upper left', title='Riparian cropland\n(ac / sq mi)',
                   extra=(), outside=False, ncol=1):
    h = [Patch(facecolor=c, edgecolor='none', label=l)
         for c, l in zip(SEQ, DENS_LABELS)] + list(extra)
    kw = (dict(loc='upper left', bbox_to_anchor=(1.01, 1.0)) if outside
          else dict(loc=loc))
    return frame(ax.legend(handles=h, frameon=True, fontsize=FS['legend'],
                           title=title, title_fontsize=FS['legend_title'],
                           borderpad=0.8, labelspacing=0.5, handlelength=1.5,
                           ncol=ncol, **kw))


# %% E1 - the statewide density surface
rule('FIGURES')

_arr, _ext = cell_raster(grid, 'rip_per_sqmi')
fig, ax = plt.subplots(figsize=map_figsize(14, ANALYSIS_GEOM.bounds, pad_in=2.6))
bare(ax)
ax.imshow(np.ma.masked_invalid(_arr), origin='lower', extent=_ext,
          cmap=DENS_CMAP, norm=DENS_NORM, interpolation='nearest')
counties.boundary.plot(ax=ax, color=MUTED, linewidth=0.5, zorder=3)
ax.set_title(f'Figure 1.  Riparian Cropland Density, {STATE_TITLE}',
             fontsize=FS['fig_title'], loc='left', pad=12)
ax.set_xlabel(f'{CELL_MI:g}-mile cells, {BAND_MI:g}-mile neighbourhood · '
              f'cropland within {RIP_BUFFER_FT:.0f} ft of water · '
              f'state average {AVG_STATE:.1f} ac/sq mi',
              color=INK2, fontsize=FS['caption'])
_UNANALYSED = ANALYSIS_AC < 0.999 * STATE_AC
legend_headroom(ax, 0.26)
density_legend(ax, loc='upper left',
               extra=([Patch(facecolor=SURFACE, edgecolor=MUTED, lw=0.6,
                             label='Not analysed (no source data)')]
                      if _UNANALYSED else []))
ax.text(0.995, 0.01,
        f'{ANALYSIS_AC / 640:,.0f} sq mi analysed · '
        f'{grid["rip_ac"].sum():,.0f} ac riparian cropland',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=FS['caption'],
        color=INK2)
F1 = os.path.join(S_FIG, 'F1_state_density.png')
fig.savefig(F1, dpi=FIG_DPI, bbox_inches='tight')
show_close(fig)
print(f'  {os.path.basename(F1)}')


# %% E2 - the statewide hotspots
_hotarr = np.full_like(_arr, np.nan)
_h_sub = grid['hot'].to_numpy()
_ha, _he = cell_raster(grid.assign(_h=_h_sub.astype(float)), '_h')
_ha = np.where(_ha > 0, 1.0, np.nan)

fig, ax = plt.subplots(figsize=map_figsize(14, ANALYSIS_GEOM.bounds, pad_in=2.6))
bare(ax)
ax.imshow(np.ma.masked_invalid(_arr), origin='lower', extent=_ext,
          cmap=ListedColormap([NEUTRAL] * len(SEQ)), norm=DENS_NORM,
          interpolation='nearest')
ax.imshow(np.ma.masked_invalid(_ha), origin='lower', extent=_he,
          cmap=ListedColormap([HOT]), vmin=0, vmax=1, interpolation='nearest',
          zorder=2)
counties.boundary.plot(ax=ax, color=MUTED, linewidth=0.5, zorder=3)
if len(zones):
    zones.boundary.plot(ax=ax, color=ZONE_INK, linewidth=0.7, zorder=4)
    _top = man[man['n_zones'] > 0].head(12)
    for _, _r in _top.iterrows():
        _cg2 = counties.loc[counties['FIPS'] == _r['FIPS'], 'geometry'].iloc[0]
        _p = _cg2.representative_point()
        ax.annotate(_r['NAME'], (_p.x, _p.y), fontsize=FS['annot'], color=INK,
                    ha='center', va='center', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.2', fc=SURFACE,
                              ec=MUTED, lw=0.4, alpha=0.92))
ax.set_title(f'Figure 2.  High-Density Areas — neighbourhood density at least '
             f'{HOT_MULTIPLE:g}× the {BASIS} average and at least '
             f'{fmt(FLOOR)} ac / sq mi',
             fontsize=FS['fig_title'], loc='left', pad=12)
ax.set_xlabel(f'{int(grid["hot"].sum()):,} qualifying cells · {len(zones):,} zones · '
              f'{man["zone_sqmi"].sum():,.0f} sq mi · '
              f'{100 * grid.loc[grid["hot"], "rip_ac"].sum() / max(grid["rip_ac"].sum(), 1e-9):.0f}% '
              f'of the state\'s riparian cropland',
              color=INK2, fontsize=FS['caption'])
legend_headroom(ax, 0.26)
frame(ax.legend(handles=[
    Patch(facecolor=HOT, edgecolor='none', label='High-density cell'),
    Line2D([0], [0], color=ZONE_INK, lw=1.4, label='Zone boundary'),
    Patch(facecolor=NEUTRAL, edgecolor='none', label='Analysed, below threshold')]
    + ([Patch(facecolor=SURFACE, edgecolor=MUTED, lw=0.6,
              label='Not analysed (no source data)')] if _UNANALYSED else []),
    loc='upper left', frameon=True, fontsize=FS['legend'], borderpad=0.8,
    labelspacing=0.5, title='Classification',
    title_fontsize=FS['legend_title']))
F2 = os.path.join(S_FIG, 'F2_state_hotspots.png')
fig.savefig(F2, dpi=FIG_DPI, bbox_inches='tight')
show_close(fig)
print(f'  {os.path.basename(F2)}')


# %% E3 - the floor sweep, i.e. where the absolute cut came from
# TERMINOLOGY. On the validation figures the class is called a HIGH-DENSITY
# AREA rather than a "hotspot" or "ground called hot". The column in every
# table is still `hot`, because renaming data to match a caption is how a
# figure and its source stop agreeing - but a figure that goes in front of a
# programme reviewer should say what was measured, and what was measured is
# riparian-cropland density above a stated threshold.
fig, axes = plt.subplots(1, 2, figsize=(15, 6.0), layout='constrained')
_s = floor_sweep.sort_values('floor_ac_per_sqmi')
_x = _s['floor_ac_per_sqmi'].to_numpy()

ax = axes[0]
# two series on one axis: the ramp's darkest step and ink. Two steps of the
# same ramp would be the wrong choice here - a ramp encodes magnitude, and
# these two lines are different measures, not more or less of one thing.
ax.plot(_x, _s['pct_resource_hot'], color=SEQ[5], lw=2.4,
        marker='o', ms=6, label='Riparian cropland captured')
ax.plot(_x, _s['pct_area_hot'], color=INK2, lw=2.4, marker='s', ms=5.5,
        ls='--', label='Land area designated high-density')
ax.set_xlabel('Absolute density floor (riparian cropland ac / sq mi)')
ax.set_ylabel('Percent of the analysed state')
# The floor markers are annotated ABOVE the axes (see below), so the panel
# titles are padded clear of them. Fixed pad rather than a measured one: the
# two panels then share a baseline whatever the data does.
_E3_PAD = 40 if (KNEE_OK and abs(KNEE_FLOOR - FLOOR) > 0.02 * max(FLOOR, 1.0)) \
    else 26
ax.set_title('A.  Resource Captured Against Area Designated',
             fontsize=FS['panel_title'], loc='left', pad=_E3_PAD)
ax.grid(axis='y', color=NEUTRAL, lw=1.0)
ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
_E3_H, _E3_L = ax.get_legend_handles_labels()

ax = axes[1]
ax.plot(_x, _s['counties_with_zones'], color=SEQ[5], lw=2.4, marker='o', ms=6)
ax.set_xlabel('Absolute density floor (riparian cropland ac / sq mi)')
ax.set_ylabel('Counties retaining a high-density cell')
ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
ax.set_title('B.  Counties Retaining at Least One High-Density Cell',
             fontsize=FS['panel_title'], loc='left', pad=_E3_PAD)
ax.grid(axis='y', color=NEUTRAL, lw=1.0)
ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
# The floor marker is anchored to the BOTTOM of each panel and carries an
# opaque background. Anchored to the top it landed on whichever curve happened
# to be high there, which on this data was the one the reader is looking at.
for _a in axes:
    # The plateau BEHIND the floor line, so the reader sees the range the rule
    # regards as equivalent before the single number drawn inside it. Light
    # grey rather than a tint of the ramp: the ramp means magnitude everywhere
    # else in this run, and this is an interval, not a quantity.
    if (FLOOR_MODE == 'knee' and KNEE_OK
            and all(np.isfinite(_v) for _v in KNEE_PLATEAU)
            and KNEE_PLATEAU[1] > KNEE_PLATEAU[0]):
        _a.axvspan(KNEE_PLATEAU[0], KNEE_PLATEAU[1], color=NEUTRAL,
                   ec='none', zorder=0)
    _a.axvline(FLOOR, color=INK2, lw=1.4, ls=':')
    _a.annotate(f'Adopted floor  {fmt(FLOOR)} ac / sq mi',
                xy=(FLOOR, 1.0), xycoords=('data', 'axes fraction'),
                xytext=(0, 5), textcoords='offset points',
                fontsize=FS['annot'], color=INK2, ha='center', va='bottom',
                annotation_clip=False)
    # The knee the data itself picks, drawn whether or not it was adopted -
    # a reviewer should be able to see the gap between the recommendation and
    # the decision without opening a CSV. Only drawn when it is somewhere
    # else, so the usual case is one clean line rather than two on top of
    # each other.
    if KNEE_OK and abs(KNEE_FLOOR - FLOOR) > 0.02 * max(FLOOR, 1.0):
        _a.axvline(KNEE_FLOOR, color=SEQ[5], lw=1.4, ls='--')
        _a.annotate(f'Knee ({FLOOR_KNEE_RULE})  {fmt(KNEE_FLOOR)} ac / sq mi',
                    xy=(KNEE_FLOOR, 1.0), xycoords=('data', 'axes fraction'),
                    xytext=(0, 19), textcoords='offset points',
                    fontsize=FS['annot'], color=SEQ[5], ha='center',
                    va='bottom', annotation_clip=False)
# The legend goes OUTSIDE both panels. Inside panel A it sat on the descending
# capture curve at every floor value that matters; there is no "safe" corner
# on a chart whose whole point is two lines crossing.
if (FLOOR_MODE == 'knee' and KNEE_OK
        and all(np.isfinite(_v) for _v in KNEE_PLATEAU)
        and KNEE_PLATEAU[1] > KNEE_PLATEAU[0]):
    _E3_H = list(_E3_H) + [Patch(facecolor=NEUTRAL, edgecolor='none')]
    _E3_L = list(_E3_L) + [
        f'{FLOOR_KNEE_RULE} plateau  {fmt(KNEE_PLATEAU[0])}\u2013'
        f'{fmt(KNEE_PLATEAU[1])} ac / sq mi']
frame(fig.legend(_E3_H, _E3_L, loc='outside lower center',
                 ncol=min(3, len(_E3_H)),
                 frameon=True, fontsize=FS['legend'], borderpad=0.8,
                 columnspacing=2.4, handlelength=2.6))
fig.suptitle('Figure 3.  Derivation of the Absolute Density Floor',
             fontsize=FS['fig_title'], fontweight='bold', color=INK, x=0.012,
             ha='left')
F3 = os.path.join(S_FIG, 'F3_floor_sweep.png')
fig.savefig(F3, dpi=FIG_DPI, bbox_inches='tight')
show_close(fig)
print(f'  {os.path.basename(F3)}')


# %% E3b - Figure 3b: where each rule would put the floor
# Figure 3 shows ONE floor. This shows the comparison behind it: every rule's
# pick on the capture curve, and the curve each argued from.
#
# Two panels, not one: the capture curve is percent against percent, the net
# benefit is a normalised index. A dual axis would let them share a panel.
# That is exactly why they do not.
#
# No commentary on the figure. The verdict - whether the rules agree, and what
# each one keeps - is in the run log and in floor_knee_rules.csv, where it can
# be argued with. A figure that states a conclusion in its own subtitle is
# harder to disagree with than it should be.
_SHOW = tuple(r for r in KNEE_RULES
              if r in tuple(globals().get('KNEE_RULES_SHOWN', KNEE_RULES)))

if KNEE_OK and len(knee_rules):
    _kc = knee_curve[np.isfinite(knee_curve['concentration'].to_numpy(float))]
    _kc = _kc[_kc['hot_cells'] > 0].reset_index(drop=True)
    _kr_shown = knee_rules[knee_rules['rule'].isin(_SHOW)]

    def _row_at(pctl):
        """the swept row nearest a rule's pick"""
        return _kc.iloc[int((_kc['pctl'] - pctl).abs().idxmin())]

    # rules that chose the same floor share one label rather than stacking two
    # on top of each other
    _groups = {}
    for _, _r in _kr_shown.iterrows():
        _groups.setdefault(round(float(_r['floor_ac_per_sqmi']), 2),
                           []).append(_r['rule'])

    # ONE PANEL PER RULE, not one axis with three lines on it. Each rule
    # normalises its own gain and cost over its own range, so the curves share
    # no unit and their heights are not comparable - on a single axis the
    # tallest line reads as the best rule, which is exactly the thing it does
    # not mean. Separate panels make that misreading impossible.
    # 'marginal' has no *_net column - it argues from an exchange rate in
    # ac/sq mi against a threshold in the same units, so its panel carries its
    # own y-axis rather than borrowing the normalised one. Separate panels
    # already make heights non-comparable, which is what lets it sit here.
    _net = [r for r in _SHOW
            if f'{r}_net' in _kc.columns
            or (r == 'marginal' and 'marginal_ac_per_sqmi' in _kc.columns)]
    fig, _axs = plt.subplots(
        1, 1 + len(_net), figsize=(6.4 + 4.1 * len(_net), 5.9),
        width_ratios=[1.55] + [1.0] * len(_net), layout='constrained')
    _axs = np.atleast_1d(_axs)      # one panel comes back bare, not in an array
    axa, _axn = _axs[0], list(_axs[1:])

    # --- A. the capture curve --------------------------------------------
    ax = axa
    ax.plot(_kc['pct_area_hot'], _kc['pct_resource_hot'], color=MUTED,
            lw=2.0, marker='o', ms=3.6, mec='none', zorder=2)
    _x0, _y0 = _kc['pct_area_hot'].iloc[-1], _kc['pct_resource_hot'].iloc[-1]
    _x1, _y1 = _kc['pct_area_hot'].iloc[0], _kc['pct_resource_hot'].iloc[0]
    ax.plot([_x0, _x1], [_y0, _y1], color=MUTED, lw=1.2, ls=':', zorder=1)

    # The picks often sit within a whisker of each other on this curve - that
    # IS the finding when they agree - so the labels are stacked down the
    # panel with a leader line to each marker, rather than written next to the
    # points where they become one smear.
    _pts = []
    for _fl, _rs in sorted(_groups.items()):
        _rr = _row_at(_kr_shown.loc[_kr_shown['rule'] == _rs[0],
                                    'pctl'].iloc[0])
        _c = RULE_COLORS[_rs[0]]
        ax.plot([_rr['pct_area_hot']], [_rr['pct_resource_hot']], 'o', ms=11,
                color=_c, mec=SURFACE, mew=2.0, zorder=5)
        _pts.append((_rs, _fl, _c,
                     (_rr['pct_area_hot'], _rr['pct_resource_hot'])))
    _ly = 0.34 + 0.085 * (len(_pts) - 1)
    for _rs, _fl, _c, _xy in _pts:
        # the join is done OUTSIDE the f-string: before Python 3.12 an
        # f-string expression may not contain a backslash, and Colab is not
        # guaranteed to be 3.12
        _lab = ' \u00b7 '.join(_rs)
        _pr = knee_rules.loc[knee_rules['rule'] == _rs[0]].iloc[0]
        _band = ('' if not (np.isfinite(float(_pr['plateau_lo']))
                            and np.isfinite(float(_pr['plateau_hi'])))
                 else f'  [{fmt(float(_pr["plateau_lo"]))}\u2013'
                      f'{fmt(float(_pr["plateau_hi"]))}]')
        ax.annotate(f'{_lab}   {fmt(_fl)} ac / sq mi{_band}', xy=_xy,
                    xycoords='data', xytext=(0.05, _ly),
                    textcoords=ax.transAxes, fontsize=FS['annot'], color=_c,
                    ha='left', va='center', zorder=6,
                    arrowprops=dict(arrowstyle='-', color=_c, lw=0.9,
                                    shrinkA=2, shrinkB=9, alpha=0.75))
        _ly -= 0.085
    _fr = _kc.iloc[int((_kc['floor_ac_per_sqmi'] - FLOOR).abs().idxmin())]
    ax.plot([_fr['pct_area_hot']], [_fr['pct_resource_hot']], 'o', ms=17,
            mfc='none', mec=INK, mew=1.6, zorder=4)
    ax.set_xlabel('Land area designated high-density (% of the analysed state)')
    ax.set_ylabel('Riparian cropland captured (%)')
    ax.set_title('A.  Capture Curve', fontsize=FS['panel_title'], loc='left')

    # --- B, C, D. one rule each ------------------------------------------
    _first_net = next((i for i, r in enumerate(_net)
                       if f'{r}_net' in _kc.columns), None)
    for _k, (_r, ax) in enumerate(zip(_net, _axn)):
        _prow = knee_rules.loc[knee_rules['rule'] == _r].iloc[0]
        _plo, _phi = float(_prow['plateau_lo']), float(_prow['plateau_hi'])
        if np.isfinite(_plo) and np.isfinite(_phi) and _phi > _plo:
            ax.axvspan(_plo, _phi, color=NEUTRAL, ec='none', zorder=0)
        _pk = _row_at(_prow['pctl'])
        if _r == 'marginal':
            # the raw exchange rate thin behind the smoothed one, so the
            # reader can see what the rolling median removed
            ax.plot(_kc['floor_ac_per_sqmi'],
                    _kc['marginal_raw_ac_per_sqmi'].to_numpy(float),
                    color=RULE_COLORS[_r], lw=0.9, alpha=0.40, zorder=2)
            _yv = _kc['marginal_ac_per_sqmi'].to_numpy(float)
            ax.plot(_kc['floor_ac_per_sqmi'], _yv, color=RULE_COLORS[_r],
                    lw=2.2, zorder=3)
            ax.axhline(AVG_STATE, color=INK2, lw=1.2, ls='--', zorder=1)
            _yp = float(_kc.loc[_kc.index[int((_kc['pctl'] -
                        _prow['pctl']).abs().argmin())],
                        'marginal_ac_per_sqmi'])
            ax.plot([_pk['floor_ac_per_sqmi']], [_yp], 'o', ms=10,
                    color=RULE_COLORS[_r], mec=SURFACE, mew=2.0, zorder=5)
            ax.set_ylabel('Ground the next step drops (ac / sq mi)')
        else:
            _col = f'{_r}_net'
            ax.plot(_kc['floor_ac_per_sqmi'], _kc[_col].to_numpy(float),
                    color=RULE_COLORS[_r], lw=2.2, zorder=3)
            ax.plot([_pk['floor_ac_per_sqmi']], [_pk[_col]], 'o', ms=10,
                    color=RULE_COLORS[_r], mec=SURFACE, mew=2.0, zorder=5)
            ax.axhline(0, color=MUTED, lw=1.0, zorder=1)
            if _k == _first_net:
                ax.set_ylabel('Net benefit  (gain \u2212 cost, normalised '
                              'within each rule)')
        ax.axvline(FLOOR, color=INK, lw=1.4, ls=':', zorder=4)
        ax.set_xlabel('Absolute density floor (ac / sq mi)')
        _ttl = f'{chr(66 + _k)}.  {_r}'
        if np.isfinite(_plo) and np.isfinite(_phi):
            _ttl += f'   [{fmt(_plo)}\u2013{fmt(_phi)}]'
        ax.set_title(_ttl, fontsize=FS['panel_title'],
                     loc='left', color=RULE_COLORS[_r])

    for _a in [axa] + _axn:
        _a.grid(axis='y', color=NEUTRAL, lw=1.0)
        _a.set_axisbelow(True)
        for _sp in ('top', 'right'):
            _a.spines[_sp].set_visible(False)

    def _rule_label(r):
        _rw = knee_rules.loc[knee_rules['rule'] == r].iloc[0]
        _t = f'{r}  \u2192  {fmt(float(_rw["floor_ac_per_sqmi"]))} ac / sq mi'
        if (np.isfinite(float(_rw['plateau_lo']))
                and np.isfinite(float(_rw['plateau_hi']))):
            _t += (f'  [{fmt(float(_rw["plateau_lo"]))}\u2013'
                   f'{fmt(float(_rw["plateau_hi"]))}]')
        if r == FLOOR_KNEE_RULE and FLOOR_MODE == 'knee':
            _t += '   [adopted]'
        return _t

    _H = [Line2D([], [], color=RULE_COLORS[r], lw=2.6, marker='o', ms=9,
                 mec=SURFACE, mew=1.6, label=_rule_label(r))
          for r in _SHOW]
    _H += [Line2D([], [], color=INK, lw=1.4, ls=':', marker='o', ms=11,
                  mfc='none', mec=INK, mew=1.6,
                  label=f'Floor used: {fmt(FLOOR)} ac / sq mi'),
           Patch(facecolor=NEUTRAL, edgecolor='none',
                 label=f'[plateau]: within {FLOOR_KNEE_PLATEAU_TOL:.0%} of '
                       f"that rule's best"),
           Line2D([], [], color=MUTED, lw=2.0, marker='o', ms=5, mec='none',
                  label='Sweep'),
           Line2D([], [], color=MUTED, lw=1.2, ls=':', label='Chord')]
    if 'marginal' in _net:
        _H += [Line2D([], [], color=INK2, lw=1.2, ls='--',
                      label=f'State average  {fmt(AVG_STATE)} ac / sq mi')]
    frame(fig.legend(handles=_H, loc='outside lower center',
                     ncol=min(3, len(_H)), frameon=True, fontsize=FS['legend'],
                     borderpad=0.8, columnspacing=2.4, handlelength=2.6,
                     title='Where each rule would put the floor',
                     title_fontsize=FS['legend_title']))
    fig.suptitle('Figure 3b.  Where Each Rule Would Put the Absolute Floor',
                 fontsize=FS['fig_title'], fontweight='bold', color=INK,
                 x=0.012, ha='left')
    F3B = os.path.join(S_FIG, 'F3b_floor_rules.png')
    fig.savefig(F3B, dpi=FIG_DPI, bbox_inches='tight')
    show_close(fig)
    print(f'  {os.path.basename(F3B)}')
elif FLOOR_MODE != 'off':
    print('  F3b skipped - the sweep has no knee to compare rules on')

# %% E4 - the county ranking
_rk = man[man['rip_crop_ac'] > 0].sort_values('zone_rip_ac',
                                              ascending=False).head(25)
if len(_rk):
    _rk = _rk.iloc[::-1]
    _y = np.arange(len(_rk))
    fig, axes = plt.subplots(1, 2, figsize=(15, max(6.0, 0.46 * len(_rk))),
                             sharey=True, layout='constrained')
    ax = axes[0]
    ax.barh(_y, _rk['zone_rip_ac'], height=0.66, color=SEQ[4],
            edgecolor='none')
    ax.set_yticks(_y)
    ax.set_yticklabels(_rk['NAME'], fontsize=FS['tick'], color=INK)
    ax.set_xlabel('Riparian cropland within a high-density area (acres)')
    ax.set_title('A.  Absolute Resource Within High-Density Areas',
                 fontsize=FS['panel_title'], loc='left')
    for _i, _v in zip(_y, _rk['zone_rip_ac']):
        ax.annotate(f'{_v:,.0f}', (_v, _i), xytext=(5, 0),
                    textcoords='offset points', va='center',
                    fontsize=FS['annot'], color=INK2)

    ax = axes[1]
    ax.barh(_y, _rk['pct_resource_hot'].fillna(0), height=0.66, color=SEQ[4],
            edgecolor='none')
    ax.set_xlabel("Percent of the county's riparian cropland so classified")
    ax.set_title('B.  Degree of Concentration Within the County',
                 fontsize=FS['panel_title'], loc='left')
    for _i, _v in zip(_y, _rk['pct_resource_hot'].fillna(0)):
        ax.annotate(f'{_v:.0f}%', (_v, _i), xytext=(5, 0),
                    textcoords='offset points', va='center',
                    fontsize=FS['annot'], color=INK2)
    for ax in axes:
        ax.grid(axis='x', color=NEUTRAL, lw=1.0)
        ax.set_axisbelow(True)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        ax.margins(x=0.14)
    fig.suptitle(f'Figure 4.  {STATE_NAME} {COUNTY_WORDS.title()} Ranked by '
                 f'Riparian Cropland in '
                 f'High-Density Areas\n'
                 f'Reference basis: {BASIS} · absolute floor '
                 f'{fmt(FLOOR)} ac / sq mi · ratio threshold {HOT_MULTIPLE:g}×',
                 fontsize=FS['fig_title'], fontweight='bold', color=INK,
                 x=0.012, ha='left')
    F4 = os.path.join(S_FIG, 'F4_county_ranking.png')
    fig.savefig(F4, dpi=FIG_DPI, bbox_inches='tight')
    show_close(fig)
    print(f'  {os.path.basename(F4)}')


# %% E5 - the same surface, three references
_BASIS_TITLE = {'county': 'A.  County Reference',
                'state': 'B.  Statewide Reference',
                'region': 'C.  Regional Reference (HUC8)'}
fig, axes = plt.subplots(1, 3, figsize=map_figsize(
    17, ANALYSIS_GEOM.bounds, pad_in=2.6, n_panels=3), layout='constrained')
for ax, _b in zip(axes, ('county', 'state', 'region')):
    bare(ax)
    _r = grid[f'avg_{_b}'].to_numpy(float)
    _hb = ((grid['rip_per_sqmi'].to_numpy() >= HOT_MULTIPLE * _r) &
           (grid['own_per_sqmi'].to_numpy() >= _r) &
           (grid['rip_per_sqmi'].to_numpy() >= FLOOR))
    _a2, _e2 = cell_raster(grid.assign(_h=_hb.astype(float)), '_h')
    ax.imshow(np.ma.masked_invalid(_arr), origin='lower', extent=_ext,
              cmap=ListedColormap([NEUTRAL] * len(SEQ)), norm=DENS_NORM,
              interpolation='nearest')
    ax.imshow(np.ma.masked_invalid(np.where(_a2 > 0, 1.0, np.nan)),
              origin='lower', extent=_e2, cmap=ListedColormap([HOT]),
              vmin=0, vmax=1, interpolation='nearest', zorder=2)
    counties.boundary.plot(ax=ax, color=MUTED, linewidth=0.4, zorder=3)
    _row = basis_cmp[(basis_cmp['basis'] == _b) &
                     (basis_cmp['floor'] != 'no floor')].iloc[0]
    ax.set_title(_BASIS_TITLE[_b] + ('   [adopted]' if _b == BASIS else ''),
                 fontsize=FS['panel_title'], loc='left')
    ax.set_xlabel(f'{_row["hot_sqmi"]:,.0f} sq mi designated · '
                  f'{_row["pct_resource_hot"]:.0f}% of the resource · '
                  f'{_row["concentration"]:.1f}× concentration\n'
                  f'{int(_row["counties_with_hot"])} counties with '
                  f'high-density ground'
                  + ('' if _row['floor_binds']
                     else '\n(floor inert on this basis — see C6)'),
                  fontsize=FS['caption'], color=INK2)
fig.suptitle(f'Figure 5.  Basis area reference comparison — reference area '
             f'decision visualization\n'
             f'All areas use same {fmt(FLOOR)} ac / sq mi absolute floor',
             fontsize=FS['fig_title'], fontweight='bold', color=INK,
             x=0.012, ha='left')
F5 = os.path.join(S_FIG, 'F5_basis_comparison.png')
fig.savefig(F5, dpi=FIG_DPI, bbox_inches='tight')
show_close(fig)
print(f'  {os.path.basename(F5)}')


# %% E6 - one sheet per county: what was measured, and where it concentrates
# Step 2. The county packet figure, in the two-panel form the notebook makes
# for a single county:
#
#   1. WHAT WAS MEASURED - the cropland, the water it sits next to, and the
#      overlap between them. This is the panel a biologist argues with.
#   2. WHERE IT CONCENTRATES - the density surface as a ratio against the
#      chosen reference, with the zones on top.
#
# Panel 2 keeps the one thing the single-county version cannot have: zones are
# drawn WHOLE, and the part outside the county is hatched rather than cut off,
# so a zone continuing into the next county is visible and its number matches
# zone_county_split.csv.
#
# These two panels use the NOTEBOOK's colours, not the statewide palette. This
# is the sheet that goes in a packet beside the Christian figures people have
# already been shown; the state-level maps (F1, F2, F5) are a different
# audience and keep their own scheme.
def _draw(gdf, bbox, ax, **kw):
    """Slice one layer to the panel, generalise it, draw it. Returns the row
    count so a legend entry is only added for a layer that is actually there."""
    if gdf is None or not len(gdf):
        return 0
    g = gdf.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
    if not len(g):
        return 0
    g = g[['geometry']].copy()
    g['geometry'] = g.geometry.simplify(DRAW_SIMPLIFY_M)
    g = g[~g.geometry.is_empty & g.geometry.notna()]
    if len(g):
        g.plot(ax=ax, **kw)
    return len(g)


def county_figure(fips, name, out_png):
    sub = grid[grid['FIPS'] == fips]
    if not len(sub):
        return False
    r = man.loc[man['FIPS'] == fips].iloc[0]
    cext = counties.loc[counties['FIPS'] == fips, 'extent'].iloc[0]
    cg = counties.loc[counties['FIPS'] == fips, 'geometry'].iloc[0]
    xmin, ymin, xmax, ymax = cext.bounds
    pad = 0.03 * max(xmax - xmin, ymax - ymin)
    bb = (xmin - pad, ymin - pad, xmax + pad, ymax + pad)

    fig, axes = plt.subplots(1, 2, figsize=map_figsize(
        18, bb, pad_in=2.8, n_panels=2), layout='constrained')

    # --- PANEL 1: cropland, water, and the overlap ---------------------------
    ax = axes[0]
    bare(ax)
    gpd.GeoSeries([cext], crs=WORKING_CRS).plot(
        ax=ax, facecolor='#f7f7f5', edgecolor='#666666', linewidth=0.8)
    leg = []
    if _draw(csb, bb, ax, facecolor=CROP_C, edgecolor='none', alpha=0.65,
             zorder=2):
        leg.append(Patch(facecolor=CROP_C, label='cropland (CSB fields)'))
    if _draw(wetlands, bb, ax, facecolor=WET_C, edgecolor='none', alpha=0.55,
             zorder=3):
        leg.append(Patch(facecolor=WET_C, alpha=0.55, label='NWI wetlands'))
    if _draw(nhd_wb, bb, ax, facecolor=WB_C, edgecolor='none', zorder=4):
        leg.append(Patch(facecolor=WB_C, label='NHD waterbodies'))
    if _draw(nhd_area, bb, ax, facecolor=RIV_C, edgecolor='none', zorder=4):
        leg.append(Patch(facecolor=RIV_C, label='NHD rivers (area)'))
    if _draw(nhd_flow, bb, ax, color=FLOW_C, linewidth=0.35, zorder=5):
        leg.append(Line2D([], [], color=FLOW_C, lw=1.2, label='NHD flowlines'))
    if _draw(rip_crop, bb, ax, facecolor=OVER_C, edgecolor='none', zorder=6):
        leg.append(Patch(facecolor=OVER_C,
                         label=f'OVERLAP ({r["rip_crop_ac"]:,.0f} ac)'))
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3])
    legend_headroom(ax, 0.18)
    frame(ax.legend(handles=leg, fontsize=FS['legend'], loc='upper left',
                    framealpha=0.96, borderpad=0.8, labelspacing=0.5,
                    title='Measured layers',
                    title_fontsize=FS['legend_title']))
    ax.set_title(f'A.  Cropland Within {RIP_BUFFER_FT:.0f} ft of Water',
                 fontsize=FS['panel_title'], loc='left')
    ax.set_xlabel(f'{r["analysed_sqmi"]:,.0f} sq mi analysed '
                  f'({r["extent_basis"]}) · {r["rip_crop_ac"]:,.0f} ac of '
                  f'cropland on water', fontsize=FS['caption'], color=INK2)

    # --- PANEL 2: the density surface and the zones --------------------------
    ax = axes[1]
    bare(ax)
    # `vs_cut` is per CELL, so a context cell from the next county is drawn
    # against ITS own threshold rather than this county's - which is what
    # makes a zone that crosses the line look continuous.
    ctx = grid.cx[bb[0]:bb[2], bb[1]:bb[3]]
    _a, _e = cell_raster(ctx, 'vs_cut')
    ax.imshow(np.ma.masked_invalid(_a), origin='lower', extent=_e,
              cmap=RAMP_NAME, vmin=0, vmax=CUT_VMAX,
              interpolation='nearest', alpha=0.35, zorder=1)
    _a2, _e2 = cell_raster(sub, 'vs_cut')
    im = ax.imshow(np.ma.masked_invalid(_a2), origin='lower', extent=_e2,
                   cmap=RAMP_NAME, vmin=0, vmax=CUT_VMAX,
                   interpolation='nearest', zorder=2)
    _draw(nhd_flow, bb, ax, color='#4a4a4a', linewidth=0.25, alpha=0.55,
          zorder=3)
    counties.boundary.plot(ax=ax, color=MUTED, linewidth=0.6, zorder=4)
    gpd.GeoSeries([cg], crs=WORKING_CRS).boundary.plot(
        ax=ax, color=INK, linewidth=1.4, zorder=5)

    cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02,
                      extend='max' if _CUT_SAT > 0 else 'neither')
    cb.set_label(f'riparian cropland density \u00f7 the cut this ground had '
                 f'to clear\nblack line = 1.0, the cut itself '
                 f'({fmt(r["threshold_applied"])} ac / sq mi here)',
                 fontsize=FS['caption'], color=INK2)
    cb.ax.tick_params(labelsize=FS['tick'])
    cb.ax.plot([0, 1], [1.0] * 2, color=ZONE_INK, lw=1.6,
               transform=cb.ax.get_yaxis_transform(), clip_on=False)

    _shared = 0
    if len(zones):
        ids = split.loc[split['FIPS'] == fips, 'zone_id'].unique()
        zz = zones.loc[zones.index.isin(ids)]
        for zid, zr in zz.iterrows():
            out = zr.geometry.difference(cg)
            if not out.is_empty:
                _shared += 1
                gpd.GeoSeries([out], crs=WORKING_CRS).plot(
                    ax=ax, facecolor='none', edgecolor=ZONE_INK, linewidth=1.0,
                    hatch='///', alpha=0.55, zorder=6)
        if len(zz):
            zz.boundary.plot(ax=ax, color=ZONE_INK, linewidth=2.2, zorder=7)
            for zid, zr in zz.iterrows():
                c = zr.geometry.intersection(cg)
                c = (c if not c.is_empty else zr.geometry).representative_point()
                ax.annotate(str(zid), (c.x, c.y), fontsize=FS['panel_title'],
                            fontweight='bold', ha='center', va='center',
                            zorder=8,
                            bbox=dict(boxstyle='circle,pad=.28', fc='white',
                                      ec=ZONE_INK))
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3])
    _h2 = [Line2D([], [], color=ZONE_INK, lw=2.4, label='High-density zone')]
    if _shared:
        _h2.append(Patch(facecolor='none', edgecolor=ZONE_INK, hatch='///',
                         label='Zone continues past the county line'))
    legend_headroom(ax, 0.18)
    frame(ax.legend(handles=_h2, fontsize=FS['legend'], loc='upper left',
                    framealpha=0.96, borderpad=0.8, labelspacing=0.5,
                    title='High-density areas',
                    title_fontsize=FS['legend_title']))
    ax.set_title(f'B.  High-Density Areas — {BAND_MI:g}-mile density at least '
                 f'{HOT_MULTIPLE:g}× the {BASIS} average '
                 f'({int(r["n_zones"])} zones)',
                 fontsize=FS['panel_title'], loc='left')
    _bits = [f'{BASIS} average {r["reference_used"]:.1f} ac/sq mi \u00b7 '
             f'threshold {fmt(r["threshold_applied"])} ac/sq mi '
             f'({HOT_MULTIPLE:g}x, floor {fmt(FLOOR)})']
    # A county whose reference is so low that the ratio test is trivial is
    # judged by the FLOOR, and saying which cut bound it is the difference
    # between a sheet that reads as a bug and one that reads as a finding.
    if fips in _FLOOR_LIMITED:
        _mx = float(sub['rip_per_sqmi'].max())
        _bits.append(f'the binding cut here is the {fmt(FLOOR)} FLOOR, not '
                     f'the {HOT_MULTIPLE:g}x ratio '
                     f'({HOT_MULTIPLE:g} \u00d7 {r["reference_used"]:.1f} = '
                     f'{HOT_MULTIPLE * r["reference_used"]:.1f}); this '
                     f'county\'s densest neighbourhood is {fmt(_mx)} ac/sq mi')
    if r['n_zones']:
        _bits.append(f'{r["zone_sqmi"]:,.0f} sq mi of zone holding '
                     f'{r["zone_rip_ac"]:,.0f} ac '
                     f'({r["pct_resource_hot"]:.0f}% of the county\'s '
                     f'riparian cropland)')
    else:
        _bits.append(f'no zones - {r["status"].replace("_", " ")}')
    ax.set_xlabel('\n'.join(_bits), fontsize=FS['caption'], color=INK2)

    fig.suptitle(f'{name} {COUNTY_WORD}, {STATE_NAME}', fontsize=FS['fig_title'],
                 fontweight='bold', color=INK, x=0.012, ha='left')
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches='tight')
    show_close(fig)
    return True


if MAKE_COUNTY_FIGURES:
    _todo = man[man['n_cells'] > 0]
    if COUNTY_FIG_MAX:
        _todo = _todo.head(COUNTY_FIG_MAX)
    _ok = 0
    for _, _r in _todo.iterrows():
        _d = os.path.join(C_DIR, _r['county_id'], 'figures')
        os.makedirs(_d, exist_ok=True)
        _png = os.path.join(_d, f'{_r["county_id"]}_overlap_and_hotspots.png')
        try:
            _ok += bool(county_figure(_r['FIPS'], _r['NAME'], _png))
        except Exception as _e:
            man.loc[man['FIPS'] == _r['FIPS'], 'error'] = f'figure: {_e}'
            print(f'  {_r["NAME"]}: figure failed - {_e}')
    print(f'  {_ok} county sheets under counties/*/figures/ '
          f'(overlap + hotspots)')
    man.to_csv(MANIFEST_CSV, index=False)


# %% [markdown]
# ## Pass C — the parcel work, for the counties that have parcel data
#
# Passes A and B never touch a parcel. That is the whole point of them: the
# zone map, the density surface and the county ranking cost nothing in parcel
# licensing, which is what makes "buy Union and Ohio, not all 120" an argument
# with numbers behind it.
#
# This pass is the other side of that argument. For every county whose parcel
# file is already on disk — the `has_parcels` column of the manifest — it runs
# the same parcel analysis the single-county notebook runs, against the
# **statewide** surface and the **statewide** zones, and writes the result into
# that county's folder beside the zone map. Same fields, same file names, same
# meanings as the Christian County packet people have already been shown:
#
# | file | what is in it |
# |---|---|
# | `outreach_parcels_ranked.csv` | one row per candidate parcel, in rank order: owner, address, watershed, zone, rank, score, acreage, frontage, cover, water-feature ids, neighbours |
# | `riparian_cover_by_parcel_cdl.csv` | one row per parcel × CDL class — the unrolled version |
# | `huc12_ranked.csv` | the same ranking at subwatershed scale |
# | `zone_*.csv` | per-zone parcels, acres, owners, residency, owner class, concentration |
# | `parcels_enriched.gpkg` | every column above, as geometry |
#
# **Two things differ from the notebook, and both are improvements.**
# The neighbourhood density each parcel inherits comes from the statewide
# lattice, so a parcel three miles from the county line is no longer measured
# over a truncated band. And `zone_id` refers to the statewide zone, so a
# parcel in a zone that crosses into the next county carries the id that
# `zone_county_split.csv` uses — the biologist next door is looking at the same
# numbered zone.
#
# Everything here is skipped, with a printed reason, for a county with no
# parcel file. Nothing above this point depends on it.

# %% P0 - the helpers Pass C needs
# Every switch, threshold and reference table Pass C uses is in fl_config.py
# section 6 and section 9: RUN_PASS_C, PARCEL_COUNTIES, KEY, MIN_RIP_AC,
# W_PLACE, TOP_N_CONTACTS, FRONTAGE_STEP_M, ADJ_TOL_M, WANT_CONTACT,
# OWNER_PAT, ID_SOURCES, CDL_NAMES, CDL_GROUPS.
_CDL_GRP = {c: g for g, codes in CDL_GROUPS.items() for c in codes}

_joined = (lambda s: ','.join(sorted({str(x).strip() for x in s.dropna()
                                      if str(x).strip() not in ('', 'nan')})))


def _open_parcels(path):
    """Read a parcel file whatever shape it arrived in.

    Regrid ships a county as a ZIPPED GeoPackage. GDAL reads straight out of a
    zip through its virtual filesystem, so a download never has to be
    extracted onto Drive first; if reading in place fails the zip is extracted
    once into the cache folder and read from there."""
    low = str(path).lower()
    if low.endswith('.zip'):
        try:
            with zipfile.ZipFile(path) as z:
                names = z.namelist()
        except Exception as e:
            raise RuntimeError(f'{os.path.basename(path)} is not a readable '
                               f'zip: {e}')
        gdb = sorted({n[:n.lower().index('.gdb') + 4] for n in names
                      if '.gdb/' in n.lower()})
        hits = ([n for n in names if n.lower().endswith('.gpkg')] or gdb
                or [n for n in names if n.lower().endswith('.shp')])
        if not hits:
            raise RuntimeError(f'no dataset inside {os.path.basename(path)}')
        inner = sorted(hits, key=len)[0]
        try:
            return gpd.read_file(f'/vsizip/{path}/{inner}')
        except Exception:
            dest = os.path.join(os.path.dirname(CACHE), '_unzipped',
                                os.path.splitext(os.path.basename(path))[0])
            if not os.path.isdir(dest) or not os.listdir(dest):
                os.makedirs(dest, exist_ok=True)
                with zipfile.ZipFile(path) as z:
                    z.extractall(dest)
            for root, dirs, files in os.walk(dest):
                for n in sorted(dirs) + sorted(files):
                    if n.lower().endswith(('.gpkg', '.gdb', '.shp')):
                        return gpd.read_file(os.path.join(root, n))
            raise RuntimeError(f'nothing readable extracted from {path}')
    return gpd.read_file(path)


def water_edge_segments(bbox, step=None):
    """Every source's water edge, chopped into short straight SEGMENTS.

    WHY SEGMENTS, AND WHY ONLY THE WATER IS EVER BUFFERED
      Frontage is "how much water edge does this field run alongside".
      Proximity is symmetric, so buffering the cropland and hitting the edge
      gives the same answer as buffering the edge and hitting the cropland;
      switching sides changes nothing on its own. What it does change is that
      a SEGMENT can be buffered with FLAT ends. A round buffer reaches the
      full distance past the end of whatever it was built from, so edge lying
      beyond the end of a field still counts. Buffering each short piece with
      square ends confines it to its own length, and a piece past the end of
      the field finds no cropland in its own rectangle and drops out.

      The cropland is never buffered, dissolved or altered anywhere in this
      pipeline - not for acres, not for frontage.

    The answer is quantised to `step`: up to one step of error at each end of
    a run, i.e. +/- 10 m, against the 30 m over-run a round buffer produced.
    10 m is well inside NHD's own positional accuracy."""
    from shapely import segmentize, get_coordinates
    step = step or FRONTAGE_STEP_M
    out = {}
    for name in USE_SOURCES:
        spec = WATER_SOURCES[name]
        src = globals().get(spec['var'])
        if src is None:
            continue
        g = src.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        g = g[(~g.geometry.is_empty) &
              (g.geometry.area + g.geometry.length > 0)]
        if name == 'wetland' and WETLAND_EXCLUDE and len(g):
            g = drop_nwi_types(g)
        if not len(g):
            continue
        edge = (g.geometry if spec['kind'] == 'line' else g.geometry.boundary)
        edge = gpd.GeoSeries(edge.values,
                             crs=WORKING_CRS).explode(index_parts=False)
        edge = edge[edge.geom_type == 'LineString']
        if not len(edge):
            continue
        dense = segmentize(edge.values, step)
        xy, idx = get_coordinates(dense, return_index=True)
        keep = idx[:-1] == idx[1:]
        a, b = xy[:-1][keep], xy[1:][keep]
        segs = shapely.linestrings(
            np.stack([np.stack([a[:, 0], b[:, 0]], axis=1),
                      np.stack([a[:, 1], b[:, 1]], axis=1)], axis=2))
        out[name] = gpd.GeoDataFrame({'seg_id': np.arange(len(segs))},
                                     geometry=segs, crs=WORKING_CRS)
    return out


def _place_names(fips):
    """Place names inside this county, for the in-county mail test.

    DERIVED from the cached county-subdivision layer, never hand-typed. A town
    missing from a typed list silently reclassifies everyone who mails there
    as out-of-county, and that percentage goes in a briefing."""
    names, srcs = set(), []
    if cousub is not None and len(cousub) and 'COUNTYFP' in cousub.columns:
        _c = cousub[cousub['COUNTYFP'].astype(str).str.zfill(3) == fips]
        _n = next((c for c in _c.columns if c.upper() == 'NAME'), None)
        if _n and len(_c):
            names |= set(_c[_n].astype(str).str.upper().str.strip())
            srcs.append(f'cousub ({len(_c)})')
    _hit = sorted(glob.glob(os.path.join(RAW_DIR, '**', 'tl_*place*.shp'),
                            recursive=True))
    if _hit:
        try:
            _g = gpd.read_file(_hit[0]).to_crs(WORKING_CRS)
            _n = next((c for c in _g.columns if c.upper() == 'NAME'), None)
            _cg = counties.loc[counties['FIPS'] == fips, 'geometry'].iloc[0]
            _g = _g[_g.geometry.intersects(_cg)]
            if _n and len(_g):
                names |= set(_g[_n].astype(str).str.upper().str.strip())
                srcs.append(f'TIGER places ({len(_g)})')
        except Exception:
            pass
    return {n for n in names if n and n != 'NAN'}, (srcs or ['none found'])


print(f'Pass C: RUN_PASS_C = {RUN_PASS_C} | PARCEL_COUNTIES = '
      f'{PARCEL_COUNTIES} | floor {MIN_RIP_AC:g} ac | W_PLACE {W_PLACE:g}')


# %% P1 - the parcel analysis for one county
# Step 2. One function, so a county either produces the whole packet or
# produces none of it and says why. It reads the statewide objects (grid,
# zones, split, rip_crop, csb, nhd_*) and writes only into this county's
# folder.
def pass_c(fips, name, county_id, row):
    out_dir = os.path.join(C_DIR, county_id, 'output')
    fig_dir = os.path.join(C_DIR, county_id, 'figures')
    val_dir = os.path.join(C_DIR, county_id, 'validation')
    for _d in (out_dir, fig_dir, val_dir):
        os.makedirs(_d, exist_ok=True)
    res = {'FIPS': fips, 'NAME': name}

    cext = counties.loc[counties['FIPS'] == fips, 'extent'].iloc[0]
    cg = counties.loc[counties['FIPS'] == fips, 'geometry'].iloc[0]

    # --- 2a. read and clip -------------------------------------------------
    cl = _open_parcels(row['parcel_path'])
    if cl.crs is None:
        cl = cl.set_crs(WORKING_CRS)
    cl = cl.to_crs(WORKING_CRS)
    if KEY not in cl.columns:
        raise KeyError(
            f'the parcel file has no {KEY!r} column. Columns are '
            f'{[c for c in cl.columns if c != "geometry"][:20]}. Set KEY in '
            f'P0 to whatever this source calls its parcel id.')
    _pre = len(cl)
    cl = clip_to(cl, cext, f'{name} parcels')
    cl = cl.drop_duplicates(subset=KEY).reset_index(drop=True)
    cl['parcel_ac'] = cl.geometry.area / AC
    print(f'    parcels        {len(cl):>8,} in the analysed extent '
          f'(file held {_pre:,})')
    if not len(cl):
        raise ValueError('no parcels inside this county\'s analysed extent - '
                         'the file is for a different county')
    pts = gpd.GeoDataFrame(cl[[KEY]],
                           geometry=cl.geometry.representative_point().values,
                           crs=WORKING_CRS)
    bx = cl.total_bounds

    # --- 2b. riparian cropland on each parcel ------------------------------
    # ONE overlay, with the source flags carried through, so the per-source
    # acres come out of the same intersection rather than four more.
    _ON = [f'on_{n}' for n in USE_SOURCES if f'on_{n}' in rip_crop.columns]
    _cols = _ON + (['CDL'] if 'CDL' in rip_crop.columns else []) + ['geometry']
    _rc = rip_crop.cx[bx[0]:bx[2], bx[1]:bx[3]][_cols]
    pc = gpd.overlay(cl[[KEY, 'geometry']], _rc, how='intersection',
                     keep_geom_type=True)
    pc['ac'] = pc.geometry.area / AC
    cl['rip_crop_ac'] = cl[KEY].map(pc.groupby(KEY)['ac'].sum()).fillna(0.0)
    for _n in USE_SOURCES:
        _c = f'on_{_n}'
        cl[f'rip_ac_{_n}'] = (
            cl[KEY].map(pc[pc[_c]].groupby(KEY)['ac'].sum()).fillna(0.0)
            if _c in pc.columns else 0.0)
    # A piece inside two corridors counts in both columns, so the per-source
    # acres can sum past rip_crop_ac. rip_crop_ac (the union) is the total.
    _sc = [f'rip_ac_{n}' for n in USE_SOURCES]
    cl['water_source'] = (cl[_sc].idxmax(axis=1)
                          .str.replace('rip_ac_', '', regex=False)
                          .where(cl['rip_crop_ac'] > 0, 'none'))
    print(f'    riparian crop  {cl["rip_crop_ac"].sum():>8,.0f} ac on '
          f'{int((cl["rip_crop_ac"] > 0).sum()):,} parcels')

    # --- 2c. the two whole-parcel resource totals --------------------------
    _cs = csb.cx[bx[0]:bx[2], bx[1]:bx[3]]
    if len(_cs):
        _cx = gpd.overlay(cl[[KEY, 'geometry']], _cs[['geometry']],
                          how='intersection', keep_geom_type=True)
        cl['csb_overlap_ac'] = cl[KEY].map(
            (_cx.geometry.area / AC).groupby(_cx[KEY]).sum()).fillna(0.0)
    else:
        cl['csb_overlap_ac'] = 0.0
    if nhd_flow is not None and len(nhd_flow):
        _nf = nhd_flow.cx[bx[0]:bx[2], bx[1]:bx[3]]
        _sx = gpd.overlay(cl[[KEY, 'geometry']], _nf[['geometry']],
                          how='intersection', keep_geom_type=False)
        _sx = _sx[_sx.geometry.geom_type.isin(['LineString',
                                               'MultiLineString'])]
        cl['nhd_stream_mi'] = cl[KEY].map(
            (_sx.geometry.length / MI_M).groupby(_sx[KEY]).sum()).fillna(0.0)
    else:
        cl['nhd_stream_mi'] = 0.0

    # --- 2d. frontage, measured along the water edge -----------------------
    _edges = water_edge_segments(bx)
    _d_fr = RIP_BUFFER_FT * 0.3048
    _M_FT = 1.0 / 0.3048
    _frcols = []
    for _n in USE_SOURCES:
        _col = f'frontage_ft_{_n}'
        _frcols.append(_col)
        cl[_col] = 0.0
        _flag = f'on_{_n}'
        if _n not in _edges or _flag not in pc.columns:
            continue
        _seg = pc[pc[_flag]][[KEY, 'geometry']]
        if not len(_seg):
            continue
        _eb = _edges[_n].copy()
        _eb['len_ft'] = _eb.geometry.length * _M_FT
        _eb['geometry'] = _eb.geometry.buffer(_d_fr, cap_style='flat')
        _hit = gpd.sjoin(_eb, _seg, predicate='intersects', how='inner')
        if not len(_hit):
            continue
        # one parcel can own several cropland slivers against the same piece
        # of edge; without this the piece is counted once per sliver
        _hit = _hit.drop_duplicates(subset=[KEY, 'seg_id'])
        cl[_col] = cl[KEY].map(
            _hit.groupby(KEY)['len_ft'].sum()).fillna(0.0).round(0)
    cl['frontage_ft_total'] = cl[_frcols].sum(axis=1).round(0)
    print(f'    frontage       {cl["frontage_ft_total"].sum() / 5280:>8,.0f} mi '
          f'measured along the water edge ({FRONTAGE_STEP_M:g} m steps, '
          f'flat caps)')

    # --- 2e. what is actually growing in the buffer ------------------------
    percdl = None
    if 'CDL' in pc.columns:
        _cd = pc.copy()
        _cd['CDL'] = pd.to_numeric(_cd['CDL'], errors='coerce').astype('Int64')
        _cd['grp'] = [_CDL_GRP.get(int(c), 'row_crop') if pd.notna(c)
                      else 'unknown' for c in _cd['CDL']]
        percdl = (_cd.groupby([KEY, 'CDL'], observed=True)['ac'].sum()
                  .round(3).reset_index().rename(columns={'ac': 'buffer_ac'}))
        percdl['cdl_name'] = percdl['CDL'].map(
            lambda c: CDL_NAMES.get(int(c), f'CDL {c}') if pd.notna(c)
            else 'unknown')
        percdl['grp'] = [_CDL_GRP.get(int(c), 'row_crop') if pd.notna(c)
                         else 'unknown' for c in percdl['CDL']]
        percdl = percdl.sort_values([KEY, 'buffer_ac'],
                                    ascending=[True, False])
        percdl.to_csv(os.path.join(out_dir,
                                   'riparian_cover_by_parcel_cdl.csv'),
                      index=False)
        cl['rip_cdl_acres'] = cl[KEY].map(
            percdl.groupby(KEY).apply(
                lambda d: ','.join(f'{n}:{a:g}' for n, a
                                   in zip(d['cdl_name'], d['buffer_ac'])),
                include_groups=False)).fillna('')
        cl['rip_n_cdl'] = cl[KEY].map(
            percdl.groupby(KEY).size()).fillna(0).astype(int)
        _gc = []
        for _g in GRP_ORDER:
            _c = f'rip_{_g}_ac'
            _gc.append(_c)
            cl[_c] = cl[KEY].map(
                _cd[_cd['grp'] == _g].groupby(KEY)['ac'].sum()
            ).fillna(0.0).round(2)
        cl['rip_cover'] = (cl[_gc].idxmax(axis=1)
                           .str.replace('rip_', '', regex=False)
                           .str.replace('_ac', '', regex=False)
                           .where(cl[_gc].max(axis=1) > 0, 'none'))
        _dom = percdl.drop_duplicates(KEY).set_index(KEY)
        cl['rip_cdl_dominant'] = cl[KEY].map(_dom['CDL'])
        cl['rip_cdl_dominant_name'] = cl[KEY].map(_dom['cdl_name'])
        _nc = sum(cl[f'rip_{g}_ac'].sum() for g in GRP_ORDER
                  if g != 'row_crop')
        print(f'    buffer cover   row crop '
              f'{cl["rip_row_crop_ac"].sum():,.0f} ac · other '
              f'{_nc:,.0f} ac (already-vegetated, idle, grass, water, built)')
    else:
        print('    buffer cover   SKIPPED - no CDL column on the overlay')

    # --- 2f. the size-neutral parcel metrics -------------------------------
    cl['rip_crop_per_ac'] = (cl['rip_crop_ac'] /
                             cl['parcel_ac'].replace(0, np.nan)).clip(0, 1)
    cl['rip_crop_frac'] = (cl['rip_crop_ac'] /
                           cl['csb_overlap_ac'].replace(0, np.nan)
                           ).fillna(0).clip(0, 1)

    # --- 2g. where each parcel sits on the STATEWIDE surface ---------------
    # This is the line that makes Pass C different from running the notebook
    # county by county: the cell a parcel inherits its neighbourhood density
    # from was measured over a full 3-mile band, not one cut at the county
    # line. `zone_id` is the STATEWIDE zone id, so it matches
    # zone_county_split.csv and the county packet's zones.gpkg.
    _gsub = grid.reset_index()
    _gsub = _gsub.cx[bx[0] - BAND_M:bx[2] + BAND_M,
                     bx[1] - BAND_M:bx[3] + BAND_M]
    _cell_of = gpd.sjoin(
        pts, _gsub[['cell', 'rip_per_sqmi', 'ratio', 'ratio_county',
                    'reference', 'p_state', 'hot', 'geometry']],
        predicate='within', how='left')
    _cell_of = _cell_of.loc[~_cell_of[KEY].duplicated(keep='first')] \
                       .set_index(KEY)
    cl['cell'] = cl[KEY].map(_cell_of['cell'])
    cl['nbhd_rip_per_sqmi'] = cl[KEY].map(
        _cell_of['rip_per_sqmi']).fillna(0.0)
    cl['nbhd_ratio'] = cl[KEY].map(_cell_of['ratio']).fillna(0.0)
    cl['nbhd_vs_county'] = cl[KEY].map(_cell_of['ratio_county']).fillna(0.0)
    cl['cell_reference'] = cl[KEY].map(_cell_of['reference'])
    cl['cell_p_state'] = cl[KEY].map(_cell_of['p_state'])
    cl['cell_is_hot'] = cl[KEY].map(_cell_of['hot']).fillna(False).astype(bool)

    if len(zones):
        _zj = gpd.sjoin(pts, gpd.GeoDataFrame(
            zones[[]].reset_index()[['zone_id']],
            geometry=np.asarray(zones.geometry.values, dtype=object),
            crs=WORKING_CRS), predicate='within', how='left')
        _zj = _zj.loc[~_zj[KEY].duplicated(keep='first')].set_index(KEY)
        cl['zone_id'] = cl[KEY].map(_zj['zone_id'])
    else:
        cl['zone_id'] = pd.NA
    cl['in_hotspot'] = cl['zone_id'].notna()
    print(f'    in a zone      {int(cl["in_hotspot"].sum()):>8,} parcels '
          f'({cl["in_hotspot"].mean():.0%}), holding '
          f'{cl.loc[cl["in_hotspot"], "rip_crop_ac"].sum():,.0f} ac')

    # --- 2h. HUC12 label, through one interior point -----------------------
    # A parcel straddling a subwatershed boundary intersects two polygons, so
    # an `intersects` join returns it twice and every acre summed afterwards
    # is inflated. A point cannot fall in two polygons.
    cl['huc12'] = pd.NA
    cl['huc12_name'] = pd.NA
    if huc12 is not None and len(huc12):
        _hc = next((c for c in huc12.columns if c.lower() == 'huc12'), None)
        _hn = next((c for c in huc12.columns
                    if c.lower() in ('name', 'huc12_name')), None)
        if _hc:
            _h12 = huc12[huc12[_hc].astype(str).str.len() == 12]
            _h12 = _h12.cx[bx[0]:bx[2], bx[1]:bx[3]]
            _keep = [_hc] + ([_hn] if _hn else [])
            if len(_h12):
                _hj = gpd.sjoin(pts, _h12[_keep + ['geometry']],
                                predicate='within', how='left')
                _hj = _hj.loc[~_hj[KEY].duplicated(keep='first')] \
                         .set_index(KEY)
                cl['huc12'] = cl[KEY].map(_hj[_hc]).astype('string')
                if _hn:
                    cl['huc12_name'] = cl[KEY].map(_hj[_hn]).astype('string')

    # --- 2i. the materiality floor and the ranking -------------------------
    # Both terms are DENSITIES and both are converted to percentile ranks
    # first, so the weight between them means something and a parcel is never
    # a priority for being large.
    _fl = pd.DataFrame([{
        'min_rip_ac': f,
        'parcels': int((cl['rip_crop_ac'] > f).sum()),
        'rip_crop_ac': round(cl.loc[cl['rip_crop_ac'] > f,
                                    'rip_crop_ac'].sum()),
        'pct_of_rip': round(100 * cl.loc[cl['rip_crop_ac'] > f,
                                         'rip_crop_ac'].sum() /
                            max(cl['rip_crop_ac'].sum(), 1e-9)),
        'median_frontage_ft': (
            round(cl.loc[cl['rip_crop_ac'] > f, 'frontage_ft_total'].median())
            if (cl['rip_crop_ac'] > f).any() else None),
        'median_parcel_ac': (
            round(cl.loc[cl['rip_crop_ac'] > f, 'parcel_ac'].median(), 1)
            if (cl['rip_crop_ac'] > f).any() else None),
    } for f in FLOOR_SWEEP_PARCEL]).set_index('min_rip_ac')
    _fl.to_csv(os.path.join(val_dir, 'parcel_floor_sweep.csv'))

    cand = cl[cl['rip_crop_ac'] > MIN_RIP_AC].copy()
    if not len(cand):
        raise ValueError(f'no parcels above MIN_RIP_AC = {MIN_RIP_AC:g} - '
                         f'lower it, or widen RIP_BUFFER_FT')
    cand['own_rank'] = cand['rip_crop_per_ac'].rank(pct=True).round(3)
    cand['place_rank'] = cand['nbhd_rip_per_sqmi'].rank(pct=True).round(3)
    cand['score'] = (W_PLACE * cand['place_rank'] +
                     (1 - W_PLACE) * cand['own_rank']).round(3)
    _own_col = next((c for c in ['owner', 'owner1', 'ownername']
                     if c in cand.columns), None)
    cand['has_contact'] = (
        (cand[_own_col].notna() &
         cand[_own_col].astype(str).str.strip().ne('') &
         cand[_own_col].astype(str).str.lower().ne('none'))
        if _own_col else False)
    # in_hotspot is deliberately NOT a sort key. The zone boundary moves with
    # the grid origin, so it is carried as a COLUMN and does not decide who
    # gets called first.
    cand = cand.sort_values(['has_contact', 'score'], ascending=[False, False])
    cand['rank'] = np.arange(1, len(cand) + 1)
    res['n_candidates'] = len(cand)

    # size-neutrality, checked rather than asserted
    _pre_f = cl[cl['rip_crop_ac'] > 0]
    res['rho_density_vs_size_prefloor'] = round(float(
        _pre_f['rip_crop_per_ac'].rank(pct=True)
        .corr(_pre_f['parcel_ac'], method='spearman')), 2) if len(_pre_f) > 2 else np.nan
    res['rho_score_vs_size'] = round(float(
        cand['score'].corr(cand['parcel_ac'], method='spearman')), 2)

    # --- 2j. the outreach table --------------------------------------------
    _low = {c.lower(): c for c in cand.columns}
    _found, _missing = {}, []
    for _k, _opts in WANT_CONTACT.items():
        _h = next((_low[o] for o in _opts if o in _low), None)
        (_found.setdefault(_k, _h) if _h else _missing.append(_k))
    contacts = pd.DataFrame({'rank': cand['rank'].values,
                             'score': cand['score'].values})
    for _k, _src in _found.items():
        contacts[_k] = cand[_src].values
    for _c in ['huc12', 'huc12_name', 'zone_id', 'in_hotspot', 'cell',
               'cell_is_hot', 'cell_p_state', 'nbhd_rip_per_sqmi',
               'nbhd_ratio', 'nbhd_vs_county', 'cell_reference',
               'rip_crop_ac', 'rip_crop_per_ac', 'rip_crop_frac', 'parcel_ac',
               'csb_overlap_ac', 'nhd_stream_mi', 'water_source',
               'frontage_ft_total', 'rip_cover', 'rip_cdl_dominant',
               'rip_cdl_dominant_name', 'rip_cdl_acres', 'rip_n_cdl',
               'own_rank', 'place_rank', 'has_contact'] \
        + [f'rip_ac_{n}' for n in USE_SOURCES] \
        + [f'frontage_ft_{n}' for n in USE_SOURCES] \
        + [f'rip_{g}_ac' for g in GRP_ORDER]:
        if _c in cand.columns:
            contacts[_c] = cand[_c].values
    contacts['basis'] = BASIS
    contacts['floor_ac_per_sqmi'] = round(FLOOR, 2)
    contacts[KEY] = cand[KEY].values
    _ll = gpd.GeoSeries(cand.geometry.representative_point().values,
                        crs=WORKING_CRS).to_crs('EPSG:4326')
    contacts['lon'] = _ll.x.round(6).values
    contacts['lat'] = _ll.y.round(6).values

    # --- 2k. the field-rep columns -----------------------------------------
    CTY_PLACES, _src_places = _place_names(fips)
    if 'mail_state' in contacts.columns:
        # fillna BEFORE astype: pandas 2.x astype(str) preserves NA, and
        # NaN.ne('KY') is True, which is how unknown-owner rows get called
        # absentee.
        _st = contacts['mail_state'].fillna('').astype(str).str.strip() \
                                    .str.upper().replace(
                                        {'NAN': '', 'NONE': '', 'NA': ''})
        _ci = (contacts['mail_city'].fillna('').astype(str).str.strip()
               .str.upper() if 'mail_city' in contacts.columns
               else pd.Series('', index=contacts.index))
        contacts['contact_out_of_state'] = _st.ne(HOME_STATE) & _st.ne('')
        contacts['contact_in_state_not_county'] = (
            _st.eq(HOME_STATE) & ~_ci.isin(CTY_PLACES))
    else:
        contacts['contact_out_of_state'] = False
        contacts['contact_in_state_not_county'] = False

    _own = (contacts['owner'].fillna('').astype(str).str.upper()
            if 'owner' in contacts.columns
            else pd.Series('', index=contacts.index))
    _use = (contacts['land_use'].fillna('').astype(str).str.upper()
            if 'land_use' in contacts.columns
            else pd.Series('', index=contacts.index))
    _cls = pd.Series('private', index=contacts.index)
    for _k in ['institutional', 'state_local', 'federal']:   # federal wins
        _cls[_own.str.contains(OWNER_PAT[_k], regex=True, na=False)] = _k
    _blank = _own.str.strip().isin(['', 'NAN', 'NONE', 'NA'])
    _cls[_blank & _use.str.contains(USE_PAT, regex=True, na=False)] = 'review'
    _cls[_blank & ~_use.str.contains(USE_PAT, regex=True,
                                     na=False)] = 'unknown_owner'
    contacts['owner_class'] = _cls

    # water-feature ids, at the SAME distance the acres were measured at
    _pg = cand[[KEY, 'geometry']]
    _pbx = _pg.total_bounds
    _d_id = RIP_BUFFER_FT * 0.3048
    for _n, _cnt, _var, _pref, _fb in ID_SOURCES:
        _src = globals().get(_var)
        contacts[f'{_n}_ids'] = ''
        contacts[_cnt] = 0
        if _src is None or not len(_src):
            continue
        _idc = next((c for c in _pref if c in _src.columns),
                    next((c for c in _fb if c in _src.columns), None))
        if _idc is None:
            continue
        _g = _src.cx[_pbx[0]:_pbx[2], _pbx[1]:_pbx[3]][[_idc, 'geometry']] \
                 .copy()
        if not len(_g):
            continue
        _g['geometry'] = _g.geometry.buffer(_d_id)
        _hj = gpd.sjoin(_pg, _g, predicate='intersects', how='inner')
        if not len(_hj):
            continue
        _agg = _hj.groupby(KEY)[_idc].agg([_joined, 'nunique'])
        _agg.columns = [f'{_n}_ids', _cnt]
        contacts[f'{_n}_ids'] = contacts[KEY].map(
            _agg[f'{_n}_ids']).fillna('')
        contacts[_cnt] = contacts[KEY].map(_agg[_cnt]).fillna(0).astype(int)

    # ADJACENCY - can this be worked as a block, or is it one door?
    _pbuf = cl[[KEY, 'geometry']].copy()
    _pbuf['geometry'] = _pbuf.geometry.buffer(ADJ_TOL_M)
    _on_list = set(cand[KEY])
    _rip_any = set(cl.loc[cl['rip_crop_ac'] > 0, KEY])
    _pairs = gpd.sjoin(_pg, _pbuf, predicate='intersects', how='inner',
                       lsuffix='L', rsuffix='R')
    _pairs = _pairs[_pairs[f'{KEY}_L'] != _pairs[f'{KEY}_R']]
    _hit = _pairs[_pairs[f'{KEY}_R'].isin(_on_list)]
    contacts['n_adj_parcels_on_list'] = contacts[KEY].map(
        _hit.groupby(f'{KEY}_L')[f'{KEY}_R'].nunique()).fillna(0).astype(int)
    contacts['adj_parcel_on_list'] = contacts['n_adj_parcels_on_list'] > 0
    contacts['adj_parcel_ids'] = contacts[KEY].map(
        _hit.groupby(f'{KEY}_L')[f'{KEY}_R'].agg(_joined)).fillna('')
    _off = _pairs[_pairs[f'{KEY}_R'].isin(_rip_any - _on_list)]
    contacts['n_adj_riparian_offlist'] = contacts[KEY].map(
        _off.groupby(f'{KEY}_L')[f'{KEY}_R'].nunique()).fillna(0).astype(int)

    contacts.to_csv(os.path.join(out_dir, 'outreach_parcels_ranked.csv'),
                    index=False)
    # The short list is a SLICE of the file above, which carries `rank`,
    # so TOP_N_CONTACTS = None removes a FILE and no information:
    # head(250) of outreach_parcels_ranked.csv is the same 250 rows, and
    # it cannot go stale beside a re-run of the full list. Set an
    # integer in fl_config.py to get the slice written again.
    if TOP_N_CONTACTS:
        contacts.head(int(TOP_N_CONTACTS)).to_csv(
            os.path.join(out_dir,
                         f'outreach_top{int(TOP_N_CONTACTS)}.csv'),
            index=False)
    print(f'    outreach       {len(contacts):>8,} rows x '
          f'{contacts.shape[1]} columns'
          + (f'   (contact fields missing: {", ".join(_missing)})'
             if _missing else ''))
    print(f'    adjacency      '
          f'{int(contacts["adj_parcel_on_list"].sum()):>8,} parcels '
          f'({contacts["adj_parcel_on_list"].mean():.0%}) touch another on '
          f'this list')
    res['n_parcels'] = len(cl)
    res['parcel_rip_ac'] = round(float(cl['rip_crop_ac'].sum()))
    res['parcels_in_zone'] = int(cl['in_hotspot'].sum())
    res['pct_contactable'] = round(100 * float(contacts['has_contact'].mean()))
    res['owners'] = int(contacts[_found['owner']].nunique()) \
        if 'owner' in _found and _found['owner'] in contacts.columns else \
        (int(contacts['owner'].nunique()) if 'owner' in contacts.columns else 0)
    res['pct_out_of_state'] = round(
        100 * float(contacts['contact_out_of_state'].mean()))
    res['frontage_mi'] = round(float(cl['frontage_ft_total'].sum()) / 5280, 1)

    # --- 2l. HUC12 rollup ---------------------------------------------------
    if cl['huc12'].notna().any() and huc12 is not None:
        _hc = next((c for c in huc12.columns if c.lower() == 'huc12'), None)
        _hs = (cl[cl['huc12'].notna()].groupby('huc12')
               .agg(n_parcels=(KEY, 'size'),
                    rip_crop_ac=('rip_crop_ac', 'sum'),
                    crop_ac=('csb_overlap_ac', 'sum'),
                    stream_mi=('nhd_stream_mi', 'sum'),
                    frontage_ft=('frontage_ft_total', 'sum'),
                    in_zone=('in_hotspot', 'sum')).reset_index())
        _hg = huc12[huc12[_hc].astype(str).str.len() == 12].copy()
        _hg['huc12'] = _hg[_hc].astype(str)
        _hg['geometry'] = _hg.geometry.intersection(cext)
        _hg = _hg[_hg.geometry.area > 0]
        _hg['analysed_sqmi'] = _hg.geometry.area / SQMI
        _hs = _hs.merge(_hg[['huc12', 'analysed_sqmi']], on='huc12',
                        how='inner')
        _hs = _hs[_hs['analysed_sqmi'] >= 1.0]
        _hs['rip_per_sqmi'] = (_hs['rip_crop_ac'] /
                               _hs['analysed_sqmi']).round(2)
        _hs['vs_state'] = (_hs['rip_per_sqmi'] / AVG_STATE).round(2)
        _hs = _hs.sort_values('rip_per_sqmi', ascending=False)
        _hs.round(2).to_csv(os.path.join(out_dir, 'huc12_ranked.csv'),
                            index=False)

    # --- 2m. the zone tables ------------------------------------------------
    zc = contacts.copy()
    zc['zone_lab'] = zc['zone_id'].map(
        lambda z: f'zone {int(z)}' if pd.notna(z) else 'outside any zone')
    _ZORD = ([f'zone {int(z)}' for z in sorted(
        pd.to_numeric(contacts['zone_id'], errors='coerce').dropna().unique())]
        + ['outside any zone'])
    _IN = zc['zone_lab'].ne('outside any zone')
    NZ = len(zc)
    TOTAC = float(zc['rip_crop_ac'].sum())

    _rows = []
    for _lab in _ZORD:
        d = zc[zc['zone_lab'] == _lab]
        if not len(d):
            continue
        r = {'zone': _lab, 'parcels': len(d),
             'pct_parcels': round(len(d) / NZ, 4),
             'rip_crop_ac': round(float(d['rip_crop_ac'].sum()), 1),
             'pct_acres': round(float(d['rip_crop_ac'].sum()) /
                                max(TOTAC, 1e-9), 4),
             'median_rip_ac': round(float(d['rip_crop_ac'].median()), 2),
             'median_parcel_ac': round(float(d['parcel_ac'].median()), 1),
             'median_frontage_ft': round(float(
                 d['frontage_ft_total'].median())),
             'pct_contactable': round(100 * float(
                 d['has_contact'].astype(bool).mean()))}
        if 'owner' in d.columns:
            r['owners'] = int(d['owner'].nunique())
            # the leverage number: a zone where one conversation covers six
            # tracts is worth more than its acreage suggests
            r['parcels_per_owner'] = round(len(d) /
                                           max(d['owner'].nunique(), 1), 1)
        _rows.append(r)
    ztbl = pd.DataFrame(_rows)
    if len(zones) and len(ztbl):
        _zg2 = zones[[c for c in ['zone_ac', 'rip_per_sqmi', 'vs_state',
                                  'n_counties', 'counties', 'stream_mi',
                                  'n_stream_segments', 'named_waters']
                      if c in zones.columns]].copy()
        _zg2.index = [f'zone {int(i)}' for i in _zg2.index]
        ztbl = ztbl.merge(_zg2.reset_index().rename(columns={'index': 'zone'}),
                          on='zone', how='left')
        if 'zone_ac' in ztbl.columns:
            ztbl['zone_sqmi'] = (ztbl['zone_ac'] / 640).round(1)
            ztbl = ztbl.drop(columns='zone_ac')
    ztbl.to_csv(os.path.join(out_dir, 'zone_counts.csv'), index=False)

    if 'water_source' in zc.columns:
        _xn = pd.crosstab(zc['water_source'], zc['zone_lab']).reindex(
            columns=_ZORD, fill_value=0)
        _xn['TOTAL'] = _xn.sum(axis=1)
        _xn.loc['TOTAL'] = _xn.sum()
        _xn.to_csv(os.path.join(out_dir, 'zone_by_water_parcels.csv'))
        _xa = zc.pivot_table(index='water_source', columns='zone_lab',
                             values='rip_crop_ac', aggfunc='sum',
                             fill_value=0.0).reindex(columns=_ZORD,
                                                     fill_value=0.0)
        _xa['TOTAL'] = _xa.sum(axis=1)
        _xa.loc['TOTAL'] = _xa.sum()
        _xa.round(2).to_csv(os.path.join(out_dir, 'zone_by_water_acres.csv'))

    # PRESENCE, not dominance: a parcel is counted under every water type it
    # touches, so these rows do NOT sum to the list total. pct_of_type_in_zone
    # is the question - do the zones favour one kind of water over another?
    _pres = []
    for _n in USE_SOURCES:
        _c = f'rip_ac_{_n}'
        if _c not in zc.columns:
            continue
        m = zc[_c] > 0
        _pres.append({'water_type': _n, 'parcels_any': int(m.sum()),
                      'pct_of_list': round(float(m.mean()), 4),
                      'riparian_ac': round(float(zc.loc[m, _c].sum()), 1),
                      'in_zone': int((m & _IN).sum()),
                      'pct_of_type_in_zone': round(
                          100 * (m & _IN).sum() / max(m.sum(), 1))})
    if _pres:
        pd.DataFrame(_pres).to_csv(
            os.path.join(out_dir, 'zone_water_presence.csv'), index=False)

    ZRES = ['In county', 'In state, out of county', 'Out of state',
            'No owner record']
    _res_s = pd.Series('In county', index=zc.index)
    _res_s[zc['contact_in_state_not_county'].astype(bool)] = ZRES[1]
    _res_s[zc['contact_out_of_state'].astype(bool)] = ZRES[2]
    _res_s[~zc['has_contact'].astype(bool)] = ZRES[3]
    zc['residency'] = _res_s
    _r3 = (zc.groupby('residency')
           .agg(parcels=('rip_crop_ac', 'size'), acres=('rip_crop_ac', 'sum'))
           .reindex(ZRES).fillna(0))
    _r3['pct_parcels'] = (_r3['parcels'] / NZ).round(4)
    _r3['in_zone'] = [int((_IN & _res_s.eq(k)).sum()) for k in ZRES]
    _r3['pct_of_group_in_zone'] = (100 * _r3['in_zone'] /
                                   _r3['parcels'].replace(0, np.nan)).round(0)
    _r3.round(2).to_csv(os.path.join(out_dir, 'zone_residency.csv'))

    _r4 = (zc.groupby('owner_class')
           .agg(parcels=('rip_crop_ac', 'size'), acres=('rip_crop_ac', 'sum')))
    _r4['pct_parcels'] = (_r4['parcels'] / NZ).round(4)
    _r4['in_zone'] = zc[_IN].groupby(
        zc.loc[_IN, 'owner_class']).size().reindex(_r4.index).fillna(0).astype(int)
    _r4.round(2).to_csv(os.path.join(out_dir, 'zone_owner_class.csv'))

    _s10 = zc['rip_crop_ac'].sort_values(ascending=False).reset_index(drop=True)
    _conc = [{'top_n': k, 'acres': round(float(_s10.head(k).sum()), 1),
              'pct_of_acres': round(float(_s10.head(k).sum()) /
                                    max(TOTAC, 1e-9), 4),
              'min_ac_at_n': round(float(_s10.iloc[k - 1]), 2)}
             for k in [10, 25, 50, 100, 250, 500, NZ] if k <= NZ]
    pd.DataFrame(_conc).to_csv(
        os.path.join(out_dir, 'zone_concentration.csv'), index=False)

    if 'owner' in zc.columns and len(zones):
        _zg3 = zc[_IN].groupby('zone_id')
        _zt = _zg3.agg(parcels=(KEY, 'size'),
                       rip_crop_ac=('rip_crop_ac', 'sum'),
                       frontage_ft=('frontage_ft_total', 'sum'))
        _zt['owners'] = _zg3['owner'].nunique()
        _zt['parcels_per_owner'] = (_zt['parcels'] /
                                    _zt['owners'].replace(0, np.nan)).round(1)
        _zt['pct_contactable'] = (_zg3['has_contact']
                                  .apply(lambda s: s.astype(bool).mean())
                                  * 100).round(0)
        _zt = _zt.join(zones[[c for c in ['rip_per_sqmi', 'vs_state',
                                          'stream_mi', 'n_stream_segments']
                              if c in zones.columns]])
        _zt.round(2).to_csv(os.path.join(out_dir, 'zone_ownership.csv'))

    # IN-ZONE vs OUT-OF-ZONE. The claim the method makes is that it selects
    # PLACES, not parcels, so neighbourhood density should separate the two
    # groups and parcel SIZE should not. This is the evidence for it.
    _cmpz = []
    for _c, _lab in [('parcel_ac', 'parcel acres'),
                     ('rip_crop_ac', 'riparian cropland ac'),
                     ('rip_crop_per_ac', 'riparian share of parcel'),
                     ('frontage_ft_total', 'measured water frontage (ft)'),
                     ('nbhd_ratio', f'neighbourhood density (x {BASIS})')]:
        if _c not in zc.columns:
            continue
        v = pd.to_numeric(zc[_c], errors='coerce')
        _cmpz.append({'measure': _lab,
                      'median_in_zone': round(float(v[_IN].median()), 3),
                      'median_outside': round(float(v[~_IN].median()), 3),
                      'ratio': round(float(v[_IN].median()) /
                                     max(float(v[~_IN].median()), 1e-9), 2)})
    cmp_t = pd.DataFrame(_cmpz)
    cmp_t.to_csv(os.path.join(val_dir, 'zone_in_vs_out.csv'), index=False)
    _sz = cmp_t[cmp_t['measure'] == 'parcel acres']
    if len(_sz):
        res['size_ratio_in_vs_out'] = float(_sz.iloc[0]['ratio'])

    # PARCEL CROSS-CHECK. The surface never saw a parcel, so parcels-per-
    # square-mile is a genuinely independent measure of the same
    # concentration. High agreement is the check.
    _ck = []
    if cl['cell'].notna().any():
        for _flr, _lab in [(0.0, 'any riparian cropland'),
                           (MIN_RIP_AC, f'> {MIN_RIP_AC:g} ac')]:
            _cnt = (cl[cl['rip_crop_ac'] > _flr].groupby('cell').size())
            _sub = grid.loc[grid.index.isin(cl['cell'].dropna().unique())]
            _dens = (_cnt.reindex(_sub.index).fillna(0) /
                     _sub['land_ac'].replace(0, np.nan) * 640)
            _ck.append({'filter': _lab,
                        'rho_vs_neighbourhood': round(float(
                            _dens.corr(_sub['rip_per_sqmi'],
                                       method='spearman')), 3),
                        'rho_vs_own_cell': round(float(
                            _dens.corr(_sub['own_per_sqmi'],
                                       method='spearman')), 3)})
        pd.DataFrame(_ck).to_csv(
            os.path.join(val_dir, 'parcel_crosscheck.csv'), index=False)
        if _ck:
            res['rho_parcels_vs_density'] = _ck[0]['rho_vs_neighbourhood']
            print(f'    cross-check    parcels per sq mi vs the surface '
                  f'density: rho = {_ck[0]["rho_vs_neighbourhood"]:+.2f} '
                  f'(independent of it)')

    # --- 2n. the parcel layer ----------------------------------------------
    _gout = cand.merge(
        contacts[[KEY] + [c for c in contacts.columns
                          if c not in cand.columns and c != KEY]],
        on=KEY, how='left')
    _gout = gpd.GeoDataFrame(_gout, geometry='geometry', crs=WORKING_CRS)
    gpkg_safe(_gout).to_file(os.path.join(out_dir, 'parcels_enriched.gpkg'),
                             layer='ranked_parcels', driver='GPKG')
    gpkg_safe(cl).to_file(os.path.join(out_dir, 'parcels_enriched.gpkg'),
                          layer='all_parcels', driver='GPKG')

    if MAKE_PARCEL_FIGURES:
        parcel_figures(fips, name, county_id, cl, cand, contacts, zc, ztbl,
                       cmp_t, _IN, fig_dir, cext, cg)
        # Guarded on its own: the cover maps are the last thing built and the
        # least load-bearing, and a county whose CDL codes surprise them
        # should still keep its outreach list and its zone map.
        try:
            cover_figures(fips, name, county_id, cl, fig_dir, out_dir,
                          cext, cg)
        except Exception as _e:
            print(f'    cover maps    FAILED: {type(_e).__name__}: {_e}')
    return res


# %% P2 - the parcel figures
# Step 3. Four figures, in the same serif and at the same sizes as the state
# validation set, so a county packet and the state deck read as one document.
def parcel_figures(fips, name, county_id, cl, cand, contacts, zc, ztbl,
                   cmp_t, _IN, fig_dir, cext, cg):
    NZ = len(zc)
    TOTAC = float(zc['rip_crop_ac'].sum())
    CAP = (f'{name} {COUNTY_WORD}, {STATE_NAME} · {NZ:,} parcels above the '
           f'{MIN_RIP_AC:g} ac floor · {TOTAC:,.0f} ac riparian cropland · '
           f'basis {BASIS}, floor {fmt(FLOOR)} ac / sq mi')
    _pct = FuncFormatter(lambda v, p: f'{v:.0f}%')
    _num = FuncFormatter(lambda v, p: f'{v:,.0f}')

    def _bare2(ax, xgrid=True):
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color(MUTED)
        if xgrid:
            ax.set_axisbelow(True)
            ax.xaxis.grid(True, color=NEUTRAL, lw=1.0)
        ax.tick_params(length=0)

    def _save(fig, stem):
        p = os.path.join(fig_dir, f'{county_id}_{stem}.png')
        fig.savefig(p, dpi=FIG_DPI, bbox_inches='tight')
        show_close(fig)

    # --- MAP: the parcels, on the zones -------------------------------------
    xmin, ymin, xmax, ymax = cext.bounds
    pad = 0.03 * max(xmax - xmin, ymax - ymin)
    bb = (xmin - pad, ymin - pad, xmax + pad, ymax + pad)
    fig, ax = plt.subplots(figsize=map_figsize(12, bb, pad_in=3.0))
    bare(ax)
    gpd.GeoSeries([cext], crs=WORKING_CRS).plot(
        ax=ax, facecolor='#f7f7f5', edgecolor='#666666', linewidth=0.8)
    _ctx = grid.cx[bb[0]:bb[2], bb[1]:bb[3]]
    if len(_ctx):
        # the same scale as the county sheet: 1.0 is the cut, everywhere
        _a, _e = cell_raster(_ctx, 'vs_cut')
        ax.imshow(np.ma.masked_invalid(_a), origin='lower', extent=_e,
                  cmap=RAMP_NAME, vmin=0, vmax=CUT_VMAX, alpha=0.30,
                  interpolation='nearest', zorder=1)
    cl.boundary.plot(ax=ax, color='#c9c9c4', linewidth=0.22, zorder=2)
    _leg = [Line2D([], [], color='#c9c9c4', lw=1.0,
                   label=f'Parcel boundary ({len(cl):,})')]
    _c0 = cand[~cand['in_hotspot']]
    _c1 = cand[cand['in_hotspot']]
    if len(_c0):
        _c0.plot(ax=ax, facecolor=SEQ[2], edgecolor='none', zorder=3)
        _leg.append(Patch(facecolor=SEQ[2],
                          label=f'Candidate parcel, outside a zone '
                                f'({len(_c0):,})'))
    if len(_c1):
        _c1.plot(ax=ax, facecolor=HOT, edgecolor='none', zorder=4)
        _leg.append(Patch(facecolor=HOT,
                          label=f'Candidate parcel, inside a zone '
                                f'({len(_c1):,})'))
    if len(zones):
        _ids = split.loc[split['FIPS'] == fips, 'zone_id'].unique()
        _zz = zones.loc[zones.index.isin(_ids)]
        if len(_zz):
            _zz.boundary.plot(ax=ax, color='white', linewidth=4.0, zorder=5)
            _zz.boundary.plot(ax=ax, color=ZONE_INK, linewidth=2.2, zorder=6)
            _leg.append(Line2D([], [], color=ZONE_INK, lw=2.2,
                               label=f'High-density zone ({len(_zz)})'))
            for _zid, _zr in _zz.iterrows():
                _p = _zr.geometry.intersection(cg)
                _p = (_p if not _p.is_empty
                      else _zr.geometry).representative_point()
                ax.annotate(str(_zid), (_p.x, _p.y),
                            fontsize=FS['panel_title'], fontweight='bold',
                            ha='center', va='center', zorder=8,
                            bbox=dict(boxstyle='circle,pad=.28', fc='white',
                                      ec=ZONE_INK))
    ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[1], bb[3])
    legend_headroom(ax, 0.20)
    frame(ax.legend(handles=_leg, loc='upper left', fontsize=FS['legend'],
                    borderpad=0.8, labelspacing=0.5,
                    title='Ranked parcels', title_fontsize=FS['legend_title']))
    ax.set_title(f'{name} {COUNTY_WORD} — Ranked Parcels on the '
                 f'High-Density Surface',
                 fontsize=FS['fig_title'], loc='left', pad=12)
    ax.set_xlabel(CAP, fontsize=FS['caption'], color=INK2)
    _save(fig, 'parcels_ranked')

    # --- FIG 1: by zone ------------------------------------------------------
    if len(ztbl):
        _ZL = list(ztbl['zone'])
        _f = ztbl.set_index('zone').reindex(_ZL)
        fig, ax = plt.subplots(figsize=(12, 0.75 * len(_ZL) + 3.4))
        y = np.arange(len(_ZL))[::-1]
        h = 0.36
        b1 = ax.barh(y + h / 2 + .02, _f['pct_parcels'] * 100, h,
                     color=SEQ[2], label='Percent of parcels', zorder=3)
        b2 = ax.barh(y - h / 2 - .02, _f['pct_acres'] * 100, h,
                     color=SEQ[5], label='Percent of riparian acres', zorder=3)
        for _b, _n in zip(b1, _f['parcels']):
            ax.text(_b.get_width() + .8,
                    _b.get_y() + _b.get_height() / 2,
                    f'{_b.get_width():.1f}%   ({int(_n):,} parcels)',
                    va='center', fontsize=FS['annot'], color=INK2)
        for _b, _a in zip(b2, _f['rip_crop_ac']):
            ax.text(_b.get_width() + .8,
                    _b.get_y() + _b.get_height() / 2,
                    f'{_b.get_width():.1f}%   ({_a:,.0f} ac)',
                    va='center', fontsize=FS['annot'], color=INK2)
        ax.set_yticks(y)
        ax.set_yticklabels(_ZL, fontsize=FS['tick'], color=INK)
        ax.set_xlim(0, max(_f['pct_parcels'].max(),
                           _f['pct_acres'].max()) * 150)
        ax.xaxis.set_major_formatter(_pct)
        _bare2(ax)
        frame(ax.legend(loc='lower right', fontsize=FS['legend'],
                        borderpad=0.8))
        ax.set_title('Riparian Cropland by High-Density Zone',
                     fontsize=FS['fig_title'], loc='left', pad=14)
        ax.set_xlabel(CAP + '\nEach parcel sits in exactly one zone or in '
                            'none, so both series sum to 100 percent.',
                      fontsize=FS['caption'], color=INK2)
        _save(fig, 'zfig1_by_zone')

    # --- FIG 2: water type, inside vs outside --------------------------------
    _SRC = [s for s in USE_SOURCES if f'rip_ac_{s}' in zc.columns]
    if _SRC:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(15, 6.0),
                                       layout='constrained')
        _o = sorted(_SRC, key=lambda s: zc.loc[_IN, f'rip_ac_{s}'].sum())
        yy = np.arange(len(_o))
        vin = [float(zc.loc[_IN, f'rip_ac_{s}'].sum()) for s in _o]
        vout = [float(zc.loc[~_IN, f'rip_ac_{s}'].sum()) for s in _o]
        axa.barh(yy + .19, vin, .36, color=SEQ[5], label='Inside a zone',
                 zorder=3)
        axa.barh(yy - .19, vout, .36, color=MUTED, label='Outside', zorder=3)
        _mx = max(max(vin), max(vout), 1.0)
        for i, (a, b) in enumerate(zip(vin, vout)):
            axa.text(a + _mx * .02, i + .19, f'{a:,.0f}', va='center',
                     fontsize=FS['annot'], color=INK2)
            axa.text(b + _mx * .02, i - .19, f'{b:,.0f}', va='center',
                     fontsize=FS['annot'], color=INK2)
        axa.set_xlim(0, _mx * 1.32)
        axa.xaxis.set_major_formatter(_num)
        axa.set_xlabel('Riparian cropland (acres)')
        axa.set_title('A.  Acreage by Water Type', fontsize=FS['panel_title'],
                      loc='left')
        frame(axa.legend(loc='lower right', fontsize=FS['legend'],
                         borderpad=0.8))
        vsh = [100 * v / max(v + w, 1e-9) for v, w in zip(vin, vout)]
        axb.barh(yy, vsh, .58, color=SEQ[5], zorder=3)
        _ov = 100 * float(zc.loc[_IN, 'rip_crop_ac'].sum()) / max(TOTAC, 1e-9)
        axb.axvline(_ov, color=INK2, lw=1.4, ls=':', zorder=4)
        axb.annotate(f'All water combined: {_ov:.0f}%', (_ov, len(_o) - 0.35),
                     xytext=(6, 0), textcoords='offset points',
                     fontsize=FS['annot'], color=INK2, va='center')
        for i, v in enumerate(vsh):
            axb.text(v + 1.6, i, f'{v:.0f}%', va='center',
                     fontsize=FS['annot'], color=INK2)
        axb.set_xlim(0, max(max(vsh) * 1.32, _ov * 1.5, 10))
        axb.xaxis.set_major_formatter(_pct)
        axb.set_xlabel('Share of that water type inside a zone')
        axb.set_title('B.  Zone Water Type Distribution',
                      fontsize=FS['panel_title'], loc='left')
        for _a in (axa, axb):
            _a.set_yticks(yy)
            _a.set_yticklabels([s.capitalize() for s in _o],
                               fontsize=FS['tick'], color=INK)
            _bare2(_a)
        fig.suptitle('Water Type Inside and Outside the High-Density Zones',
                     fontsize=FS['fig_title'], fontweight='bold', color=INK,
                     x=0.012, ha='left')
        _save(fig, 'zfig2_water_type')

    # --- FIG 3: ownership ----------------------------------------------------
    if 'residency' in zc.columns:
        ZRES = ['In county', 'In state, out of county', 'Out of state',
                'No owner record']
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(15, 5.6),
                                       layout='constrained')
        g = (zc.groupby('residency')['rip_crop_ac'].agg(['size', 'sum'])
             .reindex(ZRES).fillna(0))
        cols = [SEQ[1], SEQ[3], SEQ[5], INK2]
        yy = np.arange(len(ZRES))[::-1]
        axa.barh(yy, g['size'] / max(NZ, 1) * 100, .58, color=cols, zorder=3)
        for i, v in enumerate(g['size']):
            axa.text(v / max(NZ, 1) * 100 + 1.6, yy[i],
                     f'{v / max(NZ, 1) * 100:.1f}%   ({int(v):,})',
                     va='center', fontsize=FS['annot'], color=INK2)
        axa.set_title('A.  All Candidate Parcels',
                      fontsize=FS['panel_title'], loc='left')
        _gz = (zc[_IN].groupby('residency')['rip_crop_ac'].size()
               .reindex(ZRES).fillna(0))
        _nin = max(int(_IN.sum()), 1)
        axb.barh(yy, _gz / _nin * 100, .58, color=cols, zorder=3)
        for i, v in enumerate(_gz):
            axb.text(v / _nin * 100 + 1.6, yy[i],
                     f'{v / _nin * 100:.1f}%   ({int(v):,})', va='center',
                     fontsize=FS['annot'], color=INK2)
        axb.set_title(f'B.  Inside a High-Density Zone ({_nin:,} parcels)',
                      fontsize=FS['panel_title'], loc='left')
        for _a in (axa, axb):
            _a.set_xlim(0, 92)
            _a.set_yticks(yy)
            _a.set_yticklabels(ZRES, fontsize=FS['tick'], color=INK)
            _a.xaxis.set_major_formatter(_pct)
            _bare2(_a)
        axb.set_yticklabels([])
        fig.suptitle('Owner Residency Distribution', fontsize=FS['fig_title'], fontweight='bold',
                     color=INK, x=0.012, ha='left')
        _save(fig, 'zfig3_ownership')

    # --- FIG 4: distribution -------------------------------------------------
    _pos = zc.loc[zc['rip_crop_ac'] > 0, 'rip_crop_ac']
    if len(_pos) > 5:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(15, 6.0),
                                       layout='constrained')
        bins = np.logspace(np.log10(max(float(_pos.min()), 1e-3)),
                           np.log10(float(_pos.max())), 30)
        for _d, _c, _lb in [
                (zc.loc[~_IN & (zc['rip_crop_ac'] > 0), 'rip_crop_ac'],
                 INK2, 'Outside a zone'),
                (zc.loc[_IN & (zc['rip_crop_ac'] > 0), 'rip_crop_ac'],
                 SEQ[5], 'Inside a zone')]:
            if len(_d):
                axa.hist(_d, bins=bins, histtype='step', linewidth=2.4,
                         color=_c, weights=np.ones(len(_d)) / len(_d) * 100,
                         label=_lb, zorder=3)
        axa.set_xscale('log')
        axa.axvline(float(_pos.median()), color=INK2, lw=1.4, ls=':', zorder=5)
        axa.set_xlabel('Riparian cropland per parcel (acres, log scale)')
        axa.set_ylabel('Percent of that group')
        axa.set_title('A.  Distribution of the Resource Across Parcels',
                      fontsize=FS['panel_title'], loc='left')
        frame(axa.legend(loc='upper left', fontsize=FS['legend'],
                         borderpad=0.8))
        axa.xaxis.set_major_formatter(FuncFormatter(
            lambda v, p: f'{v:g}' if v >= 0.1 else f'{v:.2f}'))
        _bare2(axa, xgrid=False)
        axa.yaxis.grid(True, color=NEUTRAL, lw=1.0)
        axa.set_axisbelow(True)

        _s10 = zc['rip_crop_ac'].sort_values(ascending=False) \
                                .reset_index(drop=True)
        cum = _s10.cumsum() / max(TOTAC, 1e-9) * 100
        rank = np.arange(1, len(_s10) + 1)
        axb.plot(rank, cum, color=SEQ[5], lw=2.6, zorder=3)
        axb.fill_between(rank, cum, color=SEQ[5], alpha=.10, zorder=2)
        for k in (100, 250, 500):
            if k < len(_s10):
                axb.plot([k], [cum.iloc[k - 1]], 'o', ms=7, color=INK2,
                         mec=SURFACE, mew=1.6, zorder=5)
                axb.annotate(f'Top {k}  →  {cum.iloc[k - 1]:.0f}%',
                             (k, cum.iloc[k - 1]), xytext=(8, -12),
                             textcoords='offset points',
                             fontsize=FS['annot'], color=INK2, ha='left')
        axb.set_xlim(0, len(_s10)); axb.set_ylim(0, 103)
        axb.set_xlabel('Parcels, ranked by riparian cropland acreage')
        axb.set_ylabel('Cumulative share of riparian acres')
        axb.yaxis.set_major_formatter(_pct)
        axb.xaxis.set_major_formatter(_num)
        axb.set_title('B.  Resource Concentration',
                      fontsize=FS['panel_title'], loc='left')
        _bare2(axb, xgrid=False)
        axb.yaxis.grid(True, color=NEUTRAL, lw=1.0)
        axb.set_axisbelow(True)
        fig.suptitle('Resource Distribution Across '
                     'Parcels', fontsize=FS['fig_title'], fontweight='bold',
                     color=INK, x=0.012, ha='left')
        _save(fig, 'zfig4_distribution')


# %% P2b - what is actually growing there, at its own resolution
# The parcel table reports the crop mix as TOTALS. These two maps report it as
# GROUND - and they colour each overlap polygon by ITS OWN class, not by the
# dominant crop of the parcel it sits in. Colouring a whole parcel by its most
# common crop paints hundreds of acres to describe a strip a few hundred feet
# wide, and it mislabels real ground: a grass buffer inside a corn parcel is
# drawn as corn, which is exactly the ground where a row-crop practice does not
# apply.
#
# Ported from the single-county notebook. What is different here: the corridor
# is cut from the STATEWIDE rip_crop, and a pull-out window on a zone that
# crosses the county line is clipped to THIS county's share of it.

# Families, and one flat colour each. Six of the nine are reserved neutrals -
# grey, water blue - because "developed", "open water" and "other" are context,
# not crops, and giving them a hue would put them in competition with the
# classes the map is actually about.
_FAM_COLOR = {
    'row crop':      '#D9A404',   # amber   - corn, sorghum, cotton, ...
    'soybeans':      '#E06A5E',   # coral   - split out: the single largest
    'wheat / dbl':   '#6A3FC4',   # violet  - wheat and every double-crop
    'grass / hay':   '#2E7D32',   # green   - pasture, hay, alfalfa, sod
    'forest / wet':  '#00A3B8',   # teal    - the NLCD-derived natural cover
    'fallow / idle': '#6E6E6E',
    'developed':     '#9E9E9E',
    'open water':    '#7FA9D0',
    'other':         '#C8C8C8',
}
_FAM_ORDER = ['row crop', 'soybeans', 'wheat / dbl', 'grass / hay',
              'forest / wet', 'fallow / idle', 'open water', 'developed',
              'other']
_F_SOY = {5}
_F_WHEAT = {21, 22, 23, 24, 25, 26, 27, 28, 205,
            225, 226, 228, 236, 237, 238, 239, 240, 241, 254}
_F_GRASS = {36, 37, 58, 59, 60, 62, 176}
_F_FALLOW = {61}
_F_NATURAL = {63, 64, 65, 87, 88, 112, 141, 142, 143, 152, 171, 190, 195}
_F_WATER = {83, 111}
_F_DEV = {121, 122, 123, 124}
_F_BARE = {131}


def _cover_family(c):
    """CDL/NLCD code -> family. Anything unlisted is row crop, which is what
    CDL's long tail of minor crops actually is."""
    if pd.isna(c):
        return 'other'
    c = int(c)
    if c in _F_SOY:
        return 'soybeans'
    if c in _F_WHEAT:
        return 'wheat / dbl'
    if c in _F_GRASS:
        return 'grass / hay'
    if c in _F_FALLOW:
        return 'fallow / idle'
    if c in _F_NATURAL:
        return 'forest / wet'
    if c in _F_WATER:
        return 'open water'
    if c in _F_DEV:
        return 'developed'
    if c in _F_BARE:
        return 'other'
    return 'row crop'


def cover_figures(fips, name, county_id, cl, fig_dir, out_dir, cext, cg):
    """The two cover maps for one county, plus the legend as a CSV."""
    if 'CDL' not in rip_crop.columns:
        print('    cover maps    SKIPPED - no CDL column on rip_crop')
        return

    # These can be moved into fl_config.py whenever they need arguing about;
    # globals().get keeps this working with a config that has never heard of
    # them.
    N_ZOOM = int(globals().get('COVER_N_ZOOM', 4))
    ZOOM_PAD = float(globals().get('COVER_ZOOM_PAD', 0.18))
    MIN_LEG_AC = float(globals().get('COVER_MIN_LEGEND_AC', 1.0))
    HUC_LABELS = bool(globals().get('COVER_HUC_LABELS', True))

    _bx = cext.bounds
    cov = rip_crop.cx[_bx[0]:_bx[2], _bx[1]:_bx[3]]
    cov = clip_to(cov, cext, f'{name} cover') if len(cov) else cov
    if not len(cov):
        print('    cover maps    SKIPPED - no riparian cropland in this county')
        return
    cov = cov.copy()
    cov['ac'] = cov.geometry.area / AC
    cov['CDL'] = pd.to_numeric(cov['CDL'], errors='coerce').astype('Int64')
    cov['cdl_name'] = [CDL_NAMES.get(int(c), f'CDL {int(c)}') if pd.notna(c)
                       else 'unclassified' for c in cov['CDL']]
    cov['family'] = [_cover_family(c) for c in cov['CDL']]

    # THE TABLE THE MAP HAS TO AGREE WITH. Built first, and written out, so the
    # legend is checkable without reading colours off a PNG.
    tab = (cov.groupby(['family', 'cdl_name'], dropna=False)['ac'].sum()
              .rename('ac').reset_index())
    tab['pct'] = (100 * tab['ac'] / max(tab['ac'].sum(), 1e-9)).round(2)
    tab['family'] = pd.Categorical(tab['family'], _FAM_ORDER, ordered=True)
    tab = tab.sort_values(['family', 'ac'], ascending=[True, False])
    tab.round(2).to_csv(os.path.join(out_dir, 'riparian_cover_true_class.csv'),
                        index=False)
    _famtot = tab.groupby('family', observed=True)['ac'].sum()
    _TOT = float(max(_famtot.sum(), 1e-9))

    # this county's zones, from the statewide set
    _ids = (split.loc[split['FIPS'] == fips, 'zone_id'].unique()
            if len(zones) and len(split) else [])
    zz = zones.loc[zones.index.isin(_ids)] if len(zones) else zones

    # PULL-OUT WINDOWS. A zone that crosses the county line is cut to this
    # county's share before the window is taken: the panel is for the
    # biologist who works HERE, and centring it on ground in the next county
    # would be a map of someone else's problem.
    if len(zz):
        _zs = zz.sort_values('rip_crop_ac', ascending=False).head(N_ZOOM)
        wins = []
        for _zid, _zr in _zs.iterrows():
            _here = _zr.geometry.intersection(cg)
            wins.append((int(_zid),
                         (_zr.geometry if _here.is_empty else _here).bounds,
                         _zr))
    else:
        _cs = (grid[grid['FIPS'] == fips]
               .sort_values('rip_per_sqmi', ascending=False).head(N_ZOOM))
        wins = [(None, _r.geometry.bounds, _r) for _, _r in _cs.iterrows()]

    # HUC12s, clipped to the county
    _hucs = None
    if huc12 is not None and len(huc12):
        _hc = next((c for c in huc12.columns if c.lower() == 'huc12'), None)
        if _hc:
            _h = huc12[huc12[_hc].astype(str).str.len() == 12]
            _h = _h.cx[_bx[0]:_bx[2], _bx[1]:_bx[3]]
            if len(_h):
                _hg = gpd.GeoDataFrame(
                    {'huc12': _h[_hc].astype(str).to_numpy()},
                    geometry=_h.geometry.intersection(cext).to_numpy(),
                    crs=WORKING_CRS)
                _hg = _hg[~_hg.geometry.is_empty & (_hg.geometry.area > 0)]
                _hucs = _hg if len(_hg) else None

    _ext_gs = gpd.GeoSeries([cext], crs=WORKING_CRS)
    _flow = (nhd_flow.cx[_bx[0]:_bx[2], _bx[1]:_bx[3]]
             if nhd_flow is not None and len(nhd_flow) else None)
    xmin, ymin, xmax, ymax = _bx
    _COUNTY_AC = float(cov['ac'].sum())

    def _square(b, pad):
        """A square window, so the pull-outs share one scale and a long thin
        zone is not drawn at a different magnification from a round one."""
        x0, y0, x1, y1 = b
        cx_, cy_ = (x0 + x1) / 2, (y0 + y1) / 2
        h = max(x1 - x0, y1 - y0, CELL_M) * (1 + pad) / 2
        return cx_ - h, cy_ - h, cx_ + h, cy_ + h

    def _scalebar(ax, x0, x1, y0):
        w = x1 - x0
        mi = 1.0 if w > 3 * MI_M else 0.5
        bx_, by_ = x0 + 0.06 * w, y0 + 0.055 * w
        ax.plot([bx_, bx_ + mi * MI_M], [by_, by_], color='#111111', lw=3,
                solid_capstyle='butt', zorder=9)
        ax.text(bx_ + mi * MI_M / 2, by_ + 0.018 * w, f'{mi:g} mi',
                ha='center', va='bottom', fontsize=8, color='#111111',
                zorder=9)

    def _draw_cover(ax, sub):
        """One flat fill per family, drawn in _FAM_ORDER so the neutrals go
        down first and a 2 ac corn sliver is never buried under a wetland.
        No hatch and no edge: on a 200 ft ribbon an outline is most of the
        mark, and a darker rim reads as a different colour."""
        for fam in _FAM_ORDER:
            f = sub[sub['family'] == fam]
            if len(f):
                f.plot(ax=ax, facecolor=_FAM_COLOR[fam], edgecolor='none',
                       linewidth=0, zorder=4)

    def _legend_handles(parcel_line=False):
        hs = []
        for fam in _FAM_ORDER:
            rows = tab[(tab['family'] == fam) & (tab['ac'] >= MIN_LEG_AC)]
            if not len(rows):
                continue
            hs.append(Patch(
                facecolor=_FAM_COLOR[fam], edgecolor='none',
                label=f'{fam}  \u2014  {_famtot.get(fam, 0):,.0f} ac '
                      f'({100 * _famtot.get(fam, 0) / _TOT:.1f}%)'))
            for _, r in rows.iterrows():
                hs.append(Patch(facecolor='none', edgecolor='none',
                                label=f"      {r['cdl_name']}  "
                                      f"{r['ac']:,.0f} ac"))
        _tiny = tab[tab['ac'] < MIN_LEG_AC]
        if len(_tiny):
            hs.append(Patch(facecolor='none', edgecolor='none',
                            label=f'      + {len(_tiny)} classes under '
                                  f'{MIN_LEG_AC:g} ac '
                                  f'({_tiny["ac"].sum():,.1f} ac)'))
        if _flow is not None and len(_flow):
            hs.append(Line2D([0], [0], color='#cfe0ee', lw=1.4,
                             label='NHD flowline'))
        if parcel_line:
            hs.append(Line2D([0], [0], color='#b9b9b4', lw=1.2,
                             label='parcel boundary (zoom panels)'))
        return hs

    # Wrapped to three short lines. As one long line it ran past the legend
    # column and collided with panel A's title.
    LEG_TITLE = ('Cover in the riparian corridor\n'
                 'colour = family\n'
                 'indented = its classes')
    _SUB = (f'{_COUNTY_AC:,.0f} ac inside {RIP_BUFFER_FT:.0f} ft of water')

    # ---- FIGURE: the county, plus the pull-out panels ---------------------
    fig = plt.figure(figsize=(23, 11))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.55, 1.55, 0.82, 1.18, 1.18],
                          hspace=0.10, wspace=0.05)
    axo = fig.add_subplot(gs[:, 0:2])
    _ext_gs.plot(ax=axo, facecolor='#fcfcfb', edgecolor='#666666',
                 linewidth=0.9)
    if _flow is not None and len(_flow):
        _flow.plot(ax=axo, color='#cfe0ee', linewidth=0.35, zorder=2)
    _draw_cover(axo, cov)
    if len(zz):
        zz.boundary.plot(ax=axo, color='#111111', linewidth=1.8, zorder=6)
    for _i, (_zid, _b, _r) in enumerate(wins):
        wx0, wy0, wx1, wy1 = _square(_b, ZOOM_PAD)
        axo.add_patch(plt.Rectangle(
            (wx0, wy0), wx1 - wx0, wy1 - wy0, fill=False, edgecolor='#111111',
            linewidth=1.3, linestyle=(0, (2, 2)), zorder=7))
        axo.annotate(chr(65 + _i), (wx1, wy1), fontsize=12, fontweight='bold',
                     ha='center', va='center', zorder=8,
                     bbox=dict(boxstyle='square,pad=.3', fc='white',
                               ec='#111111', ls=':', lw=1.1))
    axo.set_xlim(xmin, xmax); axo.set_ylim(ymin, ymax); axo.set_axis_off()
    _scalebar(axo, xmin, xmax, ymin)
    axo.set_title(f'{name} {COUNTY_WORD} \u2014 riparian corridor cover, every '
                  f'CDL / NLCD class at its own resolution\n'
                  f'{_SUB} \u00b7 black outlines = high-density zones \u00b7 '
                  f'dotted boxes = the panels at right',
                  fontsize=12.5, loc='left')

    axl = fig.add_subplot(gs[:, 2]); axl.set_axis_off()
    axl.legend(handles=_legend_handles(parcel_line=True), loc='upper left',
               bbox_to_anchor=(-0.02, 1.0), fontsize=7.8, framealpha=0,
               borderpad=.4, labelspacing=.34, handlelength=1.6,
               handleheight=1.1, title=LEG_TITLE, title_fontsize=8.6,
               alignment='left')

    _cells = [gs[0, 3], gs[0, 4], gs[1, 3], gs[1, 4]]
    for _i, (_zid, _b, _r) in enumerate(wins[:len(_cells)]):
        ax = fig.add_subplot(_cells[_i])
        wx0, wy0, wx1, wy1 = _square(_b, ZOOM_PAD)
        _bb = shp_box(wx0, wy0, wx1, wy1)
        _ext_gs.plot(ax=ax, facecolor='#fcfcfb', edgecolor='none')
        _pw = cl.cx[wx0:wx1, wy0:wy1]
        if len(_pw):
            _pw.boundary.plot(ax=ax, color='#b9b9b4', linewidth=.45, zorder=3)
        if _flow is not None:
            _fw = _flow.cx[wx0:wx1, wy0:wy1]
            if len(_fw):
                _fw.plot(ax=ax, color='#9ecae1', linewidth=.8, zorder=2)
        _draw_cover(ax, cov.cx[wx0:wx1, wy0:wy1])
        if len(zz):
            _zw = zz.cx[wx0:wx1, wy0:wy1]
            if len(_zw):
                gpd.GeoSeries(_zw.geometry.intersection(_bb).to_numpy(),
                              crs=WORKING_CRS).boundary.plot(
                    ax=ax, color='#111111', linewidth=1.8, zorder=6)
        # the panel's own frame, dotted, to tie it back to its box on the
        # overview - the only cue that says which window this is
        for _sp in ax.spines.values():
            _sp.set_visible(True); _sp.set_linestyle((0, (2, 2)))
            _sp.set_linewidth(1.1); _sp.set_color('#111111')
        ax.set_xlim(wx0, wx1); ax.set_ylim(wy0, wy1)
        ax.set_xticks([]); ax.set_yticks([])
        # geopandas names the axes after the CRS ('Easting [metre]'), which on
        # a tickless pull-out is a label for something the reader cannot read
        ax.set_xlabel(''); ax.set_ylabel('')
        _lab = (f'{chr(65 + _i)} \u00b7 zone {_zid} \u00b7 '
                f'{_r["rip_crop_ac"]:,.0f} ac \u00b7 '
                f'{_r["vs_state"]:.1f}x state'
                if _zid is not None else
                f'{chr(65 + _i)} \u00b7 densest cell \u00b7 '
                f'{_r["rip_per_sqmi"]:,.0f} ac/sq mi')
        ax.set_title(_lab, fontsize=9, loc='left', pad=3)
        _scalebar(ax, wx0, wx1, wy0)

    _P1 = os.path.join(fig_dir, f'{county_id}_riparian_cover_true_class.png')
    fig.savefig(_P1, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
    show_close(fig)

    # ---- FIGURE: the cover on the parcel and HUC12 frame ------------------
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.7, 1.0], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    _ext_gs.plot(ax=ax, facecolor='#fcfcfb', edgecolor='#666666',
                 linewidth=0.9)
    # EVERY cropland parcel, not just the candidates, so the map shows the
    # real fabric rather than a pre-filtered subset.
    cl.plot(ax=ax, facecolor='#f0f0ee', edgecolor='#c9c9c4', linewidth=.22,
            zorder=2)
    if _flow is not None and len(_flow):
        _flow.plot(ax=ax, color='#cfe0ee', linewidth=.35, zorder=3)
    if _hucs is not None:
        _hucs.boundary.plot(ax=ax, color='#111111', linewidth=1.5, zorder=5)
    _draw_cover(ax, cov)          # corridor colours last, on top of both frames
    # ZONES DOTTED here, solid on the figure above. Two different jobs: there
    # the zone is the subject, here it is one of three frames competing for
    # the same ink, and a dotted line reads as "boundary of a class" rather
    # than as another watershed edge.
    if len(zz):
        zz.boundary.plot(ax=ax, color=ZONE_INK, linewidth=1.8,
                         linestyle=(0, (3, 2)), zorder=7)
    if HUC_LABELS and _hucs is not None:
        for _, _hr in _hucs.iterrows():
            _c = _hr.geometry.representative_point()
            ax.annotate(str(_hr['huc12'])[-6:], (_c.x, _c.y), fontsize=6.5,
                        color='#333333', ha='center', va='center', zorder=8,
                        bbox=dict(boxstyle='round,pad=.15', fc='white',
                                  ec='none', alpha=.72))
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_axis_off()
    _scalebar(ax, xmin, xmax, ymin)
    _nh = len(_hucs) if _hucs is not None else 0
    ax.set_title(f'{name} {COUNTY_WORD} \u2014 riparian corridor cover on the '
                 f'parcel and HUC12 frame\n{_SUB} \u00b7 {len(cl):,} parcels '
                 f'(grey) \u00b7 {_nh} HUC12 subwatersheds (black) \u00b7 '
                 f'{len(zz)} high-density zones (dotted)',
                 fontsize=12.5, loc='left')

    _h3 = _legend_handles()
    _h3.append(Patch(facecolor='#f0f0ee', edgecolor='#c9c9c4',
                     label=f'parcel ({len(cl):,})'))
    if _nh:
        _h3.append(Line2D([0], [0], color='#111111', lw=1.5,
                          label=f'HUC12 subwatershed ({_nh})'))
    if len(zz):
        _h3.append(Line2D([0], [0], color=ZONE_INK, lw=1.8,
                          ls=(0, (3, 2)),
                          label=f'high-density zone ({len(zz)})'))
    axl = fig.add_subplot(gs[0, 1]); axl.set_axis_off()
    axl.legend(handles=_h3, loc='upper left', bbox_to_anchor=(-0.04, 1.0),
               fontsize=8, framealpha=0, borderpad=.4, labelspacing=.36,
               handlelength=1.6, handleheight=1.1, title=LEG_TITLE,
               title_fontsize=9, alignment='left')
    _P2 = os.path.join(fig_dir,
                       f'{county_id}_riparian_cover_parcels_huc12.png')
    fig.savefig(_P2, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
    show_close(fig)

    # WHAT THE DOMINANT-CLASS SHORTCUT WOULD HAVE COST. `rip_n_cdl` is on the
    # parcel table already, so this is free - and it is the number that says
    # whether the true-class map was worth drawing on THIS county.
    _mixed = ''
    if 'rip_n_cdl' in cl.columns:
        _m = cl['rip_n_cdl'] > 0
        if int(_m.sum()):
            _mx = (cl.loc[_m, 'rip_n_cdl'] > 1)
            _mixed = (f'; {int(_mx.sum()):,} of {int(_m.sum()):,} parcels '
                      f'({_mx.mean():.0%}) hold more than one cover class')
    print(f'    cover maps    {len(tab)} classes in '
          f'{int(tab["family"].nunique())} families over '
          f'{_COUNTY_AC:,.0f} ac{_mixed}')
    print(f'                  {os.path.basename(_P1)}, '
          f'{os.path.basename(_P2)}, riparian_cover_true_class.csv')


# %% P3 - run it, county by county
# Step 4. One county at a time, and a failure in one is recorded in the
# manifest rather than ending the run. A 120-county pass that dies on county
# 40 because one parcel file has an odd schema is the failure mode this
# structure exists to prevent.
PARCEL_COLS = ['n_parcels', 'n_candidates', 'parcel_rip_ac', 'parcels_in_zone',
               'pct_parcels_in_zone', 'owners', 'pct_contactable',
               'pct_out_of_state', 'frontage_mi', 'rho_parcels_vs_density',
               'rho_score_vs_size', 'size_ratio_in_vs_out']
for _c in PARCEL_COLS:
    if _c not in man.columns:
        man[_c] = np.nan

if not RUN_PASS_C:
    print('\nRUN_PASS_C = False - the parcel work is skipped. The zone map, '
          'the county\nsheets and every state figure above are unaffected.')
elif not man['has_parcels'].any():
    rule('PASS C - PARCELS')
    print('  No county in the manifest has a parcel file, so there is nothing '
          'to run.')
    print(f'  That is the expected state for most of {STATE_NAME} and it is '
          f'the point of the\n  three-pass split: Passes A and B produced the '
          f'zone map for all {len(man)} counties\n  without a single parcel '
          f'licence.')
    print(f'  The buy-list is rank_by_zone_resource in {MANIFEST_CSV}:')
    print('    ' + ', '.join(man.sort_values('rank_by_zone_resource')
                             .head(10)['NAME']))
else:
    rule('PASS C - PARCELS')
    _todo = man[man['has_parcels']].copy()
    if PARCEL_COUNTIES != 'auto':
        _want = {str(c).zfill(3) for c in PARCEL_COUNTIES}
        _todo = _todo[_todo['FIPS'].isin(_want)]
    print(f'  {len(_todo)} county/counties with a parcel file, of '
          f'{len(man)} in the manifest')
    print(f'  every other county stops at the zone footprint - which is the '
          f'cost argument:\n  Passes A and B cost nothing in parcel '
          f'licensing.\n')
    import time as _tp
    _t0p = _tp.time()
    _done, _failed = 0, 0
    for _, _r in _todo.iterrows():
        print(f'  {_r["NAME"]} {COUNTY_WORD} ({_r["FIPS"]}) - '
              f'{os.path.basename(str(_r["parcel_path"]))}')
        try:
            _res = pass_c(_r['FIPS'], _r['NAME'], _r['county_id'], _r)
            _m = man['FIPS'] == _r['FIPS']
            for _k, _v in _res.items():
                if _k in PARCEL_COLS:
                    man.loc[_m, _k] = _v
            if _res.get('n_parcels'):
                man.loc[_m, 'pct_parcels_in_zone'] = round(
                    100 * _res['parcels_in_zone'] / _res['n_parcels'], 1)
            _done += 1
            print(f'    -> {os.path.join(C_DIR, _r["county_id"])}')
        except Exception as _e:
            _failed += 1
            man.loc[man['FIPS'] == _r['FIPS'], 'error'] = f'pass C: {_e}'
            print(f'    FAILED: {type(_e).__name__}: {_e}')
            print(f'    The zone outputs for this county are unaffected; only '
                  f'the parcel\n    files are missing. The reason is in the '
                  f'manifest `error` column.')
        print()
    man.to_csv(MANIFEST_CSV, index=False)
    man.to_csv(os.path.join(S_OUT, 'county_summary.csv'), index=False)
    print(f'  {_done} county packet(s) built'
          + (f', {_failed} failed' if _failed else '')
          + f' in {_tp.time() - _t0p:.0f}s')
    if _done:
        _p = man[man['n_parcels'].notna()]
        print(f'\n  {int(_p["n_parcels"].sum()):,} parcels examined; '
              f'{int(_p["parcels_in_zone"].sum()):,} sit inside a '
              f'high-density zone,\n  and '
              f'{int(_p["n_candidates"].sum()):,} clear the {MIN_RIP_AC:g} ac '
              f'materiality floor and reach the outreach list.')
        print(f'  Those are different counts of different things: the zone '
              f'test is about WHERE\n  a parcel sits, the floor is about how '
              f'much riparian cropland is ON it.')
        print(_p[['NAME', 'n_parcels', 'n_candidates', 'parcel_rip_ac',
                  'parcels_in_zone', 'pct_parcels_in_zone', 'owners',
                  'pct_contactable', 'pct_out_of_state', 'frontage_mi']]
              .to_string(index=False))


# %% [markdown]
# ## What was produced
#
# Every county gets a folder. What is in it depends on one thing — whether
# that county's parcel file is on disk — and the split is deliberate:
#
# | | every county | counties with parcel data |
# |---|---|---|
# | `output/` | `cells.gpkg`, `zones.gpkg`, `summary.json` | + `outreach_parcels_ranked.csv`, `riparian_cover_by_parcel_cdl.csv`, `huc12_ranked.csv`, `zone_*.csv`, `parcels_enriched.gpkg` |
# | `figures/` | the two-panel county sheet | + the ranked-parcel map and four briefing figures |
# | `validation/` | — | `zone_in_vs_out.csv`, `parcel_crosscheck.csv`, `parcel_floor_sweep.csv` |
#
# **That asymmetry is the finding, not a gap.** Passes A and B produced a zone
# map for all 120 counties without a single parcel licence, because the
# hotspot method never looks at a parcel. Pass C is the only part that costs
# money per county, and `rank_by_zone_resource` in the manifest is the order
# to spend it in: it ranks counties by how much riparian cropland already sits
# inside a zone, measured without any parcel data at all.
#
# The county packets are produced by THIS file now, against the statewide
# surface, so a parcel's neighbourhood density is measured over a full 3-mile
# band and its `zone_id` matches `zone_county_split.csv` — the biologist in
# the next county is looking at the same numbered zone. Running
# `parcel_priority_pipeline.ipynb` for a single county still works and is
# still the place for the grid-origin sensitivity test, which is a
# per-county question.

# %% F1 - closing summary
rule('SUMMARY')

print(f'  basis                 : {BASIS}'
      + ('  (the notebook\'s own rule)' if BASIS == 'county' and FLOOR == 0
         else ''))
if REPORT_MASK:
    print('  reporting mask        : priority area - the surface is '
          'statewide, every\n                          number below is over '
          'PA ground and is NOT comparable\n                          with a '
          'statewide run')
print(f'  ratio                 : {HOT_MULTIPLE:g}x')
print(f'  floor                 : {FLOOR:.2f} ac/sq mi ({FLOOR_MODE})')
print(f'  analysed              : {ANALYSIS_AC / 640:,.0f} sq mi across '
      f'{int((man["n_cells"] > 0).sum())} counties')
print(f'  riparian cropland     : {grid["rip_ac"].sum():,.0f} ac')
print(f'  hot ground            : {grid.loc[grid["hot"], "land_ac"].sum() / 640:,.0f} '
      f'sq mi holding '
      f'{100 * grid.loc[grid["hot"], "rip_ac"].sum() / max(grid["rip_ac"].sum(), 1e-9):.0f}% '
      f'of it')
print(f'  zones                 : {len(zones):,}'
      + (f', {int((zones["n_counties"] > 1).sum()):,} crossing a county line'
         if len(zones) else ''))
print(f'  counties with zones   : {int((man["n_zones"] > 0).sum())} of {len(man)}')
print(f'  counties with parcels : {int(man["has_parcels"].sum())}')
if 'n_candidates' in man.columns and man['n_candidates'].notna().any():
    print(f'  Pass C packets        : '
          f'{int(man["n_candidates"].notna().sum())} built, '
          f'{int(man["n_candidates"].sum()):,} candidate parcels, '
          f'{int(man["parcels_in_zone"].sum()):,} inside a zone')

_nocov = man[man['status'] == 'no_source_data']
if len(_nocov):
    print(f'\n  {len(_nocov)} counties could not be analysed - the cached '
          f'layers do not reach them.')
    print(f'  They are in the manifest with status = no_source_data, and they '
          f'are NOT in the\n  state or region averages. See the note in B1 for '
          f'the one-time fix.')

_none = man[man['status'].str.startswith('no_zones')]
if len(_none):
    print(f'\n  {len(_none)} analysed counties produced no zones. That is a '
          f'result, not a failure -\n  it is where the programme should not '
          f'spend a biologist\'s time:')
    print('  ' + ', '.join(_none['NAME'].head(20)) +
          (' ...' if len(_none) > 20 else ''))

print(f'\n  files:')
print(f'    {MANIFEST_CSV}')
print(f'    {SURFACE_GPKG}  (grid, grid_scored, rip_crop, zones)')
print(f'    {os.path.join(S_OUT, "county_summary.csv")}')
print(f'    {os.path.join(S_OUT, "zone_summary.csv")}      one row per zone: '
      f'acres, density, per-source acres,\n'
      f'                                       stream miles, NHD counts, '
      f'named waters, neighbours,\n'
      f'                                       HUC12 membership, cover by '
      f'CDL group')
print(f'    {os.path.join(S_OUT, "zone_cover_by_cdl.csv")} one row per zone '
      f'x CDL class - the unrolled cover')
print(f'    {os.path.join(S_OUT, "zone_huc12.csv")}        one row per zone '
      f'x HUC12 - acres and share')
print(f'    {os.path.join(S_OUT, "zone_county_split.csv")}  the same, split '
      f'per county')
print(f'    {os.path.join(S_VAL, "floor_sweep.csv")}')
print(f'    {os.path.join(S_VAL, "basis_comparison.csv")}')
print(f'    {os.path.join(S_VAL, "band_truncation.csv")}')
print(f'    {os.path.join(S_VAL, "equivalence_checks.txt")}')
print(f'    {S_FIG}/F1..F5*.png')
print(f'    {C_DIR}/<FIPS>_<Name>/output|figures|validation/')

# Step 2. If the run was built on scratch disk, copy it somewhere durable NOW.
# Printing the command and trusting someone to notice it at the end of a
# six-hour run is the weakest link in any pipeline. Does nothing when
# SCRATCH_DIR is None, which is the normal local case.
copy_to_final()
