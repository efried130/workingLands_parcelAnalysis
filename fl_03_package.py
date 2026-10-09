# %% [markdown]
# # fl_03_package.py — the stakeholder handoff
#
# `state_surface.gpkg` is a working file that carries four
# layers, two of which share the same geometry and only one of which has a
# verdict on it (`grid` is the Pass A density surface with no `hot` column;
# `grid_scored` is the same cells after Pass B). 
#
# This script reads both of those files IN PLACE and writes one package
# GeoPackage that a person who has never seen the project can open and read.
#
#
#     python fl_03_package.py          # or paste as one Colab cell
#
# Reads   S_OUT/state_surface.gpkg   (zones, grid_scored)
#         CACHE                      (counties, huc12, watersheds, pa_geo,
#                                     nhd_flow, nhd_wb, nhd_area, wetlands)
#         S_OUT/zone_summary.csv     (optional, for the dictionary's coverage)
# Writes  WORK_DIR/_state/package/<st>_riparian_hotspots.gpkg
#                                  /README.txt
#                                  /data_dictionary.csv
#                                  /overview.png
#                                  /<st>_riparian_hotspots_package.zip

# %% P1 - setup
from fl_common import *            # noqa: F401,F403 - and fl_config through it

import csv
import datetime
import glob
import json
import shutil
import sqlite3
import textwrap
import xml.etree.ElementTree as ET
import zipfile

# ---- the two knobs -------------------------------------------------------
ADD_PRIORITY_RANK = True   # the presentational ordering described above
CONTEXT_BUFFER_MI = BAND_MI  # streams/wetlands are clipped to this distance
#                              from the zones. BAND_MI (3 mi) is the same
#                              neighbourhood the analysis reasons about, so
#                              the context shown is the context used.

PKG_DIR = os.path.join(WORK_DIR, '_state', 'package')
PKG_NAME = f'{HOME_STATE.lower()}_riparian_hotspots.gpkg'
PKG = os.path.join(PKG_DIR, PKG_NAME)
os.makedirs(PKG_DIR, exist_ok=True)

rule(f'{STATE_NAME}')
print(f'  surface  {SURFACE_GPKG}')
print(f'  cache    {CACHE}')
print(f'  package  {PKG}')

if not os.path.exists(SURFACE_GPKG):
    raise SystemExit(f'\n  {SURFACE_GPKG} is not there. Run '
                     f'fl_02_statewide.py first.')

# A GeoPackage is a SQLite file and this one has been copied across a FUSE
# mount. Ask it whether it is intact BEFORE spending minutes reading it, so a
# truncated copy gives one sentence instead of a pyogrio traceback.
for _p, _what in ((SURFACE_GPKG, 'the surface'), (CACHE, 'the cache')):
    if not os.path.exists(_p):
        continue
    _c = sqlite3.connect(f'file:{_p}?mode=ro', uri=True)
    _q = _c.execute('PRAGMA quick_check(1)').fetchone()[0]
    _c.close()
    if _q != 'ok':
        raise SystemExit(f'\n  {_what} is corrupt: {str(_q)[:90]}\n'
                         f'  {_p}\n  Run tools/colab_surface_check.py.')
print('  integrity  both files pass PRAGMA quick_check')

if os.path.exists(PKG):
    os.remove(PKG)          # a package is rebuilt whole, never appended to


# %% P2 - helpers
def _first_col(frame, *names):
    """The first column whose name matches, ignoring case. Cached layers come
    from different agencies and spell the same field four ways."""
    low = {str(c).lower(): c for c in frame.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def pkg_layers():
    """Layer names in the package, sorted.

    pyogrio's list_layers returns a 2-D numpy array, not a list of tuples, so
    `sorted(list_layers(p))` raises "truth value of an array ... is
    ambiguous". Unpacking the rows is fine; sorting them is not. One helper
    so that mistake is made in one place and already made."""
    # layer_styles is a GeoPackage housekeeping table, not data. pyogrio
    # lists it like any other table, and leaving it in puts styleQML rows in
    # the data dictionary and the layer list a stakeholder reads.
    return sorted(str(n) for n, _ in list_layers(PKG)
                  if str(n) not in ('layer_styles', 'qgis_projects'))


def put(gdf, layer, keep=None, note=''):
    """Write one layer and say what went in. Empty frames are SKIPPED rather
    than written as an empty layer - a GIS showing an empty layer reads as a
    data problem, and a layer that is simply absent reads as 'not included'."""
    if gdf is None or not len(gdf):
        print(f'  {layer:<22} SKIPPED - nothing to write')
        return None
    g = gdf.copy()
    if keep:
        cols = [c for c in keep if c in g.columns]
        g = g[cols + ['geometry']] if 'geometry' not in cols else g[cols]
    g = g.to_crs(WORKING_CRS)
    gpkg_safe(g).to_file(PKG, layer=layer, driver='GPKG')
    print(f'  {layer:<22} {len(g):>9,} features · {len(g.columns) - 1:>2} '
          f'fields{("  " + note) if note else ""}')
    return g


# %% P3 - the analysis layers, copied unchanged
rule('RESULT')

zones = gpd.read_file(SURFACE_GPKG, layer='zones').to_crs(WORKING_CRS)
grid_scored = gpd.read_file(SURFACE_GPKG,
                            layer='grid_scored').to_crs(WORKING_CRS)
print(f'  read {len(zones):,} zones and {len(grid_scored):,} scored cells')

# the hot cells, and NOTHING else from the grid. The full surface is a
# different conversation: a stakeholder asked where to go, and 40,000 cells of
# which 3,200 are hot answers that question by burying it.
hot = grid_scored[grid_scored['hot'].astype(bool)].copy() \
    if 'hot' in grid_scored.columns else grid_scored.iloc[:0].copy()
print(f'  of those cells, {len(hot):,} are hot '
      f'({100 * len(hot) / max(len(grid_scored), 1):.1f}%)')

if ADD_PRIORITY_RANK and 'rip_crop_ac' in zones.columns:
    # PRESENTATIONAL ONLY - see the header. A dense ordering of a column the
    # analysis already produced, so a GIS can colour 150 polygons by one
    # field without the viewer having to pick breaks.
    _ac = pd.to_numeric(zones['rip_crop_ac'], errors='coerce')
    # method='first', not 'min': equal-acreage zones are common (two zones of
    # the same few cells) and 'min' would hand three zones rank 15 and leave
    # 16 and 17 unused, which reads as a bug in a attribute table. Ties break
    # by row order, which is zone_id order - arbitrary but stable, so the
    # same surface always produces the same ranks.
    # A PLAIN int64, not a nullable Int64: gpkg_safe converts nullable dtypes
    # to object on the way out, and the field would land in the GeoPackage as
    # TEXT - sorting 10 before 9 in every GIS that opens it.
    _r = _ac.rank(ascending=False, method='first')
    zones['priority_rank'] = _r.fillna(len(zones) + 1).astype('int64')
    zones['priority_tier'] = np.where(
        _r <= 0.20 * len(zones), '1 - highest acreage',
        np.where(_r <= 0.50 * len(zones), '2 - high acreage',
                 '3 - moderate acreage'))
    print(f'  added priority_rank / priority_tier (ordering of rip_crop_ac, '
          f'presentational)')

put(zones, 'hotspot_zones', note='<- the headline layer')
put(hot, 'hotspot_cells',
    keep=['cell', 'county', 'rip_ac', 'land_ac', 'rip_per_sqmi',
          'own_per_sqmi', 'nbhd_rip_ac', 'nbhd_land_ac', 'band_cover',
          'reference', 'threshold_applied', 'ratio', 'vs_cut', 'p_state',
          'hot', 'geometry'],
    note='<- the cells the zones are made of')


# %% P4 - the reference geography, from the cache
rule('CONTEXT')

# Every context layer comes from the preprocessing cache. Without it the
# package is still valid - it just has no reference geography, which is worth
# one clear sentence rather than nine "SKIPPED" lines.
if not os.path.exists(CACHE):
    print(f'  the cache is not on this machine:\n    {CACHE}\n'
          f'  The package will carry the result layers and NO context. Run '
          f'fl_01_preprocess.py,\n  or point CACHE at the copy in the project '
          f'folder, then re-run this.')

counties = load_cached('counties', required=False) \
    if os.path.exists(CACHE) else None
if counties is not None and len(counties):
    # the per-county result, so clicking a county answers "does my county
    # have any?" without opening a CSV
    _man = os.path.join(S_OUT, 'county_summary.csv')
    if os.path.exists(_man):
        _m = pd.read_csv(_man, dtype={'FIPS': str})
        _keep = [c for c in ('FIPS', 'NAME', 'status', 'n_zones', 'zone_ac',
                             'rip_crop_ac', 'avg_county', 'floor_applied',
                             'threshold_applied') if c in _m.columns]
        # the boundaries already carry NAME; a second copy lands as
        # NAME_summary and reads like a data error in the attribute table
        _keep = [c for c in _keep
                 if c == 'FIPS' or _first_col(counties, c) is None]
        _fc = _first_col(counties, 'FIPS', 'GEOID', 'COUNTYFP')
        if _fc and 'FIPS' in _keep:
            counties = counties.merge(_m[_keep], left_on=_fc,
                                      right_on='FIPS', how='left',
                                      suffixes=('', '_summary'))
            print(f'  joined {len(_keep) - 1} county_summary field(s) onto '
                  f'the county boundaries')
    put(counties, 'counties', note='<- with the per-county result joined')

huc12 = (load_cached('huc12', required=False)
         if os.path.exists(CACHE) else None)
if huc12 is not None and len(huc12):
    _id = _first_col(huc12, 'huc12', 'HUC12', 'huc_12')
    _nm = _first_col(huc12, 'name', 'NAME', 'hu_12_name')
    put(huc12, 'huc12_subwatersheds',
        keep=[c for c in (_id, _nm) if c] + ['geometry'])

watersheds = (load_cached('watersheds', required=False)
              if os.path.exists(CACHE) else None)
put(watersheds, 'priority_huc8', note='<- the HUC8 priority watersheds')

pa_geo = (load_cached('pa_geo', required=False)
          if os.path.exists(CACHE) else None)
put(pa_geo, 'priority_area', note='<- the ACF priority area, clipped to state')

# ---- streams and wetlands, clipped to the zones' working range -----------
# The whole state is 316,000 stream segments and 193,000 wetland polygons.
# Nobody opening a hotspot map needs the streams 80 miles from the nearest
# zone, and shipping them makes the package too large to email. CLIPPED, not
# sampled: every segment within CONTEXT_BUFFER_MI of a zone is here whole.
if len(zones):
    _near = shapely.union_all(
        shapely.buffer(np.asarray(zones.geometry.values, dtype=object),
                       CONTEXT_BUFFER_MI * MI_M))
    _mask = gpd.GeoSeries([_near], crs=WORKING_CRS)
    print(f'\n  clipping context to {CONTEXT_BUFFER_MI:g} mi around the '
          f'zones ({_near.area / AC / 640:,.0f} sq mi)')
    for _name, _layer, _note in (('nhd_flow', 'streams_near_zones',
                                  'streams within range'),
                                 ('nhd_area', 'rivers_near_zones',
                                  'mapped river surfaces'),
                                 ('nhd_wb', 'waterbodies_near_zones',
                                  'lakes and ponds'),
                                 ('wetlands', 'wetlands_near_zones',
                                  'NWI wetlands')):
        _g = (load_cached(_name, required=False)
              if os.path.exists(CACHE) else None)
        if _g is None or not len(_g):
            print(f'  {_layer:<22} SKIPPED - not in the cache')
            continue
        _sub = _g.iloc[shapely.STRtree(
            np.asarray(_g.geometry.values, dtype=object)).query(
                _near, predicate='intersects')]
        put(gpd.clip(_sub, _mask), _layer, note=f'<- {_note}')


# %% P5 - symbology, embedded where it is read automatically
rule('SYMBOLOGY')

# A SEQUENTIAL RAMP IS ONE HUE, LIGHT TO DARK. The quantity being shown -
# riparian cropland density - has a magnitude and no meaningful midpoint, so a
# diverging or rainbow ramp would invent structure that is not in the data.
# One hue, four steps, darkest where there is most.
_ORANGE = ['254,230,206', '253,174,107', '230,85,13', '127,39,4']
_OUTLINE = '64,64,64'


def _fill_symbol(name, rgb, outline=_OUTLINE, outline_w='0.26', alpha='1'):
    s = ET.Element('symbol', {'type': 'fill', 'name': name, 'alpha': alpha,
                              'clip_to_extent': '1', 'force_rhr': '0'})
    lay = ET.SubElement(s, 'layer', {'class': 'SimpleFill', 'enabled': '1',
                                     'pass': '0', 'locked': '0'})
    opt = ET.SubElement(lay, 'Option', {'type': 'Map'})
    for k, v in (('color', f'{rgb},255'), ('style', 'solid'),
                 ('outline_color', f'{outline},255'),
                 ('outline_style', 'solid'), ('outline_width', outline_w),
                 ('outline_width_unit', 'MM'), ('joinstyle', 'bevel'),
                 ('offset', '0,0'), ('offset_unit', 'MM')):
        ET.SubElement(opt, 'Option', {'name': k, 'type': 'QString',
                                      'value': v})
    return s


def _line_symbol(name, rgb, width='0.4'):
    s = ET.Element('symbol', {'type': 'line', 'name': name, 'alpha': '1',
                              'clip_to_extent': '1', 'force_rhr': '0'})
    lay = ET.SubElement(s, 'layer', {'class': 'SimpleLine', 'enabled': '1',
                                     'pass': '0', 'locked': '0'})
    opt = ET.SubElement(lay, 'Option', {'type': 'Map'})
    for k, v in (('line_color', f'{rgb},255'), ('line_width', width),
                 ('line_width_unit', 'MM'), ('line_style', 'solid'),
                 ('capstyle', 'round'), ('joinstyle', 'round')):
        ET.SubElement(opt, 'Option', {'name': k, 'type': 'QString',
                                      'value': v})
    return s


def graduated_qml(attr, breaks, colours, labels):
    """QML for a graduated fill. breaks is n+1 edges for n classes."""
    q = ET.Element('qgis', {'version': '3.28.0',
                            'styleCategories': 'Symbology'})
    r = ET.SubElement(q, 'renderer-v2',
                      {'type': 'graduatedSymbol', 'attr': attr,
                       'graduatedMethod': 'GraduatedColor',
                       'symbollevels': '0', 'forceraster': '0',
                       'enableorderby': '0'})
    rng = ET.SubElement(r, 'ranges')
    for i in range(len(colours)):
        ET.SubElement(rng, 'range',
                      {'lower': f'{breaks[i]:.6f}',
                       'upper': f'{breaks[i + 1]:.6f}',
                       'symbol': str(i), 'label': labels[i],
                       'render': 'true'})
    sym = ET.SubElement(r, 'symbols')
    for i, c in enumerate(colours):
        sym.append(_fill_symbol(str(i), c))
    return ('<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
            + ET.tostring(q, encoding='unicode'))


def single_qml(symbol):
    q = ET.Element('qgis', {'version': '3.28.0',
                            'styleCategories': 'Symbology'})
    r = ET.SubElement(q, 'renderer-v2', {'type': 'singleSymbol',
                                         'symbollevels': '0',
                                         'forceraster': '0',
                                         'enableorderby': '0'})
    s = ET.SubElement(r, 'symbols')
    s.append(symbol)
    return ('<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
            + ET.tostring(q, encoding='unicode'))


def _breaks(series, n=4):
    """Quantile edges, widened so the top and bottom classes cannot clip a
    value. Quantiles rather than equal intervals because density is skewed -
    equal intervals would put 140 of 150 zones in the lowest class."""
    v = pd.to_numeric(series, errors='coerce').dropna()
    if not len(v):
        return list(range(n + 1))
    e = [float(v.quantile(i / n)) for i in range(n + 1)]
    e[0], e[-1] = float(v.min()) - 1e-6, float(v.max()) + 1e-6
    for i in range(1, n):                     # strictly increasing
        e[i] = max(e[i], e[i - 1] + 1e-6)
    return e


STYLES = []
if len(zones) and 'rip_per_sqmi' in zones.columns:
    _b = _breaks(zones['rip_per_sqmi'])
    STYLES.append(('hotspot_zones', 'Riparian cropland density',
                   graduated_qml('rip_per_sqmi', _b, _ORANGE,
                                 [f'{_b[i]:,.0f} - {_b[i + 1]:,.0f} '
                                  f'ac/sq mi' for i in range(4)])))
if len(hot) and 'ratio' in hot.columns:
    _b = _breaks(hot['ratio'])
    STYLES.append(('hotspot_cells', 'How far above the local reference',
                   graduated_qml('ratio', _b, _ORANGE,
                                 [f'{_b[i]:,.1f}x - {_b[i + 1]:,.1f}x '
                                  f'the reference' for i in range(4)])))
STYLES += [
    ('counties', 'Boundaries only',
     single_qml(_fill_symbol('0', '0,0,0', outline='120,120,120',
                             outline_w='0.2', alpha='0'))),
    ('huc12_subwatersheds', 'Boundaries only',
     single_qml(_fill_symbol('0', '0,0,0', outline='140,160,180',
                             outline_w='0.16', alpha='0'))),
    ('priority_area', 'Priority area outline',
     single_qml(_fill_symbol('0', '0,0,0', outline='40,90,140',
                             outline_w='0.6', alpha='0'))),
    ('priority_huc8', 'Priority watersheds',
     single_qml(_fill_symbol('0', '0,0,0', outline='70,130,170',
                             outline_w='0.4', alpha='0'))),
    ('streams_near_zones', 'Streams',
     single_qml(_line_symbol('0', '70,130,180', '0.3'))),
]

# QGIS reads a `layer_styles` table inside the GeoPackage on load and applies
# any row with useAsDefault = 1. This is the ONLY way to ship symbology
# inside a single .gpkg; a sidecar .qml is only read for file-based layers.
_have = set(pkg_layers())
_con = sqlite3.connect(PKG)
_con.execute("""
CREATE TABLE IF NOT EXISTS layer_styles (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    f_table_catalog TEXT, f_table_schema TEXT, f_table_name TEXT,
    f_geometry_column TEXT, styleName TEXT, styleQML TEXT, styleSLD TEXT,
    useAsDefault BOOLEAN, description TEXT, owner TEXT, ui TEXT,
    update_time DATETIME DEFAULT CURRENT_TIMESTAMP)""")
_n = 0
for _layer, _label, _qml in STYLES:
    if _layer not in _have:
        continue
    ET.fromstring(_qml.split('\n', 1)[1])      # refuse to store broken XML
    _con.execute(
        'INSERT INTO layer_styles (f_table_catalog, f_table_schema, '
        'f_table_name, f_geometry_column, styleName, styleQML, styleSLD, '
        'useAsDefault, description, owner, ui) '
        "VALUES ('', '', ?, 'geom', ?, ?, '', 1, ?, '', '')",
        (_layer, _label, _qml, f'{STATE_NAME} riparian hotspots — {_label}'))
    _n += 1
_con.commit()
_con.close()
print(f'  embedded {_n} QGIS style(s) in the package (layer_styles, '
      f'useAsDefault)')
print(f'  ArcGIS Pro ignores these - the README gives the field, breaks and '
      f'colours to set by hand')


# %% P6 - the data dictionary
rule('DATA DICTIONARY')

LAYER_NOTES = {
    # (layer, field) wins over the by-name table below. The county layer
    # reuses names like zone_ac and rip_crop_ac with a COUNTY meaning, and
    # inheriting the zone wording there would be actively wrong.
    ('counties', 'FIPS'): 'County FIPS code.',
    ('counties', 'NAME'): 'County name.',
    ('counties', 'status'): 'What the run found here: ok, no_zones (analysed, '
                            'nothing cleared the bar - a result, not a gap), '
                            'or no_source_data (the cached layers do not '
                            'reach this county).',
    ('counties', 'n_zones'): 'Hotspot zones whose largest share is in this '
                             'county.',
    ('counties', 'zone_ac'): 'Total hotspot-zone acres in this county.',
    ('counties', 'rip_crop_ac'): 'Riparian cropland acres in this county, '
                                 'all of it, in or out of a zone.',
    ('counties', 'avg_county'): 'This county\'s average riparian cropland '
                                'density, acres per square mile - the '
                                'reference its own cells were judged against '
                                'under the county basis.',
    ('counties', 'floor_applied'): 'Whether this county was limited by the '
                                   'absolute floor rather than by the ratio.',
    ('counties', 'threshold_applied'): 'The density a cell in this county '
                                       'had to clear.',
    ('huc12_subwatersheds', 'huc12'): 'HUC12 subwatershed code.',
    ('huc12_subwatersheds', 'name'): 'Subwatershed name.',
}

FIELD_NOTES = {
    'zone_id': 'Zone identifier. Joins to zone_summary.csv, zone_huc12.csv '
               'and zone_cover_by_cdl.csv.',
    'primary_FIPS': 'County FIPS holding the largest share of the zone.',
    'primary_county': 'Name of that county.',
    'counties': 'Every county the zone touches, semicolon separated.',
    'n_counties': 'How many counties the zone spans.',
    'n_cells': 'One-mile cells dissolved into this zone.',
    'zone_ac': 'Total ground area of the zone, acres.',
    'pct_inside': 'Share of the zone inside the analysis extent.',
    'rip_crop_ac': 'RIPARIAN CROPLAND acres in the zone - cropland within '
                   f'{RIP_BUFFER_FT:.0f} ft of mapped water. The quantity '
                   'the whole analysis is about.',
    'rip_per_sqmi': 'Riparian cropland acres per square mile of the zone. '
                    'The density figure.',
    'vs_state': 'This zone\'s density as a multiple of the statewide figure.',
    'basis': 'What each cell was compared against: county, state or region '
             '(HUC8).',
    'report_mask': 'Whether figures were restricted to a priority area '
                   '("none" if not).',
    'reference_used': 'The local reference density the zone\'s cells were '
                      'judged against, acres per square mile.',
    'floor_ac_per_sqmi': 'The absolute minimum density for the whole run. No '
                         'cell is hot below it however it compares locally.',
    'threshold_applied': 'The density a cell actually had to clear: the '
                         'greater of the ratio cut and the floor.',
    'median_ratio': 'Median of the zone\'s cells\' density / reference.',
    'max_ratio': 'The highest such multiple in the zone.',
    'n_stream_segments': 'NHD flowline segments intersecting the zone.',
    'stream_mi': 'Total mapped stream miles in the zone.',
    'n_rivers': 'NHD mapped river surfaces (wide channels) in the zone.',
    'n_waterbodies': 'Lakes and ponds in the zone.',
    'n_wetlands': 'NWI wetland polygons in the zone.',
    'named_waters': 'Named streams and waterbodies, for talking to a '
                    'landowner about "their" creek.',
    'n_huc12': 'HUC12 subwatersheds the zone intersects.',
    'primary_huc12': 'The HUC12 holding most of the zone.',
    'primary_huc12_name': 'Its name.',
    'pct_in_primary_huc12': 'Share of the zone in that subwatershed.',
    'huc12_ids': f'Up to {ZONE_HUC12_TOP} HUC12 codes, semicolon separated.',
    'rip_cover': 'Plain-language cover mix of the riparian cropland.',
    'rip_cdl_dominant': 'The single most common CDL code.',
    'rip_cdl_dominant_name': 'Its crop name.',
    'rip_n_cdl': 'How many distinct CDL classes appear.',
    'rip_cdl_acres': 'Acres per CDL class, semicolon separated.',
    'rip_seq_class': f'What the fields have BEEN over the last '
                     f'{CROP_SEQ_YEARS} CSB years - continuous row crop, '
                     f'row crop to grass, continuous pasture, fallow in '
                     f'sequence and so on. GRAZING IS INFERRED FROM A GRASS '
                     f'CLASS, NEVER OBSERVED.',
    'n_zones_within_band': f'Other zones within {ZONE_NEIGHBOR_MI:g} miles - '
                           f'how clustered the work is.',
    'nearest_zone_mi': 'Distance to the nearest other zone, miles.',
    'priority_rank': 'PRESENTATIONAL. Rank by rip_crop_ac, 1 = most acres. '
                     'Not an analytical output; an ordering of a column the '
                     'analysis already produced, added for symbology.',
    'priority_tier': 'PRESENTATIONAL. Top 20% / next 30% / remainder of '
                     'priority_rank. Same caveat.',
    'cell': 'Cell identifier within the statewide one-mile lattice.',
    'county': 'County FIPS the cell falls in.',
    'rip_ac': 'Riparian cropland acres in this cell alone.',
    'land_ac': 'Ground area of the cell inside the extent.',
    'own_per_sqmi': 'This cell\'s own density, ignoring its neighbours.',
    'nbhd_rip_ac': f'Riparian cropland acres within {BAND_MI:g} miles.',
    'nbhd_land_ac': f'Ground acres within {BAND_MI:g} miles.',
    'band_cover': 'Share of that neighbourhood that is real ground rather '
                  'than beyond the extent. Low near the state edge.',
    'reference': 'The reference density this cell was compared against.',
    'ratio': 'Neighbourhood density / reference. The comparison itself.',
    'vs_cut': 'Density minus threshold_applied. Positive means hot.',
    'p_state': 'Percentile of this cell\'s density statewide.',
    'hot': 'True for every cell in this layer, by construction.',
}

def _fields_of(layer):
    """Field names, without reading the geometry. read_info is a pyogrio
    thing and fl_common falls back to fiona on a machine without it, so this
    degrades to a one-row read rather than failing the dictionary."""
    try:
        return [str(f) for f in read_info(PKG, layer=layer).get('fields', [])]
    except Exception:
        try:
            _one = gpd.read_file(PKG, layer=layer, rows=1)
            return [c for c in _one.columns if c != 'geometry']
        except Exception:
            return []


rows = []
for _layer in pkg_layers():
    for _f in _fields_of(_layer):
        rows.append({'layer': _layer, 'field': _f,
                     'meaning': LAYER_NOTES.get(
                         (_layer, _f),
                         FIELD_NOTES.get(
                             _f, 'From the source dataset, carried through '
                                 'unchanged.'))})
DICT_CSV = os.path.join(PKG_DIR, 'data_dictionary.csv')
pd.DataFrame(rows).to_csv(DICT_CSV, index=False)
_undoc = sum(1 for r in rows if r['meaning'].startswith('From the source'))
print(f'  {len(rows):,} fields across {len({r["layer"] for r in rows})} '
      f'layers -> {os.path.basename(DICT_CSV)}')
print(f'  {len(rows) - _undoc} described here, {_undoc} carried through from '
      f'the source agencies')


# %% P7 - the README a person reads first
_surf = {}
if os.path.exists(SURFACE_JSON):
    with open(SURFACE_JSON) as _f:
        _surf = json.load(_f)

_tot_ac = float(pd.to_numeric(zones.get('rip_crop_ac'),
                              errors='coerce').sum()) if len(zones) else 0.0

# THE FLOOR IS NOT A CONFIG VALUE. FLOOR_MODE is, but the floor itself is
# DERIVED inside fl_02 and does not exist in this namespace - reading it from
# the config would print the mode's name where a number belongs. Every zone
# row carries the floor the run actually used, which is the authoritative
# copy, so take it from there and say so plainly if it is missing.
_fl = pd.to_numeric(zones.get('floor_ac_per_sqmi'), errors='coerce').dropna() \
    if len(zones) else pd.Series(dtype=float)
FLOOR_TXT = (f'{_fl.iloc[0]:,.2f} acres per square mile' if len(_fl)
             else 'see floor_ac_per_sqmi in the zones layer')
BASIS_TXT = {'county': 'county average', 'state': 'statewide average',
             'region': 'HUC8 regional average'}.get(BASIS, f'{BASIS} average')
RULE_TXT = (f'{FLOOR_MODE} ({FLOOR_KNEE_RULE} rule)' if FLOOR_MODE == 'knee'
            else str(FLOOR_MODE))

_P1 = textwrap.fill(
    f"Where cropland sits closest to water, concentrated enough that a visit "
    f"is worth someone's day. {len(zones):,} hotspot zones across "
    f"{STATE_NAME}, holding {_tot_ac:,.0f} acres of riparian cropland in "
    f"total.", 74, initial_indent='  ', subsequent_indent='  ')
_P2 = textwrap.fill(
    f'"Riparian cropland" means cropland within {RIP_BUFFER_FT:.0f} feet of '
    f"water mapped by the NHD - streams, rivers, lakes and NWI wetlands. "
    f"Cropland comes from USDA's Crop Sequence Boundaries, which are field "
    f"boundaries rather than pixels.", 74,
    initial_indent='  ', subsequent_indent='  ')
_READ = f"""{STATE_NAME.upper()} RIPARIAN CROPLAND HOTSPOTS
{'=' * (len(STATE_NAME) + 28)}
Prepared {datetime.date.today().isoformat()} · {PKG_NAME}

WHAT THIS SHOWS
{_P1}

{_P2}

HOW A ZONE WAS CHOSEN
  The state is cut into {CELL_MI:g}-mile cells. Each cell gets the riparian
  cropland density of its {BAND_MI:g}-mile neighbourhood - a field rep's
  working range, not the cell alone. A cell is HOT when all three hold:

    1. its neighbourhood density is at least {HOT_MULTIPLE:g}x the local
       reference (the {BASIS_TXT});
    2. its own cell density is at least that reference too, so a cell cannot
       ride on its neighbours alone;
    3. its density clears an absolute floor of
       {FLOOR_TXT},
       so being {HOT_MULTIPLE:g}x a quiet neighbourhood does not on its own
       qualify.

  Touching hot cells are dissolved into one zone. Zones are therefore made of
  merged one-mile squares and look blocky - that is honest to the method
  rather than a drafting shortcut. Do not read a zone edge as a property line.

WHAT IS IN THE FILE
  hotspot_zones          the result. One row per zone, ~40 attributes.
  hotspot_cells          the one-mile cells the zones are made of, with each
                         cell's density, reference and ratio.
  counties               county boundaries with the per-county result joined
                         (how many zones, how many acres, or why none).
  huc12_subwatersheds    subwatershed boundaries.
  priority_huc8          the HUC8 priority watersheds.
  priority_area          the priority-area boundary, clipped to the state.
  *_near_zones           streams, rivers, waterbodies and wetlands within
                         {CONTEXT_BUFFER_MI:g} miles of a zone - clipped to
                         keep the file sendable. Nothing is thinned or
                         sampled: everything in range is here whole.

  data_dictionary.csv    every field in every layer, in plain language.
  overview.png           what the map should look like once styled.

OPENING IT
  QGIS      Layer > Add Layer > Add Vector Layer, pick the .gpkg, add all.
            The symbology is embedded and applies on load. Drag
            hotspot_zones to the top if the draw order needs it.

  ArcGIS Pro  Map > Add Data, pick the .gpkg, add the layers. ArcGIS does not
            read QGIS styles, so set them once:
              hotspot_zones  Symbology > Graduated Colors, field
                             rip_per_sqmi, 4 quantile classes, a single-hue
                             sequential ramp (light to dark - Oranges works),
                             outline dark grey 0.5 pt.
              hotspot_cells  the same on the field `ratio`.
              everything else  hollow fill, grey outline.

  Coordinates are EPSG:5070 (NAD83 Conus Albers), an equal-area projection -
  acres measured in it are correct. Both packages reproject on the fly, so
  this will line up with anything else you add.

TWO CAVEATS
  GRAZING IS INFERRED, NEVER OBSERVED.

  A ZONE IS A PLACE TO LOOK, NOT A LIST OF PARCELS. These footprints are
  parcel-invariant by design - they reference only cropland geometry and the
  state's extent, never ownership. Parcel-level work happens separately and
  only where a county parcel file exists.

PROVENANCE
  surface built   {_surf.get('built_at', 'see state_surface.json')}
  reference basis {BASIS}
  ratio cut       {HOT_MULTIPLE:g}x
  floor mode      {RULE_TXT}
  floor applied   {FLOOR_TXT}
  cells / band    {CELL_MI:g} mi cell · {BAND_MI:g} mi neighbourhood
  crop sequence   newest {CROP_SEQ_YEARS} CSB years

  Method and code: see the project README.
"""
READ_TXT = os.path.join(PKG_DIR, 'README.txt')
with open(READ_TXT, 'w') as _f:
    _f.write(_READ)
print(f'\n  wrote {os.path.basename(READ_TXT)} '
      f'({len(_READ.splitlines())} lines)')


# %% P8 - an overview PNG, so they can see what it should look like
try:
    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.set_facecolor(SURFACE)
    fig.patch.set_facecolor(SURFACE)
    if counties is not None and len(counties):
        counties.boundary.plot(ax=ax, color=MUTED, linewidth=0.4, zorder=1)
    if pa_geo is not None and len(pa_geo):
        pa_geo.boundary.plot(ax=ax, color='#285a8c', linewidth=1.2, zorder=2)
    if len(zones):
        _v = pd.to_numeric(zones['rip_per_sqmi'], errors='coerce')
        zones.assign(_v=_v).plot(ax=ax, column='_v', cmap='Oranges',
                                 scheme=None, linewidth=0.3,
                                 edgecolor='#404040', legend=True, zorder=3,
                                 legend_kwds={'label': 'riparian cropland '
                                                       'ac / sq mi',
                                              'shrink': 0.55})
    ax.set_axis_off()
    ax.set_title(f'{STATE_NAME} — riparian cropland hotspots\n'
                 f'{len(zones):,} zones · {_tot_ac:,.0f} ac riparian '
                 f'cropland · {BASIS} basis · {HOT_MULTIPLE:g}x cut',
                 fontsize=FS['fig_title'], loc='left', pad=14)
    PNG = os.path.join(PKG_DIR, 'overview.png')
    fig.savefig(PNG, dpi=150, bbox_inches='tight', facecolor=SURFACE)
    show_close(fig)
    print(f'  wrote {os.path.basename(PNG)}')
except Exception as _e:
    PNG = None
    print(f'  overview.png skipped ({type(_e).__name__}: {_e})')


# %% P9 - one zip to send, and the integrity check on the way out
rule('PACKAGE')

_c = sqlite3.connect(f'file:{PKG}?mode=ro', uri=True)
assert _c.execute('PRAGMA quick_check(1)').fetchone()[0] == 'ok', \
    'the package GeoPackage did not come out intact - rebuild it'
_c.close()

ZIP = os.path.join(PKG_DIR, f'{HOME_STATE.lower()}_riparian_hotspots'
                            f'_package.zip')
with zipfile.ZipFile(ZIP, 'w', zipfile.ZIP_DEFLATED) as _z:
    for _f in [PKG, READ_TXT, DICT_CSV] + ([PNG] if PNG else []):
        _z.write(_f, os.path.basename(_f))

print(f'  layers    {", ".join(pkg_layers())}')
print(f'  gpkg      {os.path.getsize(PKG) / 1e6:,.1f} MB')
print(f'  zip       {os.path.getsize(ZIP) / 1e6:,.1f} MB  <- send this one')
print(f'\n  {PKG_DIR}')
for _f in sorted(os.listdir(PKG_DIR)):
    print(f'    {_f:<42} {os.path.getsize(os.path.join(PKG_DIR, _f)) / 1e6:>8,.1f} MB')

copy_to_final()
