# ---
# jupytext:
#   text_representation:
#     extension: .py
#     format_name: percent
# ---

# %% [markdown]
# # `fl_common.py` — the shared machinery
#
# Imports, file-finding, the loader, and the geometry helpers both scripts use.
#
# **Nothing in this file makes a decision.** Every number and path it uses
# comes from `fl_config.py`. If you find yourself wanting to edit a value here,
# the value belongs in the config instead.

# %% M0 - libraries
import os
import re
import glob
import json
import hashlib
import zipfile
import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import shapely
from shapely import union_all, set_precision, make_valid
from shapely.geometry import box as shp_box
from scipy.spatial import cKDTree
from scipy import ndimage

from fl_config import *            # noqa: F401,F403 - every switch lives there

try:
    from pyogrio import list_layers, read_info, read_dataframe
    _HAVE_PYOGRIO = True
except ImportError:                 # fiona fallback: slower, no pushdown
    _HAVE_PYOGRIO = False
    from fiona import listlayers as _fl

    def list_layers(p):
        return [(n, None) for n in _fl(p)]

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 220)

# Running as a plain script with no display, `plt.show()` either blocks or
# warns. Every figure here is saved to disk anyway, so the file-only backend is
# the correct default; a notebook keeps its inline backend.
_INTERACTIVE = (COLAB or hasattr(__builtins__, '__IPYTHON__')
                or 'ipykernel' in str(matplotlib.get_backend()).lower())
if not _INTERACTIVE:
    try:
        matplotlib.use('Agg')
    except Exception:
        pass


def setup_matplotlib():
    """Apply the config's type and colour tokens to every figure in the run."""
    plt.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE,
        'font.family': 'serif', 'font.serif': SERIF_STACK,
        'mathtext.fontset': 'dejavuserif',
        'font.size': FS['base'],
        'axes.edgecolor': MUTED, 'axes.labelcolor': INK2,
        'axes.labelsize': FS['axis_label'],
        'xtick.color': INK2, 'ytick.color': INK2,
        'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'],
        'legend.fontsize': FS['legend'],
        'legend.title_fontsize': FS['legend_title'],
        'axes.titlecolor': INK, 'axes.titlesize': FS['panel_title'],
        'axes.titleweight': 'bold',
    })


def show_close(fig):
    """Every figure ends here: shown if SHOW_FIGS, and always closed.

    Closing matters more than it looks - a 120-county run makes several
    hundred figures, and matplotlib keeps every one of them alive until it is
    closed. That is how a long run turns into a memory problem that looks like
    a geometry problem."""
    if SHOW_FIGS:
        try:
            plt.show()
        except Exception:
            pass
    plt.close(fig)


# %% M1 - printing
def say(msg, indent=0):
    print(' ' * indent + str(msg))


def rule(title):
    """A banner, so a long run stays readable in a terminal."""
    print(f'\n{"=" * 74}\n{title}\n{"=" * 74}')


def fmt(v, width=1):
    """Enough decimals to be true. A floor of 1.12 must not print as '1'."""
    return (f'{v:,.0f}' if abs(v) >= 10
            else f'{v:,.{width + 1}f}'.rstrip('0').rstrip('.'))


def safe_name(s):
    """A file-system-safe county name."""
    return re.sub(r'[^A-Za-z0-9]+', '_', str(s)).strip('_') or 'Unknown'


# %% M2 - finding the input files
# Downloads arrive as folders as often as as files, with download IDs in the
# names and inconsistent capitalisation, so a single literal glob is the one
# thing that reliably does NOT find them.
def find_all(pattern, root=None):
    """Every path matching `pattern`, tried four ways: exact at the top level,
    exact anywhere beneath, then both again case-insensitively. The first form
    that matches anything wins."""
    def nocase(p):
        return ''.join(f'[{c.lower()}{c.upper()}]' if c.isalpha() else c
                       for c in p)

    root = root or RAW_DIR
    seen, out = set(), []
    for pat in (pattern, os.path.join('**', pattern),
                nocase(pattern), os.path.join('**', nocase(pattern))):
        for h in glob.glob(os.path.join(root, pat), recursive=True):
            if h not in seen:
                seen.add(h)
                out.append(h)
        if out:
            break
    return sorted(out)


def find_one(pattern, root=None, required=True, label=''):
    """Resolve a glob to exactly one path.

    Source files carry download IDs, so a hard-coded path goes stale the moment
    someone re-downloads. Matching a pattern and saying which file won is the
    difference between a clear message and a silent wrong-file run."""
    hits = find_all(pattern, root)
    if len(hits) == 1:
        return hits[0]
    if not hits:
        if required:
            raise FileNotFoundError(
                f'{label or pattern}: nothing matched {pattern!r} in '
                f'{root or RAW_DIR}')
        return None
    print(f'  {label or pattern}: {len(hits)} matched, using the newest '
          f'({os.path.basename(max(hits, key=os.path.getmtime))})')
    return max(hits, key=os.path.getmtime)


def resolve_dataset(spec, label=''):
    """First pattern in a registry row that matches anything wins. An absolute
    path that exists is taken as-is, which is how you point at a file outside
    the project folder."""
    pats = spec['pattern']
    pats = [pats] if isinstance(pats, str) else list(pats)
    for p in pats:
        if os.path.isabs(p) and os.path.exists(p):
            return p
        hit = find_one(p, required=False, label=label)
        if hit:
            return hit
    return None


def show_raw_dir(root=None, max_rows=60):
    """List what is actually in raw_data, two levels deep.

    Printed whenever a required pattern misses. A pattern that does not match
    is almost always a naming guess rather than a missing download, and the
    only thing that settles which is seeing the folder."""
    root = root or RAW_DIR
    if not os.path.isdir(root):
        print(f'  {root} does not exist')
        return
    rows = []
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if os.path.isdir(p) and not name.lower().endswith(('.gdb', '.gpkg')):
            kids = sorted(os.listdir(p))[:6]
            rows.append(f'  {name}/')
            rows += [f'      {k}' for k in kids]
            if len(os.listdir(p)) > 6:
                rows.append(f'      ... {len(os.listdir(p)) - 6} more')
        else:
            rows.append(f'  {name}')
    print(f'\nwhat is actually in {root}:')
    for r in rows[:max_rows]:
        print(r)
    if len(rows) > max_rows:
        print(f'  ... {len(rows) - max_rows} more entries')


def check_raw_dir():
    """raw_data is an INPUT, so it is checked and never created.

    Creating it was a real trap: a PROJECT_DIR with a typo in it produced an
    empty folder and then twelve "pattern matched nothing" rows blaming the
    registry for a path problem."""
    if os.path.isdir(RAW_DIR) and os.listdir(RAW_DIR):
        return
    why = 'does not exist' if not os.path.isdir(RAW_DIR) else 'is EMPTY'
    near = os.path.abspath(RAW_DIR)
    while near != os.path.dirname(near) and not os.path.isdir(near):
        near = os.path.dirname(near)
    print(f'\n{"=" * 74}\nraw_data {why}\n{"=" * 74}')
    print(f'  looked in : {RAW_DIR}')
    print(f'  deepest folder that does exist: {near}')
    try:
        kids = sorted(os.listdir(near))[:25]
        print(f'  it contains: {", ".join(kids) if kids else "(nothing)"}')
    except Exception:
        pass
    raise FileNotFoundError(
        f'{RAW_DIR} {why}. This is a PATH problem, not a pattern problem - '
        f'fix PROJECT_DIR in fl_config.py, or put the downloads in that '
        f'folder. Nothing can run until it holds the source files.')


DATA_EXT = ('.gpkg', '.shp', '.gdb', '.geojson', '.json', '.tif', '.tiff',
            '.img', '.vrt', '.csv')


def open_path(path, verbose=True):
    """The path GDAL should actually open.

    Downloads often arrive zipped. GDAL reads straight out of a zip through its
    virtual filesystem, so nothing has to be extracted first. This finds the
    dataset inside and returns the `/vsizip/...` path; anything that is not a
    zip comes back unchanged."""
    if not str(path).lower().endswith('.zip'):
        return path
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
    except Exception as e:
        print(f'  could not read {os.path.basename(path)} as a zip: {e}')
        return path
    # a .gdb is a directory, so match on the folder rather than a file entry
    gdb = sorted({n[:n.lower().index('.gdb') + 4] for n in names
                  if '.gdb/' in n.lower()})
    hits = ([n for n in names if n.lower().endswith('.gpkg')]
            or gdb
            or [n for n in names if n.lower().endswith('.shp')]
            or [n for n in names if n.lower().endswith(DATA_EXT)
                and not n.lower().endswith(('.dbf', '.shx', '.prj', '.cpg'))])
    if not hits:
        print(f'  {os.path.basename(path)}: no dataset inside '
              f'({len(names)} entries, e.g. {names[:3]})')
        return path
    inner = sorted(hits, key=len)[0]
    if verbose:
        print(f'  {os.path.basename(path)} -> reading {inner!r} in place')
    return f'/vsizip/{path}/{inner}'


def stage_zip(path, dest=None, verbose=True):
    """Extract a zip once and return the dataset path inside it.

    The fallback for when reading in place fails - some GeoPackages inside a
    zip need random access the virtual filesystem cannot give cheaply."""
    dest = dest or os.path.join(CLIPPED_DIR, '_unzipped',
                                os.path.splitext(os.path.basename(path))[0])
    if not os.path.isdir(dest) or not os.listdir(dest):
        os.makedirs(dest, exist_ok=True)
        with zipfile.ZipFile(path) as z:
            z.extractall(dest)
        if verbose:
            print(f'  extracted {os.path.basename(path)} -> {dest}')
    for root, dirs, files in os.walk(dest):
        for n in sorted(dirs) + sorted(files):
            if n.lower().endswith(('.gpkg', '.gdb', '.shp')):
                return os.path.join(root, n)
    return dest


# %% M3 - reading, normalising, caching
def resolve_cols(available, wanted):
    """Case-insensitive column matching.

    NHD High Resolution ships lowercase field names; CSB and NWI use mixed
    case. pyogrio's `columns=` is case-sensitive, so asking for 'FTYPE' in a
    layer that calls it 'ftype' silently returns a frame without it."""
    lut = {c.lower(): c for c in available}
    return [lut[w.lower()] for w in wanted if w.lower() in lut]


def newest_spec(spec):
    """Decode a registry row's `newest` into (alias, pattern, n_years).

    ('CDL', pattern)      -> keep the newest year only
    ('CDL', pattern, 3)   -> keep the newest three

    THE ONLY PLACE THAT DECODES THIS TUPLE. It grew a third element for
    the crop sequence and normalise() was still unpacking two, which
    raised on the first cached csb load and stopped the whole run. Four
    callers indexing a tuple by hand is how that happens; one accessor
    is how it stops. Returns None when the row has no `newest`."""
    n = (spec or {}).get('newest')
    if not n:
        return None
    return (n[0], n[1], int(n[2]) if len(n) > 2 else 1)


def normalise(g, name, verbose=True):
    """Give the newest year-stamped column its stable alias (CDL2025 -> CDL).

    THIS RUNS ON THE CACHE PATH TOO, and that is the whole point. A cache
    written before the rename existed holds the RAW column name. Load it
    without this and `csb` arrives with no `CDL` column at all - which does not
    fail, it silently skips every crop table and map. A normalisation that runs
    on only one of two load paths is a trap."""
    spec = DATASETS.get(name)
    _ns = newest_spec(spec)
    if g is None or _ns is None:
        return g
    want, pat = _ns[0], _ns[1]
    if want in g.columns:
        return g
    hits = cdl_year_cols(g, pat)
    if hits:
        if verbose:
            print(f'  {name:<12} using {hits[-1]!r} as {want!r}'
                  + (f' (also present: {", ".join(hits[:-1])})'
                     if len(hits) > 1 else ''))
        # A COPY, NOT A RENAME. The crop sequence needs the years under
        # their own names, and `CDL` is an alias for the newest one that
        # everything written before the sequence existed still reads.
        # Renaming gave the newest year away to get the alias.
        g = g.copy()
        g[want] = g[hits[-1]]
        return g
    if verbose:
        print(f'  {name:<12} NO {want} COLUMN and nothing matching {pat!r}. '
              f'Crop tables will be skipped.')
    return g


def cdl_year_cols(frame, pat=r'CDL\d{4}'):
    """Year-stamped CDL columns on a frame, OLDEST FIRST.

    One place that knows the naming, so the preprocessing read, Pass A's
    carry-through and the sequence classifier cannot disagree about which
    columns are years. Sorted, so [-1] is always the newest and
    [-N:] is always the newest N."""
    if frame is None:
        return []
    return sorted(c for c in frame.columns
                  if re.fullmatch(pat, str(c), re.I))


def cache_layers(path=None):
    """The set of layer names already in the cache."""
    path = path or CACHE
    if not os.path.exists(path):
        return set()
    try:
        return {n for n, _ in list_layers(path)}
    except Exception:
        return set()


def gpkg_safe(gdf):
    """pandas nullable dtypes (string, Int64, boolean) have no GeoPackage
    equivalent and raise on write. Convert them to object with real None."""
    out = gdf.copy()
    for c in out.columns:
        if c != 'geometry' and str(out[c].dtype) in ('string', 'Int64',
                                                     'boolean'):
            out[c] = out[c].astype(object).where(out[c].notna(), None)
    return out


# %% M3b - provenance: has the source changed since the cache was built?
# WHY THIS EXISTS. `RUN_PREPROCESS = 'auto'` means "build only what is
# missing", and without a record of WHAT each cached layer was built from, a
# re-downloaded source file is not missing - it is simply never read again.
# The cache goes stale silently, and every number downstream is computed from
# last month's data while the new file sits in raw_data/ untouched.
#
# The manifest is a small JSON sidecar next to the cache. If it is absent or
# unreadable every layer is treated as stale, which is the safe direction.
CACHE_MANIFEST = os.path.splitext(CACHE)[0] + '.cache.json'

# A SIDECAR THAT HAS TO TRAVEL WITH THE CACHE. The GeoPackage moves
# between scratch and the project folder, so the manifest has to be
# findable from either side. It is not: a cache restored from the
# project folder without its sidecar has no record of what built it,
# every layer reads as stale, and a run that should have rebuilt one
# layer rebuilds all nine. Two halves to the fix - the end of fl_01
# copies the sidecar with the GeoPackage, and a READ falls back to the
# sidecar beside the other candidate cache location, which covers a
# GeoPackage someone copied by hand on its own.
CACHE_MANIFEST_CANDIDATES = [
    os.path.splitext(os.path.join(_d, CACHE_NAME))[0] + '.cache.json'
    for _d in (CLIPPED_DIR, FINAL_CLIPPED_DIR)]
_MANIFEST_SAID = set()


def source_fingerprint(path):
    """What a source file looked like when we read it.

    A File Geodatabase is a DIRECTORY, so its own mtime says nothing about the
    tables inside it - walk it and take the newest member and the total size.
    Size alone is not enough (an edit can preserve it) and mtime alone is not
    enough (a copy resets it), so both are recorded."""
    if path is None or not os.path.exists(path):
        return None
    if os.path.isdir(path):
        mt, sz, n = 0.0, 0, 0
        for root, _d, files in os.walk(path):
            for f in files:
                try:
                    st = os.stat(os.path.join(root, f))
                except OSError:
                    continue
                mt = max(mt, st.st_mtime)
                sz += st.st_size
                n += 1
        return {'path': path, 'mtime': round(mt, 3), 'size': sz, 'members': n}
    st = os.stat(path)
    return {'path': path, 'mtime': round(st.st_mtime, 3), 'size': st.st_size}


def read_cache_manifest(verbose=True):
    """Provenance for the cached layers: the sidecar beside the cache in
    use, or failing that the one beside the other candidate copy.

    THE FALLBACK IS NOT A CONVENIENCE. Without it, reading the cache
    from a different folder than the one it was written in silently
    discards every layer's provenance and rebuilds the lot. Says so
    once when it falls back, because reading a record that sits next to
    a DIFFERENT copy of the cache is worth one line of output."""
    _seen = []
    for _p in [CACHE_MANIFEST] + list(CACHE_MANIFEST_CANDIDATES):
        if _p in _seen:
            continue
        _seen.append(_p)
        try:
            with open(_p) as f:
                d = json.load(f)
        except Exception:
            continue
        if not d:
            continue
        if _p != CACHE_MANIFEST and verbose and _p not in _MANIFEST_SAID:
            _MANIFEST_SAID.add(_p)
            print('  manifest     no record beside the cache in use - '
                  'reading ' + _p)
        return d
    return {}


def write_cache_manifest(d):
    tmp = CACHE_MANIFEST + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(d, f, indent=2, sort_keys=True)
    os.replace(tmp, CACHE_MANIFEST)      # atomic: never a half-written manifest


def stale_reason(name, spec, manifest=None):
    """Why this cached layer should be rebuilt, or None to keep it.

    TWO DIFFERENT FACTS THAT USED TO PRINT THE SAME SENTENCE. A source
    file that really changed and a cache with no record of what built it
    both mean "rebuild", but only the first is news about raw_data/. The
    second happens every time a cache is restored without its sidecar,
    and calling that a changed source sends you looking for a download
    that never happened. Returns the reason so the caller can print it,
    and names the size change when there is one - a re-download that
    arrives at a different size is worth seeing, and one that arrives at
    the same size but a new mtime is worth distinguishing from it."""
    man = read_cache_manifest() if manifest is None else manifest
    was = man.get(name, {}).get('source')
    now = source_fingerprint(spec.get('path'))
    if now is None:                      # source gone: keep what we have
        return None
    if not was:
        return 'no record of what this layer was built from'
    if was.get('path') != now.get('path'):
        _old = os.path.basename(str(was.get('path')))
        return f'now built from a different file (was {_old})'
    if was.get('size') != now.get('size'):
        _a, _b = was.get('size') or 0, now.get('size') or 0
        _f = ((lambda n: f'{n / 1e6:,.1f} MB') if max(_a, _b) >= 1e6
              else (lambda n: f'{n:,} bytes'))
        return f'the source file changed size ({_f(_a)} -> {_f(_b)})'
    if was.get('mtime') != now.get('mtime'):
        return 'the source file was modified (same size, new timestamp)'
    return None


def layer_is_stale(name, spec, manifest=None):
    """True when the cached layer should be rebuilt. stale_reason() is
    the same question with the answer attached; this stays for callers
    that only need the boolean."""
    return stale_reason(name, spec, manifest) is not None


# %% M3c - county names, for any state
_COUNTY_SUFFIX = re.compile(
    r'\s+(County|Parish|Borough|Census Area|City and Borough|Municipality|'
    r'Municipio|District|Island|city)$', re.I)


def census_county_names(verbose=True):
    """FIPS -> name for STATE_FIPS, from the Census national county list.

    The last automatic step in the naming cascade, and the only one that
    touches the network. It is fetched at most ONCE: the answer is cached as
    county_names_<FIPS>.json beside the GeoPackage, and a machine with no
    network gets an empty dict and a printed line rather than an exception.
    Nothing measured depends on it - these are labels and folder names."""
    path = os.path.join(os.path.dirname(CACHE),
                        f'county_names_{STATE_FIPS}.json')
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except Exception:
            pass
    if not COUNTY_NAME_FETCH:
        return {}
    try:
        import urllib.request
        with urllib.request.urlopen(COUNTY_NAME_URL, timeout=20) as r:
            txt = r.read().decode('latin-1')
    except Exception as e:
        if verbose:
            print(f'  county names: could not reach the Census list '
                  f'({type(e).__name__}) - falling back to FIPS.\n'
                  f'                Fill COUNTY_NAMES in fl_config.py §9 to '
                  f'name them offline.')
        return {}
    # STATE|STATEFP|COUNTYFP|COUNTYNS|COUNTYNAME|CLASSFP|FUNCSTAT
    out = {}
    for line in txt.splitlines()[1:]:
        f = line.split('|')
        if len(f) > 4 and f[1].strip() == STATE_FIPS:
            out[f[2].strip().zfill(3)] = _COUNTY_SUFFIX.sub('', f[4].strip())
    if out:
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                json.dump(out, f, indent=1, sort_keys=True)
            if verbose:
                print(f'  county names: fetched {len(out)} for {STATE_NAME} '
                      f'and cached them as {os.path.basename(path)}')
        except OSError:
            pass
    return out


def county_name_field(gdf):
    """The column in a boundary file that holds a COUNTY's name, if any.

    TIGER's cousub files carry COUNTYFP but the NAME column is the
    SUBDIVISION's, so this deliberately does not match a bare 'NAME'. Some
    state-published equivalents do carry the county, under a name that says
    so."""
    if gdf is None or not len(gdf):
        return None
    for c in gdf.columns:
        k = str(c).strip().lower()
        if k in ('countyname', 'county_name', 'cnty_name', 'county', 'cnty',
                 'co_name', 'namelsadco', 'county_nm'):
            return c
    return None


def load_cached(name, required=True, verbose=True):
    """Read one layer back out of the cache, reprojected and normalised."""
    if name not in cache_layers():
        if required:
            raise FileNotFoundError(
                f'layer {name!r} is not in {CACHE}.\n'
                f'  Run fl_01_preprocess.py first - the statewide script '
                f'reads what\n  preprocessing leaves in clipped_layers/ and '
                f'never opens a source file.')
        return None
    g = gpd.read_file(CACHE, layer=name)
    if g.crs is None:
        g = g.set_crs(WORKING_CRS)
    if str(g.crs) != WORKING_CRS:
        g = g.to_crs(WORKING_CRS)
    if verbose:
        print(f'  {name:<12} {len(g):>9,} features  (from cache)')
    return normalise(g, name, verbose=verbose)


def resolve_layer(path, spec, verbose=False):
    """Which layer inside this file does the registry mean?

    NWI names its wetlands layer after the state (KY_Wetlands, TN_Wetlands);
    the WBD is WBDHU12 nationally and whatever a state agency called its
    export locally. So a registry entry's `layer` may be a LIST of spellings,
    and the first one the file actually has wins. A file with exactly one
    layer needs no name at all, which covers every shapefile and most
    single-purpose GeoPackages.

    Returns the layer name to hand the driver, or None for 'the only one'."""
    want = spec.get('layer')
    cands = ([want] if isinstance(want, str) else list(want or []))
    if not cands:
        return None
    if path is None or not str(path).lower().endswith(('.gpkg', '.gdb')):
        return cands[0]
    try:
        have = [n for n, _ in list_layers(open_path(path, verbose=False))]
    except Exception:
        return cands[0]
    low = {n.lower(): n for n in have}
    for c in cands:                      # exact, then case-insensitive
        if c in have:
            return c
        if c.lower() in low:
            return low[c.lower()]
    for c in cands:                      # then 'contains', for a prefixed copy
        for n in have:
            if c.lower() in n.lower():
                if verbose:
                    print(f'  {"":<12} layer {c!r} -> {n!r}')
                return n
    if len(have) == 1:
        if verbose:
            print(f'  {"":<12} no layer named {cands[0]!r}; the file has one '
                  f'layer, {have[0]!r} - using it')
        return have[0]
    return cands[0]                      # let read_source report it properly


def read_source(name, bbox=None, verbose=True):
    """Read one registered dataset from its ORIGINAL file.

    Order matters and is the same every time: push the row filter and the
    bounding box down to the driver, read, reproject to WORKING_CRS, keep
    columns. The pushdown is not an optimisation - it is what makes a national
    file readable at all.

    `bbox` is in WORKING_CRS and is reprojected to the source's own CRS before
    being handed to the driver. A bounding box compared across coordinate
    systems silently returns nothing at all."""
    spec = DATASETS[name]
    path = spec.get('path')
    if not path:
        if spec['required']:
            raise FileNotFoundError(f'{name}: no path resolved')
        return None
    if not _HAVE_PYOGRIO:
        g = gpd.read_file(path, layer=resolve_layer(path, spec))
        if g.crs is None:
            g = g.set_crs(WORKING_CRS)
        return normalise(g.to_crs(WORKING_CRS), name, verbose)

    path = open_path(path, verbose=verbose)
    layer = resolve_layer(spec.get('path'), spec, verbose=verbose)
    if str(path).startswith('/vsizip/'):
        try:
            read_info(path, layer=layer) if layer else read_info(path)
        except Exception:
            print(f'  {name:<12} reading in place failed - extracting instead')
            path = stage_zip(spec['path'])
    try:
        info = read_info(path, layer=layer) if layer else read_info(path)
    except Exception as e:
        # a file that is present but whose LAYER name is wrong lands here
        if not spec['required']:
            print(f'  {name:<12} unreadable ({type(e).__name__}); '
                  f'{spec["degrades"]}')
            return None
        raise RuntimeError(
            f'{name}: could not open {os.path.basename(path)}'
            + (f' layer {layer!r}' if layer else '') + f' - {e}\n'
            f'  layers actually in the file: '
            f'{[n for n, _ in list_layers(path)]}') from None

    cols = (resolve_cols(info['fields'], spec['wanted'])
            if spec.get('wanted') else None)
    newest_src = None
    if spec.get('newest'):
        # ('CDL', pattern) keeps one year; a third element keeps that many
        # of the NEWEST years, which is what the crop sequence reads. The
        # years are read under their own names and `CDL` is added as an
        # alias for the newest, so nothing that predates the sequence
        # notices the difference.
        _as, _pat, _keep_n = newest_spec(spec)
        hits = sorted(f for f in info['fields'] if re.fullmatch(_pat, f, re.I))
        if hits:
            _take = hits[-max(_keep_n, 1):]
            newest_src = _take[-1]
            cols = (cols or []) + _take
            if verbose:
                print(f'  {name:<12} newest {_as} field: {newest_src}'
                      + (f'  (+{len(_take) - 1} earlier year(s): '
                         f'{", ".join(_take[:-1])})'
                         if len(_take) > 1 else ''))
        elif verbose:
            print(f'  {name:<12} no field matching {_pat!r} - {_as} omitted')

    kw = {}
    if layer:
        kw['layer'] = layer
    if bbox is not None:
        kw['bbox'] = tuple(gpd.GeoSeries([shp_box(*bbox)], crs=WORKING_CRS)
                           .to_crs(info['crs']).total_bounds)
    if cols is not None:
        kw['columns'] = cols
    if spec.get('where'):
        kw['where'] = spec['where']

    try:
        g = read_dataframe(path, force_2d=True, **kw)
    except Exception as e:
        if not spec['required']:
            print(f'  {name:<12} unreadable ({e}); {spec["degrades"]}')
            return None
        raise

    if g.crs is None:
        g = g.set_crs(info['crs'])
    if str(g.crs) != WORKING_CRS:
        g = g.to_crs(WORKING_CRS)
    if newest_src:
        g[newest_spec(spec)[0]] = g[newest_src]   # alias, not a rename
    if spec.get('wanted'):
        lost = [c for c in spec['wanted']
                if c.lower() not in {x.lower() for x in g.columns}]
        if lost and verbose:
            print(f'  {name:<12} columns not in the file: {lost}')
    if verbose:
        f = []
        if spec.get('where'):
            f.append(spec['where'])
        if bbox is not None:
            f.append('bbox')
        print(f'  {name:<12} {len(g):>9,} features'
              + (f'  (pushdown: {"; ".join(f)})' if f else ''))
    return g


# %% M4 - geometry
DIM = {'Point': 0, 'MultiPoint': 0, 'LineString': 1, 'MultiLineString': 1,
       'LinearRing': 1, 'Polygon': 2, 'MultiPolygon': 2}


def keep_dimension(clipped, like=None, label=''):
    """Drop the slivers a clip leaves behind.

    Cutting a polygon layer at a boundary can leave a line or a point where a
    feature only grazed the edge. Those are real geometries with zero area, and
    they make later overlays fail outright with "mixed geometry types"."""
    if clipped is None or not len(clipped):
        return clipped
    want = (max(DIM.get(t, 2) for t in like.geometry.geom_type.unique())
            if like is not None and len(like)
            else max(DIM.get(t, 2)
                     for t in clipped.geometry.geom_type.unique()))
    ok = clipped.geometry.geom_type.map(lambda t: DIM.get(t, -1) == want)
    if (~ok).any():
        print(f'  {label:<12} dropped {int((~ok).sum())} zero-'
              f'{"area" if want == 2 else "length"} sliver(s) left by the clip')
    return clipped[ok].copy()


def repair(gdf):
    """Fix invalid geometries rather than letting a later overlay fail."""
    bad = ~gdf.geometry.is_valid
    if bad.any():
        print(f'  {"":<12} repairing {int(bad.sum()):,} invalid geometries')
        gdf.loc[bad, 'geometry'] = gdf.loc[bad, 'geometry'].apply(make_valid)
    return gdf


def area_ac(geom_array):
    """Areas in acres for a numpy array of geometries."""
    return shapely.area(np.asarray(geom_array, dtype=object)) / AC


def _subdivide(geom, max_coords=None, _depth=0):
    """Recursive bounding-box bisection: one big polygon -> many small ones
    that tile it exactly.

    An overlay indexes the frame on its RIGHT, and a single dissolved geometry
    is an index with ONE entry in it - so every left-hand row is tested against
    the whole thing. Measured at 34.7 s dissolved against 1.7 s subdivided, for
    identical acres. Depth-capped so a pathological sliver cannot run away."""
    max_coords = max_coords or SUBDIV_MAX_COORDS
    if shapely.get_num_coordinates(geom) <= max_coords or _depth >= 14:
        return [geom]
    x0, y0, x1, y1 = geom.bounds
    if (x1 - x0) >= (y1 - y0):
        xm = (x0 + x1) / 2
        halves = (shp_box(x0, y0, xm, y1), shp_box(xm, y0, x1, y1))
    else:
        ym = (y0 + y1) / 2
        halves = (shp_box(x0, y0, x1, ym), shp_box(x0, ym, x1, y1))
    out = []
    for h in halves:
        c = geom.intersection(h)
        if c.is_empty:
            continue
        gs = (c.geoms if c.geom_type.startswith('Multi')
              or c.geom_type == 'GeometryCollection' else [c])
        for g in gs:
            if g.geom_type in ('Polygon', 'MultiPolygon') and not g.is_empty:
                out += _subdivide(g, max_coords, _depth + 1)
    return out


def subdivided_gdf(geom, label='corridor', verbose=False):
    import time as _t
    t0 = _t.time()
    parts = (list(geom.geoms) if geom.geom_type.startswith('Multi')
             else [geom])
    pieces = []
    for p in parts:
        pieces += _subdivide(p)
    g = gpd.GeoDataFrame(geometry=pieces, crs=WORKING_CRS)
    if verbose:
        print(f'  {label}: {len(parts):,} component(s), '
              f'{shapely.get_num_coordinates(geom):,} vertices -> '
              f'{len(g):,} tiles in {_t.time() - t0:.1f}s')
    return g


def fast_overlay(left, big_geom, label='corridor', verbose=False):
    """`left` x one big geometry, via subdivision. Same result as overlaying
    against the whole thing, minus the wait."""
    return gpd.overlay(left,
                       subdivided_gdf(big_geom, label, verbose)[['geometry']],
                       how='intersection', keep_geom_type=True)


def clip_to(gdf, geom, label='extent', verbose=False):
    """Clip to one polygon, but INTERSECT only the rows that cross its
    boundary.

    A row completely inside is already its own answer, and on this data that is
    the overwhelming majority - intersecting them anyway means rebuilding
    thousands of geometries against an outline for nothing. `shapely.prepare`
    builds the index that makes `contains_properly` cheap."""
    shapely.prepare(geom)
    g = np.asarray(gdf.geometry.values, dtype=object)
    inside = shapely.contains_properly(geom, g)
    touches = shapely.intersects(geom, g) & ~inside
    out = gdf.loc[inside].copy()
    if touches.any():
        edge = gdf.loc[touches].copy()
        edge['geometry'] = shapely.intersection(
            np.asarray(edge.geometry.values, dtype=object), geom)
        out = pd.concat([out, edge], ignore_index=True)
    if verbose:
        print(f'  {label} clip: {int(inside.sum()):,} inside, '
              f'{int(touches.sum()):,} intersected, '
              f'{len(gdf) - int(inside.sum()) - int(touches.sum()):,} dropped')
    return out[~out.geometry.is_empty].reset_index(drop=True)


def corridor_one(gdf, kind, interior, ft):
    """One water source's corridor.

    `set_precision` on every union is not optional: without it, overlapping
    buffers from adjacent reaches raise 'ring edge missing'."""
    d = ft * 0.3048
    outer = union_all(set_precision(gdf.geometry.buffer(d).values, 0.01),
                      grid_size=0.01)
    if kind == 'line' or interior:
        return outer
    inner = union_all(set_precision(gdf.geometry.values, 0.01), grid_size=0.01)
    return outer.difference(inner)          # ring only: drop the open water


def drop_nwi_types(g, verbose=False):
    """Remove NWI systems that duplicate NHD. Prefers the WETLAND_TYPE text
    field; falls back to the first letter of ATTRIBUTE, which is the NWI system
    code (R = Riverine, L = Lacustrine)."""
    if 'WETLAND_TYPE' in g.columns:
        keep = ~g['WETLAND_TYPE'].astype(str).str.strip().isin(WETLAND_EXCLUDE)
    elif 'ATTRIBUTE' in g.columns:
        letters = {'Riverine': 'R', 'Lake': 'L'}
        bad = {letters[t] for t in WETLAND_EXCLUDE if t in letters}
        keep = ~g['ATTRIBUTE'].astype(str).str[:1].isin(bad)
    else:
        return g
    if verbose and (~keep).any():
        print(f'  wetland    excluding {int((~keep).sum()):,} NWI polygons '
              f'that duplicate NHD')
    return g[keep]


def band_sums(values, groups_flat, starts):
    """Sum `values` over pre-flattened neighbour lists. Equivalent to
    [values[i].sum() for i in neighbour_lists], ~40x faster on a state."""
    return np.add.reduceat(values[groups_flat], starts)


def copy_to_final(extra_files=()):
    """Copy WORK_DIR to FINAL_DIR once, sequentially, and verify it landed.

    Only does anything when SCRATCH_DIR is set. Printing the command and
    trusting someone to notice it at the end of a six-hour run is the weakest
    link in the whole pipeline.

    VERIFIED FILE BY FILE, FROM THE SOURCE SIDE. Comparing the two folders'
    total size was wrong in a way that cried wolf on every Colab run: the
    destination legitimately holds things the source never had. The per-tile
    checkpoints go straight to the durable folder when TILES_ON_FINAL is set,
    and a copy MERGES rather than replaces, so anything an earlier run left is
    still sitting there. Both make the destination bigger, and the old check
    called that a mismatch.

    What matters is whether every file THIS run copied arrived at the right
    size. Extra files in the destination are reported, not failed."""
    if WORK_DIR == FINAL_DIR:
        return
    import shutil
    import time as _t

    def files_under(p):
        """{relative path: size} for everything under p."""
        out = {}
        for r, _d, fs in os.walk(p):
            for x in fs:
                fp = os.path.join(r, x)
                try:
                    out[os.path.relpath(fp, p)] = os.path.getsize(fp)
                except OSError:
                    pass
        return out

    jobs = [(WORK_DIR, FINAL_DIR, 'outputs')] + list(extra_files)
    rule('COPYING TO THE PROJECT FOLDER')
    print('  The equivalent commands, if you would rather do it by hand:')
    for src, dst, _w in jobs:
        print(f'    cp -r "{src}" "{os.path.dirname(dst)}/"')
    if not COPY_TO_FINAL_AT_END:
        print('\n  COPY_TO_FINAL_AT_END = False - not copying. Run the above.')
        return

    for src, dst, what in jobs:
        isdir = os.path.isdir(src)
        want = (files_under(src) if isdir
                else {os.path.basename(src): os.path.getsize(src)})
        mb = sum(want.values()) / 1e6
        print(f'\n  {what}: {len(want):,} files, {mb:,.0f} MB -> {dst}',
              flush=True)
        t0 = _t.time()
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if isdir:
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)
        except Exception as e:
            print(f'    COPY FAILED: {e}\n    Run the command above by hand.')
            continue

        got = (files_under(dst) if isdir
               else {os.path.basename(dst): os.path.getsize(dst)})
        missing = [f for f in want if f not in got]
        wrong = [f for f in want if f in got and got[f] != want[f]]
        extra = [f for f in got if f not in want]
        el = _t.time() - t0
        if missing or wrong:
            print(f'    INCOMPLETE in {el:,.0f}s - {len(missing):,} file(s) '
                  f'missing, {len(wrong):,} the wrong size'
                  f'  <- re-copy by hand and check')
            for f in (missing + wrong)[:6]:
                print(f'      {f}')
            if len(missing) + len(wrong) > 6:
                print(f'      ... {len(missing) + len(wrong) - 6} more')
        else:
            print(f'    done in {el:,.0f}s - all {len(want):,} files verified '
                  f'({mb:,.0f} MB)')
        if extra:
            print(f'    the destination also holds {len(extra):,} file(s) '
                  f'({sum(got[f] for f in extra) / 1e6:,.0f} MB) this run did '
                  f'not write:\n    tile checkpoints'
                  + (' (TILES_ON_FINAL = True keeps them there)'
                     if TILES_ON_FINAL else '')
                  + ' and whatever an earlier run left.\n    Nothing was '
                    'deleted - a copy merges rather than replaces.')

    # a GeoPackage is the one thing worth opening again after the trip
    surf = os.path.join(FINAL_DIR, '_state', 'output', 'state_surface.gpkg')
    if os.path.exists(surf):
        try:
            ls = sorted(n for n, _ in list_layers(surf))
            want_l = {'grid', 'grid_scored', 'rip_crop'}
            print(f'\n  CHECK state_surface.gpkg: {", ".join(ls) or "(none)"}')
            if not want_l <= set(ls):
                print(f'    MISSING {", ".join(sorted(want_l - set(ls)))} - the '
                      f'copy did not land intact.\n    Delete it and copy '
                      f'again; the original is still in {WORK_DIR}.')
        except Exception as e:
            print(f'\n  CHECK state_surface.gpkg: UNREADABLE ({e})')


print(f'fl_common ready | geopandas {gpd.__version__} | '
      f'shapely {shapely.__version__} | '
      f'{"pyogrio" if _HAVE_PYOGRIO else "fiona (no pushdown - slower)"}')
