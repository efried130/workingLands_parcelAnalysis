# ---
# jupytext:
#   text_representation:
#     extension: .py
#     format_name: percent
# ---

# %% [markdown]
# # `fl_config.py` — every decision this pipeline makes, in one file
#
# **Nothing in this file reads data, computes anything, or imports geopandas.**
# It is paths, switches, thresholds and the dataset registry, and it is the
# only file you should have to edit. `fl_01_preprocess.py` and
# `fl_02_statewide.py` both start with `from fl_config import *`, so a value
# set here is the value the whole run uses.
#
# Read it top to bottom once. The sections are ordered by how often you will
# touch them:
#
# | § | what | how often you edit it |
# |---|---|---|
# | 0 | which state, and where you are running | once, per state |
# | 1 | where your files are | once, on a new machine |
# | 2 | the run switches | every run |
# | 3 | the study area | when you change state or county set |
# | 4 | the threshold | when you argue about the method |
# | 5 | Pass A cost and memory | rarely |
# | 6 | Pass C, the parcel work | when you buy parcel data |
# | 7 | figures | when someone asks for a different look |
# | 8 | the dataset registry — INPUT PATHS | when a download is renamed |
# | 9 | reference tables | almost never |

import os
import sys

# %% [markdown]
# ## 0 — Which state, and where you are running it
#
# **Two settings.** `STATE_FIPS` is the only thing that makes this a Kentucky
# run rather than a Tennessee one; `COLAB` is the only thing that makes it a
# Colab run rather than a local one. Everything downstream — the state's name,
# its postal abbreviation, whether its counties are called parishes, the cache
# file's name, the download patterns, every figure title — is derived from
# those two.
#
# There is no Kentucky anywhere below this section.

# %% C0a - the state
# Two digits, quoted. '21' Kentucky · '47' Tennessee · '01' Alabama · '13'
# Georgia. FL_STATE_FIPS overrides it without editing this file.
STATE_FIPS = os.environ.get('FL_STATE_FIPS') or '21'

# FIPS -> (name, postal code). Reference data: the only reason it is here is so
# that setting STATE_FIPS is genuinely the only edit.
STATE_TABLE = {
    '01': ('Alabama', 'AL'),        '02': ('Alaska', 'AK'),
    '04': ('Arizona', 'AZ'),        '05': ('Arkansas', 'AR'),
    '06': ('California', 'CA'),     '08': ('Colorado', 'CO'),
    '09': ('Connecticut', 'CT'),    '10': ('Delaware', 'DE'),
    '11': ('District of Columbia', 'DC'),
    '12': ('Florida', 'FL'),        '13': ('Georgia', 'GA'),
    '15': ('Hawaii', 'HI'),         '16': ('Idaho', 'ID'),
    '17': ('Illinois', 'IL'),       '18': ('Indiana', 'IN'),
    '19': ('Iowa', 'IA'),           '20': ('Kansas', 'KS'),
    '21': ('Kentucky', 'KY'),       '22': ('Louisiana', 'LA'),
    '23': ('Maine', 'ME'),          '24': ('Maryland', 'MD'),
    '25': ('Massachusetts', 'MA'),  '26': ('Michigan', 'MI'),
    '27': ('Minnesota', 'MN'),      '28': ('Mississippi', 'MS'),
    '29': ('Missouri', 'MO'),       '30': ('Montana', 'MT'),
    '31': ('Nebraska', 'NE'),       '32': ('Nevada', 'NV'),
    '33': ('New Hampshire', 'NH'),  '34': ('New Jersey', 'NJ'),
    '35': ('New Mexico', 'NM'),     '36': ('New York', 'NY'),
    '37': ('North Carolina', 'NC'), '38': ('North Dakota', 'ND'),
    '39': ('Ohio', 'OH'),           '40': ('Oklahoma', 'OK'),
    '41': ('Oregon', 'OR'),         '42': ('Pennsylvania', 'PA'),
    '44': ('Rhode Island', 'RI'),   '45': ('South Carolina', 'SC'),
    '46': ('South Dakota', 'SD'),   '47': ('Tennessee', 'TN'),
    '48': ('Texas', 'TX'),          '49': ('Utah', 'UT'),
    '50': ('Vermont', 'VT'),        '51': ('Virginia', 'VA'),
    '53': ('Washington', 'WA'),     '54': ('West Virginia', 'WV'),
    '55': ('Wisconsin', 'WI'),      '56': ('Wyoming', 'WY'),
    '60': ('American Samoa', 'AS'), '66': ('Guam', 'GU'),
    '69': ('Northern Mariana Islands', 'MP'),
    '72': ('Puerto Rico', 'PR'),    '78': ('U.S. Virgin Islands', 'VI'),
}
assert STATE_FIPS in STATE_TABLE, (
    f'STATE_FIPS {STATE_FIPS!r} is not a state FIPS code. It is two digits, '
    f'quoted, from: {", ".join(sorted(STATE_TABLE))}')

STATE_NAME, HOME_STATE = STATE_TABLE[STATE_FIPS]   # 'Kentucky', 'KY'

# What this state calls a county. Nothing but labels and folder names.
COUNTY_WORD = {'22': 'Parish', '02': 'Borough'}.get(STATE_FIPS, 'County')
COUNTY_WORDS = {'County': 'counties', 'Parish': 'parishes',
                'Borough': 'boroughs'}[COUNTY_WORD]

# The long form used in figure titles. Four states call themselves a
# commonwealth and two of them are in this programme's footprint.
STATE_TITLE = ('Commonwealth of ' + STATE_NAME
               if STATE_FIPS in ('21', '42', '25', '51') else
               'State of ' + STATE_NAME)

# %% C0b - Colab
# True   force the Colab path: mount Drive, build on /content, copy back
# False  force the local path
# 'auto' decide by whether this is actually running in Colab   <- the default
#
# On Colab this is not a convenience. Google Drive is a FUSE mount, and GDAL
# updates a GeoPackage in place through SQLite - thousands of small random
# writes. A big .gpkg written straight to Drive can come back with layers
# missing AND NOT RAISE. The Colab path therefore builds on /content (real
# local disk) and copies the finished folder to Drive once, while the per-tile
# checkpoints - one-shot creates, which Drive does handle - stay on Drive so a
# runtime reset costs one tile rather than the pass.
COLAB = 'auto'

if os.environ.get('FL_COLAB'):           # FL_COLAB=1 / FL_COLAB=0 overrides
    COLAB = os.environ['FL_COLAB'] not in ('0', 'false', 'False', '')
if COLAB == 'auto':
    # importable, not merely /content existing - other hosted notebooks have a
    # /content too, and none of them want a Drive mount.
    import importlib.util as _ilu
    COLAB = ('google.colab' in sys.modules
             or _ilu.find_spec('google.colab') is not None)

# Where your project folder lives ON DRIVE. Only read when COLAB is true.
COLAB_DRIVE_DIR = '/content/drive/MyDrive/FergusonLynch/Parcel_priority_maps'
COLAB_MOUNT = True        # False if you mounted Drive yourself in a cell above
COLAB_SCRATCH = '/content'

if COLAB and COLAB_MOUNT and not os.path.isdir('/content/drive/MyDrive'):
    try:
        from google.colab import drive as _gdrive
        _gdrive.mount('/content/drive')
    except Exception as _e:                      # not Colab after all, or refused
        print(f'  Drive not mounted ({type(_e).__name__}: {_e}) - '
              f'carrying on with local paths')

# %% [markdown]
# ## 1 — Where your files are
#
# **This is the only path you have to get right.** Everything else is derived
# from it. Set the environment variable `FL_PROJECT_DIR` to override without
# editing the file, which is how the test fixture and any CI run do it.
#
# The folder must already contain `raw_data/` with your downloads in it.
# `clipped_layers/` and `outputs_statewide/` are created for you.

# %% C1 - the project folder
# Order of precedence: the environment variable, then the Drive folder if this
# is Colab, then a folder in your home directory.
PROJECT_DIR = (os.environ.get('FL_PROJECT_DIR')
               or (f'{COLAB_DRIVE_DIR}_{HOME_STATE}' if COLAB else
                   os.path.expanduser(
                       f'~/FergusonLynch/Parcel_priority_maps_{HOME_STATE}')))

RAW_DIR = os.path.join(PROJECT_DIR, 'raw_data')              # INPUTS, never created
OUT_DIR = os.path.join(PROJECT_DIR, 'outputs_statewide')     # everything this makes

# The one cache file, named for the state so two states can share a disk.
# Preprocessing writes it; the statewide run reads it and never opens a source
# file again. Rename it and you get a second, independent cache - which is how
# you keep a statewide build and a priority-area build side by side without
# either overwriting the other.
CACHE_NAME = f'{HOME_STATE.lower()}_state_layers.gpkg'

# SCRATCH_DIR - build the big GeoPackages on a fast local disk, then copy the
# finished folder to PROJECT_DIR once.
#
#   None            build in place under PROJECT_DIR      <- right on a local machine
#   '/content'      Colab: local disk, because Drive is a FUSE mount
#   '/scratch/...'  a cluster, or a laptop whose project folder is a network share
#
# WHY THIS SWITCH EXISTS AT ALL. GDAL updates a GeoPackage in place through
# SQLite - thousands of small random-access writes. A local SSD does not care.
# A FUSE-mounted network drive (Google Drive, some SMB shares) handles that
# worst of anything it does, and a large .gpkg written straight to one can come
# back with most of its layers missing WITHOUT RAISING. If PROJECT_DIR is on
# ordinary local disk, leave this None.
#
# On Colab this is set for you (§0, COLAB_SCRATCH) because Drive is exactly
# that FUSE mount.
SCRATCH_DIR = (os.environ.get('FL_SCRATCH_DIR')
               or (COLAB_SCRATCH if COLAB and PROJECT_DIR.startswith('/content/drive')
                   else None))

# Where the heavy writing happens, and where the finished copy lands.
WORK_DIR  = (os.path.join(SCRATCH_DIR, 'outputs_statewide') if SCRATCH_DIR
             else OUT_DIR)
FINAL_DIR = OUT_DIR
COPY_TO_FINAL_AT_END = True      # only does anything when SCRATCH_DIR is set

# The cache is WRITTEN to scratch when there is one, but READ from wherever it
# actually is - scratch first, then the project folder. That is what lets you
# build on a fast disk this session and re-read last week's copy the next.
CLIPPED_DIR = (os.path.join(SCRATCH_DIR, 'clipped_layers') if SCRATCH_DIR
               else os.path.join(PROJECT_DIR, 'clipped_layers'))
FINAL_CLIPPED_DIR = os.path.join(PROJECT_DIR, 'clipped_layers')
CACHE = next((_p for _p in (os.path.join(CLIPPED_DIR, CACHE_NAME),
                            os.path.join(FINAL_CLIPPED_DIR, CACHE_NAME))
              if os.path.exists(_p)),
             os.path.join(CLIPPED_DIR, CACHE_NAME))

# Pass A writes one checkpoint per tile. Those are one-shot creates, which even
# a FUSE mount does reliably, so they go to the DURABLE folder: a machine that
# reboots four hours into a six-hour pass then costs you one tile, not the pass.
TILES_ON_FINAL = bool(SCRATCH_DIR)

# %% [markdown]
# ## 2 — The run switches
#
# The two expensive steps, and whether to do them again.

# %% C2 - what to run
# PREPROCESSING - cut the national/statewide downloads into the cache.
#   True     rebuild every cached layer
#   False    build nothing; fail loudly if the cache is missing a layer
#   'auto'   build only the layers that are not cached yet   <- the usual setting
RUN_PREPROCESS = 'auto'

# PASS A - the statewide density surface. This is the expensive one: on
# Kentucky it is 89 tiles and about five hours. Everything downstream of it -
# the threshold, the zones, the counties, every figure - re-runs in SECONDS off
# the cached surface.
#
#   True     rebuild the surface
#   False    load it; fail loudly if absent
#   'auto'   build it only if state_surface.gpkg is missing   <- the usual setting
#
# LEAVE THIS ON 'auto'. Set it to True only when something UPSTREAM of the
# surface changes: the cache, RIP_BUFFER_FT, CELL_MI, BAND_MI, TILE_KM.
# Changing BASIS, the floor, a colour or a figure does NOT require it.
RUN_SURFACE = 'auto'

# The surface on disk carries a SIGNATURE of what it was built from: the cache,
# the extent rule, the county list, the analysed acreage, the cell and band
# size, the buffer, the tile size and the water sources. With the guard on, a
# run whose settings no longer match that signature does not quietly reuse the
# surface -
#   RUN_SURFACE = 'auto'   rebuilds it, and says which setting changed
#   RUN_SURFACE = False    stops, because you asked for a reuse that is unsafe
# This is what stops `FL_COUNTIES=047,225` - the two-county test the README
# recommends - from leaving a two-county surface behind that the next statewide
# run silently reports the whole state out of.
#
# False disables the check entirely: the mismatch is printed and the surface is
# used anyway. Only do that if you know exactly why the signature differs.
SURFACE_GUARD = True

# PASS C - the parcel work, for the counties that have a parcel file on disk.
# False keeps the zone-only deliverable, which is what 119 of 120 Kentucky
# counties get anyway.
RUN_PASS_C = True

# Pass A checkpointing. RESUME picks up a stopped run where it left off;
# KEEP_TILES decides whether the checkpoints survive a successful build.
#
# KEEP_TILES = True is insurance against a crash in a LATER stage: the surface
# exists but if you had to rebuild it you would re-tile from zero. The price is
# a few hundred MB.
RESUME_TILES = True
KEEP_TILES   = True

WRITE_COUNTY_FILES = True        # one folder per county (Pass B)

# %% [markdown]
# ## 3 — The study area
#
# Which counties, and what each county's extent is measured on. **The state
# itself is §0** — `STATE_FIPS`, and the name, postal code and county word that
# follow from it.

# %% C3 - the study area
# 'all', or a list of 3-digit FIPS strings. A short list is how you test the
# whole pipeline in minutes before committing to a full run.
#   COUNTIES = ['047', '225']
#
# FL_COUNTIES=047,225 in the environment overrides it without editing the file,
# which is the least error-prone way to do a quick test run: you cannot forget
# to set it back afterwards.
COUNTIES = ([c.strip().zfill(3) for c in os.environ['FL_COUNTIES'].split(',')]
            if os.environ.get('FL_COUNTIES') else 'all')

# WHAT PASS A MEASURES. This sets each county's extent - i.e. what
# `avg_county` divides by - and therefore the ground the density surface is
# built over, so CHANGING IT REBUILDS PASS A.
#
# IT IS NOT THE SWITCH FOR "ANALYSE ONLY THE PRIORITY AREA". That is
# REPORT_WITHIN_PA, below: applied in Pass B, costs seconds, and leaves every
# cell's 3-mile neighbourhood measured over its real surroundings instead of
# a band cut at the PA boundary.
#
#   'whole_county'  every county measured edge to edge. The right answer with
#                   a statewide cache, because it is what makes avg_county
#                   comparable across all 120 - and it is what
#                   REPORT_WITHIN_PA requires.               <- default
#   'pa_only'       only the counties overlapping the priority area, each cut
#                   to it. Pass A then costs what the PA costs rather than
#                   what the state costs, which is worth it for exactly one
#                   case: a cold start on a new state where you will never
#                   want the statewide surface. The price is that the 3-mile
#                   band is truncated at the PA boundary - the county-line
#                   defect the statewide surface exists to remove - and that
#                   wanting the state later means paying for Pass A twice.
#
# 'pa_where_available' is LEGACY and no longer documented in the README. It
# measured PA counties on county n PA and every other county edge to edge,
# which produces a manifest whose rows are not comparable on avg_county; it
# existed only because there was no way to keep a statewide surface and
# report on a subset of it. REPORT_WITHIN_PA does that properly now. The
# value still RUNS and the assert in section 10 still accepts it, so a
# surface built under it loads and reports a mismatch rather than raising.
#
# Two counties on different bases are NOT comparable on avg_county. Compare
# `threshold_applied`, which is in absolute acres per square mile either way.
COUNTY_EXTENT_MODE = 'whole_county'

# A county overlapping the priority area by less than this is treated as
# outside it. Two Kentucky counties only touch the boundary, and clipping one
# to a sliver would hand it a county average computed over a square mile.
MIN_PA_SHARE = 0.05

# THE PRIORITY-AREA LAYER IS ENABLED BY BEING THERE. There is no on/off
# switch for it, because it is a required=False dataset like every other one
# in the registry: a file in raw_data/ matching the 'pa' pattern is used, and
# a missing one degrades with a printed reason. To turn it off, take the file
# out of raw_data/ and rebuild pa_geo (RUN_PREPROCESS = 'auto').
#
# It does three jobs, in increasing order of consequence:
#
#   1. it says which counties are in the framework, for REPORTING -
#      `pa_share` and `in_pa`. Always live, decides nothing.
#   2. it can decide where results are REPORTED - REPORT_WITHIN_PA, below.
#      Applied in Pass B, so it costs seconds.
#   3. it can decide what Pass A MEASURES - COUNTY_EXTENT_MODE, above.
#      That one rebuilds the surface.

# REPORT INSIDE THE PRIORITY AREA WITHOUT REBUILDING THE SURFACE.
#
# COUNTY_EXTENT_MODE, above, decides what Pass A MEASURES. A pa_* value
# there moves analysis_ac and the county list, so the surface signature no
# longer matches and Pass A rebuilds - hours - and what comes back is
# WORSE, because the 1-mile lattice is then laid only over the priority
# area and every edge cell's 3-mile band is measured over a truncated
# piece of ground. That is the county-line defect the statewide surface
# exists to remove, reintroduced against a programmatic boundary.
#
# This switch decides where the results are REPORTED. It is applied in
# Pass B (cell B6), after the surface has loaded, so it costs seconds and
# never rebuilds anything. Every reference, the floor sweep, the zones,
# the manifest and every figure are then computed on priority-area ground
# - while rip_per_sqmi and the nbhd_* columns stay exactly as Pass A
# measured them statewide, so a cell just inside the boundary still sees
# its real surroundings. `band_in_report` records, per cell, how much of
# its 3-mile band lies inside the mask.
#
# THE USUAL COMBINATION IS 'whole_county' PLUS True. Counties under
# MIN_PA_SHARE are dropped rather than clipped to a sliver.
#
# A masked run's averages and derived floor are NOT comparable with a
# statewide run's - they are computed over different ground, which is the
# point - so `report_mask` is stamped on zone_summary.csv and the county
# manifest and the two cannot be read as one table.
REPORT_WITHIN_PA = False

# A county whose cached layers do not reach this much of its own extent cannot
# be analysed. It is reported as `no_source_data` and kept OUT of the state and
# region averages, so missing data cannot drag a real average down. A county
# with no zones and a county nobody looked at are different sentences.
MIN_SOURCE_COVER = 0.95

# %% [markdown]
# ## 4 — The threshold
#
# A cell is **hot** when it clears a RATIO against its reference AND an
# absolute FLOOR in acres per square mile, and holds enough itself to be a
# concentration rather than a neighbour of one.
#
# ```
# hot  =  (neighbourhood density >= HOT_MULTIPLE x reference)
#     AND (own-cell density      >= reference)
#     AND (neighbourhood density >= FLOOR)
# ```
#
# `BASIS = 'county'` with `FLOOR_MODE = 'off'` reproduces the single-county
# notebook's rule exactly. That is how the two are kept comparable.

# HOW MANY YEARS OF CSB TO READ, for the crop-sequence label.
#
# CSB polygons are drawn so the crop sequence is CONSTANT within each one,
# and every year of the release rides on the SAME polygon as its own
# integer column (CDL2018 ... CDL2025 in the 2018-2025 release). So this
# costs columns, not geometry: no extra corridor, lattice or overlay, and
# the run-time cost is a groupby.
#
# 3 answers "what has this field been lately" - the question a field rep
# asks - without reaching back to land use that has since changed. Raise it
# for a longer rotation history; 1 turns the sequence off entirely and
# restores the single-year behaviour.
#
# Changing this makes the cached csb layer short of columns, which
# fl_01_preprocess.py notices on its own and fixes by rebuilding that ONE
# layer. It does NOT require a Pass A rebuild: an existing surface gets the
# years joined on by _fid and verified against the CDL column it already
# carries (see the crop-sequence block below B5b).
CROP_SEQ_YEARS = 3

# %% C4 - the rule
CELL_MI        = 1.0     # grid cell
BAND_MI        = 3.0     # neighbourhood radius - a field rep's working range
HOT_MULTIPLE   = 1.5     # hot = neighbourhood density >= this x the reference
MIN_ZONE_CELLS = 2       # ignore single isolated cells

# WHICH AVERAGE THE RATIO IS TAKEN AGAINST.
#   'county'  each county against itself - the notebook's rule. Answers "where
#             in this county?", and zones do not travel between counties.
#                                                                 <- default
#   'state'   one average for the state. Answers "which counties deserve a
#             biologist?", and flattens real regional differences.
#   'region'  one average per HUC8. The defensible middle: the Jackson Purchase
#             and the Bluegrass have different riparian-cropland baselines, and
#             comparing a Purchase county to the state mean makes it look hot
#             for reasons that are not about where anyone should stand.
# All three ratios are computed and carried as columns whatever this says, so
# switching the basis is a REPORTING decision, not a rebuild.
BASIS = 'county'

# THE ABSOLUTE CUT, in riparian-crop acres per square mile.
#   'percentile'  FLOOR = FLOOR_PCTL of the statewide cell-density distribution.
#                 Derived from the data, reported, and swept.
#   'knee'        FLOOR = the knee of the sweep, found by the rule below
#                                                                 <- default
#   'absolute'    FLOOR = FLOOR_AC_PER_SQMI, set by hand
#   'off'         no floor - pure ratio, i.e. the notebook's behaviour
#
# THE KNEE IS ALWAYS COMPUTED AND ALWAYS REPORTED, whatever this is set to, so
# a run on 'percentile' still tells you what the data would have chosen and how
# far 75 is from it. Only 'knee' actually adopts it.
#
# WHAT THE FLOOR IS FOR: a ratio alone always finds a top tail, so a county
# with almost no riparian cropland manufactures "hotspots" out of three fields
# on a creek. The floor is what makes a zone in Ballard and a zone in Christian
# mean the same thing.
#
# WHEN THE FLOOR DOES NOTHING: on a basis whose reference is high enough,
# HOT_MULTIPLE x reference already exceeds the floor everywhere, so the floor
# removes no cells. On the STATE basis that is arithmetic - one reference
# applies everywhere, so the ratio cut is the single number
# HOT_MULTIPLE x AVG_STATE, and the floor only ever binds when FLOOR_PCTL puts
# it ABOVE HOT_MULTIPLE x the state average. The run says so, and the
# basis_comparison table carries `cells_cut_by_floor` and `floor_binds`.
FLOOR_MODE        = 'knee'        # adopt the knee found by FLOOR_KNEE_RULE
FLOOR_PCTL        = 75.0          # still the fallback if no knee is found
FLOOR_AC_PER_SQMI = 20.0
FLOOR_SWEEP_PCTLS = [0, 25, 50, 60, 70, 75, 80, 85, 90, 95]   # the REPORTED table

# FINDING THE KNEE
# ----------------
# Raising the floor buys CONCENTRATION - the hot ground holds a larger share of
# the resource than it does of the area - and costs COUNTIES, which drop out of
# the deliverable one by one. The knee is where that trade stops paying.
#
# Raising the floor always buys CONCENTRATION - the designated ground holds a
# larger share of the resource than it does of the area - and always costs
# CAPTURE, because some riparian cropland ends up outside the map. What differs
# between the rules below is what they count as the cost.
#
# Three of them share one shape: put gain and cost on a 0-1 scale over the
# swept range and take the percentile that maximises
#
#       net(p) = gain(p) - weight x cost(p)
#       gain(p) = (concentration(p) - min) / (max - min)
#
# Normalising is what makes "concentration rose 0.4x" and "nine counties left"
# comparable at all; without it the two are in different units and the phrase
# "when the trade stops paying" has no arithmetic behind it. net is 0 at both
# ends by construction, so the knee is the high point of a curve pinned at
# zero. Ties go to the LOWEST percentile, because the floor should be the
# smallest cut that earns its keep.
# WHICH RULE `FLOOR_MODE = 'knee'` ADOPTS. All four are computed and printed
# on every run whatever this says, with what each one keeps, so the choice is
# made against numbers rather than in the abstract:
#
#   'resource'  cost = riparian cropland pushed outside the map.  <- default
#               Monotone, smooth, and the number a reviewer asks about.
#   'counties'  cost = counties that end up with no zone at all. Protects
#               statewide delivery; coarse, because with few counties one of
#               them is a large fraction of the whole range.
#   'curve'     NO PARAMETERS. Normalise resource captured and area designated
#               to 0-1 and take the point furthest above the chord - the
#               classic knee of the capture curve.
#   'marginal'  The exchange rate, in ac/sq mi. Keep raising the floor while
#               the ground it throws away is WORSE than the state average;
#               stop when it would start costing better-than-average ground.
#
# ZONE COUNT IS DELIBERATELY NOT A COST. Raising the floor can knock out one
# bridging cell and SPLIT a zone in two, so n_zones moves both ways as the
# floor rises; an argmax over that finds where a bridge broke, not a knee. It
# is reported in the sweep table and never optimised against.
FLOOR_KNEE_RULE = 'marginal'

# Which rules appear on FIGURE 3b. Every rule in KNEE_RULES is still computed,
# still printed in the run log and still written to floor_knee_rules.csv -
# this only decides what is drawn.
#
# 'marginal' is drawn again, and its panel deliberately carries a DIFFERENT
# y-unit from the others: an exchange rate in ac/sq mi against the state
# average, not a normalised index. That is the reason it was once dropped and
# the reason it earns its place back - it is the only rule whose y-axis is in
# the units the argument is actually conducted in, so it is the only one whose
# pick can be checked without trusting a normalisation.
KNEE_RULES_SHOWN = ('resource', 'counties', 'curve', 'marginal')

FLOOR_KNEE_WEIGHT = 1.0      # 'resource' and 'counties' only. > 1 protects the
                             # cost side, < 1 chases concentration
FLOOR_KNEE_MARGINAL_X = 1.0  # 'marginal' only: the multiple of the state
                             # average that the discarded ground may reach
FLOOR_KNEE_MARGINAL_SMOOTH = 3  # 'marginal' only: steps in the centred rolling
                             # MEDIAN taken before the crossing is tested. The
                             # raw exchange rate is a ratio of two small
                             # differences and spikes wherever a step happens
                             # to drop very little ground; without smoothing a
                             # single spike decides the floor. 1 = raw.
FLOOR_KNEE_GRID   = 2.5      # percentile step for the search (the reported
                             # table stays at FLOOR_SWEEP_PCTLS)

# HOW SHARP IS THE PICK. An argmax returns one number whatever the curve looks
# like, and these curves are flat: on the Kentucky county-basis run the
# 'resource' net beat its own runner-up by 1.3% of its range and reversed
# direction seven times across the span it was choosing within. Reporting
# 3.69 there implies a precision the curve does not have.
#
# So every rule also reports its PLATEAU: the contiguous run of swept floors
# whose net is within this fraction of the rule's own net RANGE of the best
# one. It changes no adopted value - the pick is still the argmax - it only
# says how much of the axis the rule regards as equivalent to it. Where two
# rules' plateaus overlap is a far better argument for a floor than either
# rule's argmax.
#
# 'marginal' is not an argmax but a crossing, so its range is the bracket
# between the last floor under the threshold and the first one over it.
FLOOR_KNEE_PLATEAU_TOL = 0.05
FLOOR_KNEE_MAX_PCTL = 95.0   # never recommend a floor above this percentile
# A knee is only meaningful if concentration actually moves across the sweep.
# Below this it is flat, the search is noise, and the run says so and falls
# back to FLOOR_PCTL.
FLOOR_KNEE_MIN_RANGE = 0.05  # in concentration units (e.g. 2.10x -> 2.15x)

# A HUC8 mostly outside the analysed ground has a truncated denominator, and
# one with very few cells is too small a sample to be its own reference. Either
# way those cells fall back to the state average, and WHICH reference each cell
# actually used is recorded on the cell (`region_ref`).
REGION_MIN_COVER = 0.50
REGION_MIN_CELLS = 30

# A zone is a union of whole 1-mile cells, so it laps a little way over the
# boundary of the county it belongs to. Less than this next door is that
# overhang, not a share of the zone.
#
# ONE WHOLE CELL (640 ac), not a fraction of one. At 64 ac - a TENTH of a cell
# - a zone lapping a few hundred metres across a line gave the county next
# door a zone row, a packet and a line in the summary, without that county
# holding a single hot cell of its own. That is how `counties with zones`
# came out ABOVE `counties_with_hot_cells` in the sweep: the manifest credits
# a county by the zone's GEOMETRY, the sweep counts only counties owning a
# cell that survived. Requiring a full cell makes the two answer the same
# question.
MIN_ZONE_COUNTY_AC = 640.0

# The same argument one level down. A zone laps into the SUBWATERSHED
# next door as well, and a HUC12 is 20-40 sq mi where a county is ~300 -
# so one whole cell of a HUC12 is a real share of it, not an overhang,
# and the cut is smaller. Raise it if `n_huc12` reads high.
MIN_ZONE_HUC12_AC = 64.0

# How many entries the digest columns on zone_summary.csv carry before
# they stop. Both have a full, unrolled table beside them
# (zone_cover_by_cdl.csv, zone_huc12.csv), so this is about how wide a
# spreadsheet cell should be, not about what was measured.
ZONE_CDL_TOP = 8          # classes listed in rip_cdl_acres
ZONE_HUC12_TOP = 8        # codes listed in huc12_ids

# Zone adjacency. Zones are disjoint by construction, so "adjacent" has to mean
# "within working distance". None = use BAND_MI, which is the honest answer:
# two zones one neighbourhood apart are one trip.
ZONE_NEIGHBOR_MI = None

# %% [markdown]
# ## 5 — Pass A: what gets measured, and what it costs
#
# The corridor is every stream, river, waterbody and wetland buffered by
# `RIP_BUFFER_FT` and **unioned**, so ground that is both streamside and
# wetland is counted once. Statewide that is one geometry with millions of
# vertices, so it is built one tile at a time with a halo.

# %% C5 - the corridor and the tiling
RIP_BUFFER_FT = 100.0    # THE adjacency distance, outward from the water's edge

# Each source, and how its buffer is measured.
#   line     100 ft from the CENTERLINE, so land per bank is 100 ft minus half
#            the channel width. Every other source measures from an edge.
#   polygon  interior False = ring only, i.e. the open water itself is dropped
WATER_SOURCES = {
    'stream':    {'var': 'nhd_flow', 'kind': 'line',    'interior': None},
    'river':     {'var': 'nhd_area', 'kind': 'polygon', 'interior': False},
    'waterbody': {'var': 'nhd_wb',   'kind': 'polygon', 'interior': False},
    'wetland':   {'var': 'wetlands', 'kind': 'polygon', 'interior': False},
}
USE_SOURCES = ['stream', 'river', 'waterbody', 'wetland']

# NWI systems that duplicate NHD. Left in, a riverine NWI polygon and the NHD
# flowline inside it both contribute the same corridor.
WETLAND_EXCLUDE = ['Riverine', 'Lake']

# MEMORY AND TIME. The halo only has to exceed the buffer width; a mile is
# fifty times it. Tiling changes how many PIECES a field is split into, never
# how many acres it has.
#
# Smaller TILE_KM = less peak memory, more tiles, slightly more total time.
# 40 km on Kentucky is 89 tiles. The dominant per-tile cost is the water
# buffer-and-union, so eastern Kentucky's headwater density makes late tiles
# slower than early ones - a rising ETA is expected, not a problem.
TILE_KM  = 40.0
HALO_MI  = 1.0
SUBDIV_MAX_COORDS = 512  # vertices per overlay piece; near the flat part of the trade

# %% [markdown]
# ## 6 — Pass C: the parcel work
#
# Runs only for counties whose parcel file is already on disk. Everything here
# governs the outreach list; none of it touches the zones.

# %% C6 - parcels
#   'auto'  every county in the manifest with has_parcels = True  <- default
#   a list  ['047', '225'] runs those counties only
PARCEL_COUNTIES = 'auto'

KEY = 'll_uuid'          # the parcel identifier; one row per parcel

# MATERIALITY FLOOR. Below this a parcel is not worth a visit. This is a floor
# on the RESOURCE, not a size rank - a 500-acre parcel with 1 riparian acre
# also fails - but it IS size-selective at the bottom, because a parcel cannot
# hold more riparian cropland than it has land.
#
# LEAVE IT FIXED ACROSS COUNTIES. The argument is the opposite of the one for
# zones, on purpose: for zones you want local contrast plus a global floor,
# because the consumer is a map; for parcels you want one global rule, because
# the consumer is a contact list a human acts on, and two lists built to
# different cut-offs cannot be compared. Each county's own sweep is written to
# validation/parcel_floor_sweep.csv so you can see where it would have landed.
MIN_RIP_AC = 2.0

# THE SCORE. Two terms, both DENSITIES, both converted to percentile ranks
# first so the weight between them means something:
#
#   score = W_PLACE x (how dense the surrounding 3 miles are)
#         + (1 - W_PLACE) x (how much of this parcel is riparian cropland)
#
# Weighted toward PLACE because the point is clusters a biologist can work as
# one trip. Acres are reported everywhere and scored nowhere.
W_PLACE = 0.75

# None = do not write outreach_top<N>.csv at all. That file was only ever
# head(N) of outreach_parcels_ranked.csv, which is already in rank order
# and carries the `rank` column - so turning it off costs a file and no
# information, and stops a stale short list outliving a re-run. Set an
# integer to get it back.
TOP_N_CONTACTS     = None   # was 250
FRONTAGE_STEP_M    = 10.0   # water edge is chopped into pieces this long
ADJ_TOL_M          = 1.0    # parcel fabrics rarely share exact edges
FLOOR_SWEEP_PARCEL = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

# %% [markdown]
# ## 7 — Figures
#
# One sequential family everywhere: **YlOrRd**, light = little, dark = a lot.
# Nothing maps a magnitude with any other ramp, so a colour means the same
# thing on the state map, the comparison panels and all 120 county sheets.
# Nothing is coloured by rank, so a county does not change colour when another
# county is added.

# %% C7 - figures
MAKE_COUNTY_FIGURES = True
MAKE_PARCEL_FIGURES = True
COUNTY_FIG_MAX = None     # None = every county with cells; a number caps it
FIG_DPI = 200

# Draw each figure into the notebook as well as saving it. Every figure is
# written to disk either way; this is only about seeing them go by. On Colab
# that is usually what you want, so it follows COLAB unless you say otherwise -
# but a 120-county run puts 120 maps in the notebook, which is worth knowing
# before you set it True by hand.
SHOW_FIGS = bool(COLAB)

# The ramp, and the three inks. HOT is the top step of the ramp rather than a
# separate accent colour.
SEQ = ['#ffffb2', '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026']
RAMP_NAME = 'YlOrRd'      # the continuous version of the same ramp
HOT = SEQ[-1]             # "hot" is the TOP STEP of the ramp, not a new colour
# The four knee rules are CATEGORICAL - four different arguments, not more or
# less of one thing - so they get four hues rather than four steps of the ramp,
# and the order matters: it is validated for colour-vision deficiency pairwise
# in this sequence (worst adjacent pair dE 20.9 protan, 25.1 normal), and
# swapping blue and purple to adjacent positions fails it. Every rule is also
# directly labelled on the figure, so identity never rests on colour alone.
RULE_COLORS = {'resource': '#bd0026', 'counties': '#1464a8',
               'curve': '#5c8a1e', 'marginal': '#6a3d9a'}

SURFACE = '#fcfcfb'
INK     = '#0b0b0b'
INK2    = '#52514e'
MUTED   = '#b5b4ac'
NEUTRAL = '#f0efec'
ZONE_INK = '#111111'

# Panel A of the county sheet names four layers - cropland, wetlands,
# waterbodies, rivers. That is an IDENTITY job, not a magnitude one, so it
# keeps four distinct hues. A ramp would be the wrong tool there.
CROP_C, WET_C, WB_C   = '#e8e2c8', '#74c476', '#08519c'
RIV_C, FLOW_C, OVER_C = '#2171b5', '#3182bd', '#d94801'

# TYPE. Serif at report sizes, because these go in front of a State Biologist
# and eventually a programme reviewer. matplotlib always ships DejaVu Serif;
# the others are tried first and fall back silently.
SERIF_STACK = ['Charter', 'Georgia', 'Times New Roman', 'STIXGeneral',
               'DejaVu Serif', 'serif']
FS = {           # one place to change every type size in the run
    'base': 12, 'axis_label': 12.5, 'tick': 11.5, 'legend': 11,
    'legend_title': 11.5, 'panel_title': 14, 'fig_title': 17,
    'caption': 11, 'annot': 10.5,
}

# The county sheets' ratio colour bar runs 0 to RATIO_HEADROOM x HOT_MULTIPLE,
# so the black threshold line always sits at the same height on the bar and
# raising HOT_MULTIPLE re-scales every map with it. Anchoring on a percentile
# of the data instead was tried and is worse: one heavy-tailed county drags the
# maximum up and squashes the range everyone else lives in.
RATIO_HEADROOM = 2.67

# Generalisation for DRAWING ONLY. At county scale 25 m is invisible, and it is
# what keeps 120 sheets to minutes instead of an hour. No MEASUREMENT anywhere
# in this pipeline uses it.
DRAW_SIMPLIFY_M = 25.0

# %% [markdown]
# ## 8 — The dataset registry: your input paths
#
# **Every file this pipeline reads is described here and nowhere else.** To
# swap a source, change one row. To add one, add a row.
#
# `pattern` is a glob, or a list of globs tried in order. It is matched at the
# top of `raw_data/`, then anywhere beneath it, then both again
# case-insensitively — downloads arrive as folders as often as as files, with
# download IDs in the names, so an exact literal path is the one thing that
# reliably goes stale. An **absolute path** is also accepted, which is how you
# point at a file that lives outside the project folder.

# %% C8 - where every input comes from
#   layer     the layer inside a GeoPackage or File Geodatabase; None otherwise
#   wanted    columns to keep, matched case-insensitively (NHD ships lowercase,
#             CSB and NWI do not). None keeps everything.
#   newest    ('CDL', r'CDL\d{4}') = keep the highest-numbered matching field,
#             renamed to 'CDL'. Source files add a new year's column each
#             release; naming one freezes the analysis to that year.
#   where     an OGR SQL row filter pushed down to the driver, so a national
#             file never fully enters memory
#   stage     'boundary' is read first and defines the study area;
#             'layer' is cut to it
#   required  False = carry on and say what was lost
#   degrades  plain-English consequence of the file being absent
DATASETS = {

    # ---- boundaries: these define the study area -------------------------
    # Counties are DERIVED from the county-subdivision file by dissolving on
    # COUNTYFP, so there is one less download to keep in step. cousub also
    # carries the town names Pass C needs for the owner-residency split.
    'cousub': dict(
        pattern=[f'tl_*_{STATE_FIPS}_cousub*.shp', 'tl_*_cousub*.shp',
                 '*cousub*.shp'],
        layer=None,
        wanted=['STATEFP', 'COUNTYFP', 'NAME', 'NAMELSAD', 'GEOID'],
        stage='boundary', required=True,
        degrades='nothing runs - the counties are dissolved from this'),

    # The priority area. Supplied as a multi-state file; preprocessing trims it
    # to the state before anything else touches it. REPORTING ONLY while
    # COUNTY_EXTENT_MODE = 'whole_county'.
    'pa': dict(
        pattern=['ACF_PA_*.shp', '*_PA_*.shp', '*priority*area*.shp',
                 '*priority*area*.gpkg'],
        layer=None,
        wanted=None,
        stage='boundary', required=False,
        degrades='pa_share and in_pa are 0 for every county; nothing else changes'),

    'watersheds': dict(
        pattern=['Priority_HUC8*.shp', '*Priority_HUC8*.shp'],
        layer=None,
        wanted=None,
        stage='boundary', required=False,
        degrades='the priority-HUC8 context layer is missing from the maps'),

    # ---- the resource layers ---------------------------------------------
    # CSB is a NATIONAL File Geodatabase. The `where` is what keeps that
    # tractable: the state filter is pushed down to the driver, so only this
    # state's fields are ever read.
    'csb': dict(
        pattern=['CSB*.gdb', '*CSB*.gdb', 'CSB*.gpkg'],
        layer='national1825',
        wanted=['CSBID', 'CSBACRES', 'STATEFIPS', 'CNTY', 'CNTYFIPS'],
        newest=('CDL', r'CDL\d{4}', CROP_SEQ_YEARS),
        where=f"STATEFIPS = '{STATE_FIPS}'",
        stage='layer', required=True,
        degrades='nothing runs - this is the cropland'),

    # One NHD High Resolution GeoPackage holds all three water layers.
    'nhd_flow': dict(
        pattern=['NHD_H_*_State_GPKG.gpkg', 'NHD*.gpkg'],
        layer='NHDFlowline',
        wanted=['permanent_identifier', 'gnis_name', 'ftype', 'fcode',
                'lengthkm', 'reachcode'],
        # 460 stream/river, 558 artificial path, 336 canal/ditch. Everything
        # else NHD calls a flowline (pipelines, coastlines) is not a stream.
        where='ftype IN (460, 558, 336)',
        stage='layer', required=True,
        degrades='nothing runs - these are the streams'),

    'nhd_wb': dict(
        pattern=['NHD_H_*_State_GPKG.gpkg', 'NHD*.gpkg'],
        layer='NHDWaterbody',
        wanted=['permanent_identifier', 'gnis_name', 'ftype', 'fcode',
                'areasqkm'],
        stage='layer', required=False,
        degrades='lakes and ponds contribute no corridor; waterbody ids omitted'),

    'nhd_area': dict(
        pattern=['NHD_H_*_State_GPKG.gpkg', 'NHD*.gpkg'],
        layer='NHDArea',
        wanted=['permanent_identifier', 'gnis_name', 'ftype', 'fcode'],
        stage='layer', required=False,
        degrades='WIDE RIVERS CONTRIBUTE NOTHING - a flowline buffered 100 ft '
                 'from the centre of a 300 ft river is entirely inside the '
                 'channel, so the biggest rivers in the state vanish'),

    # NWI ships one GeoPackage per state, and the layer inside it is named
    # for the state: KY_Wetlands, TN_Wetlands. `layer` may be a LIST - the
    # first name that matches a layer actually in the file wins, and a file
    # with only one layer needs no name at all.
    'wetlands': dict(
        pattern=[f'{HOME_STATE}_*wetlands*.gpkg', '*wetlands*.gpkg',
                 '*Wetlands*.gpkg', '*wetlands*.shp', '*Wetlands*.shp'],
        layer=[f'{HOME_STATE}_Wetlands', f'{HOME_STATE}_Wetlands_Project',
               'Wetlands', 'NWI_Wetlands'],
        wanted=['ATTRIBUTE', 'WETLAND_TYPE', 'ACRES'],
        stage='layer', required=False,
        degrades='wetlands contribute no corridor; wetland ids omitted'),

    # The national WBD ships as WBD_<FIPS>_HU2_GPKG with a WBDHU12 layer; a
    # state agency's own copy is usually named after its catalogue entry. Both
    # spellings are listed, and an unlisted one is reported with the layers the
    # file actually has, so the fix is one line.
    'huc12': dict(
        pattern=['*Hydrologic_Units*.gpkg', '*huc12*.gpkg', '*HUC12*.gpkg',
                 f'WBD_{STATE_FIPS}_*.gpkg', 'WBD*.gpkg', '*huc12*.shp',
                 '*HUC12*.shp', '*WBDHU12*.shp'],
        layer=['WBDHU12', 'Hydrologic_Units_81012', 'HUC12', 'huc12'],
        wanted=['huc12', 'name'],
        stage='layer', required=True,
        degrades='nothing runs - the region basis and the HUC12 rollup need these'),
}

# PARCELS ARE NOT IN THE REGISTRY, ON PURPOSE. Pass C finds them per county by
# matching the county NAME against the file names under raw_data/, because
# parcels are bought one county at a time and the set changes between runs. A
# file called `ky_christian.gpkg` or `christian_parcels.zip` anywhere under
# raw_data/ is found for Christian County. `has_parcels` and `parcel_path` in
# the manifest record what was matched.
PARCEL_EXTENSIONS = ('.gpkg', '.gpkg.zip', '.zip', '.shp')

# %% [markdown]
# ## 9 — Reference tables and units
#
# Facts, not decisions. `WANT_CONTACT` and `OWNER_PAT` are the two you might
# genuinely edit — the first when a parcel vendor renames a field, the second
# when a public owner in your state is not being caught.

# %% C9 - units and CRS
# EPSG:5070 is equal-area for the continental US, so an area measured in it is
# a real area. Every layer is reprojected into it on load and never leaves it.
WORKING_CRS = 'EPSG:5070'
AC   = 4046.8564224        # square metres in an acre
SQMI = 2589988.110336      # square metres in a square mile
MI_M = 1609.344            # metres in a mile

# Preprocessing generalisation. Applied ONCE, at cache time, and before the
# spatial test rather than after - a spatial test costs vertices on both sides,
# so simplifying after the test pays the full cost and saves nothing.
#   'select' keeps a feature WHOLE if any part of it is in the study area
#   'clip'   cuts it at the boundary
# 'select' is faster and the better input: the acres the analysis reports come
# from its own clip to the analysed ground, so this affects the cache, not the
# answer.
CLIP_MODE       = 'select'
SIMPLIFY_M      = 1
MASK_SIMPLIFY_M = 30

# %% C9b - contact fields, owner screens, crop classes
# Parcel vendor field names vary by vintage, so they are DETECTED, never
# assumed. First match wins. Add a spelling to a row if your extract uses one.
WANT_CONTACT = {
    'owner':        ['owner', 'owner1', 'ownername', 'own1'],
    'owner2':       ['owner2', 'ownername2'],
    'mail_address': ['mailadd', 'mail_address', 'mailaddress', 'maddress'],
    'mail_city':    ['mail_city', 'mailcity', 'mcity'],
    'mail_state':   ['mail_state2', 'mail_state', 'mailstate', 'mstate'],
    'mail_zip':     ['mail_zip', 'mailzip', 'mzip'],
    'site_address': ['address', 'saddress', 'situs_address',
                     'original_address', 'physaddress'],
    'site_city':    ['scity', 'city', 'situs_city'],
    'site_zip':     ['szip', 'zip', 'situs_zip'],
    'parcel_id':    ['parcelnumb', 'parcel_id', 'pin', 'apn'],
    'land_use':     ['usedesc', 'usecode', 'landuse'],
    'county':       ['county', 'cntyname'],
}

# Owner-name screens. THESE READ NAME TEXT, so they find a parcel owned by an
# agency that says so and CANNOT find one whose owner field is blank. A missing
# record is a data-acquisition backlog item, not evidence of public ownership,
# which is why blanks are classed `unknown_owner` and KEPT on the list.
OWNER_PAT = {
    'federal': r'\b(?:UNITED STATES|U\.?S\.?A?\b|FEDERAL|ARMY|DEPT OF DEFENSE|'
               r'DEPARTMENT OF DEFENSE|FORT CAMPBELL|CORPS OF ENGINEERS|'
               r'USDA|FOREST SERVICE|FISH (?:AND|&) WILDLIFE|NATIONAL PARK)\b',
    'state_local': r'\b(?:COMMONWEALTH|STATE OF|KENTUCKY (?:DEPT|DEPARTMENT|'
                   r'TRANSPORTATION)|KY DEPT|COUNTY OF|CITY OF|FISCAL COURT|'
                   r'HOUSING AUTHORITY|WATER DISTRICT|SEWER|UTILITY|'
                   r'BOARD OF EDUCATION|SCHOOL (?:DIST|BOARD)|PARK)\b',
    'institutional': r'\b(?:CHURCH|BAPTIST|METHODIST|PRESBYTERIAN|CATHOLIC|'
                     r'CEMETERY|UNIVERSITY|COLLEGE|HOSPITAL|FOUNDATION|'
                     r'CONSERVANCY|LAND TRUST|DUCKS UNLIMITED)\b',
}
# land_use is a WEAK signal, so it only speaks for rows with NO owner name at
# all, and it flags them 'review' - never 'public'.
USE_PAT = r'\b(?:EXEMPT|PUBLIC|GOVERNMENT|MILITARY|RIGHT[ -]OF[ -]WAY)\b'

# Which water-feature id column each source actually carries. The FALLBACK list
# is types or names rather than unique ids, and the run prints which was used -
# a count of distinct wetland TYPES is not a count of wetlands.
ID_SOURCES = [
    ('stream',    'n_stream_segments', 'nhd_flow',
     ['permanent_identifier', 'reachcode'], ['gnis_name']),
    ('river',     'n_rivers', 'nhd_area',
     ['permanent_identifier'], ['gnis_name']),
    ('waterbody', 'n_waterbodies', 'nhd_wb',
     ['permanent_identifier'], ['gnis_name']),
    ('wetland',   'n_wetlands', 'wetlands',
     ['GLOBALID', 'OBJECTID', 'WETLAND_ID'], ['ATTRIBUTE', 'WETLAND_TYPE']),
]

# CDL classes rolled up for the buffer-cover table. Anything not listed is
# treated as row crop.
CDL_GROUPS = {
    'pasture':   {62, 176},
    'hay':       {36, 37, 58, 59, 60},
    'fallow':    {61},
    'vegetated': {63, 64, 65, 87, 88, 112, 131, 141, 142, 143, 152, 171,
                  190, 195},
    'water':     {83, 111},
    'developed': {121, 122, 123, 124},
}
GRP_ORDER = ['row_crop', 'pasture', 'hay', 'fallow', 'vegetated', 'water',
             'developed']

CDL_NAMES = {
    1: 'Corn', 2: 'Cotton', 3: 'Rice', 4: 'Sorghum', 5: 'Soybeans',
    6: 'Sunflower', 10: 'Peanuts', 11: 'Tobacco', 12: 'Sweet Corn',
    21: 'Barley', 22: 'Durum Wheat', 23: 'Spring Wheat', 24: 'Winter Wheat',
    25: 'Other Small Grains', 26: 'Dbl WinWht/Soy', 27: 'Rye', 28: 'Oats',
    29: 'Millet', 31: 'Canola', 33: 'Safflower', 36: 'Alfalfa',
    37: 'Other Hay', 41: 'Sugarbeets', 42: 'Dry Beans', 43: 'Potatoes',
    44: 'Other Crops', 53: 'Peas', 59: 'Sod/Grass Seed', 61: 'Fallow/Idle',
    205: 'Triticale', 216: 'Peppers', 221: 'Strawberries',
    225: 'Dbl WinWht/Corn', 226: 'Dbl Oats/Corn', 228: 'Dbl Triticale/Corn',
    229: 'Pumpkins', 236: 'Dbl WinWht/Sorghum', 237: 'Dbl Barley/Corn',
    238: 'Dbl WinWht/Cotton', 239: 'Dbl Soy/Cotton', 240: 'Dbl Soy/Oats',
    241: 'Dbl Corn/Soy', 254: 'Dbl Barley/Soy', 111: 'Open Water',
    121: 'Developed/Open', 122: 'Developed/Low', 123: 'Developed/Med',
    124: 'Developed/High', 131: 'Barren', 141: 'Deciduous Forest',
    142: 'Evergreen Forest', 143: 'Mixed Forest', 152: 'Shrubland',
    176: 'Grass/Pasture', 190: 'Woody Wetlands', 195: 'Herbaceous Wetlands',
}

# %% C9b - county names
# WHERE COUNTY NAMES COME FROM, in order. Nothing here is state-specific, and
# nothing but labels, folder names and the parcel-file search reads them.
#
#   1. the county-subdivision file (cousub), if it carries a county-name field
#   2. CSB's CNTY, which is already in the cache and costs nothing
#   3. a TIGER county file in raw_data/  (tl_*_county*.shp)  if you downloaded one
#   4. the Census national county list, fetched once and cached beside the
#      GeoPackage as county_names_<FIPS>.json          (COUNTY_NAME_FETCH)
#   5. COUNTY_NAMES below, if you filled it in
#   6. 'County 047', so a nameless county still gets a row and a folder
#
# TIGER's cousub files carry COUNTYFP but usually NOT the county's name - the
# NAME column is the subdivision's. Step 1 uses it only where a file really
# does have a county-name field; it is step 2 that names the counties in
# practice, because CSB is a required input and reaches every county in the
# cache.
COUNTY_NAME_FETCH = 'auto'   # 'auto' fetch once if steps 1-3 left gaps and the
                             #        machine has a network; False never fetch
COUNTY_NAME_URL = ('https://www2.census.gov/geo/docs/reference/codes2020/'
                   'national_county2020.txt')

# Your own overrides, or a whole table for an offline machine:
#   COUNTY_NAMES = {'047': 'Christian', '225': 'Union'}
COUNTY_NAMES = {}

# %% C10 - derived paths and self-checks
# Nothing below is a decision. It turns the switches above into the folders the
# two scripts write into, and it fails HERE rather than four hours in.
HALO_M   = HALO_MI * MI_M
CELL_M   = CELL_MI * MI_M
BAND_M   = BAND_MI * MI_M
ZONE_NEIGHBOR_MI = BAND_MI if ZONE_NEIGHBOR_MI is None else ZONE_NEIGHBOR_MI

S_OUT = os.path.join(WORK_DIR, '_state', 'output')
S_VAL = os.path.join(WORK_DIR, '_state', 'validation')
S_FIG = os.path.join(WORK_DIR, '_state', 'figures')
C_DIR = os.path.join(WORK_DIR, 'counties')
P_VAL = os.path.join(WORK_DIR, '_preprocess', 'validation')
P_FIG = os.path.join(WORK_DIR, '_preprocess', 'figures')

MANIFEST_CSV = os.path.join(WORK_DIR, 'county_manifest.csv')
SURFACE_GPKG = os.path.join(S_OUT, 'state_surface.gpkg')
SURFACE_JSON = os.path.join(S_OUT, 'state_surface.json')
TILE_DIR = os.path.join(FINAL_DIR if TILES_ON_FINAL else WORK_DIR,
                        '_state', 'output', '_tiles')
TILE_SIG = os.path.join(TILE_DIR, '_signature.json')

_MODES = (True, False, 'auto')
assert RUN_PREPROCESS in _MODES, f'RUN_PREPROCESS must be one of {_MODES}'
assert RUN_SURFACE in _MODES, f'RUN_SURFACE must be one of {_MODES}'
assert BASIS in ('county', 'state', 'region'), f'bad BASIS: {BASIS!r}'
assert FLOOR_MODE in ('percentile', 'knee', 'absolute', 'off'), \
    f'bad FLOOR_MODE: {FLOOR_MODE!r}'
assert FLOOR_KNEE_RULE in ('resource', 'counties', 'curve', 'marginal'), \
    f'bad FLOOR_KNEE_RULE: {FLOOR_KNEE_RULE!r}'
assert FLOOR_KNEE_MARGINAL_SMOOTH >= 1, 'FLOOR_KNEE_MARGINAL_SMOOTH >= 1'
assert 0 < FLOOR_KNEE_GRID <= 10, 'FLOOR_KNEE_GRID is a percentile step'
assert 0 < FLOOR_KNEE_PLATEAU_TOL < 1, \
    'FLOOR_KNEE_PLATEAU_TOL is a fraction of each rule\'s net range'
assert 0 < FLOOR_KNEE_MAX_PCTL <= 100
# 'pa_where_available' is legacy (section 3) and is kept in this tuple ON
# PURPOSE: dropping it would turn a surface built under it from a reported
# signature mismatch into an exception on load.
assert COUNTY_EXTENT_MODE in ('pa_where_available', 'pa_only', 'whole_county')
assert isinstance(REPORT_WITHIN_PA, bool), 'REPORT_WITHIN_PA is True or False'
# Whether there is a priority area to mask WITH cannot be known here - this
# file never reads data. B6 raises on it, where the answer is in hand and the
# message can say which file to put where.
assert CLIP_MODE in ('select', 'clip')
assert 0.0 <= W_PLACE <= 1.0, 'W_PLACE is a weight between 0 and 1'
assert set(USE_SOURCES) <= set(WATER_SOURCES), \
    'USE_SOURCES names a source that is not in WATER_SOURCES'


def ensure_dirs():
    """Create the output folders. Called by each script, never by import, so
    reading the config cannot leave folders behind on someone's disk."""
    for _p in (CLIPPED_DIR, S_OUT, S_VAL, S_FIG, C_DIR, P_VAL, P_FIG,
               TILE_DIR):
        os.makedirs(_p, exist_ok=True)


def describe():
    """One screen that says what this run is. Printed by both scripts, so the
    top of any log answers "what settings produced this?"."""
    print(f'{"=" * 74}\nCONFIGURATION\n{"=" * 74}')
    print(f'  running on   {"Google Colab" if COLAB else "a local machine"}'
          + ('  (Drive mounted, building on local disk)'
             if COLAB and SCRATCH_DIR else ''))
    print(f'  project      {PROJECT_DIR}')
    print(f'  raw inputs   {RAW_DIR}'
          + ('' if os.path.isdir(RAW_DIR) else '   <- DOES NOT EXIST'))
    print(f'  cache        {CACHE}'
          + ('' if os.path.exists(CACHE) else '   (not built yet)'))
    print(f'  outputs      {WORK_DIR}'
          + (f'\n               -> copied at the end to {FINAL_DIR}'
             if WORK_DIR != FINAL_DIR else ''))
    print(f'  switches     RUN_PREPROCESS={RUN_PREPROCESS!r}  '
          f'RUN_SURFACE={RUN_SURFACE!r}  RUN_PASS_C={RUN_PASS_C}  '
          f'SURFACE_GUARD={SURFACE_GUARD}')
    print(f'  study area   state {STATE_FIPS} ({STATE_NAME}, {HOME_STATE}) · '
          f'{COUNTY_WORDS} {COUNTIES} · extent {COUNTY_EXTENT_MODE}'
          + ('\n               REPORTED INSIDE THE PRIORITY AREA '
             '(REPORT_WITHIN_PA) - the surface stays statewide'
             if REPORT_WITHIN_PA else ''))
    print(f'  surface      {CELL_MI:g}-mile cells · {BAND_MI:g}-mile band · '
          f'{RIP_BUFFER_FT:.0f} ft buffer · {TILE_KM:g} km tiles')
    print(f'  threshold    basis {BASIS} · ratio {HOT_MULTIPLE:g}x · floor '
          f'{FLOOR_MODE}'
          + (f' ({FLOOR_PCTL:g}th pctl)' if FLOOR_MODE == 'percentile'
             else f' (knee, {FLOOR_KNEE_RULE!r} rule)' if FLOOR_MODE == 'knee'
             else f' ({FLOOR_AC_PER_SQMI:g} ac/sq mi)'
             if FLOOR_MODE == 'absolute' else ''))
    print(f'  parcels      floor {MIN_RIP_AC:g} ac · W_PLACE {W_PLACE:g} · '
          + (f'top {TOP_N_CONTACTS}' if TOP_N_CONTACTS
             else 'ranked list only'))
    print(f'  CRS          {WORKING_CRS} (equal-area)')


if __name__ == '__main__':
    describe()
