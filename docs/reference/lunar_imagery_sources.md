# Building lunar terrain and imagery

The optional Python tools in `scripts/dev/terrain/` prepare a local terrain
bundle for the viewer and the lunar terrain loader. Ordinary simulations do
not need these tools, Python packages, or network downloads.

Run commands from the repository root. Use a separate Python environment if
NumPy and Pillow are not already available:

```text
python3 -m venv .venv-terrain
.venv-terrain/bin/python -m pip install -r scripts/dev/terrain/requirements.txt
```

The commands below use `python3`; substitute that environment's Python when
needed. On Windows its executable is `.venv-terrain\Scripts\python.exe`.
Python 3.9 or later is required.

## Start with elevation data

This example requests LOLA elevation data around the configured Apollo 11 site
and skips the optional NAC raster and imagery:

```text
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --name apollo11 --out output/terrain/moon/apollo11 --no-nac --no-imagery
```

`--lola-half-deg` controls the local elevation window. A second, coarser LOLA
window covers the configured descent corridor. Its extent depends on
`--root-zoom`, `--uprange-km`, `--approach-azimuth` and the altitude profile.
`--profile` accepts a JSON list of `[distance_to_go_km, height_above_site_km]`
pairs. The built-in profile is a visualization coverage template, not a new
validation of Apollo flight reconstruction.

To request imagery as well, omit `--no-imagery`. To request the registered
Apollo 11 NAC DTM, omit `--no-nac`. NAC download reads the complete source TIFF,
which is larger than a local subset. LOLA and the tiled PDS imagery use bounded
HTTP range reads. Source availability and rate limits can stop a download;
partial caches are retained for later reuse.

The supplied NAC registry covers Apollo 11. Another site needs its own
verified product entry or `--no-nac`. The reader supports the explicitly
validated lunar sphere and north-up equirectangular GeoTIFF layout; an
unsupported layout is rejected rather than interpreted approximately.

## Files and height convention

`site.json` lists the DEMs and optional imagery index. Each DEM has a JSON
header and a `.f32` binary array. Values are **elevations in metres above a
1,737,400 m sphere**, not planet-centred radii. The coordinates are
planetocentric latitude and east-positive longitude. Columns run west to east;
rows run north to south. Float32 values are little-endian and row-major.

Header latitude and longitude limits are the outer edges of the cells. Samples
lie at cell centres. Regional longitude bounds may extend beyond the usual
longitude interval to keep a region crossing the seam continuous. Use the
same reference radius when connecting the grid to guidance, touchdown, or
plume models.

The coarse LOLA producer averages source columns and a bounded number of evenly
spaced rows inside each output cell. This is a stratified approximation, not a
full two-dimensional area mean. The NAC stride preserves the coordinates of
the selected source centres. Invalid or missing elevations are rejected;
the producer does not replace them with invented terrain.

## Reusing and rebuilding results

DEM reuse requires the same request, producer version, datum, shape and binary
checksum. A mismatched bundle must be regenerated in a fresh output directory.
`--reuse DIR` can borrow matching DEMs and source caches from another bundle.
Exported DEM files are copied into the new bundle so it remains self-contained;
raw download caches may be linked to avoid duplicating large source rasters.

Imagery is stored in `imagery/tiles.json` and JPEG tiles. Reuse checks the build
request, including the track, source registry, registration corrections,
texture and image settings, and checks every referenced tile's checksum.
Changing those inputs requires `--retile` or a fresh output directory.
`--retile` builds a replacement in a temporary directory and keeps the previous
imagery if building fails. It does not delete the reusable source caches.

The stub tool, `stub_quadtree_site.py`, is for **older nested-square bundles**
with `imagery/imagery.json`, or for a development fixture made from a local
global texture. It does not convert the current producer's quadtree again.
It refuses to overwrite existing nonempty tile output or use the source site
as its output directory. Both producers require a root image. If the stub
reports that its source is too coarse, enlarge `--root-deg`, reduce
`--tile-px`, or supply a finer source covering the root.

## Registration and resolution

`check_site_registration.py` compares archive imagery with a DEM hillshade,
optional Moon Trek imagery, and a bright-object heuristic near the supplied
site coordinate. For example:

```text
python3 scripts/dev/terrain/check_site_registration.py --site 0.67416 23.47314 --dem output/terrain/moon/apollo11/dem_nac.json --cache output/terrain/moon/apollo11/pds_cache
```

It can fetch missing archive data. Reported shifts are metres north and east,
with separate latitude and longitude spacings. Correlation peaks and the
bright-object heuristic are diagnostics, not proof that one product is an
absolute position reference. Recheck registration after changing a product,
site, or datum before treating the imagery as aligned.

The inherited `LAYER_REGISTRATION` constants describe the original Apollo 11
study. They are retained as configured inputs, not independently remeasured
acceptance results for this extraction. The deepest imagery level may sample
more finely than the source resolves. Viewer texture sharpness is not evidence
of terrain accuracy, and image tone/detail processing does not improve the
physical DEM.

## Sources and attribution

The producer's named registries record product URLs and source metadata:

- [LOLA elevation convention](https://pds-geosciences.wustl.edu/missions/lro/lola_faq.htm): each integer sample is multiplied by 0.5 to obtain metres above the 1,737,400 m sphere.
- [LOLA GDR archive](https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/): the `LDEM_128` cylindrical elevation raster.
- [NASA PDS LROC archive](https://pds.mcp.nasa.gov/data/store/img/lunar_reconnaissance_orbiter/pds4/lroc/): the Apollo 11 NAC DTM and `NAC_PHO_E010N0230_M175124932R` browse imagery.
- [NASA Moon Trek](https://trek.nasa.gov/moon/): the configured WAC, Kaguya and Apollo 11 imagery layers.
- [LROC terms and attribution](https://lroc.im-ldi.com/about/terms): consult the applicable product terms before redistributing an imagery bundle. Preserve the producer's product URLs, credits and citations.

The Apollo 11 coordinate is inherited from the original study's use of
Wagner et al. (2017), *Icarus* 283, 92–103,
[doi:10.1016/j.icarus.2016.05.011](https://doi.org/10.1016/j.icarus.2016.05.011).
LROC imagery metadata retains the NASA/GSFC/Arizona State University credit and
Robinson et al. (2010), *Space Science Reviews* 150, 81–124.

## What is tested

Run the deterministic Python tests without downloads:

```text
python3 -m unittest discover -s test/scripts -p test_terrain_pipeline.py -v
```

They exercise synthetic raster geometry, cell-centre preservation, longitude
seams, HTTP range rejection, cache identity, failed-rebuild preservation,
registration distance units, and legacy-stub output. A bounded real 4×4 LOLA
subset was also read from the current default WUSTL endpoint with TLS
verification enabled. Its 16 elevations match the earlier MIT-mirror sample
and checksum. This is a source-subset check, not a complete live CLI run.

Live NAC download, current source registration and browser appearance require
separate acceptance with real assets. Passing the synthetic tests does not
establish those results. No new live product measurements or physical
calibration are claimed by this guide.
