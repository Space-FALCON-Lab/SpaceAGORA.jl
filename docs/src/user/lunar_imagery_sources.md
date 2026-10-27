# Building lunar terrain and imagery

The optional Python tools in `scripts/dev/terrain/` build a local terrain
bundle for the viewer and the lunar terrain loader. Ordinary simulations need
neither these tools nor their Python packages or downloads.

Run the commands from the repository root. The tools need Python 3.9 or later
with NumPy and Pillow; if those are not installed, use a separate environment:

```text
python3 -m venv .venv-terrain
.venv-terrain/bin/python -m pip install -r scripts/dev/terrain/requirements.txt
```

The commands below call `python3`; substitute that environment's Python when
needed. On Windows its executable is `.venv-terrain\Scripts\python.exe`.

## Start with elevation data

This command requests LOLA elevation data around the configured Apollo 11
site and skips the optional NAC terrain model and the imagery:

```text
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --name apollo11 --out output/terrain/moon/apollo11 --no-nac --no-imagery
```

`--lola-half-deg` sets the half-width of the local elevation window in
degrees. A second, coarser LOLA window covers the whole imagery root, which is
sized to contain the descent corridor, so its extent depends on
`--root-zoom`, `--uprange-km`, `--approach-azimuth` and the altitude profile.
`--profile` takes a JSON file holding a list of
`[distance_to_go_km, height_above_site_km]` pairs; the tool rescales the
distances so that the longest equals `--uprange-km`. Without `--profile` it
uses the built-in Apollo 11 descent profile. That profile is a template for
visualization coverage, not a new validation of the Apollo flight
reconstruction.

To request imagery as well, leave out `--no-imagery`; to request the
registered Apollo 11 NAC DTM, leave out `--no-nac`. The NAC download fetches
the complete source TIFF, which is larger than the local subset it keeps.
LOLA and the tiled PDS imagery are read with bounded HTTP range requests.
Source availability and rate limits can stop a download; whatever was already
cached is kept for the next run.

The NAC registry in the tool covers only Apollo 11. Another site needs its
own verified product entry there. Without one, pass `--no-nac`, or give the
site its own `--name`, which the tool reports as unregistered and builds
without a NAC terrain model. The NAC reader accepts only the layout it was
written for: a north-up equirectangular GeoTIFF on the lunar sphere with a
2 m grid. It rejects any other layout rather than interpret it approximately.

## Files and height convention

`site.json` lists the DEMs, finest first, and the imagery index when imagery
was built. Each DEM (`dem_nac`, `dem_lola` and `dem_lola_wide`) is a JSON
header plus a `.f32` binary array. The values are elevations in meters above
a 1,737,400 m sphere, not radii from the Moon's center. Coordinates are
planetocentric latitude and east-positive longitude. Columns run west to east
and rows north to south; the values are little-endian Float32, stored
row-major.

The header's latitude and longitude limits are the outer edges of the cells,
and the samples lie at cell centers. A regional window that crosses the
longitude seam keeps continuous bounds, so they may fall outside the usual
longitude interval. Use the same reference radius when you connect the grid to
guidance, touchdown or plume models.

The coarse LOLA window averages the source columns and up to four evenly
spaced source rows inside each output cell. This is a stratified
approximation, not a full two-dimensional area mean. The NAC stride keeps the
coordinates of the source cell centers it selects. Invalid or missing
elevations are rejected; the tool does not fill them with invented terrain.

## Reusing and rebuilding results

A DEM is reused only when the request, producer version, datum, shape and
binary checksum all match. A bundle that does not match must be regenerated in
a fresh output directory. `--reuse DIR` borrows matching DEMs and source
caches from another bundle. Borrowed DEM files are copied, so the new bundle
stays self-contained; raw download caches may be linked instead, to avoid
duplicating large source rasters. A borrowed DEM pair is validated before it
is copied, and a mismatched source leaves the fresh output directory
untouched, so a retry without `--reuse` can still build the DEM.

Moon Trek tiles are reused only if they decode completely at the expected
size, and archive tiles only if they have the byte count that the raster index
declares. Invalid entries are skipped and fetched again into the current
output's cache, without modifying a borrowed cache. Cache writes publish
complete files atomically. A marker file records a Trek HTTP 404; empty files
left by older versions of the tool are fetched again, because they could be
interrupted writes. An invalid image response raises an error instead of
leaving imagery silently missing.

Imagery goes into `imagery/tiles.json` and JPEG tiles. Reuse checks the whole
build request, including the track, the source registry, the registration
corrections, the texture and the image settings, and it checks every
referenced tile's checksum. If any of those inputs change, pass `--retile` or
use a fresh output directory. `--retile` builds the replacement in a
temporary directory and keeps the previous imagery if the build fails. It does
not delete the reusable source caches.

The stub tool, `stub_quadtree_site.py`, is for older bundles made of nested
squares (`imagery/imagery.json`), or for a development fixture cut from a
local global texture. It does not re-cut a bundle that the current tool has
already written as a quadtree. It refuses to write into existing nonempty tile
output or to use the source site as its output directory. Both tools require
imagery for the root tile. If the stub reports that its source is too coarse,
enlarge `--root-deg`, reduce `--tile-px`, or supply a finer source that covers
the root.

## Registration and resolution

`check_site_registration.py` compares archive imagery with a hillshade of a
DEM, optionally with Moon Trek imagery, and with a bright-object heuristic
near the given site coordinate. The example below needs `dem_nac.json`, which
the LOLA-only command above does not create. For the configured Apollo 11
site, first rerun `fetch_moon_site.py` without `--no-nac` to produce the NAC
DEM (you can keep `--no-imagery`), then run:

```text
python3 scripts/dev/terrain/check_site_registration.py --site 0.67416 23.47314 --dem output/terrain/moon/apollo11/dem_nac.json --cache output/terrain/moon/apollo11/pds_cache
```

The tool can fetch archive data that is not cached yet. It reports shifts in
meters north and east, computed with separate latitude and longitude
spacings. Correlation peaks and the bright-object heuristic are diagnostics;
they do not show that either product is an absolute position reference. Check
the registration again after changing a product, the site or the datum,
before you treat the imagery as aligned.

The `LAYER_REGISTRATION` offsets in `fetch_moon_site.py` come from the
original Apollo 11 study. They are configured inputs, not results remeasured
for this repository. The deepest imagery level may be sampled more finely than
its source resolves. A sharp texture in the viewer is not evidence of terrain
accuracy, and the tone and detail processing applied to the imagery does not
improve the physical DEM.

## Sources and attribution

The tool's source registries record the product URLs and source metadata:

- [LOLA elevation convention](https://pds-geosciences.wustl.edu/missions/lro/lola_faq.htm): each integer sample times 0.5 gives meters above the 1,737,400 m sphere.
- [LOLA GDR archive](https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/): the `LDEM_128` cylindrical elevation raster.
- [NASA PDS cloud archive](https://pds.mcp.nasa.gov/), LROC PDS4 products under `data/store/img/lunar_reconnaissance_orbiter/pds4/lroc/` (the directory itself serves no listing): the Apollo 11 NAC DTM and the `NAC_PHO_E010N0230_M175124932R` browse imagery.
- [NASA Moon Trek](https://trek.nasa.gov/moon/): the configured WAC, Kaguya and Apollo 11 imagery layers.
- [LROC terms and attribution](https://lroc.im-ldi.com/about/terms): read the terms for each product before you redistribute an imagery bundle, and keep the product URLs, credits and citations that the tool records.

The Apollo 11 coordinate comes from the original study, which used
Wagner et al. (2017), *Icarus* 283, 92–103,
[doi:10.1016/j.icarus.2016.05.011](https://doi.org/10.1016/j.icarus.2016.05.011).
The LROC imagery metadata keeps the NASA/GSFC/Arizona State University credit
and the citation Robinson et al. (2010), *Space Science Reviews* 150, 81–124.

## What is tested

Run the deterministic Python tests, which download nothing:

```text
python3 -m unittest discover -s test/scripts -p test_terrain_pipeline.py -v
```

They cover synthetic raster geometry, cell-center preservation, longitude
seams, rejection of bad HTTP range responses, cache identity, preservation of
the previous output when a rebuild fails, registration distance units and
legacy-stub output. One test runs the LOLA-only command end to end against a
small synthetic raster in place of the network.

Because nothing is downloaded, the tests do not cover a live NAC download, the
registration of the current source products, or how the imagery looks in the
viewer. Those need a separate check with the real assets.
