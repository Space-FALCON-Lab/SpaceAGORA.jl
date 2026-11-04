#!/usr/bin/env python3
"""Build lunar DEM grids and optional imagery for SpaceAGORA.

Requires Python 3.9+, NumPy and Pillow. Run with --help for the options and see
`docs/reference/lunar_imagery_sources.md` for the datum, cache and source contract.
DEM samples are elevations above the 1737400 m sphere on north-to-south,
west-to-east cell-centred grids with explicit outer edges. LOLA uses bounded
HTTP ranges; coarse rows are stratified samples, not full-area averages.
NAC DTM download is optional (--no-nac); --no-imagery skips image requests.
Reusing imagery requires identical build inputs and unchanged tile checksums.
--retile stages a replacement and preserves the previous imagery on failure.
"""
import argparse, hashlib, io, json, math, os, pathlib, shutil, struct, tempfile, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from PIL import Image

MOON_RADIUS_M = 1737400.0
M_PER_DEG = math.pi / 180.0 * MOON_RADIUS_M
LOLA_URL = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_128.img"
LOLA_PPD, LOLA_LINES, LOLA_SAMPLES, LOLA_RECORD_BYTES = 128, 23040, 46080, 92160
NAC_DTM = {
    # site name: (TIF url, label geometry). Add sites here as their NAC DTMs are found on the PDS.
    "apollo11": {
        "url": "https://pds.mcp.nasa.gov/data/store/img/lunar_reconnaissance_orbiter/pds4/lroc/lro-l-lroc-5-rdr/LROLRC_2001/DATA/SDP/NAC_DTM/APOLLO11/NAC_DTM_APOLLO11.TIF",
        "lines": 13978, "samples": 2111, "res_px_per_deg": 15161.67521207, "center_lat_deg": 1.0, "center_lon_deg": 180.0,
        "line_offset": 18747.5, "sample_offset": 2374376.5, "lat_min": 0.314609, "lat_max": 1.2365388, "lon_min": 23.3722758, "lon_max": 23.5115296,
        "nodata": -3.4028226550889045e+38,
    },
}
TREK = "https://trek.nasa.gov/tiles/Moon/EQ/{layer}/1.0.0/default/default028mm/{z}/{y}/{x}.{ext}"
# Moon Trek raster layers, coarse to fine. bbox is (lon_min, lat_min, lon_max, lat_max) from the
# Trek product index and zmax the deepest zoom the service serves; a request outside either is
# skipped rather than sent, and a tile the service does not have comes back 404 and lets the next
# coarser layer through.
TREK_LAYERS = {
    "wac":    {"layer": "LRO_WAC_Mosaic_Global_303ppd_v02", "ext": "jpg", "zmax": 8,  "bbox": (-180.0, -90.0, 180.0, 90.0)},
    "kaguya": {"layer": "Kaguya_TCortho_Mosaic_Global_4096ppd", "ext": "png", "zmax": 10, "bbox": (-180.0, -90.0, 180.0, 90.0)},
    "a11_60": {"layer": "A11_60x60km.eq", "ext": "png", "zmax": 15, "bbox": (22.4806483, -0.3193414, 24.459338, 1.6593483)},
    "a11_nac": {"layer": "LRO_NAC_Apollo11_Mosaic_p", "ext": "png", "zmax": 15, "bbox": (23.4484977, 0.1464484, 23.5396839, 1.1149051)},
    "a11_26cm": {"layer": "apollo11_26cm_mosaic_byte_geo_1_2_highContrast", "ext": "png", "zmax": 15, "bbox": (23.4484963, 0.1464532, 23.5396925, 1.1149077)},
}
# ---- PDS archive rasters -----------------------------------------------------------------
# Moon Trek caps every layer it serves (see docs/reference/lunar_imagery_sources.md): its own
# Apollo 11 mosaic is held at about 0.26 m/px and served at 0.651, so the deepest Trek level of
# the quadtree is 2.5x coarser than the mosaic behind it. The archive has no such cap. These are
# the products read straight out of the PDS, by HTTP range requests into the raster, for the
# levels Trek cannot reach. Only PDS-archived rasters belong here: the LROC terms
# (https://lroc.im-ldi.com/about/terms) put the PDS products in the public domain and keep the
# copyright on everything served through the project's own non-PDS interfaces.
ARCHIVE_SOURCES = {
    # The map-projected browse raster of LROC NAC frame M175124932R, imaged 2011-11-05 from a
    # 24 km low-altitude pass - the finest frame of the 110 the Orbital Data Explorer lists over
    # this site. 8-bit, uncompressed, internally tiled 256 x 256, 915 MB, so a window of it is a
    # few contiguous byte ranges. Bounds and pixel size are cross-checks read from the PDS4 label
    # of the companion .IMG; the reader takes its geometry from the TIFF itself.
    "a11_pho_r": {
        "product": "NAC_PHO_E010N0230_M175124932R",
        "url": "https://pds.mcp.nasa.gov/data/store/img/lunar_reconnaissance_orbiter/pds4/lroc/"
               "lro-l-lroc-5-rdr/LROLRC_2001/EXTRAS/BROWSE/NAC_PHO/E010N0230/"
               "NAC_PHO_E010N0230_M175124932R.PYR.TIF",
        # west, south, east, north, from the PDS4 label of NAC_PHO_E010N0230_M175124932R.IMG
        "bbox": (23.448279815425, 0.31671962244668, 23.496770554946, 1.1142051945617),
        "label_m_per_px": 0.2302,
        "credit": "NASA/GSFC/Arizona State University",
        "citation": "LROC NAC frame M175124932R (2011-11-05, 24 km altitude), PDS LRO-L-LROC-5-RDR-V1.0, "
                    "volume LROLRC_2001, product NAC_PHO_E010N0230_M175124932R; "
                    "Robinson et al. (2010), Space Sci. Rev. 150, 81-124",
    },
}
# TIFF tags this reader needs, and the GeoTIFF keys it checks.
_TIFF_TYPE_BYTES = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 11: 4, 12: 8}
_TIFF_TYPE_FMT = {1: "B", 3: "H", 4: "I", 11: "f", 12: "d"}


# Which layers to try at each Trek zoom, finest first. The WAC global mosaic is the sharpest
# source through zoom 8 (its native 100 m/px); Kaguya's TC ortho mosaic carries the corridor to
# zoom 10 (20.8 m/px); the Apollo 11 mosaics take the last five zooms at the site itself.
ZOOM_LAYERS = {
    9: ["a11_60", "kaguya"], 10: ["a11_60", "kaguya"],
    11: ["a11_nac", "a11_60"], 12: ["a11_nac", "a11_60"], 13: ["a11_nac", "a11_60"],
    14: ["a11_26cm", "a11_nac", "a11_60"], 15: ["a11_26cm", "a11_nac", "a11_60"],
    # Trek serves nothing past zoom 15 (0.651 m/px), so zoom 16 comes from the PDS itself:
    # the browse raster of the frame Trek's own Apollo 11 mosaic was made from, read in place
    # by range requests (ARCHIVE_SOURCES). A node the archive does not cover is simply not
    # built, and the viewer inherits it from level 12.
    16: ["a11_pho_r"],
}
ZOOM_LAYERS_DEFAULT = ["wac"]
# Where each layer draws the ground, in meters north and east of where it really is. Trek serves
# these mosaics on different control, and a node takes its detail from whichever layer reaches its
# zoom, so an uncorrected offset puts the same crater in two places either side of a level change -
# a seam no gutter, UV inset or tone match can touch. Measured here, over this site's own tiles:
#   a11_26cm against a11_nac, phase correlation of their common Trek tiles, 0.0 px at zoom 14 and
#     at zoom 15 (peaks 0.85 and 0.70): the two NAC mosaics share control.
#   a11_60 against them, the same way: 28 px at zoom 14 and 56 px at zoom 15, both 36.4 m, plus
#     about 5.9 m east.
#   The NAC mosaics against the site's own PDS NAC digital terrain model (the imagery resampled
#     onto the DTM grid and correlated against its hillshade; the peak is at sun azimuth 90 deg,
#     which is the illumination of the Apollo 11 approach): 24 m south, 4 m east. The PDS raster
#     carries the LROC control that the landing coordinate itself comes from.
# The two routes agree: they put A11_60x60km 12 m north of the truth, and 36 m north of the NAC
# mosaics. Kaguya and the WAC global mosaic are left alone: a 24 m error is a third of a pixel at
# zoom 10 and a twentieth at zoom 8, so there is nothing to measure and nothing to correct.
#   The PDS browse raster the archive level reads (a11_pho_r) against all three: 0.6 m south and
#     1.9 m west of the corrected Trek NAC mosaics (imagery against imagery, the same ground at
#     Trek zoom 15, correlation peak 0.71 and the same answer to 0.02 m over boxes of 167, 250
#     and 333 m); 0.9 m north and 3.2 m west of the DTM hillshade; and it draws the lunar module
#     1.8 m south and 0.2 m west of the Wagner et al. (2017) coordinate. The three agree to
#     within about 2 m, which is the spread of the two hillshade-based routes, so the entry below
#     is the imagery-against-imagery number: it is the one measured to a fifth of a meter, and it
#     is the offset the level 12 to 13 boundary would show. Measured by
#     scripts/dev/terrain/check_site_registration.py; the same run reads the uncorrected Trek
#     mosaic 23.5 m north of the archive, which is the defect the entries below correct.
LAYER_REGISTRATION = {
    "a11_pho_r": (-0.6, -1.9),
    "a11_26cm": (-24.0, 4.0),
    "a11_nac": (-24.0, 4.0),
    "a11_60": (12.4, -1.9),
}
# Apollo 11 powered descent, as flown by scripts/dev/viewer_demos/apollo11_landing.jl: distance
# still to go along the ground track (km) against height above the landing site (km).
APOLLO11_PROFILE = [
    (480.0, 15.19), (430.6, 15.22), (382.5, 15.62), (337.2, 16.48), (294.5, 17.60), (254.7, 18.77),
    (217.8, 19.81), (184.0, 20.55), (153.3, 20.84), (125.8, 20.58), (101.6, 19.70), (80.8, 18.25),
    (63.6, 16.39), (49.1, 14.17), (36.5, 11.54), (25.6, 8.68), (16.2, 5.76), (8.3, 3.02),
    (3.3, 1.21), (1.18, 0.53), (0.33, 0.20), (0.06, 0.055), (0.0, 0.0),
]


def fetch(url, timeout=120, headers=None, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers or {})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                requested = (headers or {}).get("Range")
                if requested:
                    first, last = map(int, requested.removeprefix("bytes=").split("-"))
                    prefix = f"bytes {first}-{last}/"
                    if r.status != 206 or not r.headers.get("Content-Range", "").startswith(prefix):
                        raise ValueError("Server did not honor the exact byte range; refusing a full raster download.")
                    data = r.read(last - first + 2)
                    if len(data) != last - first + 1:
                        raise ValueError("Incomplete or oversized raster range response.")
                    return data
                return r.read()
        except ValueError:
            raise
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code == 429:
                # The archive mirrors answer a burst of range requests with an hour-long block.
                # Retrying is what earns the next one, so stop and say so.
                raise RuntimeError(f"HTTP 429 from {urllib.parse.urlsplit(url).netloc}"
                                   f" (retry-after {e.headers.get('retry-after', 'unset')} s): rate limited,"
                                   " not retrying. Re-run later; cached bytes are kept.") from e
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(1.5 * (attempt + 1))
    return None


def longitude_near(lon, center):
    """Choose the equivalent longitude closest to the regional interval centre."""
    return center + (lon-center+180.0) % 360.0 - 180.0


def longitude_box_near(west, east, center):
    shift = longitude_near(0.5*(west+east), center)-0.5*(west+east)
    return west+shift, east+shift


def validate_profile(profile):
    if (len(profile) < 2 or any(len(p) != 2 or not all(map(math.isfinite,p)) or min(p) < 0 for p in profile)
            or len({p[0] for p in profile}) != len(profile) or max(p[0] for p in profile) <= 0):
        raise ValueError("Descent profile needs finite nonnegative distinct distances, positive extent and finite heights.")


def write_grid(path_base, heights, lat_min, lat_max, lon_min, lon_max, source):
    heights = np.ascontiguousarray(np.asarray(heights, dtype="<f4"))
    if heights.ndim != 2 or 0 in heights.shape or not np.isfinite(heights).all():
        raise ValueError("DEM heights must be a nonempty finite 2-D array.")
    if not all(map(math.isfinite, (lat_min, lat_max, lon_min, lon_max))):
        raise ValueError("DEM bounds must be finite.")
    if not (-90 <= lat_min < lat_max <= 90 and 0 < lon_max - lon_min < 360):
        raise ValueError("DEM bounds must be increasing regional cell edges.")
    rows, cols = heights.shape
    with open(str(path_base) + ".f32", "wb") as f:
        f.write(heights.tobytes())
    meta = {"rows": rows, "cols": cols, "lat_min": lat_min, "lat_max": lat_max, "lon_min": lon_min, "lon_max": lon_max,
            "units": "m above the 1737.4 km sphere", "layout": "row-major, north to south, west to east, little-endian float32",
            "reference_radius_m": MOON_RADIUS_M, "source": source, "producer_version": 2,
            "sha256": hashlib.sha256(heights.tobytes()).hexdigest()}
    with open(str(path_base) + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    return meta


def borrow(out, reuse, name):
    """Materialize exported DEM files; raw source downloads may remain cache links."""
    if reuse is None:
        return False
    src = pathlib.Path(reuse)/name
    dst = pathlib.Path(out)/name
    exported = pathlib.Path(name).suffix in (".json",".f32")
    if exported and dst.is_symlink() and dst.exists():
        src = dst.resolve()
    elif dst.exists() or not src.exists():
        return dst.exists()
    dst.parent.mkdir(parents=True,exist_ok=True)
    if exported:
        # Atomic replacement also materializes links made by an older producer.
        with tempfile.NamedTemporaryFile(prefix=".dem-copy-",dir=dst.parent,delete=False) as handle:
            temporary = pathlib.Path(handle.name)
        try:
            shutil.copy2(src,temporary)
            os.replace(temporary,dst)
        finally:
            temporary.unlink(missing_ok=True)
    else:
        try:
            os.symlink(os.path.abspath(src),dst)
        except OSError:
            shutil.copy2(src,dst)
    return True


def existing(base, request=None):
    header, binary = base.with_suffix(".json"), base.with_suffix(".f32")
    if not header.exists() or not binary.exists():
        return False
    meta = json.loads(header.read_text())
    if meta.get("producer_version") != 2:
        raise ValueError(f"{header}: regenerate this legacy grid with corrected cell registration.")
    data = binary.read_bytes()
    if hashlib.sha256(data).hexdigest() != meta.get("sha256"):
        raise ValueError(f"{binary}: DEM checksum mismatch.")
    if (any(isinstance(meta.get(k),bool) or not isinstance(meta.get(k),int) or meta[k] <= 0 for k in ("rows","cols"))
            or len(data) != 4*meta["rows"]*meta["cols"] or meta.get("reference_radius_m") != MOON_RADIUS_M
            or not all(math.isfinite(meta[k]) for k in ("lat_min","lat_max","lon_min","lon_max"))
            or not -90 <= meta["lat_min"] < meta["lat_max"] <= 90
            or not 0 < meta["lon_max"]-meta["lon_min"] < 360 or not np.isfinite(np.frombuffer(data,dtype="<f4")).all()):
        raise ValueError(f"{header}: invalid cached DEM shape, datum or bounds.")
    if request is not None and meta.get("request") != request:
        raise ValueError(f"{header}: cached grid does not match this request; choose a fresh output directory.")
    return True


def reuse_grid(base, request, reuse):
    """Validate a complete DEM before copying; failed reuse leaves new outputs empty."""
    if existing(base, request):
        if reuse is not None:
            for suffix in ('.f32', '.json'):
                borrow(base.parent, reuse, base.name + suffix)
        return True
    if reuse is None:
        return False
    source = pathlib.Path(reuse) / base.name
    if not existing(source, request):
        return False
    targets = [base.with_suffix(suffix) for suffix in ('.f32', '.json')]
    if any(p.exists() or p.is_symlink() for p in targets):
        raise ValueError(f"{base}: incomplete output DEM; choose a fresh output directory.")
    with tempfile.TemporaryDirectory(prefix='.dem-reuse-', dir=base.parent) as directory:
        staged = pathlib.Path(directory) / base.name
        for suffix in ('.f32', '.json'):
            shutil.copy2(source.with_suffix(suffix), staged.with_suffix(suffix))
        existing(staged, request)
        published = []
        try:
            for target in targets:  # publish metadata last
                os.replace(staged.with_suffix(target.suffix), target)
                published.append(target)
        except OSError:
            for target in published:
                target.unlink(missing_ok=True)
            raise
    return True


def atomic_write(path, data):
    """Publish complete cache bytes; interruption cannot truncate the final path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix='.cache-', dir=path.parent, delete=False) as handle:
            temporary = pathlib.Path(handle.name)
            handle.write(data)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def lola_window(lat_max, lat_min, lon_min, lon_max, step, lines_per_row, label):
    """Mean columns and evenly spaced representative rows inside each LOLA cell.

    Longitude edges remain unwrapped, including windows across zero longitude.
    At most lines_per_row representative source rows are averaged per output row;
    this is a bounded stratified approximation, not a full two-dimensional box mean.
    """
    if not all(map(math.isfinite, (lat_max, lat_min, lon_min, lon_max))):
        raise ValueError("LOLA bounds must be finite.")
    if not (-90 <= lat_min < lat_max <= 90 and 0 < lon_max - lon_min < 360):
        raise ValueError("LOLA needs regional increasing bounds inside the latitude domain.")
    if isinstance(step, bool) or not isinstance(step, int) or not 1 <= step <= LOLA_LINES:
        raise ValueError("LOLA step must be a positive bounded integer.")
    if isinstance(lines_per_row, bool) or not isinstance(lines_per_row, int) or lines_per_row < 1:
        raise ValueError("LOLA sampled row count must be positive.")
    l0 = int(math.floor((90 - lat_max) * LOLA_PPD))
    l1 = int(math.ceil((90 - lat_min) * LOLA_PPD))
    s0 = int(math.floor(lon_min * LOLA_PPD))
    s1 = int(math.ceil(lon_max * LOLA_PPD))
    rows_out, cols_out = math.ceil((l1-l0)/step), math.ceil((s1-s0)/step)
    if rows_out*step > LOLA_LINES or cols_out*step >= LOLA_SAMPLES:
        raise ValueError("Rounded LOLA window is not regional or exceeds a pole.")
    l0 = min(l0, LOLA_LINES - rows_out*step)
    l1, s1 = l0+rows_out*step, s0+cols_out*step
    cols = s1-s0
    take = min(step, lines_per_row)
    offsets = np.rint(np.linspace(0, step-1, take)).astype(int) if take > 1 else np.array([(step-1)//2])
    grid = np.zeros((rows_out, cols_out), dtype=np.float64)

    def one(idx):
        r, k = divmod(idx, take)
        line = l0+r*step+int(offsets[k])
        remaining, column, pieces = cols, s0 % LOLA_SAMPLES, []
        while remaining:
            count = min(remaining, LOLA_SAMPLES-column)
            first = line*LOLA_RECORD_BYTES+column*2
            data = fetch(LOLA_URL, headers={"Range": f"bytes={first}-{first+count*2-1}"})
            if data is None or len(data) != count*2:
                raise RuntimeError(f"LOLA line {line}: incomplete range")
            pieces.append(np.frombuffer(data, dtype="<i2").astype(np.float64)*0.5)
            remaining -= count
            column = 0
        return r, np.concatenate(pieces).reshape(cols_out, step).mean(axis=1)
    with ThreadPoolExecutor(max_workers=4) as ex:
        for r, row in ex.map(one, range(rows_out*take)):
            grid[r] += row/take
    return grid.astype(np.float32), 90-l1/LOLA_PPD, 90-l0/LOLA_PPD, s0/LOLA_PPD, s1/LOLA_PPD


def fetch_lola(out, lat, lon, half_deg, reuse):
    base = out / "dem_lola"
    request = dict(lat=lat,lon=lon,half_deg=half_deg)
    if reuse_grid(base, request, reuse):
        print("lola: reusing", base.with_suffix(".json"))
        return json.loads(base.with_suffix(".json").read_text())
    lon360 = lon % 360.0
    grid, lat_min, lat_max, lon_min, lon_max = lola_window(
        lat + half_deg, lat - half_deg, lon360 - half_deg, lon360 + half_deg, 1, 1, "")
    meta = write_grid(base, grid, lat_min, lat_max, lon_min, lon_max,
                      f"LOLA LDEM_128 (PDS lola_gdr cylindrical, 128 ppd, 237 m/px), {LOLA_URL}")
    print("lola: wrote", meta["rows"], "x", meta["cols"], "min", float(grid.min()), "max", float(grid.max()))
    meta["request"] = request
    base.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    return meta


def fetch_lola_wide(out, root, samples, reuse):
    """The coarse height window under the whole imagery root, so the horizon has relief."""
    base = out / "dem_lola_wide"
    request = dict(root=root,samples=samples)
    if reuse_grid(base, request, reuse):
        print("lola wide: reusing", base.with_suffix(".json"))
        return json.loads(base.with_suffix(".json").read_text())
    step = max(1, int(round((root["lat_max"] - root["lat_min"]) * LOLA_PPD / max(1, samples))))
    grid, lat_min, lat_max, lon_min, lon_max = lola_window(
        root["lat_max"], root["lat_min"], root["lon_min"], root["lon_max"], step, 4, " wide")
    meta = write_grid(base, grid, lat_min, lat_max, lon_min, lon_max,
                      f"LOLA LDEM_128: {step} source columns averaged, {min(step,4)} evenly spaced rows sampled per {step}-row cell ({step * 237} m/px), {LOLA_URL}")
    print("lola wide: wrote", meta["rows"], "x", meta["cols"], "min", float(grid.min()), "max", float(grid.max()))
    meta["request"] = request
    base.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    return meta


def equirectangular_geometry(scale, tie, keydir, doubles=()):
    if len(keydir) < 4 or len(keydir) < 4+4*keydir[3]:
        raise ValueError("Incomplete GeoTIFF key directory.")
    keys = {}
    for i in range(4, 4+4*keydir[3], 4):
        key, loc, count, value = keydir[i:i+4]
        if loc == 0: keys[key] = value
        elif loc == 34736 and count == 1:
            if not 0 <= value < len(doubles):
                raise ValueError("GeoTIFF key references a missing double value.")
            keys[key] = doubles[value]
    if keys.get(1024) != 1 or keys.get(3075) != 17 or keys.get(1025) != 1:
        raise ValueError("Expected a projected equirectangular PixelIsArea GeoTIFF.")
    if keys.get(3076) != 9001 or keys.get(2054,9102) != 9102:
        raise ValueError("GeoTIFF projected lengths must explicitly use metres and angular coordinates must use degrees.")
    if keys.get(3089,keys.get(3081,0)) != 0:
        raise ValueError("Only zero-latitude equirectangular origins are supported.")
    if 2057 not in keys:
        raise ValueError("GeoTIFF must declare the reference radius.")
    radius = float(keys[2057])
    if radius != MOON_RADIUS_M or (2058 in keys and float(keys[2058]) != radius):
        raise ValueError("GeoTIFF must explicitly use the 1737400 m lunar sphere.")
    if len(scale) < 2 or len(tie) < 6 or not all(map(math.isfinite, list(scale)+list(tie))):
        raise ValueError("Invalid GeoTIFF scale or tiepoint.")
    if scale[0] <= 0 or scale[1] <= 0:
        raise ValueError("GeoTIFF pixel scale must be positive.")
    idx = dict(radius_m=radius, sx=float(scale[0]), sy=float(scale[1]),
        x0=float(tie[3]-tie[0]*scale[0]), y0=float(tie[4]+tie[1]*scale[1]),
        std_parallel_deg=float(keys.get(3078,0)), center_lon_deg=float(keys.get(3088,keys.get(3080,0))),
        false_easting_m=float(keys.get(3082,0)),false_northing_m=float(keys.get(3083,0)))
    if abs(idx['std_parallel_deg']) >= 90 or not all(map(math.isfinite,idx.values())):
        raise ValueError("Invalid equirectangular projection parameters.")
    return idx


def projected_pixel(idx, lon, lat):
    x = idx['false_easting_m']+idx['radius_m']*math.cos(math.radians(idx['std_parallel_deg']))*math.radians(lon-idx['center_lon_deg'])
    y = idx['false_northing_m']+idx['radius_m']*math.radians(lat)
    return (x-idx['x0'])/idx['sx'], (idx['y0']-y)/idx['sy']


def projected_lonlat(idx, column, row):
    lon = idx['center_lon_deg']+math.degrees((idx['x0']+column*idx['sx']-idx['false_easting_m'])/(idx['radius_m']*math.cos(math.radians(idx['std_parallel_deg']))))
    lat = math.degrees((idx['y0']-row*idx['sy']-idx['false_northing_m'])/idx['radius_m'])
    return lon, lat


def fetch_nac(out, site, lat, lon, half_deg, step_m, reuse):
    spec = NAC_DTM.get(site)
    if spec is None:
        print("nac: no NAC DTM registered for site", site)
        return None
    base = out / "dem_nac"
    request = dict(site=site,lat=lat,lon=lon,half_deg=half_deg,step_m=step_m)
    if reuse_grid(base, request, reuse):
        print("nac: reusing", base.with_suffix(".json"))
        return json.loads(base.with_suffix(".json").read_text())
    raw = out / "raw" / pathlib.Path(spec["url"]).name
    borrow(out, reuse, os.path.join("raw", raw.name))
    raw.parent.mkdir(parents=True, exist_ok=True)
    if not raw.exists():
        print("nac: downloading", spec["url"])
        data = fetch(spec["url"], timeout=900)
        if data is None:
            print("nac: not available"); return None
        raw.write_bytes(data)
        print("nac: saved", raw, len(data), "bytes")
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(raw) as im:
        if 34264 in im.tag_v2:
            raise ValueError("Rotated GeoTIFF transforms are not supported.")
        geometry = equirectangular_geometry(im.tag_v2[33550], im.tag_v2[33922], im.tag_v2[34735], im.tag_v2.get(34736, ()))
        arr = np.array(im,dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[:, :, 0]
    lines, samples = arr.shape
    if (lines, samples) != (spec["lines"], spec["samples"]):
        print("nac: raster shape", arr.shape, "differs from the label", (spec["lines"], spec["samples"]), "- using the raster")
    if not math.isfinite(step_m) or step_m <= 0 or any(abs(geometry[k]-2.0) > 1e-8 for k in ("sx","sy")):
        raise ValueError("This NAC DTM reader requires an explicit 2 m grid and a positive output stride.")
    lon = longitude_near(lon, geometry["center_lon_deg"])
    left, top = projected_pixel(geometry, lon-half_deg, lat+half_deg)
    right, bottom = projected_pixel(geometry, lon+half_deg, lat-half_deg)
    l0, l1 = max(0,math.floor(top)), min(lines,math.ceil(bottom))
    s0, s1 = max(0,math.floor(left)), min(samples,math.ceil(right))
    def lat_of_line(row): return projected_lonlat(geometry,0,row)[1]
    def lon_of_sample(column): return projected_lonlat(geometry,column,0)[0]
    if l1 <= l0 or s1 <= s0:
        raise ValueError("Requested NAC window does not overlap this product.")
    sub = arr[l0:l1, s0:s1]
    nodata = sub <= -1e30
    if nodata.any():
        sub = sub.copy(); sub[nodata] = np.nan
    stride = max(1, int(round(step_m / 2.0)))
    sub = sub[::stride, ::stride]
    # LROC NAC DTM samples are elevation in metres above the 1737400 m sphere.
    # LROC's processing guide adds that radius when making an absolute-radius DEM.
    # Never infer units from the numerical magnitude or fabricate missing elevations.
    if not np.isfinite(sub).all():
        raise ValueError("NAC window contains no-data values; choose a fully covered window.")
    # GeoTIFF tiepoints locate outer PixelIsArea edges. Preserve the centres
    # selected by striding by placing output edges half a stride from them.
    north = lat_of_line(l0 + 0.5 - stride/2)
    south = lat_of_line(l0 + 0.5 + stride*(sub.shape[0]-0.5))
    west = lon_of_sample(s0 + 0.5 - stride/2)
    east = lon_of_sample(s0 + 0.5 + stride*(sub.shape[1]-0.5))
    meta = write_grid(base, sub, south, north, west, east,
                      f"LROC NAC DTM elevation (1737400 m sphere, 2 m/px, samples every {stride} px), {spec['url']}")
    print("nac: wrote", meta["rows"], "x", meta["cols"], "min", float(np.nanmin(sub)), "max", float(np.nanmax(sub)))
    meta["request"] = request
    base.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    return meta


# ---- the descent corridor ----------------------------------------------------------------

def great_circle_point(lat, lon, azimuth_deg, distance_m):
    p1, l1 = math.radians(lat), math.radians(lon)
    az, d = math.radians(azimuth_deg), distance_m / MOON_RADIUS_M
    p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(az))
    l2 = l1 + math.atan2(math.sin(az) * math.sin(d) * math.cos(p1), math.cos(d) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), math.degrees(l2)


def descent_track(lat, lon, approach_azimuth_deg, profile, samples=4000):
    """Ground-track samples (lat, lon, height above the site) from the start of the descent in.

    The profile is resampled evenly and finely: the coverage of a level is the union of disks
    around the track samples, so a spacing wider than the finest level's reach would leave gaps
    between them along the corridor."""
    uprange_az = (approach_azimuth_deg + 180.0) % 360.0
    validate_profile(profile)
    ordered = sorted(profile, key=lambda p: p[0])
    d = [p[0] for p in ordered]; h = [p[1] for p in ordered]
    far = d[-1]
    track = []
    for k in range(max(2, samples) + 1):
        to_go = far * k / max(2, samples)
        p = great_circle_point(lat, lon, uprange_az, to_go * 1.0e3)
        track.append((p[0], p[1], float(np.interp(to_go, d, h)) * 1.0e3))
    return track


def quadtree_root(track, root_zoom):
    """The square root region: one Trek tile wide at `root_zoom`, its corner on that pixel grid,
    centered on the corridor so the whole track and the horizon around it are inside."""
    side = 180.0 / 2 ** root_zoom
    pixel = side / 256.0
    lat_c = 0.5 * (min(p[0] for p in track) + max(p[0] for p in track))
    longitudes = np.rad2deg(np.unwrap(np.deg2rad([p[1] for p in track])))
    lon_c = 0.5 * (min(longitudes) + max(longitudes))
    lon_min = round((lon_c - side / 2) / pixel) * pixel
    lat_max = round((lat_c + side / 2) / pixel) * pixel
    lat_max = min(90.0, max(-90.0 + side, lat_max))
    return {"lat_min": lat_max - side, "lat_max": lat_max, "lon_min": lon_min, "lon_max": lon_min + side}


def level_coverage(track, root, level, lod_factor):
    """Tile indices (x, y) of `level` within the LOD reach of the corridor.

    A node of side s is drawn out to about 5.3 s from the camera before its parent takes over, so
    `lod_factor` s is the radius within which this level is worth having; a camera at height h can
    only reach sqrt((F s)^2 - h^2) of ground within it."""
    side = root["lat_max"] - root["lat_min"]
    n = 2 ** level
    reach = lod_factor * side * M_PER_DEG / n
    tiles = set()
    for lat, lon, height in track:
        lon = longitude_near(lon, 0.5*(root["lon_min"]+root["lon_max"]))
        if height >= reach:
            continue
        radius = math.sqrt(reach * reach - height * height)
        dlat = radius / M_PER_DEG
        dlon = dlat / max(1e-6, math.cos(math.radians(lat)))
        x0 = int(math.floor((lon - dlon - root["lon_min"]) / side * n))
        x1 = int(math.floor((lon + dlon - root["lon_min"]) / side * n))
        y0 = int(math.floor((root["lat_max"] - (lat + dlat)) / side * n))
        y1 = int(math.floor((root["lat_max"] - (lat - dlat)) / side * n))
        for x in range(max(0, x0), min(n - 1, x1) + 1):
            for y in range(max(0, y0), min(n - 1, y1) + 1):
                tiles.add((x, y))
    return tiles


# ---- Moon Trek tiles ---------------------------------------------------------------------

class TrekTiles:
    """Moon Trek WMTS tiles, cached on disk so re-tiling costs no new downloads."""

    def __init__(self, cache_dirs, workers=6):
        self.cache_dirs = [pathlib.Path(d) for d in cache_dirs]
        self.workers = workers
        self.hits = self.misses = self.absent = 0

    def path(self, key, z, x, y, root=0):
        spec = TREK_LAYERS[key]
        x = x % (2 ** (z+1))
        return self.cache_dirs[root] / key / str(z) / str(y) / f"{x}.{spec['ext']}"

    def covers(self, key, z, lon_min, lat_min, lon_max, lat_max):
        spec = TREK_LAYERS[key]
        if z > spec["zmax"]:
            return False
        w, s, e, nn = spec["bbox"]
        lon_min, lon_max = longitude_box_near(lon_min, lon_max, 0.5*(w+e))
        return not (lon_max <= w or lon_min >= e or lat_max <= s or lat_min >= nn)

    # Empty legacy entries may be interrupted writes, so use an explicit 404 marker.
    NOT_FOUND = b'SpaceAGORA Trek cache: HTTP 404\n'

    @staticmethod
    def valid_tile(data):
        try:
            with Image.open(io.BytesIO(data)) as image:
                if image.size != (256, 256):
                    return False
                image.load()
            return True
        except (OSError, ValueError):
            return False

    def _cached(self, key, z, x, y):
        for root in range(len(self.cache_dirs)):
            p = self.path(key, z, x, y, root)
            if p.is_file():
                data = p.read_bytes()
                if data == self.NOT_FOUND:
                    return True, None
                if self.valid_tile(data):
                    return True, data
        return False, None

    def get(self, key, z, x, y):
        """The raw bytes of one Trek tile, or None when the service does not have it."""
        if not 0 <= y < 2**z:
            return None
        x %= 2**(z+1)
        found, data = self._cached(key, z, x, y)
        if found:
            self.hits += 1
            return data
        spec = TREK_LAYERS[key]
        data = fetch(TREK.format(layer=spec["layer"], z=z, y=y, x=x, ext=spec["ext"]), timeout=90)
        p = self.path(key, z, x, y)
        if data is not None and not self.valid_tile(data):
            raise RuntimeError(f"{key}/{z}/{y}/{x}: invalid Trek tile response; not cached.")
        atomic_write(p, self.NOT_FOUND if data is None else data)
        if data is None:
            self.absent += 1
        else:
            self.misses += 1
        return data

    def prefetch(self, jobs):
        jobs = [j for j in jobs if not self._cached(*j)[0]]
        if not jobs:
            return
        print(f"  trek: fetching {len(jobs)} tiles", flush=True)
        with ThreadPoolExecutor(max_workers=self.workers) as ex:
            for k, _ in enumerate(ex.map(lambda j: self.get(*j), jobs)):
                if k and k % 200 == 0:
                    print(f"    {k}/{len(jobs)}", flush=True)

    def image(self, key, z, x, y):
        data = self.get(key, z, x, y)
        if not data:
            return None
        with Image.open(io.BytesIO(data)) as image:
            return image.convert("LA")


class ArchiveRaster:
    """One internally tiled, uncompressed GeoTIFF on the PDS, read in place by HTTP range.

    The file is 915 MB, so nothing is downloaded whole: the TIFF header points at the image
    file directory, the directory gives the tile grid and the byte range of every tile, and a
    window of the raster is the handful of tile rows that cover it. Offsets are read from the
    file, never assumed; a layout this reader does not understand (compressed, striped, more
    than one band, a projection other than the equirectangular one the LROC RDRs use) is an
    error rather than a silent wrong answer.

    Rate limits are the design constraint. A burst of small range requests to the archive's
    front end earns HTTP 429 with an hour-long retry-after, so a window is fetched as one range
    request per tile row (about 0.8 MB each) and every 64 kB tile is written to the cache as it
    arrives: a re-run, or a run resumed after an interruption, costs no requests at all.

    Georeferencing comes from the file. ModelPixelScaleTag and ModelTiepointTag put raster
    (0, 0) - the outer corner of the upper-left pixel, since RasterPixelIsArea - at projected
    (x0, y0) meters, and the GeoTIFF keys give the projection: equirectangular on a sphere of
    radius R about ProjCenterLong with ProjStdParallel1, so

        x = FE + R cos(phi_std) (lon - lon_0),    y = FN + R lat

    with lon, lat in radians. The product's PDS4 label bounds are carried in ARCHIVE_SOURCES
    only as a cross-check: the constructor fails if the two disagree by more than two pixels.
    """

    def __init__(self, key, cache_dirs, verbose=True):
        self.key = key
        self.spec = ARCHIVE_SOURCES[key]
        if isinstance(cache_dirs, (str, pathlib.Path)):
            cache_dirs = [cache_dirs]
        # the first cache is written to; the rest (a --reuse site's, say) are read-only
        self.dirs = [pathlib.Path(d) / self.spec["product"] for d in cache_dirs]
        self.dir = self.dirs[0]
        self.dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        self.requests = self.hits = self.fetched = 0
        self._index = None

    # -- the file's own header -------------------------------------------------------------
    def _range(self, start, end):
        """One HTTP range request, [start, end] inclusive."""
        self.requests += 1
        data = fetch(self.spec["url"], timeout=300, headers={"Range": f"bytes={start}-{end}"})
        if data is None or len(data) != end - start + 1:
            raise RuntimeError(f"{self.key}: short read of bytes {start}-{end} "
                               f"({0 if data is None else len(data)} of {end - start + 1})")
        return data

    @staticmethod
    def _entry(blob, base, tag_map, tag):
        typ, count, value, inline = tag_map[tag]
        fmt = _TIFF_TYPE_FMT.get(typ)
        if fmt is None:
            raise RuntimeError(f"tiff tag {tag}: unsupported type {typ}")
        size = _TIFF_TYPE_BYTES[typ] * count
        if size <= 4:
            return list(struct.unpack_from("<" + fmt * count, inline, 0))
        off = value - base
        if off < 0 or off + size > len(blob):
            raise RuntimeError(f"tiff tag {tag}: value at {value} is outside the header block read")
        return list(struct.unpack_from("<" + fmt * count, blob, off))

    def index(self):
        """The raster's geometry and tile table, parsed once and cached on disk."""
        if self._index is not None:
            return self._index
        path = self.dir / "index.json"
        for d in self.dirs:
            if (d / "index.json").exists():
                self._index = json.loads((d / "index.json").read_text())
                self._index["tile_offsets"] = np.asarray(self._index["tile_offsets"], dtype=np.int64)
                try:
                    self._check_label(self._index)
                except Exception:
                    self._index = None
                    raise
                return self._index
        head = self._range(0, 15)
        if head[:4] != b"II\x2a\x00":
            raise RuntimeError(f"{self.key}: not a little-endian classic TIFF ({head[:4]!r})")
        ifd_off = struct.unpack_from("<I", head, 4)[0]
        blob = self._range(ifd_off, ifd_off + 262143)     # the directory and its out-of-line arrays
        n = struct.unpack_from("<H", blob, 0)[0]
        tag_map = {}
        for i in range(n):
            tag, typ, count, value = struct.unpack_from("<HHII", blob, 2 + 12 * i)
            tag_map[tag] = (typ, count, value, blob[2 + 12 * i + 8:2 + 12 * i + 12])
        def one(tag, default=None):
            if tag not in tag_map:
                if default is None:
                    raise RuntimeError(f"{self.key}: tiff tag {tag} is missing")
                return default
            return self._entry(blob, ifd_off, tag_map, tag)[0]
        width, height = one(256), one(257)
        for tag, want, what in ((258, 8, "bits per sample"), (259, 1, "compression"),
                                (262, 1, "photometric interpretation"), (277, 1, "samples per pixel"),
                                (284, 1, "planar configuration"), (339, 1, "sample format")):
            got = one(tag, want)
            if got != want:
                raise RuntimeError(f"{self.key}: {what} is {got}, not {want}; this reader handles "
                                   "single-band 8-bit uncompressed tiled rasters only")
        tw, th = one(322), one(323)
        offsets = np.asarray(self._entry(blob, ifd_off, tag_map, 324), dtype=np.int64)
        counts = np.asarray(self._entry(blob, ifd_off, tag_map, 325), dtype=np.int64)
        across = -(-width // tw); down = -(-height // th)
        if len(offsets) != across * down or len(counts) != len(offsets):
            raise RuntimeError(f"{self.key}: {len(offsets)} tile offsets for a {across} x {down} tile grid")
        if not np.all(counts == tw * th):
            raise RuntimeError(f"{self.key}: tiles are not {tw * th} bytes; the raster is not uncompressed")
        if 34264 in tag_map:
            raise ValueError("Rotated GeoTIFF transforms are not supported.")
        scale = self._entry(blob, ifd_off, tag_map, 33550)
        tie = self._entry(blob, ifd_off, tag_map, 33922)
        if tie[0] != 0.0 or tie[1] != 0.0:
            raise RuntimeError(f"{self.key}: the tiepoint is at raster ({tie[0]}, {tie[1]}), not (0, 0)")
        keydir = self._entry(blob, ifd_off, tag_map, 34735)
        doubles = self._entry(blob, ifd_off, tag_map, 34736) if 34736 in tag_map else []
        geometry = equirectangular_geometry(scale, tie, keydir, doubles)
        idx = {"width": int(width), "height": int(height), "tile_w": int(tw), "tile_h": int(th),
               "tiles_across": int(across), "tiles_down": int(down), "tile_bytes": int(tw * th),
               **geometry, "index_version": 2, "source_url": self.spec["url"],
               "ifd_offset": int(ifd_off), "tile_offsets": offsets.tolist()}
        self._index = idx
        try:
            self._check_label(idx)
        except Exception:
            self._index = None
            raise
        atomic_write(path, json.dumps(idx).encode())
        idx["tile_offsets"] = offsets
        if self.verbose:
            print(f"  {self.key}: {width} x {height} px, {across} x {down} tiles of {tw} x {th}, "
                  f"{self.pixel_size_m()[0]:.4f} x {self.pixel_size_m()[1]:.4f} m/px")
        return idx

    # -- geometry --------------------------------------------------------------------------
    def _x_per_rad(self, idx):
        return idx["radius_m"] * math.cos(math.radians(idx["std_parallel_deg"]))

    def pixel_size_m(self):
        idx = self.index()
        return idx["sx"], idx["sy"]

    def pixel_of(self, lon_deg, lat_deg):
        """Fractional (column, row) of a longitude and latitude, the raster's own convention."""
        idx = self.index()
        center = idx["center_lon_deg"] + math.degrees((idx["x0"]+0.5*idx["width"]*idx["sx"]-idx["false_easting_m"])/self._x_per_rad(idx))
        lon_deg = longitude_near(lon_deg, center)
        x = idx["false_easting_m"] + self._x_per_rad(idx) * math.radians(lon_deg - idx["center_lon_deg"])
        y = idx["false_northing_m"] + idx["radius_m"] * math.radians(lat_deg)
        return (x - idx["x0"]) / idx["sx"], (idx["y0"] - y) / idx["sy"]

    def lonlat_of(self, col, row):
        idx = self.index()
        lon = idx["center_lon_deg"] + math.degrees((idx["x0"] + col * idx["sx"] - idx["false_easting_m"]) / self._x_per_rad(idx))
        lat = math.degrees((idx["y0"] - row * idx["sy"] - idx["false_northing_m"]) / idx["radius_m"])
        return lon, lat

    def bounds(self):
        idx = self.index()
        w, n = self.lonlat_of(0, 0)
        e, s = self.lonlat_of(idx["width"], idx["height"])
        return w, s, e, n

    def _check_label(self, idx):
        """The label's bounds and pixel size against the file's own geotransform, to a pixel."""
        if idx.get("index_version") != 2 or idx.get("source_url") != self.spec["url"]:
            raise ValueError(f"{self.key}: regenerate the legacy or mismatched archive index cache.")
        positive = ("width","height","tile_w","tile_h","tiles_across","tiles_down","tile_bytes")
        if any(isinstance(idx.get(k),bool) or not isinstance(idx.get(k),int) or idx[k] <= 0 for k in positive):
            raise ValueError(f"{self.key}: invalid cached raster dimensions.")
        numeric = ("x0","y0","sx","sy","radius_m","std_parallel_deg","center_lon_deg","false_easting_m","false_northing_m")
        if (not all(math.isfinite(idx[k]) for k in numeric) or idx["sx"] <= 0 or idx["sy"] <= 0
                or idx["radius_m"] != MOON_RADIUS_M or abs(idx["std_parallel_deg"]) >= 90
                or len(idx["tile_offsets"]) != idx["tiles_across"]*idx["tiles_down"]
                or any(int(v) != v or v < 0 for v in idx["tile_offsets"])):
            raise ValueError(f"{self.key}: invalid cached raster datum or tile offsets.")
        nominal = self.spec.get("label_m_per_px")
        if nominal and max(abs(idx["sx"] - nominal), abs(idx["sy"] - nominal)) > 0.01 * nominal:
            raise RuntimeError(f"{self.key}: the file is on a {idx['sx']:.4f} x {idx['sy']:.4f} m grid, "
                               f"the PDS4 label's is {nominal} m; the product has changed")
        want = self.spec.get("bbox")
        if not want:
            return
        got = self.bounds()
        for got_deg, want_deg, name in zip(got, want, ("west", "south", "east", "north")):
            err = abs(got_deg - want_deg) * M_PER_DEG
            if err > 2.0 * max(idx["sx"], idx["sy"]):
                raise RuntimeError(f"{self.key}: the file's {name} edge is {got_deg:.9f}, the PDS4 label's "
                                   f"{want_deg:.9f} ({err:.2f} m apart); the geotransform convention has changed")
            if self.verbose:
                print(f"    {name} edge: file {got_deg:.9f}, label {want_deg:.9f} ({err:.2f} m)")

    def covers(self, lon_min, lat_min, lon_max, lat_max):
        w, s, e, n = self.bounds()
        lon_min, lon_max = longitude_box_near(lon_min, lon_max, 0.5*(w+e))
        return not (lon_max <= w or lon_min >= e or lat_max <= s or lat_min >= n)

    # -- tiles -----------------------------------------------------------------------------
    def _tile_path(self, row, col, root=0):
        return self.dirs[root] / str(row) / f"{col}.bin"

    def _cached(self, row, col):
        size = self.index()["tile_bytes"]
        for r in range(len(self.dirs)):
            p = self._tile_path(row, col, r)
            if p.is_file() and p.stat().st_size == size:
                return p
        return None

    def _registered(self, lon_min, lat_min, lon_max, lat_max):
        """The box to read so the ground it draws lands where it belongs (LAYER_REGISTRATION)."""
        north_m, east_m = LAYER_REGISTRATION.get(self.key, (0.0, 0.0))
        if not north_m and not east_m:
            return lon_min, lat_min, lon_max, lat_max
        dlat = north_m / M_PER_DEG
        dlon = east_m / (M_PER_DEG * max(1e-6, math.cos(math.radians(0.5 * (lat_min + lat_max)))))
        return lon_min + dlon, lat_min + dlat, lon_max + dlon, lat_max + dlat

    def _tile_span(self, lon_min, lat_min, lon_max, lat_max, margin_px=2):
        """The tile row and column range covering a latitude/longitude box, clipped to the raster."""
        idx = self.index()
        lon_min, lat_min, lon_max, lat_max = self._registered(lon_min, lat_min, lon_max, lat_max)
        c0, r0 = self.pixel_of(lon_min, lat_max)
        c1, r1 = self.pixel_of(lon_max, lat_min)
        c0 = max(0, int(math.floor(c0)) - margin_px); c1 = min(idx["width"] - 1, int(math.ceil(c1)) + margin_px)
        r0 = max(0, int(math.floor(r0)) - margin_px); r1 = min(idx["height"] - 1, int(math.ceil(r1)) + margin_px)
        if c1 < c0 or r1 < r0:
            return None
        return (r0 // idx["tile_h"], r1 // idx["tile_h"], c0 // idx["tile_w"], c1 // idx["tile_w"])

    def prefetch(self, boxes):
        """Cache every tile covering these latitude/longitude boxes, a tile row per request."""
        idx = self.index()
        spans = [s for s in (self._tile_span(*b) for b in boxes) if s is not None]
        if not spans:
            return
        r0 = min(s[0] for s in spans); r1 = max(s[1] for s in spans)
        c0 = min(s[2] for s in spans); c1 = max(s[3] for s in spans)
        missing = [(r, c) for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)
                   if self._cached(r, c) is None]
        self.hits += (r1 - r0 + 1) * (c1 - c0 + 1) - len(missing)
        if not missing:
            if self.verbose:
                print(f"  {self.key}: {(r1 - r0 + 1) * (c1 - c0 + 1)} tiles already cached")
            return
        rows = sorted({r for r, _ in missing})
        mb = len(missing) * idx["tile_bytes"] / 1e6
        print(f"  {self.key}: fetching {len(missing)} tiles ({mb:.1f} MB) in {len(rows)} range requests", flush=True)
        for r in rows:
            cols = sorted(c for rr, c in missing if rr == r)
            # one request per contiguous run of columns: the tiles of a row are consecutive in
            # the file, so a run is a single byte range
            run = []
            for c in cols + [None]:
                if run and (c is None or c != run[-1] + 1 or not self._contiguous(r, run[-1], c)):
                    self._fetch_run(r, run)
                    run = []
                if c is not None:
                    run.append(c)
            print(f"    row {r} of {rows[0]}..{rows[-1]}", flush=True)

    def _contiguous(self, row, col_a, col_b):
        idx = self.index()
        off = idx["tile_offsets"]
        a = row * idx["tiles_across"] + col_a
        return int(off[a + 1] - off[a]) == idx["tile_bytes"] and col_b == col_a + 1

    def _fetch_run(self, row, cols):
        idx = self.index()
        off = idx["tile_offsets"]; nb = idx["tile_bytes"]
        first = off[row * idx["tiles_across"] + cols[0]]
        last = off[row * idx["tiles_across"] + cols[-1]] + nb - 1
        data = self._range(int(first), int(last))
        for c in cols:
            start = int(off[row * idx["tiles_across"] + c] - first)
            p = self._tile_path(row, c)
            atomic_write(p, data[start:start + nb])
            self.fetched += 1

    def _tile(self, row, col):
        idx = self.index()
        p = self._cached(row, col)
        if p is None:
            self._fetch_run(row, [col])
            p = self._tile_path(row, col)
        else:
            self.hits += 1
        return np.frombuffer(p.read_bytes(), dtype=np.uint8).reshape(idx["tile_h"], idx["tile_w"])

    # -- windows ---------------------------------------------------------------------------
    def window(self, lon_min, lat_max, lon_max, lat_min, px):
        """A `px` by `px` area-averaged window of the raster over a latitude/longitude box.

        Returns an "LA" image whose alpha is 0 where the raster has no data (outside its own
        bounds, or the zero fill outside the frame's footprint inside them), or None when the
        box misses the raster entirely. The source grid is not a power-of-two relative of the
        quadtree's, so this is a true resample - a box average, the honest filter for a
        decimation - unlike the integer-pixel crops the Trek path makes.
        """
        idx = self.index()
        lon_min, lat_min, lon_max, lat_max = self._registered(lon_min, lat_min, lon_max, lat_max)
        fx0, fy0 = self.pixel_of(lon_min, lat_max)
        fx1, fy1 = self.pixel_of(lon_max, lat_min)
        if fx1 <= 0 or fy1 <= 0 or fx0 >= idx["width"] or fy0 >= idx["height"]:
            return None
        ix0 = int(math.floor(fx0)); iy0 = int(math.floor(fy0))
        ix1 = int(math.ceil(fx1)); iy1 = int(math.ceil(fy1))
        buf = np.zeros((iy1 - iy0, ix1 - ix0), dtype=np.uint8)
        tw, th = idx["tile_w"], idx["tile_h"]
        for r in range(max(0, iy0) // th, (min(idx["height"], iy1) - 1) // th + 1):
            for c in range(max(0, ix0) // tw, (min(idx["width"], ix1) - 1) // tw + 1):
                tile = self._tile(r, c)
                sy0 = max(iy0, r * th); sy1 = min(iy1, (r + 1) * th, idx["height"])
                sx0 = max(ix0, c * tw); sx1 = min(ix1, (c + 1) * tw, idx["width"])
                if sy1 <= sy0 or sx1 <= sx0:
                    continue
                buf[sy0 - iy0:sy1 - iy0, sx0 - ix0:sx1 - ix0] = tile[sy0 - r * th:sy1 - r * th, sx0 - c * tw:sx1 - c * tw]
        box = (fx0 - ix0, fy0 - iy0, fx1 - ix0, fy1 - iy0)
        gray = Image.fromarray(buf).resize((px, px), Image.BOX, box=box)
        cover = Image.fromarray(np.where(buf > 0, 255, 0).astype(np.uint8)).resize((px, px), Image.BOX, box=box)
        alpha = Image.fromarray((np.asarray(cover) >= 128).astype(np.uint8) * 255)
        out = Image.merge("LA", (gray, alpha))
        return out


def source_covers(trek, archives, key, z, lon_min, lat_min, lon_max, lat_max):
    """Whether the layer or archive raster `key` has anything over this box at this zoom."""
    if key in ARCHIVE_SOURCES:
        src = archives.get(key)
        return src is not None and src.covers(lon_min, lat_min, lon_max, lat_max)
    return trek.covers(key, z, lon_min, lat_min, lon_max, lat_max)


def node_pixel_origin(root, root_zoom, level, x, y):
    """The node's top-left corner in Trek's global pixel grid at zoom root_zoom + level."""
    z = root_zoom + level
    per_deg = 2 ** (z + 1) * 256 / 360.0
    px = (root["lon_min"] + 180.0) * per_deg + x * 256
    py = (90.0 - root["lat_max"]) * per_deg + y * 256
    return z, int(round(px)), int(round(py))


def layer_origin(key, z, px, py):
    """Where to read `key` so that the ground it draws lands where it belongs: its own offset,
    turned into pixels of this zoom (see LAYER_REGISTRATION)."""
    north_m, east_m = LAYER_REGISTRATION.get(key, (0.0, 0.0))
    if not north_m and not east_m:
        return px, py
    m_per_px = M_PER_DEG * 360.0 / (2 ** (z + 1) * 256)
    lat = 90.0-(py+128)*360.0/(2**(z+1)*256)
    east_m_per_px = m_per_px*max(1e-6,math.cos(math.radians(lat)))
    return px + int(round(east_m / east_m_per_px)), py - int(round(north_m / m_per_px))


def node_image(trek, archives, root, root_zoom, level, x, y, tile_px):
    """One node, composited coarse source first so the finer ones paste over. Trek layers are an
    integer pixel crop of Trek's own tiles; an archive raster is area-averaged onto the node's
    grid, since its grid is not a power-of-two relative of Trek's. Returns (image, primary
    source key) or (None, None)."""
    z, px, py = node_pixel_origin(root, root_zoom, level, x, y)
    side = (root["lat_max"] - root["lat_min"]) / 2 ** level
    lon_min = root["lon_min"] + x * side
    lat_max = root["lat_max"] - y * side
    keys = [k for k in ZOOM_LAYERS.get(z, ZOOM_LAYERS_DEFAULT)
            if source_covers(trek, archives, k, z, lon_min, lat_max - side, lon_min + side, lat_max)]
    canvas, primary = None, None
    for key in reversed(keys):                      # coarse first
        got = False
        if key in ARCHIVE_SOURCES:
            layer_img = archives[key].window(lon_min, lat_max, lon_min + side, lat_max - side, tile_px)
            if layer_img is None:
                continue
            got = True
        else:
            layer_img = Image.new("LA", (tile_px, tile_px), (0, 0))
            lpx, lpy = layer_origin(key, z, px, py)
            tx0, ty0 = lpx // 256, lpy // 256
            tx1, ty1 = (lpx + tile_px - 1) // 256, (lpy + tile_px - 1) // 256
            for ty in range(ty0, ty1 + 1):
                for tx in range(tx0, tx1 + 1):
                    im = trek.image(key, z, tx, ty)
                    if im is None:
                        continue
                    layer_img.paste(im, (tx * 256 - lpx, ty * 256 - lpy))
                    got = True
        if not got:
            continue
        if canvas is None:
            canvas = Image.new("LA", (tile_px, tile_px), (0, 255))
        canvas.paste(layer_img, (0, 0), layer_img.split()[1])
        alpha = np.asarray(layer_img.split()[1], dtype=np.float32)
        if alpha.mean() > 160:                       # this layer supplies most of the node
            primary = key
    if canvas is None:
        return None, None
    return canvas.convert("L"), primary


def match_globe_texture(img, texture_path, root):
    """Level 0's gain and offset against the page's own global texture over the same box, so the
    root of the quadtree is not a dark square in the middle of the globe it sits in."""
    if not texture_path or not os.path.exists(texture_path):
        print("  (no global texture to match level 0 against)")
        return img, (1.0, 0.0)
    Image.MAX_IMAGE_PIXELS = None
    tex = Image.open(texture_path).convert("L")
    w, h = tex.size
    box = (int(round((root["lon_min"] + 180) / 360 * w)), int(round((90 - root["lat_max"]) / 180 * h)),
           int(round((root["lon_max"] + 180) / 360 * w)), int(round((90 - root["lat_min"]) / 180 * h)))
    # Longitude wraps at the texture seam; Pillow crop would pad with black.
    x0,y0,x1,y1 = box
    if y0 < 0 or y1 > h or y1 <= y0 or x1 <= x0:
        raise ValueError("Texture-match root is outside the latitude domain or subpixel.")
    ref = np.asarray(tex,dtype=np.float32)[y0:y1][:,np.arange(x0,x1)%w]
    cur = np.asarray(img, dtype=np.float32)
    if ref.std() < 1e-3 or cur.std() < 1e-3:
        return img, (1.0, 0.0)
    gain = float(ref.std() / cur.std())
    offset = float(ref.mean() - gain * cur.mean())
    out = np.clip(cur * gain + offset, 0, 255).astype(np.uint8)
    return Image.fromarray(out), (gain, offset)


def detail_over_parent(images, level, parents, tile_px, detail_std, max_gain):
    """Rebuild each node as its ancestor's tone plus its own octave of detail.

    The low frequencies a node shares with the level above it are taken from the (already built)
    parent and only the node's own detail band is added, with an amplitude normalized per source
    mosaic. Nothing of a mosaic's exposure survives, so no level and no mosaic edge draws a box,
    a node and the ancestor a missing neighbor inherits agree by construction, and the surface
    keeps the detail each level actually resolves instead of being flattened toward the coarsest
    one - which is what matching a level's contrast to its parent does over a dozen levels."""
    half = tile_px // 2
    bands, parented = {}, {}
    for (x, y), (img, key) in images.items():
        parent = parents.get((level - 1, x // 2, y // 2))
        if parent is None:
            continue
        arr = np.asarray(img, dtype=np.float32)
        low = np.asarray(img.resize((half, half), Image.BILINEAR).resize((tile_px, tile_px), Image.BICUBIC), dtype=np.float32)
        crop = parent.crop(((x % 2) * half, (y % 2) * half, (x % 2) * half + half, (y % 2) * half + half))
        base = np.asarray(crop.resize((tile_px, tile_px), Image.BICUBIC), dtype=np.float32)
        parented[(x, y)] = (base, arr - low, key)
        bands.setdefault(key, []).append(float((arr - low).std()))
    gains = {}
    for key, stds in bands.items():
        m = float(np.mean(stds))
        gains[key] = min(max_gain, detail_std / m) if m > 1e-3 else 1.0
    out = {}
    for (x, y), (img, key) in images.items():
        if (x, y) not in parented:
            out[(x, y)] = img
            continue
        base, band, key = parented[(x, y)]
        out[(x, y)] = Image.fromarray(np.clip(base + band * gains.get(key, 1.0), 0, 255).astype(np.uint8))
    return out, gains


def imagery_request(track, root, a):
    options = ("root_zoom","max_level","tile_px","lod_factor","fine_lod_factor","detail_std", "detail_max_gain",
               "quality_ramp_level","quality_coarse","quality_fine","approach_azimuth","uprange_km")
    texture = pathlib.Path(a.match_texture) if a.match_texture else None
    request = {"root":root,"options":{k:getattr(a,k) for k in options},
        "track_sha256":hashlib.sha256(json.dumps(track,allow_nan=False).encode()).hexdigest(),
        "texture_sha256":hashlib.sha256(texture.read_bytes()).hexdigest() if texture and texture.is_file() else None,
        "trek_layers":TREK_LAYERS,"zoom_layers":ZOOM_LAYERS,"archive_sources":ARCHIVE_SOURCES,
        "registration":LAYER_REGISTRATION}
    return json.loads(json.dumps(request,sort_keys=True,allow_nan=False))


def validate_imagery_cache(idir, request):
    meta = json.loads((idir/"tiles.json").read_text())
    if meta.get("producer_version") != 2 or meta.get("request") != request:
        raise ValueError("Cached imagery does not match the build request; use --retile or a fresh output directory.")
    if not meta.get("nodes") or meta.get("root") != request["root"] or meta.get("tile_px") != request["options"]["tile_px"]:
        raise ValueError("Cached imagery has invalid nodes or geometry; rebuild with --retile.")
    for node in meta["nodes"]:
        path = (idir/node["file"]).resolve()
        if idir.resolve() not in path.parents or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != node.get("sha256"):
            raise ValueError("Cached imagery tile is missing, outside its directory, or has a checksum mismatch; rebuild with --retile.")
    return meta


def build_tiles(out, track, root, a):
    idir = out/"imagery"
    request = imagery_request(track,root,a)
    if (idir/"tiles.json").exists() and not a.retile:
        return validate_imagery_cache(idir,request)
    if idir.exists() and any(idir.iterdir()) and not a.retile:
        raise ValueError("Incomplete imagery output exists; use --retile or a fresh output directory.")
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".imagery-build-",dir=out) as temporary:
        staging = pathlib.Path(temporary)/"imagery"
        meta = _build_tiles(out,track,root,a,staging,request)
        backup = pathlib.Path(temporary)/"previous"
        had_previous = idir.exists()
        if had_previous:
            os.replace(idir,backup)
        try:
            os.replace(staging,idir)
        except Exception:
            if had_previous: os.replace(backup,idir)
            raise
    return meta


def _build_tiles(out, track, root, a, idir, request):
    meta_path = idir/"tiles.json"
    tdir = idir/"tiles"
    tdir.mkdir(parents=True)
    caches = [out / "trek_cache"] + ([pathlib.Path(a.reuse) / "trek_cache"] if a.reuse else [])
    caches[0].mkdir(parents=True, exist_ok=True)
    trek = TrekTiles(caches, workers=a.workers)
    # the archive rasters any level asks for, opened once; --reuse borrows their tile cache too
    archive_keys = sorted({k for z, keys in ZOOM_LAYERS.items() if z <= a.root_zoom + a.max_level
                           for k in keys if k in ARCHIVE_SOURCES})
    pds_caches = [out / "pds_cache"] + ([pathlib.Path(a.reuse) / "pds_cache"] if a.reuse else [])
    archives = {}
    for key in archive_keys:
        print(f"archive source {key}: {ARCHIVE_SOURCES[key]['product']}")
        archives[key] = ArchiveRaster(key, pds_caches)
        archives[key].index()
    height_max = max(p[2] for p in track)
    nodes, levels, parents = [], [], {}
    total_bytes = 0
    anchor = (1.0, 0.0)
    for level in range(a.max_level + 1):
        side_m = (root["lat_max"] - root["lat_min"]) * M_PER_DEG / 2 ** level
        # the corridor levels (whose reach exceeds the highest point of the descent, so that they
        # span the whole track) are the expensive ones and get the tighter factor
        factor = a.lod_factor if a.lod_factor * side_m >= height_max else a.fine_lod_factor
        want = sorted(level_coverage(track, root, level, factor))
        z = a.root_zoom + level
        keys = ZOOM_LAYERS.get(z, ZOOM_LAYERS_DEFAULT)
        jobs = []
        boxes = {k: [] for k in keys if k in ARCHIVE_SOURCES}
        for (x, y) in want:
            zz, px, py = node_pixel_origin(root, a.root_zoom, level, x, y)
            s = (root["lat_max"] - root["lat_min"]) / 2 ** level
            lon_min = root["lon_min"] + x * s; lat_max = root["lat_max"] - y * s
            for key in keys:
                if not source_covers(trek, archives, key, zz, lon_min, lat_max - s, lon_min + s, lat_max):
                    continue
                if key in ARCHIVE_SOURCES:
                    boxes[key].append((lon_min, lat_max - s, lon_min + s, lat_max))
                    continue
                lpx, lpy = layer_origin(key, zz, px, py)
                for ty in range(lpy // 256, (lpy + a.tile_px - 1) // 256 + 1):
                    for tx in range(lpx // 256, (lpx + a.tile_px - 1) // 256 + 1):
                        jobs.append((key, zz, tx, ty))
        source_name = "trek z=" + str(z) if not any(boxes.values()) else "archive"
        print(f"level {level} ({source_name}, {side_m / 256:.3f} m/px, reach {factor * side_m / 1e3:.1f} km): {len(want)} nodes")
        trek.prefetch(sorted(set(jobs)))
        for key, bx in boxes.items():
            if bx:
                archives[key].prefetch(bx)
        images = {}
        for (x, y) in want:
            img, key = node_image(trek, archives, root, a.root_zoom, level, x, y, a.tile_px)
            if img is not None:
                images[(x, y)] = (img, key)
        if level == 0:
            gains = {}
            for xy, (img, _key) in list(images.items()):
                images[xy], anchor = match_globe_texture(img, a.match_texture, root)
        else:
            images, gains = detail_over_parent(images, level, parents, a.tile_px, a.detail_std, a.detail_max_gain)
        ramp = max(1, a.quality_ramp_level)
        quality = int(round(a.quality_coarse + (a.quality_fine - a.quality_coarse) * min(level, ramp) / ramp))
        level_bytes = 0
        ldir = tdir / f"L{level}"
        ldir.mkdir(exist_ok=True)
        for (x, y), img in sorted(images.items()):
            name = f"L{level}/{x}_{y}.jpg"
            img.save(tdir / name, quality=quality, optimize=True)
            n = (tdir / name).stat().st_size
            level_bytes += n
            nodes.append({"level": level, "x": x, "y": y, "file": f"tiles/{name}", "m_per_px": side_m / a.tile_px,
                          "sha256":hashlib.sha256((tdir/name).read_bytes()).hexdigest()})
            parents[(level, x, y)] = img
        for key in [k for k in parents if k[0] < level - 1]:
            del parents[key]                          # only the level above is matched against
        total_bytes += level_bytes
        levels.append({"level": level, "trek_zoom": z, "m_per_px": side_m / a.tile_px, "node_side_m": side_m,
                       "lod_factor": factor, "reach_m": factor * side_m, "quality": quality,
                       "nodes": len(images), "bytes": level_bytes, "layers": keys,
                       "detail_gain": {k: round(g, 3) for k, g in gains.items()}})
        print(f"  -> {len(images)} tiles, q={quality}, {level_bytes / 1e6:.3f} MB"
              f"{'' if not gains else ' detail ' + ', '.join(f'{k} x{g:.2f}' for k, g in gains.items())}")
        if not images and level > 0:
            print("  (no imagery at this level; stopping)"); break
    used = sorted({k for lv in levels for k in lv["layers"] if k in ARCHIVE_SOURCES
                   and any(n["level"] == lv["level"] for n in nodes)})
    finest = min((n["m_per_px"] for n in nodes), default=0.0)
    archive_meta = [{"key": k, "product": ARCHIVE_SOURCES[k]["product"], "url": ARCHIVE_SOURCES[k]["url"],
                     "credit": ARCHIVE_SOURCES[k]["credit"], "citation": ARCHIVE_SOURCES[k]["citation"],
                     "grid_m_per_px": round(archives[k].pixel_size_m()[0], 4),
                     "levels": [lv["level"] for lv in levels if k in lv["layers"]
                                and any(n["level"] == lv["level"] for n in nodes)]} for k in used]
    source = "NASA Moon Trek WMTS (LROC WAC global mosaic; Kaguya TC ortho mosaic; LROC NAC Apollo 11 mosaics)"
    for m in archive_meta:
        source += f"; PDS {m['citation']}"
    if not any(n["level"] == 0 for n in nodes):
        raise ValueError("No root imagery was produced; check source coverage before using this site.")
    meta = {"producer_version":2,"request":request,"scheme": "quadtree", "root": root, "tile_px": a.tile_px,
            "layer_registration_m": {k: {"north": v[0], "east": v[1]} for k, v in LAYER_REGISTRATION.items()},
            "max_level": max((n["level"] for n in nodes), default=0), "root_zoom": a.root_zoom,
            "corridor": {"approach_azimuth_deg": a.approach_azimuth, "uprange_m": a.uprange_km * 1.0e3,
                         "lod_factor": a.lod_factor, "fine_lod_factor": a.fine_lod_factor},
            "levels": levels, "nodes": nodes, "archive_sources": archive_meta,
            "globe_texture_match": {"path": os.path.basename(a.match_texture or ""), "gain": round(anchor[0], 4), "offset": round(anchor[1], 2)},
            # What the finest level is worth, which is not the grid it is drawn on. The deepest
            # source is one LROC NAC frame whose ground sample is 0.24 m across the detector line
            # and 0.56 m along track (SCALED_PIXEL_WIDTH and _HEIGHT of M175124932R in the PDS
            # volume index), and radially averaged power spectra of the archive raster lose their
            # signal at about 0.5 m of feature size: see docs/reference/lunar_imagery_sources.md,
            # section 6. So the tiles are sampled finer than the imagery resolves, and a viewer
            # that quotes the sampling alone overstates what is there by about a factor of two.
            "resolution": {"finest_m_per_px": round(finest, 4),
                           "source_grid_m_per_px": round(archives[used[0]].pixel_size_m()[0], 4) if used else None,
                           "detector_sample_m": [0.24, 0.56] if used else None,
                           "feature_scale_m": 0.5 if used else None,
                           "note": ("tiles are sampled at %.2f m/px; the original Apollo 11 study estimated "
                                    "about 0.5 m feature scale for this source. That estimate is inherited, "
                                    "not remeasured by this producer" % finest) if used else
                                   ("tiles are sampled at %.2f m/px (Moon Trek's deepest tile matrix)" % finest)},
            "attribution": sorted({ARCHIVE_SOURCES[k]["credit"] for k in used} | {"NASA Moon Trek"}),
            "source": source}
    meta_path.write_text(json.dumps(meta,indent=1))
    print(f"tiles: {len(nodes)} nodes, max level {meta['max_level']}, {total_bytes / 1e6:.3f} MB "
          f"(trek cache: {trek.hits} hits, {trek.misses} fetched, {trek.absent} absent)")
    print(f"  base64 in the payload: about {total_bytes * 4 / 3 / 1e6:.3f} MB")
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", nargs=2, type=float, metavar=("LAT", "LON"), required=True)
    ap.add_argument("--name", default="apollo11", help="site key for the NAC DTM registry")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reuse", default=None, help="another site directory to borrow DEM files, raw downloads and the Moon Trek tile cache from")
    ap.add_argument("--lola-half-deg", type=float, default=2.5)
    ap.add_argument("--lola-wide-samples", type=int, default=288, help="rows of the coarse LOLA grid under the whole root")
    ap.add_argument("--nac-half-deg", type=float, default=0.06)
    ap.add_argument("--nac-step", type=float, default=4.0, help="NAC DTM resample step, meters")
    ap.add_argument("--approach-azimuth", type=float, default=270.0, help="heading flown at the site, degrees east of north")
    ap.add_argument("--uprange-km", type=float, default=480.0, help="where the descent starts, along the track")
    ap.add_argument("--track-samples", type=int, default=4000, help="how finely the ground track is resampled for the coverage test")
    ap.add_argument("--profile", default=None, help="JSON file of [distance_to_go_km, height_above_site_km] pairs; the Apollo 11 descent by default")
    ap.add_argument("--root-zoom", type=int, default=3, help="Moon Trek zoom whose tile size is the quadtree root (3 = 22.5 deg)")
    ap.add_argument("--max-level", type=int, default=13,
                    help="deepest quadtree level; level L is Trek zoom root_zoom + L, and levels past Trek's "
                         "own cap (zoom 15, level 12) come from the PDS archive rasters in ARCHIVE_SOURCES")
    ap.add_argument("--tile-px", type=int, default=256)
    ap.add_argument("--lod-factor", type=float, default=2.0, help="coverage radius of a level, in node sides (a full pyramid is 5.3)")
    ap.add_argument("--fine-lod-factor", type=float, default=2.5, help="the same for the levels that no longer span the corridor")
    ap.add_argument("--detail-std", type=float, default=3.2, help="amplitude each level's own detail band is normalized to, 0..255")
    ap.add_argument("--detail-max-gain", type=float, default=2.2, help="cap on amplifying a level whose source has little detail left")
    ap.add_argument("--match-texture", default=str(pathlib.Path(__file__).resolve().parents[3] / "data" / "textures" / "moon_lroc_wac_4k.jpg"),
                    help="global texture the root level is matched to, so the covered region is not a box on the globe")
    ap.add_argument("--quality-ramp-level", type=int, default=12,
                    help="level at which the JPEG quality reaches --quality-fine; deeper levels keep it")
    ap.add_argument("--quality-coarse", type=int, default=44)
    ap.add_argument("--quality-fine", type=int, default=66)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--retile", action="store_true", help="rebuild the imagery quadtree even if tiles.json is there")
    ap.add_argument("--no-imagery", action="store_true")
    ap.add_argument("--no-nac", action="store_true", help="skip the optional full NAC DTM download")
    a = ap.parse_args(argv)
    if (not all(map(math.isfinite, a.site)) or not -90 < a.site[0] < 90
            or not all(math.isfinite(v) and v > 0 for v in (a.lola_half_deg,a.nac_half_deg,a.nac_step,a.lod_factor,a.fine_lod_factor))
            or not math.isfinite(a.approach_azimuth)
            or not math.isfinite(a.detail_std) or a.detail_std < 0
            or not math.isfinite(a.detail_max_gain) or a.detail_max_gain < 0
            or a.quality_ramp_level < 1 or a.root_zoom+a.max_level > 16
            or a.site[0]-a.lola_half_deg < -90 or a.site[0]+a.lola_half_deg > 90
            or a.lola_half_deg >= 180 or a.nac_half_deg >= 180
            or not 1 <= a.workers <= 16 or not 1 <= a.lola_wide_samples <= 4096
            or not 0 <= a.root_zoom <= 16 or not 0 <= a.max_level <= 16
            or a.tile_px != 256 or not 2 <= a.track_samples <= 100000
            or not math.isfinite(a.uprange_km) or a.uprange_km <= 0
            or not 1 <= a.quality_coarse <= 95 or not 1 <= a.quality_fine <= 95):
        ap.error("Invalid site, bounds, sampling, tile, worker or quality option.")
    out = pathlib.Path(a.out)
    lat, lon = a.site
    site = {"lat_deg": lat, "lon_deg": lon, "name": a.name}
    profile = APOLLO11_PROFILE
    if a.profile:
        profile = [(float(p[0]), float(p[1])) for p in json.loads(pathlib.Path(a.profile).read_text())]
    try:
        validate_profile(profile)
    except ValueError as error:
        ap.error(str(error))
    scale = a.uprange_km / max(p[0] for p in profile)
    profile = [(d * scale, h) for d, h in profile]
    track = descent_track(lat, lon, a.approach_azimuth, profile, a.track_samples)
    root = quadtree_root(track, a.root_zoom)
    center = 0.5*(root["lon_min"]+root["lon_max"])
    if any(not (root["lat_min"] <= lat <= root["lat_max"] and root["lon_min"] <= longitude_near(lon,center) <= root["lon_max"]) for lat,lon,_ in track):
        ap.error("The root does not cover the requested track; reduce --root-zoom or --uprange-km.")
    out.mkdir(parents=True,exist_ok=True)
    print(f"site {a.name} at {lat} N {lon} E, approach {a.approach_azimuth} deg, {a.uprange_km} km uprange")
    print(f"  descent starts at {track[-1][0]:.5f} N {track[-1][1]:.5f} E, {track[-1][2] / 1e3:.2f} km up")
    print(f"  root {root['lat_min']:.6f}..{root['lat_max']:.6f} N, {root['lon_min']:.6f}..{root['lon_max']:.6f} E"
          f" ({(root['lat_max'] - root['lat_min']) * M_PER_DEG / 1e3:.1f} km square)")
    lola_wide = fetch_lola_wide(out, root, a.lola_wide_samples, a.reuse)
    lola = fetch_lola(out, lat, lon, a.lola_half_deg, a.reuse)
    nac = None if a.no_nac else fetch_nac(out, a.name, lat, lon, a.nac_half_deg, a.nac_step, a.reuse)
    tiles = None if a.no_imagery else build_tiles(out, track, root, a)
    dems = [{"name": n, **m} for n, m in (("dem_nac", nac), ("dem_lola", lola), ("dem_lola_wide", lola_wide)) if m]  # finest first
    (out/"site.json").write_text(json.dumps({"site":site,"dem":dems,"tiles":"imagery/tiles.json" if tiles else None},indent=1))
    print("done:", out / "site.json")


if __name__ == "__main__":
    main()
