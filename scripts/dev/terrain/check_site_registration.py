#!/usr/bin/env python3
"""Measure where a site's imagery actually draws the ground, against the site's own DTM.

    python3 scripts/dev/terrain/check_site_registration.py --site 0.67416 23.47314 \
        --dem data/terrain/moon/apollo11/dem_nac.json --cache output/terrain/moon/apollo11/pds_cache

This is the measurement behind `LAYER_REGISTRATION` in `fetch_moon_site.py`, for the archive
rasters that script reads by HTTP range (`ARCHIVE_SOURCES`). Imagery and terrain are fetched
from different products, and a level of the quadtree takes its detail from whichever source
reaches its zoom, so an uncorrected offset between two sources puts the same crater in two
places either side of a level change. Two independent routes:

1. **Against the digital terrain model.** The imagery is resampled onto the DTM's own grid and
   phase-correlated against a hillshade of it. The illumination that makes a hillshade look
   like the image is the illumination the image was taken under, so the correlation is scanned
   over sun azimuth and the peak reported with the shift; `--sun-azimuth` fixes it instead.
   The DTM carries the LROC control the published site coordinate comes from, so a shift here
   is the imagery's error, not the terrain's.
2. **Against the level above.** The same ground rendered from the Moon Trek layer that feeds the
   deepest Trek level of the quadtree, at that level's own pixel grid and through the same
   `LAYER_REGISTRATION` correction the tile builder applies, correlated against the archive
   window. Both are imagery of the same surface under different illumination, so this
   correlation is far better conditioned than the hillshade one, and it is the quantity the
   level boundary in the viewer actually shows.

3. **Against the published coordinate of a surface object.** The brightest pixel inside a small
   box around the site coordinate is the lunar module (a bright metal object on a dark mare
   surface); its offset from the coordinate is a second, independent reading of the same error.
   Wagner et al. (2017), Icarus 283, 92-103 place the Apollo 11 descent stage at
   0.67416 N 23.47314 E to 0.3 m, which is the `A11_LM` row of the LROC ANTHROPOGENIC_OBJECTS
   shapefile and the coordinate the landing simulation targets.

Both are reported in meters north and east, in the sense `LAYER_REGISTRATION` uses: the offset
of where the source draws the ground from where the ground is.
"""
import argparse, importlib.util, json, math, pathlib, sys
import numpy as np
from PIL import Image

sys.dont_write_bytecode = True          # no __pycache__ beside the script it borrows
_here = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("fetch_moon_site", _here / "fetch_moon_site.py")
fms = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(fms)


def load_grid(meta_path):
    meta = json.load(open(meta_path))
    data = np.fromfile(str(meta_path).replace(".json", ".f32"), dtype="<f4").reshape(meta["rows"], meta["cols"])
    return meta, data


def hillshade(heights, m_per_cell, sun_azimuth_deg, sun_elevation_deg):
    """Lambertian shading of a height grid, rows north to south."""
    dz_dy, dz_dx = np.gradient(heights.astype(np.float64), m_per_cell)
    # rows run north to south, so a positive dz_dy is a slope down to the north
    nx, ny, nz = -dz_dx, dz_dy, np.ones_like(heights, dtype=np.float64)
    norm = np.sqrt(nx * nx + ny * ny + nz * nz)
    az, el = math.radians(sun_azimuth_deg), math.radians(sun_elevation_deg)
    sx, sy, sz = math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)
    return np.clip((nx * sx + ny * sy + nz * sz) / norm, 0.0, None)


def phase_correlate(a, b):
    """Shift (rows, cols) that moves `b` onto `a`, by phase correlation with a Hann window."""
    a = np.asarray(a, dtype=np.float64); b = np.asarray(b, dtype=np.float64)
    a = (a - a.mean()) / (a.std() + 1e-12); b = (b - b.mean()) / (b.std() + 1e-12)
    win = np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    A = np.fft.rfft2(a * win); B = np.fft.rfft2(b * win)
    cross = A * np.conj(B)
    cross /= np.abs(cross) + 1e-12
    corr = np.fft.irfft2(cross, s=a.shape)
    peak = np.unravel_index(int(np.argmax(corr)), corr.shape)
    def refine(i, n, axis):
        # parabolic fit through the peak and its neighbors, on the wrapped axis
        lo = corr[(peak[0] - 1) % corr.shape[0], peak[1]] if axis == 0 else corr[peak[0], (peak[1] - 1) % corr.shape[1]]
        hi = corr[(peak[0] + 1) % corr.shape[0], peak[1]] if axis == 0 else corr[peak[0], (peak[1] + 1) % corr.shape[1]]
        mid = corr[peak]
        den = lo - 2 * mid + hi
        d = 0.5 * (lo - hi) / den if abs(den) > 1e-12 else 0.0
        s = i + max(-1.0, min(1.0, d))
        return s - n if s > n / 2 else s
    return refine(peak[0], corr.shape[0], 0), refine(peak[1], corr.shape[1], 1), float(corr[peak])


def against_dem(source, dem_meta, dem, lat, lon, box_m, azimuths, sun_elevation_deg):
    """Phase correlation of the imagery against the DTM's hillshade, over a box at the site."""
    dlat = (dem_meta["lat_max"] - dem_meta["lat_min"]) / dem_meta["rows"]
    dlon = (dem_meta["lon_max"] - dem_meta["lon_min"]) / dem_meta["cols"]
    m_per_cell = dlat * fms.M_PER_DEG
    n = int(box_m / m_per_cell) // 2 * 2
    r0 = int(round((dem_meta["lat_max"] - lat) / dlat - n / 2))
    c0 = int(round((lon - dem_meta["lon_min"]) / dlon - n / 2))
    if r0 < 0 or c0 < 0 or r0 + n > dem_meta["rows"] or c0 + n > dem_meta["cols"]:
        raise SystemExit("the correlation box does not fit inside the DTM grid; use a smaller --box")
    heights = dem[r0:r0 + n, c0:c0 + n]
    lat_max = dem_meta["lat_max"] - r0 * dlat
    lat_min = lat_max - n * dlat
    lon_min = dem_meta["lon_min"] + c0 * dlon
    lon_max = lon_min + n * dlon
    img = source.window(lon_min, lat_max, lon_max, lat_min, n)
    if img is None:
        raise SystemExit("the source does not cover the correlation box")
    gray = np.asarray(img.split()[0], dtype=np.float64)
    if np.asarray(img.split()[1]).min() == 0:
        print("  (the source has no data over part of the box)")
    print(f"  box {n} x {n} cells of {m_per_cell:.2f} m ({n * m_per_cell:.0f} m), DTM relief "
          f"{heights.max() - heights.min():.1f} m")
    best = None
    for az in azimuths:
        shade = hillshade(heights, m_per_cell, az, sun_elevation_deg)
        dr, dc, peak = phase_correlate(gray, shade)
        north, east = -dr * m_per_cell, dc * m_per_cell
        print(f"  sun azimuth {az:5.1f} deg: peak {peak:.3f}, imagery {north:+7.2f} m north, {east:+7.2f} m east")
        if best is None or peak > best[0]:
            best = (peak, az, north, east)
    return best


def against_object(source, lat, lon, box_m, dark_percentile, bright_dn, search_m):
    """Where the archive draws the lunar module, against its published coordinate.

    The descent stage is a bright metal object a few meters across casting a long shadow on a
    mare surface of otherwise gentle contrast, so it is found in two steps: the darkest pixels
    of the box are its shadow, and the saturated pixels within `search_m` of that shadow are the
    sunlit top of the stage itself. The top, not the shadow, is what sits on the coordinate; the
    shadow only says where the Sun was. The photometric browse raster is clipped at DN 255 over
    bright rocks elsewhere in the box, which is why the search is anchored on the shadow.
    """
    half = box_m / 2 / fms.M_PER_DEG
    halflon = half / math.cos(math.radians(lat))
    px = int(round(box_m / source.pixel_size_m()[0]))
    img = source.window(lon - halflon, lat + half, lon + halflon, lat - half, px)
    if img is None:
        raise SystemExit("the source does not cover the object box")
    arr = np.asarray(img.split()[0], dtype=np.float64)
    m_per_px = box_m / px
    dark_t = np.percentile(arr, dark_percentile)
    dark = arr <= dark_t
    ys, xs = np.nonzero(dark)
    w = (dark_t - arr[dark]) + 1e-6
    sy, sx = float((ys * w).sum() / w.sum()), float((xs * w).sum() / w.sum())
    yy, xx = np.mgrid[0:px, 0:px]
    near = (yy - sy) ** 2 + (xx - sx) ** 2 <= (search_m / m_per_px) ** 2
    bright = (arr >= bright_dn) & near
    if bright.sum() < 4:
        raise SystemExit("no sunlit object found next to the shadow; widen --object-box")
    ys, xs = np.nonzero(bright)
    w = arr[bright] - bright_dn + 1.0
    by, bx = float((ys * w).sum() / w.sum()), float((xs * w).sum() / w.sum())
    center = px / 2 - 0.5
    print(f"  box {box_m:.0f} m at {m_per_px:.3f} m/px: shadow {dark.sum()} px at "
          f"{(center - sy) * m_per_px:+.2f} N {(sx - center) * m_per_px:+.2f} E, "
          f"sunlit top {bright.sum()} px (DN >= {bright_dn:.0f}, median {np.median(arr):.0f})")
    return (center - by) * m_per_px, (bx - center) * m_per_px


def trek_window(cache_dirs, key, zoom, lon_min, lat_max, n, corrected=True):
    """The same ground as the tile builder reads it: an integer-pixel crop of Trek's own tiles."""
    trek = fms.TrekTiles(cache_dirs)
    per_deg = 2 ** (zoom + 1) * 256 / 360.0
    px = int(round((lon_min + 180.0) * per_deg)); py = int(round((90.0 - lat_max) * per_deg))
    if corrected:
        px, py = fms.layer_origin(key, zoom, px, py)
    out = Image.new("L", (n, n), 0)
    got = False
    for ty in range(py // 256, (py + n - 1) // 256 + 1):
        for tx in range(px // 256, (px + n - 1) // 256 + 1):
            im = trek.image(key, zoom, tx, ty)
            if im is None:
                continue
            out.paste(im.convert("L"), (tx * 256 - px, ty * 256 - py))
            got = True
    if not got:
        return None
    return np.asarray(out, dtype=np.float64)


def against_trek(source, cache_dirs, key, zoom, lat, lon, sizes):
    """The archive against the Moon Trek layer that feeds the level above it."""
    per_deg = 2 ** (zoom + 1) * 256 / 360.0
    m_per_px = fms.M_PER_DEG * 360.0 / (2 ** (zoom + 1) * 256)
    out = []
    for n in sizes:
        lon_min = lon - (n / 2) / per_deg; lat_max = lat + (n / 2) / per_deg
        lon_max = lon_min + n / per_deg; lat_min = lat_max - n / per_deg
        arch = source.window(lon_min, lat_max, lon_max, lat_min, n)
        for corrected in (True, False):
            trek = trek_window(cache_dirs, key, zoom, lon_min, lat_max, n, corrected)
            if arch is None or trek is None:
                print(f"  {n} px: no imagery"); continue
            dr, dc, peak = phase_correlate(np.asarray(arch.split()[0], dtype=np.float64), trek)
            north, east = -dr * m_per_px, dc * m_per_px
            tag = f"{key} corrected" if corrected else f"{key} as served"
            print(f"  {n} px ({n * m_per_px:.0f} m) against {tag}: peak {peak:.3f}, "
                  f"archive {north:+7.2f} m north, {east:+7.2f} m east of it")
            if corrected:
                out.append((north, east, peak))
    if not out:
        return None
    return (float(np.mean([o[0] for o in out])), float(np.mean([o[1] for o in out])),
            float(np.mean([o[2] for o in out])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", nargs=2, type=float, metavar=("LAT", "LON"), required=True)
    ap.add_argument("--dem", required=True, help="the site's finest DEM grid header, e.g. dem_nac.json")
    ap.add_argument("--source", default="a11_pho_r", help="key in fetch_moon_site.ARCHIVE_SOURCES")
    ap.add_argument("--cache", required=True, help="the archive tile cache directory (pds_cache)")
    ap.add_argument("--box", type=float, default=640.0, help="side of the hillshade correlation box, meters")
    ap.add_argument("--object-box", type=float, default=40.0, help="side of the box searched for the lander, meters")
    ap.add_argument("--object-dark-percentile", type=float, default=1.0, help="percentile taken as the lander's shadow")
    ap.add_argument("--object-bright-dn", type=float, default=250.0, help="DN taken as the lander's sunlit top")
    ap.add_argument("--object-search-m", type=float, default=5.0, help="how far from the shadow the sunlit top is looked for")
    ap.add_argument("--sun-azimuth", type=float, default=None, help="fix the hillshade azimuth instead of scanning")
    ap.add_argument("--sun-elevation", type=float, default=49.0,
                    help="solar elevation of the frame (90 - incidence; 40.98 deg incidence for M175124932R)")
    ap.add_argument("--trek-cache", action="append", default=None, help="a Moon Trek tile cache to compare against (repeatable)")
    ap.add_argument("--trek-layer", default="a11_26cm", help="the Trek layer feeding the level above the archive one")
    ap.add_argument("--trek-zoom", type=int, default=15, help="the Trek zoom of that level")
    ap.add_argument("--trek-sizes", type=int, nargs="*", default=[256, 384, 512], help="correlation sizes, Trek pixels")
    a = ap.parse_args()
    lat, lon = a.site
    source = fms.ArchiveRaster(a.source, a.cache)
    source.index()
    dem_meta, dem = load_grid(a.dem)
    reg = fms.LAYER_REGISTRATION.get(a.source, (0.0, 0.0))
    print(f"{a.source}: {fms.ARCHIVE_SOURCES[a.source]['product']}")
    print(f"  LAYER_REGISTRATION in force: {reg[0]:+.2f} m north, {reg[1]:+.2f} m east "
          "(a measurement below near zero means the correction is right)")
    print(f"DTM: {dem_meta['rows']} x {dem_meta['cols']}, {dem_meta['source'].split(',')[0]}")
    azimuths = [a.sun_azimuth] if a.sun_azimuth is not None else list(np.arange(0.0, 360.0, 15.0))
    print("1. against the DTM hillshade:")
    peak, az, north, east = against_dem(source, dem_meta, dem, lat, lon, a.box, azimuths, a.sun_elevation)
    print(f"   -> best at sun azimuth {az:.1f} deg (peak {peak:.3f}): {north:+.2f} m north, {east:+.2f} m east")
    trek = None
    if a.trek_cache:
        print("2. against the level above:")
        trek = against_trek(source, a.trek_cache, a.trek_layer, a.trek_zoom, lat, lon, a.trek_sizes)
        if trek:
            print(f"   -> mean over {len(a.trek_sizes)} boxes: {trek[0]:+.2f} m north, {trek[1]:+.2f} m east "
                  f"(peak {trek[2]:.3f})")
    print("3. against the published coordinate:")
    onorth, oeast = against_object(source, lat, lon, a.object_box, a.object_dark_percentile,
                                   a.object_bright_dn, a.object_search_m)
    print(f"   -> the lander is {onorth:+.2f} m north, {oeast:+.2f} m east of {lat} N {lon} E")
    print(f"residual with the correction in force: hillshade {north:+.2f} N {east:+.2f} E"
          + (f", level above {trek[0]:+.2f} N {trek[1]:+.2f} E" if trek else "")
          + f", lander {onorth:+.2f} N {oeast:+.2f} E")
    print(f"  (requests this run: {source.requests}, tiles already cached: {source.hits})")


if __name__ == "__main__":
    main()
