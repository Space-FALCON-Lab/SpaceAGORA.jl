#!/usr/bin/env python3
"""Fetch a lunar landing site's terrain and imagery for SpaceAGORA.

    python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --out data/terrain/moon/apollo11

Writes, under --out:
  dem_lola.json / dem_lola.f32   LOLA LDEM 128 ppd (237 m/px) heights, meters above the
                                 1737.4 km reference sphere, a +-2.5 deg window around the site,
                                 read straight out of the PDS raster with HTTP range requests
  dem_nac.json / dem_nac.f32     LROC NAC DTM heights (2 m/px, resampled to --nac-step m) over the
                                 inner window, when a NAC DTM covers the site (Apollo 11 does)
  imagery/level_k.jpg + imagery.json
                                 nested square image patches centered on the site, coarse to fine,
                                 from NASA Moon Trek WMTS (LROC WAC global mosaic down to the
                                 Apollo 11 NAC mosaics), each with its lat/lon box and meters per pixel
  imagery/level_k_albedo.jpg     the same patch with its baked-in illumination divided out (see
                                 `derive_albedo`), for a page that lights the surface itself

Every grid file is little-endian Float32, row-major, rows from north to south, columns from
west to east; its JSON carries rows, cols, lat_min, lat_max, lon_min, lon_max and the source.
Needs numpy and Pillow, and network access to pds-geosciences.wustl.edu, pds.mcp.nasa.gov and
trek.nasa.gov. Files already present are reused (delete them to refetch); `--reuse DIR` copies the
DEMs and the imagery originals of an existing site directory instead, so a second copy of a site can
be derived (new albedo images, say) without touching the network.
"""
import argparse, io, json, math, os, pathlib, shutil, struct, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from PIL import Image, ImageFilter

MOON_RADIUS_M = 1737400.0
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
# Imagery levels: half-width of the square patch (deg) and the Trek layers to try (finest first).
# Each layer entry is (layer, zoom, extension); a tile outside a layer's coverage comes back
# transparent (or 404) and the next layer fills it, so the finer local mosaics sit on the global one.
LEVELS = [
    (2.0,       [("LRO_WAC_Mosaic_Global_303ppd_v02", 8, "jpg")]),
    (0.5,       [("A11_60x60km.eq", 10, "png"), ("LRO_WAC_Mosaic_Global_303ppd_v02", 10, "jpg")]),
    (0.125,     [("A11_60x60km.eq", 12, "png"), ("LRO_WAC_Mosaic_Global_303ppd_v02", 10, "jpg")]),
    (0.03125,   [("LRO_NAC_Apollo11_Mosaic_p", 14, "png"), ("A11_60x60km.eq", 14, "png")]),
    (0.0078125, [("apollo11_26cm_mosaic_byte_geo_1_2_highContrast", 15, "png"), ("LRO_NAC_Apollo11_Mosaic_p", 15, "png"), ("A11_60x60km.eq", 15, "png")]),
]


def fetch(url, timeout=120, headers=None, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers or {})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(1.5 * (attempt + 1))
    return None


def write_grid(path_base, heights, lat_min, lat_max, lon_min, lon_max, source):
    heights = np.ascontiguousarray(heights.astype("<f4"))
    rows, cols = heights.shape
    with open(str(path_base) + ".f32", "wb") as f:
        f.write(heights.tobytes())
    meta = {"rows": rows, "cols": cols, "lat_min": lat_min, "lat_max": lat_max, "lon_min": lon_min, "lon_max": lon_max,
            "units": "m above the 1737.4 km sphere", "layout": "row-major, north to south, west to east, little-endian float32",
            "reference_radius_m": MOON_RADIUS_M, "source": source}
    with open(str(path_base) + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    return meta


def fetch_lola(out, lat, lon, half_deg):
    base = out / "dem_lola"
    if (base.with_suffix(".json")).exists() and (base.with_suffix(".f32")).exists():
        print("lola: reusing", base.with_suffix(".json"))
        return json.load(open(base.with_suffix(".json")))
    lon360 = lon % 360.0
    lat_max, lat_min = lat + half_deg, lat - half_deg
    lon_min, lon_max = lon360 - half_deg, lon360 + half_deg
    l0 = int(math.floor((90.0 - lat_max) * LOLA_PPD)); l1 = int(math.ceil((90.0 - lat_min) * LOLA_PPD))
    s0 = int(math.floor(lon_min * LOLA_PPD)); s1 = int(math.ceil(lon_max * LOLA_PPD))
    l0, l1 = max(0, l0), min(LOLA_LINES, l1); s0, s1 = max(0, s0), min(LOLA_SAMPLES, s1)
    rows, cols = l1 - l0, s1 - s0
    print(f"lola: lines {l0}..{l1} samples {s0}..{s1} ({rows} x {cols}) by range requests")
    grid = np.zeros((rows, cols), dtype=np.float32)

    def one(li):
        start = (l0 + li) * LOLA_RECORD_BYTES + s0 * 2
        data = fetch(LOLA_URL, headers={"Range": f"bytes={start}-{start + cols * 2 - 1}"})
        if data is None or len(data) != cols * 2:
            raise RuntimeError(f"lola line {l0 + li}: short read")
        return li, np.frombuffer(data, dtype="<i2").astype(np.float32) * 0.5
    with ThreadPoolExecutor(max_workers=8) as ex:
        for k, (li, row) in enumerate(ex.map(one, range(rows))):
            grid[li, :] = row
            if k % 100 == 0:
                print(f"  lola row {k}/{rows}", flush=True)
    # cell centers: line l covers 90 - l/128 .. 90 - (l+1)/128; store the edges of the window
    meta = write_grid(base, grid, 90.0 - l1 / LOLA_PPD, 90.0 - l0 / LOLA_PPD, s0 / LOLA_PPD, s1 / LOLA_PPD,
                      f"LOLA LDEM_128 (PDS lola_gdr cylindrical, 128 ppd, 237 m/px), {LOLA_URL}")
    print("lola: wrote", meta["rows"], "x", meta["cols"], "min", float(grid.min()), "max", float(grid.max()))
    return meta


def fetch_nac(out, site, lat, lon, half_deg, step_m):
    spec = NAC_DTM.get(site)
    if spec is None:
        print("nac: no NAC DTM registered for site", site)
        return None
    base = out / "dem_nac"
    if base.with_suffix(".json").exists() and base.with_suffix(".f32").exists():
        print("nac: reusing", base.with_suffix(".json"))
        return json.load(open(base.with_suffix(".json")))
    raw = out / "raw" / pathlib.Path(spec["url"]).name
    raw.parent.mkdir(parents=True, exist_ok=True)
    if not raw.exists():
        print("nac: downloading", spec["url"])
        data = fetch(spec["url"], timeout=900)
        if data is None:
            print("nac: not available"); return None
        raw.write_bytes(data)
        print("nac: saved", raw, len(data), "bytes")
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(raw)
    arr = np.array(im, dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[:, :, 0]
    lines, samples = arr.shape
    if (lines, samples) != (spec["lines"], spec["samples"]):
        print("nac: raster shape", arr.shape, "differs from the label", (spec["lines"], spec["samples"]), "- using the raster")
    res, lat0, lon0 = spec["res_px_per_deg"], spec["center_lat_deg"], spec["center_lon_deg"]
    cosl = math.cos(math.radians(lat0))
    def lat_of_line(l): return (spec["line_offset"] - l) / res
    def lon_of_sample(s): return lon0 + (s - spec["sample_offset"]) / (res * cosl)
    lat_max_w, lat_min_w = min(lat + half_deg, spec["lat_max"]), max(lat - half_deg, spec["lat_min"])
    lon_min_w, lon_max_w = max(lon - half_deg, spec["lon_min"]), min(lon + half_deg, spec["lon_max"])
    l0 = int(max(0, math.floor(spec["line_offset"] - lat_max_w * res))); l1 = int(min(lines, math.ceil(spec["line_offset"] - lat_min_w * res)))
    s0 = int(max(0, math.floor(spec["sample_offset"] + (lon_min_w - lon0) * res * cosl))); s1 = int(min(samples, math.ceil(spec["sample_offset"] + (lon_max_w - lon0) * res * cosl)))
    sub = arr[l0:l1, s0:s1]
    nodata = sub <= -1e30
    if nodata.any():
        sub = sub.copy(); sub[nodata] = np.nan
    stride = max(1, int(round(step_m / 2.0)))
    sub = sub[::stride, ::stride]
    # heights are radii in meters (PC_REAL, relative to the 1737.4 km sphere in the LROC RDR)
    if np.nanmedian(sub) > 1.0e6:
        sub = sub - MOON_RADIUS_M
    # fill gaps from neighbors so the grid interpolates everywhere
    if np.isnan(sub).any():
        mask = np.isnan(sub)
        filled = sub.copy(); filled[mask] = np.nanmean(sub)
        sub = filled
    meta = write_grid(base, sub, lat_of_line(l0 + stride * sub.shape[0]), lat_of_line(l0), lon_of_sample(s0), lon_of_sample(s0 + stride * sub.shape[1]),
                      f"LROC NAC DTM (PDS LRO-L-LROC-5-RDR, 2 m/px, resampled every {stride} px), {spec['url']}")
    print("nac: wrote", meta["rows"], "x", meta["cols"], "min", float(np.nanmin(sub)), "max", float(np.nanmax(sub)))
    return meta


def fetch_level(lat, lon, half, layers, target_px=1456):
    """Composite the square patch [lat-half, lat+half] x [lon-half, lon+half] from Trek tiles."""
    lat_min, lat_max, lon_min, lon_max = lat - half, lat + half, lon - half, lon + half
    z = layers[0][1]
    cols_z = 2 ** (z + 1); dpp = 360.0 / cols_z / 256  # degrees per pixel at the finest layer's zoom
    W = H = int(round(2 * half / dpp))
    canvas = None
    for layer, zl, ext in reversed(layers):  # coarse first, so finer mosaics paste over
        cols_l, rows_l = 2 ** (zl + 1), 2 ** zl
        dpp_l = 360.0 / cols_l / 256
        x0 = int(math.floor((lon_min + 180) / 360 * cols_l)); x1 = int(math.floor((lon_max + 180) / 360 * cols_l))
        y0 = int(math.floor((90 - lat_max) / 180 * rows_l)); y1 = int(math.floor((90 - lat_min) / 180 * rows_l))
        n = (x1 - x0 + 1) * (y1 - y0 + 1)
        print(f"  level half={half} layer={layer} z={zl} tiles={n}")
        layer_img = Image.new("LA", ((x1 - x0 + 1) * 256, (y1 - y0 + 1) * 256), (0, 0))
        def one(xy):
            x, y = xy
            data = fetch(TREK.format(layer=layer, z=zl, y=y, x=x, ext=ext))
            return x, y, data
        with ThreadPoolExecutor(max_workers=8) as ex:
            for x, y, data in ex.map(one, [(x, y) for y in range(y0, y1 + 1) for x in range(x0, x1 + 1)]):
                if data is None:
                    continue
                try:
                    tile = Image.open(io.BytesIO(data)).convert("LA")
                except Exception:
                    continue
                layer_img.paste(tile, ((x - x0) * 256, (y - y0) * 256))
        # crop to the patch and resample onto the canvas
        left = ((lon_min + 180) / 360 * cols_l - x0) * 256; top = ((90 - lat_max) / 180 * rows_l - y0) * 256
        right = left + 2 * half / dpp_l; bottom = top + 2 * half / dpp_l
        patch = layer_img.crop((int(round(left)), int(round(top)), int(round(right)), int(round(bottom)))).resize((W, H), Image.BILINEAR)
        if canvas is None:
            canvas = Image.new("L", (W, H), 90)
        canvas.paste(patch.convert("L"), (0, 0), patch.split()[1])
    return canvas, lat_min, lat_max, lon_min, lon_max, dpp


def fetch_imagery(out, lat, lon):
    idir = out / "imagery"; idir.mkdir(parents=True, exist_ok=True)
    meta_path = idir / "imagery.json"
    if meta_path.exists():
        print("imagery: reusing", meta_path)
        return json.load(open(meta_path))
    levels = []
    parent = None
    for k, (half, layers) in enumerate(LEVELS):
        img, lat_min, lat_max, lon_min, lon_max, dpp = fetch_level(lat, lon, half, layers)
        # Match each level's brightness to its parent over the area they share, so the
        # mosaics' different exposures do not draw a square around every level.
        if parent is not None:
            pimg, plat_min, plat_max, plon_min, plon_max = parent
            l = int((lon_min - plon_min) / (plon_max - plon_min) * pimg.width); r = int((lon_max - plon_min) / (plon_max - plon_min) * pimg.width)
            t = int((plat_max - lat_max) / (plat_max - plat_min) * pimg.height); b = int((plat_max - lat_min) / (plat_max - plat_min) * pimg.height)
            ref = np.asarray(pimg.crop((l, t, max(r, l + 1), max(b, t + 1))), dtype=np.float32)
            cur = np.asarray(img, dtype=np.float32)
            if ref.std() > 1e-3 and cur.std() > 1e-3:
                matched = (cur - cur.mean()) * (ref.std() / cur.std()) + ref.mean()
                img = Image.fromarray(np.clip(matched, 0, 255).astype(np.uint8), "L")
        parent = (img, lat_min, lat_max, lon_min, lon_max)
        name = f"level_{k}.jpg"
        img.save(idir / name, quality=88)
        levels.append({"file": name, "half_deg": half, "lat_min": lat_min, "lat_max": lat_max, "lon_min": lon_min, "lon_max": lon_max,
                       "width": img.width, "height": img.height, "deg_per_px": dpp, "m_per_px": dpp * math.pi / 180 * MOON_RADIUS_M,
                       "layers": [l[0] for l in layers]})
        print(f"imagery: level {k} {img.width}x{img.height} {levels[-1]['m_per_px']:.2f} m/px -> {name}")
    meta = {"site": {"lat_deg": lat, "lon_deg": lon}, "levels": levels, "source": "NASA Moon Trek WMTS (LROC WAC global mosaic; LROC NAC Apollo 11 mosaics)"}
    json.dump(meta, open(meta_path, "w"), indent=1)
    return meta


ALBEDO_KERNEL_FRAC = 0.12     # Gaussian radius as a fraction of the patch, see derive_albedo


def derive_albedo(out, meta, kernel_frac=ALBEDO_KERNEL_FRAC):
    """Write an albedo-normalized copy of every imagery level.

    A mosaic tile is radiance, not albedo: the Sun stood somewhere while the frames
    were taken and its shading is baked into the picture. A page that lights the
    surface itself has to start from albedo, or the fixed illumination fights the
    moving one. The estimate here is a flat field: blur the level with a Gaussian
    whose radius is `kernel_frac` of the patch (about 175 px of 1456 by default),
    divide the level by that low-pass and restore the level's own mean. What
    survives is the local contrast; what goes is every brightness gradient broader
    than the kernel.

    Limits, and they are real. (1) Only illumination broader than the kernel is
    removed: the shading inside a crater, and the hard shadow of its rim, are at
    the scale of the crater and stay in the image, so under a sun of the page's own
    the two shadows can disagree. (2) Genuine albedo that varies slowly -- a mare
    and highland boundary, a ray from a young crater -- is flattened along with the
    illumination, so the levels lose their large-scale albedo contrast. (3) The
    division is done on the brightness-matched levels (see `fetch_imagery`), so the
    levels still agree with each other where they overlap. A physical inversion
    (dividing by the reflectance modeled from the DEM and the mosaic's own solar
    geometry) would fix (1) and (2), and needs per-tile illumination metadata Trek
    does not serve.
    """
    idir = out / "imagery"
    levels = meta["levels"]
    for k, lvl in enumerate(levels):
        src = idir / lvl["file"]
        if not src.exists():
            continue
        name = f"{pathlib.Path(lvl['file']).stem}_albedo.jpg"
        dst = idir / name
        if not dst.exists():
            img = Image.open(src).convert("L")
            radius = max(8.0, kernel_frac * min(img.width, img.height))
            low = np.asarray(img.filter(ImageFilter.GaussianBlur(radius=radius)), dtype=np.float32)
            cur = np.asarray(img, dtype=np.float32)
            floor = max(1.0, 0.05 * float(low.mean()))
            flat = cur / np.maximum(low, floor) * float(low.mean())
            Image.fromarray(np.clip(flat, 0, 255).astype(np.uint8), "L").save(dst, quality=88)
            print(f"albedo: level {k} radius {radius:.0f} px -> {name}"
                  f" (mean {cur.mean():.1f} -> {np.clip(flat, 0, 255).mean():.1f}, std {cur.std():.1f} -> {np.clip(flat, 0, 255).std():.1f})")
        else:
            print("albedo: reusing", dst)
        lvl["albedo_file"] = name
    meta["albedo"] = {"method": "flat field: level divided by its own Gaussian low-pass",
                      "kernel_frac": kernel_frac}
    json.dump(meta, open(idir / "imagery.json", "w"), indent=1)
    return meta


SITE_COPY_GLOBS = ("dem_*.json", "dem_*.f32", "site.json")


def copy_site(src, out):
    """Copy an existing site directory's DEMs and imagery originals into `out`.

    The fetched originals (the LOLA range reads, the NAC DTM and the Trek mosaics)
    are what takes the time; everything derived from them is cheap. `--reuse` takes
    them from a directory that already has them so a second copy of the site can be
    derived offline.
    """
    src = pathlib.Path(src)
    if not src.is_dir():
        raise SystemExit(f"--reuse: {src} is not a directory")
    out.mkdir(parents=True, exist_ok=True)
    copied = 0
    for pattern in SITE_COPY_GLOBS:
        for f in sorted(src.glob(pattern)):
            target = out / f.name
            if not target.exists():
                shutil.copy2(f, target)
                copied += 1
    isrc = src / "imagery"
    if isrc.is_dir():
        idst = out / "imagery"
        idst.mkdir(parents=True, exist_ok=True)
        for f in sorted(isrc.glob("level_*.jpg")):
            if f.name.endswith("_albedo.jpg"):
                continue            # derived here, not copied
            if not (idst / f.name).exists():
                shutil.copy2(f, idst / f.name)
                copied += 1
        if (isrc / "imagery.json").exists() and not (idst / "imagery.json").exists():
            shutil.copy2(isrc / "imagery.json", idst / "imagery.json")
            copied += 1
    print(f"reuse: copied {copied} file(s) from {src}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", nargs=2, type=float, metavar=("LAT", "LON"), required=True)
    ap.add_argument("--name", default="apollo11", help="site key for the NAC DTM registry")
    ap.add_argument("--out", required=True)
    ap.add_argument("--lola-half-deg", type=float, default=2.5)
    ap.add_argument("--nac-half-deg", type=float, default=0.06)
    ap.add_argument("--nac-step", type=float, default=4.0, help="NAC DTM resample step, meters")
    ap.add_argument("--no-imagery", action="store_true")
    ap.add_argument("--reuse", default=None, help="copy the DEMs and imagery originals of this site directory instead of fetching them")
    ap.add_argument("--albedo-kernel-frac", type=float, default=ALBEDO_KERNEL_FRAC,
                    help="Gaussian radius of the flat field, as a fraction of the patch")
    ap.add_argument("--no-albedo", action="store_true", help="skip the albedo-normalized copies")
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    lat, lon = a.site
    if a.reuse:
        copy_site(a.reuse, out)
    site = {"lat_deg": lat, "lon_deg": lon, "name": a.name}
    lola = fetch_lola(out, lat, lon, a.lola_half_deg)
    nac = fetch_nac(out, a.name, lat, lon, a.nac_half_deg, a.nac_step)
    imagery = None if a.no_imagery else fetch_imagery(out, lat, lon)
    if imagery is not None and not a.no_albedo:
        imagery = derive_albedo(out, imagery, a.albedo_kernel_frac)
    dems = [{"name": n, **m} for n, m in (("dem_nac", nac), ("dem_lola", lola)) if m]  # finest first
    json.dump({"site": site, "dem": dems, "imagery": "imagery/imagery.json" if imagery else None}, open(out / "site.json", "w"), indent=1)
    print("done:", out / "site.json")


if __name__ == "__main__":
    main()
