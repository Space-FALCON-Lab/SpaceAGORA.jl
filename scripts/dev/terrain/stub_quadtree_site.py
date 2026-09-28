#!/usr/bin/env python3
"""Development stub: re-cut an existing nested-square terrain site into the
quadtree tile payload the viewer's horizon renderer expects.

Older site bundles may contain nested square images (`imagery/imagery.json`).
The current `fetch_moon_site.py` already writes a quadtree and needs no conversion. The renderer in `viewer/src/terrain.js` instead
takes a quadtree over one square root region:

    tiles = { scheme: "quadtree", root: {lat_min, lat_max, lon_min, lon_max},
              tile_px: 256, max_level: N,
              nodes: [ {level, x, y, file, m_per_px}, ... ] }

Nodes exist only where imagery was built; the renderer inherits the nearest
present ancestor's texture for the rest. This script builds such a funnel from
sources that are already on disk, so the renderer can be developed and measured
before the real tile pyramid lands:

  * the global surface texture (`data/textures/moon_lroc_wac_4k.jpg`) feeds the
    coarse levels over the whole root, which is sized to cover the descent
    corridor rather than only the site;
  * each nested square of the existing site feeds the levels whose tiles it both
    contains and resolves.

A node is written only when the finest source containing it is at least as sharp
as the node itself, so the coverage narrows level by level exactly as a real
pyramid's does.

    python3 scripts/dev/terrain/stub_quadtree_site.py \
        --site data/terrain/moon/apollo11 --out output/terrain/moon/apollo11_quadtree
"""
import argparse
import json
import math
import os
import shutil
import sys

from PIL import Image

MOON_RADIUS_M = 1737400.0
Image.MAX_IMAGE_PIXELS = None


def deg_to_m(deg, radius_m=MOON_RADIUS_M):
    return math.radians(deg) * radius_m


class Source:
    """One georeferenced equirectangular image, with its own resolution."""

    def __init__(self, path, lat_min, lat_max, lon_min, lon_max, name):
        self.path = path
        self.lat_min, self.lat_max = lat_min, lat_max
        self.lon_min, self.lon_max = lon_min, lon_max
        self.name = name
        self._image = None
        self.gain = 1.0          # brightness match against the site mosaic
        with Image.open(path) as im:
            self.width, self.height = im.size
        self.m_per_px = deg_to_m(lon_max - lon_min) / self.width

    @property
    def image(self):
        if self._image is None:
            self._image = Image.open(self.path).convert("RGB")
        return self._image

    def contains(self, lat_min, lat_max, lon_min, lon_max, eps=1e-9):
        if self.lon_max-self.lon_min >= 360-eps:
            return lat_min >= self.lat_min-eps and lat_max <= self.lat_max+eps and 0 < lon_max-lon_min <= 360
        return (lat_min >= self.lat_min - eps and lat_max <= self.lat_max + eps
                and lon_min >= self.lon_min - eps and lon_max <= self.lon_max + eps)

    def intersects(self, lat_min, lat_max, lon_min, lon_max):
        return not (lat_min >= self.lat_max or lat_max <= self.lat_min
                    or lon_min >= self.lon_max or lon_max <= self.lon_min)

    def crop(self, lat_min, lat_max, lon_min, lon_max, px):
        x0 = (lon_min - self.lon_min) / (self.lon_max - self.lon_min) * self.width
        x1 = (lon_max - self.lon_min) / (self.lon_max - self.lon_min) * self.width
        y0 = (self.lat_max - lat_max) / (self.lat_max - self.lat_min) * self.height
        y1 = (self.lat_max - lat_min) / (self.lat_max - self.lat_min) * self.height
        box = (max(0.0, x0), max(0.0, y0), min(float(self.width), x1), min(float(self.height), y1))
        image = self.image
        if self.lon_max-self.lon_min >= 360-1e-9:
            x0 %= self.width
            x1 = x0+(lon_max-lon_min)/(self.lon_max-self.lon_min)*self.width
            image = Image.new("RGB",(self.width*2,self.height))
            image.paste(self.image,(0,0));image.paste(self.image,(self.width,0))
            box = (x0,max(0.0,y0),x1,min(float(self.height),y1))
        tile = image.resize((px, px), Image.BICUBIC, box=box)
        if abs(self.gain - 1.0) > 1e-3:
            tile = tile.point(lambda v: min(255, int(v * self.gain + 0.5)))
        return tile

    def mean_over(self, lat_min, lat_max, lon_min, lon_max, px=64):
        stat = self.crop_raw(lat_min, lat_max, lon_min, lon_max, px).convert("L")
        return sum(stat.getdata()) / (px * px)

    def crop_raw(self, lat_min, lat_max, lon_min, lon_max, px):
        gain, self.gain = self.gain, 1.0
        try:
            return self.crop(lat_min, lat_max, lon_min, lon_max, px)
        finally:
            self.gain = gain


def load_sources(site_dir, global_texture, texture_lon_left):
    sources = []
    if global_texture and os.path.isfile(global_texture):
        sources.append(Source(global_texture, -90.0, 90.0, texture_lon_left, texture_lon_left + 360.0, "global"))
    imagery_json = os.path.join(site_dir, "imagery", "imagery.json")
    if os.path.isfile(imagery_json):
        with open(imagery_json) as fh:
            imagery = json.load(fh)
        for lvl in imagery["levels"]:
            path = os.path.join(site_dir, "imagery", lvl["file"])
            if not os.path.isfile(path):
                continue
            sources.append(Source(path, lvl["lat_min"], lvl["lat_max"], lvl["lon_min"], lvl["lon_max"], lvl["file"]))
    sources.sort(key=lambda s: s.m_per_px, reverse=True)   # coarse first
    # The global mosaic and the site mosaics are different products: match the
    # global one's brightness to the widest site image over the ground they
    # share, the way the site fetch script matches its own layers.
    if len(sources) > 1 and sources[0].name == "global":
        ref = sources[1]
        here = sources[0].mean_over(ref.lat_min, ref.lat_max, ref.lon_min, ref.lon_max)
        there = ref.mean_over(ref.lat_min, ref.lat_max, ref.lon_min, ref.lon_max)
        if here > 1 and there > 1:
            sources[0].gain = there / here
            print("global mosaic brightness x{:.3f} to match {}".format(sources[0].gain, ref.name))
    return sources


def best_source(sources, lat_min, lat_max, lon_min, lon_max):
    """The sharpest source that fully contains the box (sources are coarse first)."""
    found = None
    for s in sources:
        if s.contains(lat_min, lat_max, lon_min, lon_max):
            found = s
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", default="data/terrain/moon/apollo11", help="site directory written by fetch_moon_site.py")
    ap.add_argument("--out", default="output/terrain/moon/apollo11_quadtree", help="directory to write the stub site into")
    ap.add_argument("--global-texture", default="data/textures/moon_lroc_wac_4k.jpg")
    ap.add_argument("--texture-lon-left", type=float, default=-180.0)
    ap.add_argument("--root-deg", type=float, default=32.0, help="side of the square root region, degrees")
    ap.add_argument("--uprange-km", type=float, default=480.0, help="distance the vehicle starts uprange of the site")
    ap.add_argument("--azimuth-deg", type=float, default=270.0, help="approach azimuth (270 = flying west)")
    ap.add_argument("--tile-px", type=int, default=256)
    ap.add_argument("--max-level", type=int, default=16)
    ap.add_argument("--sharpness", type=float, default=1.5,
                    help="a node is written when its source is at least this much sharper than the node")
    ap.add_argument("--quality-coarse", type=int, default=72)
    ap.add_argument("--quality-fine", type=int, default=86)
    args = ap.parse_args(argv)
    if (not 0 < args.root_deg < 180 or not math.isfinite(args.root_deg)
            or not math.isfinite(args.uprange_km) or args.uprange_km <= 0
            or not math.isfinite(args.azimuth_deg) or not math.isfinite(args.texture_lon_left)
            or not 1 <= args.tile_px <= 4096 or not 0 <= args.max_level <= 16
            or not math.isfinite(args.sharpness) or args.sharpness <= 0
            or not 1 <= args.quality_coarse <= 95 or not 1 <= args.quality_fine <= 95):
        ap.error("Invalid root, tile, level, quality or corridor option.")

    site_dir = os.path.abspath(args.site)
    out_dir = os.path.abspath(args.out)
    if os.path.realpath(site_dir) == os.path.realpath(out_dir):
        raise ValueError("Stub output must differ from the source site directory.")
    with open(os.path.join(site_dir, "site.json")) as fh:
        site_meta = json.load(fh)
    if site_meta.get("tiles"):
        raise ValueError("This site already has a quadtree; use it directly instead of the legacy stub.")
    site_lat = float(site_meta["site"]["lat_deg"])
    site_lon = float(site_meta["site"]["lon_deg"])

    if not all(map(math.isfinite,(site_lat,site_lon))) or not -90 < site_lat < 90:
        raise ValueError("Invalid site coordinates.")
    # Center the root half way along the descent corridor, so the root covers the
    # ground track from powered descent initiation to the site.
    uprange_deg = math.degrees(args.uprange_km * 1000.0 / MOON_RADIUS_M)
    az = math.radians(args.azimuth_deg)
    start_lat = site_lat - uprange_deg * math.cos(az)
    start_lon = site_lon - uprange_deg * math.sin(az) / max(math.cos(math.radians(site_lat)), 1e-6)
    center_lat = 0.5 * (site_lat + start_lat)
    center_lon = 0.5 * (site_lon + start_lon)
    half = 0.5 * args.root_deg
    root = {
        "lat_min": center_lat - half, "lat_max": center_lat + half,
        "lon_min": center_lon - half, "lon_max": center_lon + half,
    }
    if root["lat_min"] < -90 or root["lat_max"] > 90:
        raise ValueError("Stub root crosses a pole; reduce the root size or use the full producer.")
    span = root["lat_max"] - root["lat_min"]

    sources = load_sources(site_dir, args.global_texture, args.texture_lon_left)
    if not sources:
        print("no imagery sources found", file=sys.stderr)
        return 1
    print("root {:.4f}..{:.4f} N, {:.4f}..{:.4f} E ({:.1f} deg, {:.0f} km)".format(
        root["lat_min"], root["lat_max"], root["lon_min"], root["lon_max"], span, deg_to_m(span) / 1000.0))
    for s in sources:
        print("  source {:<22} {:8.2f} m/px  {:.4f}..{:.4f} N {:.4f}..{:.4f} E".format(
            s.name, s.m_per_px, s.lat_min, s.lat_max, s.lon_min, s.lon_max))

    tiles_dir = os.path.join(out_dir, "tiles")
    if os.path.isdir(tiles_dir) and os.listdir(tiles_dir):
        raise ValueError("Output tiles already exist; choose a fresh output directory.")
    os.makedirs(tiles_dir, exist_ok=True)

    nodes = []
    per_level = {}
    # Breadth-first over the quadtree. A node is written when a source that
    # contains it also resolves it; the walk descends into a node whenever a
    # sharper source still overlaps it, so a level with nothing new to say (the
    # global mosaic is too coarse, the site mosaic too small) is simply left
    # empty and the renderer inherits across it.
    frontier = [(0, 0, 0)]
    level = 0
    while frontier and level <= args.max_level:
        written, descend = 0, []
        node_deg = span / (1 << level)
        node_m_per_px = deg_to_m(node_deg) / args.tile_px
        for (lv, x, y) in frontier:
            lon_min = root["lon_min"] + node_deg * x
            lon_max = lon_min + node_deg
            lat_max = root["lat_max"] - node_deg * y
            lat_min = lat_max - node_deg
            src = best_source(sources, lat_min, lat_max, lon_min, lon_max)
            if src is not None and src.m_per_px <= node_m_per_px * args.sharpness:
                quality = args.quality_coarse if level < 6 else args.quality_fine
                name = "l{}_{}_{}.jpg".format(lv, x, y)
                src.crop(lat_min, lat_max, lon_min, lon_max, args.tile_px).save(
                    os.path.join(tiles_dir, name), "JPEG", quality=quality, optimize=True)
                nodes.append({"level": lv, "x": x, "y": y, "file": name,
                              "m_per_px": round(node_m_per_px, 4), "source": src.name})
                written += 1
            if any(s.m_per_px < node_m_per_px and s.intersects(lat_min, lat_max, lon_min, lon_max) for s in sources):
                descend.append((lv, x, y))
        if written:
            per_level[level] = (written, node_m_per_px)
        if not descend:
            break
        frontier = [(level + 1, 2 * x + dx, 2 * y + dy) for (_, x, y) in descend for dx in (0, 1) for dy in (0, 1)]
        level += 1

    if not any(n["level"] == 0 and n["x"] == 0 and n["y"] == 0 for n in nodes):
        raise ValueError("No root imagery is available at the requested sampling; enlarge --root-deg, reduce --tile-px, or supply a finer covering source.")
    max_level = max((n["level"] for n in nodes), default=0)
    index = {"scheme": "quadtree", "root": root, "tile_px": args.tile_px, "max_level": max_level, "nodes": nodes}
    with open(os.path.join(tiles_dir, "tiles.json"), "w") as fh:
        json.dump(index, fh, indent=1)

    total = 0
    print("\n level  tiles   m/px      bytes")
    for lv in sorted(per_level):
        count, mpp = per_level[lv]
        size = sum(os.path.getsize(os.path.join(tiles_dir, "l{}_{}_{}.jpg".format(n["level"], n["x"], n["y"])))
                   for n in nodes if n["level"] == lv)
        total += size
        print("  {:>4}  {:>5}  {:>8.3f}  {:>9,}".format(lv, count, mpp, size))
    print("  total {:>5} tiles, {:,} bytes on disk, about {:,} bytes of base64".format(
        len(nodes), total, int(total * 4 / 3)))

    # The site directory the demo points at: the DEM grids as they are, plus the tile index.
    for d in site_meta["dem"]:
        for ext in (".json", ".f32"):
            shutil.copyfile(os.path.join(site_dir, d["name"] + ext), os.path.join(out_dir, d["name"] + ext))
    out_meta = dict(site_meta)
    out_meta["tiles"] = "tiles/tiles.json"
    out_meta.pop("imagery", None)
    with open(os.path.join(out_dir, "site.json"), "w") as fh:
        json.dump(out_meta, fh, indent=1)
    print("wrote", os.path.join(out_dir, "site.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
