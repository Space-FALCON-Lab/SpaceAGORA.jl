#!/usr/bin/env python3
"""Build the viewer's surface textures from their public-domain sources.

Writes equirectangular JPEGs into data/textures/ at the requested tiers and
prints the manifest entries they correspond to. Sources (all NASA/USGS, public
domain):

  earth  NASA Earth Observatory Blue Marble Next Generation, Dec 2004, topo+bathy
         (5400x2700 for 4k, 21600x10800 for 8k/16k)
  mars   NASA Trek WMTS, Viking MDIM 2.1 colorized mosaic (zoom 3 -> 4k, 4 -> 8k)
  venus  USGS Magellan C3-MDIR colorized topography, 6600 m/px (5761x2880 source, 4k only)
  titan  USGS Cassini ISS P19658 global mosaic, 4 km/px (4k only; left edge at 0 E)
  moon   NASA SVS 4720 CGI Moon Kit, LRO WAC colour mosaic (4k and 8k TIFs)

Usage:  python3 scripts/dev/build_viewer_textures.py [--tier 4k|8k|16k] [--body earth,...]
Requires Pillow and network access. The 16k tier is Earth only and is not
committed (its JPEG is ~25 MB); build it locally when a run needs it.
"""
import argparse, io, os, sys, urllib.request
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(REPO, "data", "textures")
TIERS = {"4k": 4096, "8k": 8192, "16k": 16384}
TREK = "https://trek.nasa.gov/tiles/Mars/EQ/Mars_Viking_MDIM21_ClrMosaic_global_232m/1.0.0/default/default028mm"
BMNG_4K = "https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/world.topo.bathy.200412.3x5400x2700.jpg"
BMNG_8K = "https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/world.topo.bathy.200412.3x21600x10800.jpg"
SVS_MOON = "https://svs.gsfc.nasa.gov/vis/a000000/a004700/a004720/lroc_color_poles_{tier}.tif"
USGS_VENUS = "https://planetarymaps.usgs.gov/mosaic/Venus_Magellan_C3-MDIR_ClrTopo_Global_Mosaic_6600m.tif"
USGS_TITAN = "https://planetarymaps.usgs.gov/mosaic/Titan_ISS_P19658_Mosaic_Global_4km.tif"


def fetch(url, timeout=600):
    print("  fetching", url, flush=True)
    return urllib.request.urlopen(url, timeout=timeout).read()


def save(im, name, width):
    im = im.convert("RGB").resize((width, width // 2), Image.LANCZOS)
    path = os.path.join(OUT, name)
    im.save(path, "JPEG", quality=85, optimize=True, progressive=True)
    print(f"  wrote {name} {im.size} {os.path.getsize(path)} bytes", flush=True)


def earth(tier):
    width = TIERS[tier]
    src = fetch(BMNG_4K if tier == "4k" else BMNG_8K)
    save(Image.open(io.BytesIO(src)), f"earth_bluemarble_ng_{tier}.jpg", width)


def mars(tier):
    zoom = {"4k": 3, "8k": 4, "16k": 5}[tier]
    nx, ny = 2 ** (zoom + 1), 2 ** zoom
    mosaic = Image.new("RGB", (256 * nx, 256 * ny))
    for y in range(ny):
        for x in range(nx):
            for attempt in range(3):
                try:
                    tile = urllib.request.urlopen(f"{TREK}/{zoom}/{y}/{x}.jpg", timeout=60).read()
                    break
                except Exception:
                    if attempt == 2:
                        raise
            mosaic.paste(Image.open(io.BytesIO(tile)), (256 * x, 256 * y))
        print(f"  mars row {y + 1}/{ny}", flush=True)
    save(mosaic, f"mars_viking_mdim21_{tier}.jpg", TIERS[tier])


def moon(tier):
    if tier == "16k":
        raise SystemExit("moon: 16k not built (SVS 16k TIF is ~200 MB); use 8k")
    src = fetch(SVS_MOON.format(tier=tier))
    save(Image.open(io.BytesIO(src)), f"moon_lroc_wac_{tier}.jpg", TIERS[tier])


def venus(tier):
    if tier != "4k":
        raise SystemExit("venus: source is 5761 px wide; only the 4k tier is built")
    save(Image.open(io.BytesIO(fetch(USGS_VENUS))), "venus_magellan_clrtopo_4k.jpg", 4096)


def titan(tier):
    if tier != "4k":
        raise SystemExit("titan: source is 4 km/px; only the 4k tier is built")
    save(Image.open(io.BytesIO(fetch(USGS_TITAN))), "titan_cassini_iss_4k.jpg", 4096)


BUILDERS = {"earth": earth, "mars": mars, "venus": venus, "titan": titan, "moon": moon}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="4k", choices=list(TIERS))
    ap.add_argument("--body", default=",".join(BUILDERS))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for body in a.body.split(","):
        print(body, a.tier, flush=True)
        BUILDERS[body](a.tier)
    print("done; add or update the matching [[texture]] entries in data/textures/manifest.toml")
