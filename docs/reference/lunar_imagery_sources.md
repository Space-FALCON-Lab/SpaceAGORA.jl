# Higher-resolution imagery for the Apollo 11 landing site viewer

Research note, September 2026. Site: the Apollo 11 lunar module descent stage at
**0.67416 N, 23.47314 E** — this is exactly the Wagner et al. (2017) coordinate for `A11_LM`
(uncertainty 0.3 m), carried in the LROC `ANTHROPOGENIC_OBJECTS` shapefile. The viewer drapes an
imagery quadtree built by `scripts/dev/terrain/fetch_moon_site.py` over this point and the 480 km
descent ground track to its east.

Every number below was measured or retrieved during this study; the reproduction recipe is given
with each. Section 8 lists what could not be verified.

## 1. Summary and recommendation

> **Implemented.** Level 13 was built from `NAC_PHO_E010N0230_M175124932R.PYR.TIF` as recommended
> below; `fetch_moon_site.py` now reads PDS rasters in place (`ARCHIVE_SOURCES`, `ArchiveRaster`)
> and `--max-level` defaults to 13. Measured outcome, on the Apollo 11 site: **46 nodes,
> 253 230 bytes of JPEG (5 505 B a tile at quality 66), 0.338 MB inlined**; the whole quadtree is
> 837 nodes and 3 618 548 bytes, **4.82 MB of base64** against the 5 MB budget, and the built page
> is 13 717 620 bytes against 13 375 049 before — both inside the 14 MB ceiling. Level 14 was not
> built. Three things the study left open are now settled, and one of its numbers was wrong:
>
> - **The `.PYR.TIF` geotransform convention is resolved** (section 9 had it as unverified with a
>   ~50 m residual, which was an arithmetic error). The file is self-describing:
>   `x = R cos(phi_std) (lon - lon_0)`, `y = R lat`, with `R`, `phi_std` = 1 deg and
>   `lon_0` = 180 deg read from its GeoTIFF keys, and `ModelTiepointTag` giving the outer corner of
>   the upper-left pixel (`RasterPixelIsArea`). Against the `.IMG` label's bounding coordinates
>   that puts the four edges within **0.15, 0.00, 0.08 and 0.17 m** — under a pixel. The reader
>   uses the file's own geometry and fails loudly if the label disagrees by more than two pixels.
>   The `.IMG` label's positive `x` is the same magnitude with the opposite sign.
> - **Registration was measured three ways** (`scripts/dev/terrain/check_site_registration.py`).
>   The archive raster sits 0.6 m south and 1.9 m west of the corrected Trek NAC mosaics
>   (imagery against imagery at Trek zoom 15, correlation peak 0.71, the same answer to 0.02 m
>   over 167, 250 and 333 m boxes); 0.9 m north and 3.2 m west of the NAC DTM's hillshade at the
>   approach sun azimuth of 90 deg (peak 0.12); and it draws the lunar module 1.8 m south and
>   0.2 m west of the Wagner coordinate. The three agree to about 2 m, which is the spread of the
>   hillshade-based routes, so `LAYER_REGISTRATION["a11_pho_r"]` takes the imagery-to-imagery
>   number, `(-0.6, -1.9)`. With it in force the level 12 to 13 residual is **+0.05 m north,
>   -0.01 m east**, and in the rendered page a 256 px window either side of the level-13 coverage
>   edge moves by less than 0.05 screen pixels while its high-frequency energy rises 2.4 times.
>   The same run re-measures the Trek defect of section 5 independently: uncorrected, Trek's
>   Apollo 11 mosaic is 23.5 m north of the archive.
> - **Rate limiting was designed around and not hit.** The tile table costs two range requests;
>   the 143 tiles (9.4 MB) covering the level-13 node set and the correlation box came down in
>   **11 requests of about 0.85 MB, in 5.5 s**, one per tile row, each tile written to
>   `pds_cache/` as it arrives so a re-run or a resumed run costs nothing. An HTTP 429 now raises
>   instead of being retried.
>
> The one number that moved: level 13 costs 0.338 MB inlined rather than the 0.449 MB estimated in
> section 7.2, because the shipped quality ramp reaches 66 rather than 68 and the pipeline's own
> detail normalization attenuates this source (gain 0.41).


Trek's zoom cap is **not** the resolution limit of the underlying data: Trek serves the Apollo 11
mosaic 2.5 times coarser than the mosaic itself, and the Planetary Data System serves the same
frames at up to 0.2 m/px with no cap. The recommendation is to add **one** quadtree level
(level 13, 0.325 m/px) sourced from the PDS LROC RDR product
`NAC_PHO_E010N0230_M175124932R.PYR.TIF` by HTTP range reads of its internal 256x256 GeoTIFF tiles,
keeping Trek for levels 0–12. That fits the 5 MB inline budget exactly; a second extra level does
not, and would in any case be past the information content of the data.

A separate defect turned up while measuring this: **Trek's two finest Apollo 11 layers are
misregistered by about 22 m in latitude**, so the viewer's sharpest imagery does not line up with
the landing point the simulation targets. Fixing that is worth more visually than the extra level
is.

## 2. Candidate sources

Pixel scales are as the archive states them. "Real" is the finest structure actually resolved, as
measured here (section 6), not the grid the product is stored on.

| Source | Product identifier | Grid | Real | Covers the site | Where it lives | Access |
|---|---|---|---|---|---|---|
| **LROC NAC, low-altitude pass** (recommended) | `NAC_PHO_E010N0230_M175124932R` (and `...L`), PDS4, LRO-L-LROC-5-RDR-V1.0, `LROLRC_2001/DATA/SDP/NAC_PHO/E010N0230/` | 0.2 m/px | ~0.5–0.6 m features | yes, 23.4483–23.4968 E, 0.3167–1.1142 N | PDS / LROC data node | HTTP range on `.IMG` (14.2 GB float) or on `.PYR.TIF` (915 MB, 8-bit, internally tiled) |
| LROC NAC controlled mosaic, high Sun | `NAC_ROI_APOLLO11HIB_E006N0235`, same volume, `DATA/BDR/NAC_ROI/APOLLO11HIB/` | 0.4 m/px | ~0.6 m | yes, 23.45–23.54 E, 0.15–1.12 N | PDS | HTTP range on `.IMG` (2.04 GB float32) |
| LROC NAC controlled mosaic, low Sun | `NAC_ROI_APOLLO11LOA_E004N0235` | 0.52 m/px | — | yes | PDS | same |
| LROC NAC controlled mosaic, high Sun (A) | `NAC_ROI_APOLLO11HIA_E010N0234` | 1.3 m/px | — | yes | PDS | same |
| LROC NAC, 60x60 km area mosaic | Trek layer `A11_60x60km.eq` | 0.5 m/px source, Trek caps at 0.651 | — | yes, 22.48–24.46 E, −0.32–1.66 N | Moon Trek | WMTS, capped at zoom 15 |
| **What the viewer uses today** | Trek layer `apollo11_26cm_mosaic_byte_geo_1_2_highContrast` ("LRO LROC Image Mosaic 26cm, Apollo 11") | ~0.26 m/px source, Trek caps at **0.651** | ~1.3–2 m | yes, but misregistered ~22 m S | Moon Trek | WMTS, capped at zoom 15 |
| LROC LunaServ WMS | layer `luna_apollo_11_high_resolution_nac_mosaic` | arbitrary, no cap found | ~0.5–0.6 m | yes | ASU, non-PDS interface | OGC WMS `GetMap`, tested to 8192 px |
| Chandrayaan-2 OHRC | image not identified | ~0.25 m GSD from 100 km, 3 km swath | unmeasured | yes (published result) | ISSDC PRADAN | portal, registration required |
| Kaguya Terrain Camera | `Kaguya_TCortho_Mosaic_Global_4096ppd` (7.40 m/px native) | Trek caps at 20.8 m/px | — | yes (near-global) | JAXA DARTS; USGS COGs on AWS | WMTS capped at zoom 10; S3 COGs |
| Chang'e-2 CCD | global DOM | 7 m; 1.5 m only over Sinus Iridum | — | 7 m yes, 1.5 m no | CNSA/CLEP | not investigated |
| Apollo 15/16/17 panoramic camera | — | 1–2 m | — | **no** (those missions did not overfly Tranquility Base) | ASU Apollo Image Archive / LPI | n/a |
| Lunar Orbiter V | site 2 high-resolution frames | ~1–2 m from film | — | yes | LPI Lunar Orbiter Photo Gallery, LOIRP | worse than NAC, not pursued |
| Apollo 11 surface Hasselblad photography | e.g. magazine AS11-40 | sub-centimeter | — | yes, but oblique and **not map-projected** | ASU Apollo Image Archive, LPI Apollo Image Atlas | see section 7 |

## 3. Is Trek's zoom cap the real limit? Yes for Trek, no for the data

**The cap is declared, not incidental.** Each Trek layer's own WMTS capabilities document names the
deepest tile matrix it will serve:

```
https://trek.nasa.gov/tiles/Moon/EQ/<layer>/1.0.0/WMTSCapabilities.xml
```

| Layer | Deepest `TileMatrix` | Finest scale served |
|---|---|---|
| `apollo11_26cm_mosaic_byte_geo_1_2_highContrast` | 15 | 0.6507 m/px |
| `LRO_NAC_Apollo11_Mosaic_p` | 15 | 0.6507 m/px |
| `A11_60x60km.eq` | 15 | 0.6507 m/px |
| `Kaguya_TCortho_Mosaic_Global_4096ppd` | 10 | 20.82 m/px |
| `LRO_WAC_Mosaic_Global_303ppd_v02` | 8 | 83.29 m/px |

At zoom 15 the matrix is 65 536 x 32 768 tiles of 256 px, so 360 deg spans 16 777 216 px:
2.1458e-5 deg/px, which at R = 1737.4 km is **0.6507 m/px**. Tile requests at zoom 16 and 17
return HTTP 404 for both Apollo 11 NAC layers, in both `.png` and `.jpg`. This confirms the earlier
probe and matches the `zmax` values already hard-coded in `fetch_moon_site.py`.

Two of the caps are below the source's own resolution:

- Trek's own catalog entry for the Apollo 11 layer (`TrekServices/ws/index/eq/searchItems?key=apollo11_26cm`)
  reports `title` "LRO LROC Image Mosaic 26cm, Apollo 11", `instrument` "Narrow Angle Camera",
  `BANDS` 1, `PIXEL_DEPTH` "Byte", `fileSize` 1 205 462 934 bytes, bbox
  23.4484963–23.5396925 E by 0.1464532–1.1149077 N. That bbox is 2765 m by 29 365 m; one byte per
  pixel over 1.205e9 pixels implies a grid of about 10 640 x 113 000, i.e. **about 0.26 m/px**.
  Trek serves it at 0.651 m/px, **2.5 times coarser than the mosaic it holds**.
- Kaguya's mosaic is `4096ppd` = 7.40 m/px native; Trek caps it at 20.8 m/px, 2.8 times coarser.
  This affects quadtree levels 6–7, not the site itself.

**What the Trek layer is.** Its bbox matches, to four decimal places, the union footprint of the
NAC frames `M175124932L` and `M175124932R` (below). The layer is that frame pair, map-projected to
a ~26 cm grid and contrast-stretched. The source is therefore directly obtainable from the PDS.

**Trek offers no bulk route.** Four plausible raster subset/export endpoints
(`ws/raster/subset`, `ws/raster/export`, `ws/download/getDownloadInfo`, `TrekWS/rest/subset/bbox`)
all returned 404, as did a WMS `GetCapabilities` under the tile path. Trek's FGDC metadata for the
layer (`TrekWS/rest/cat/metadata/stream?label=...`) has `origin`, `abstract` and `purpose` all set
to "TBD", so it is not usable as a provenance citation either.

## 4. The best orbital imagery that exists for this site

### 4.1 LROC NAC — the record holder

The Orbital Data Explorer at Washington University was queried for every LROC NAC calibrated data
record whose footprint contains the site:

```
https://oderest.rsl.wustl.edu/live2/?query=product&results=m&output=JSON
  &ihid=LRO&iid=LROC&pt=CDRNAC4&target=moon
  &minlat=0.65&maxlat=0.70&westernlon=23.45&easternlon=23.50
```

110 frames come back. Ranked by the archive's `Map_resolution`, the top of the list is:

| Product | Map resolution | Date | Incidence | Emission |
|---|---|---|---|---|
| `M175124932RC` | **0.399 m** | 2011-11-05T09:34:33Z | 40.98 deg | 1.12 deg |
| `M175124932LC` | 0.400 m | 2011-11-05T09:34:33Z | 40.94 deg | 1.67 deg |
| `M162154734RC` | 0.473 m | 2011-06-08 | 73.47 deg | 7.14 deg |
| `M177481212RC` | 0.473 m | 2011-12-02 | 69.07 deg | 6.22 deg |
| `M150361817RC` | 0.474 m | 2011-01-22 | 62.61 deg | 9.41 deg |

`M175124932` is the only frame at the site from a low-altitude pass, and it is 15 percent finer
than anything else. It is the frame LROC published as its closest look at Tranquility Base; the
LROC featured-image page for `M175124932R` gives the altitude as "just 24 km (15 miles) above the
surface".

The PDS volume index record carried in the ODE response gives the exact geometry:

```
52224, 5064, 16,   0.24,  0.56,  0.399, 1.12, 40.98, ... ,  24.03, 1761.55
       lines/samples  ^SCALED_PIXEL_WIDTH  ^HEIGHT  ^RESOLUTION        ^SPACECRAFT_ALTITUDE (km)
```

So: 52 224 lines by 5064 samples, **spacecraft altitude 24.03 km** (target center distance
1761.55 km minus R = 1737.4 km gives 24.15 km, consistent), and a **0.24 m by 0.56 m pixel**.
0.24 m is exactly 24.03 km times the NAC's 10 microradian instantaneous field of view (Robinson et al. 2010: the NACs are monochrome push-broom scanners of 5064 pixels with a 700 mm effective focal length and a 10 microradian IFOV, giving 0.5 m/px from the nominal 50 km orbit — a figure the PDS4 label of `M175124932RC` itself repeats as "the two NACs are monochrome narrow-angle linescan imagers (0.5 m/pixel)").

**This is the honest reading of "25 cm".** The frame samples the ground at 0.24 m across the
detector line but only 0.56 m along track, because the line rate was not sped up for the low pass.
A "26 cm" map product oversamples the along-track direction by a factor of about 2.3. The gain over
Trek's 0.651 m/px is therefore large in one direction and small in the other — see the measurement
in section 6, which finds exactly that anisotropy in the data.

### 4.2 The PDS products derived from it

Retrieved labels, all in `LRO-L-LROC-5-RDR-V1.0`, volume `LROLRC_2001`, host
`pds.lroc.im-ldi.com/data/` (the ASU LROC data node; the same tree is mirrored at
`pds.mcp.nasa.gov/data/store/img/lunar_reconnaissance_orbiter/pds4/lroc/lro-l-lroc-5-rdr/`, which
is the host `fetch_moon_site.py` already uses for the NAC DTM):

**`DATA/SDP/NAC_PHO/E010N0230/NAC_PHO_E010N0230_M175124932R.IMG`** — "Apollo 11 four band
map-projected NAC photometric product from M175124932R. Band 1 is I/F, Bands 2-4 are backplanes
for Phase, Emission, and Incidence."

| | |
|---|---|
| Array | 4 bands x 120 913 lines x 7352 samples, band-sequential |
| Type | `IEEE754LSBSingle` (float32) |
| Header | 29 408 bytes, detached PDS4 label |
| File size | 14 223 267 424 bytes (= 29 408 + 4 x 120 913 x 7352 x 4, exact) |
| Projection | Equirectangular, standard parallel 1.0 deg, central meridian 180 deg |
| Body | Planetocentric, a = b = c = 1737.4 km, positive east |
| Pixel | 0.2 m/px; 6.595577243361499e-06 deg/px in both axes |
| Bounds | 23.448279815425–23.496770554946 E, 0.31671962244668–1.1142051945617 N |
| Upper-left | x = 4 746 449.7 m, y = 33 786.5 m |

The site falls inside those bounds. The `L` frame of the pair is a separate product covering the
adjacent strip.

**`EXTRAS/BROWSE/NAC_PHO/E010N0230/NAC_PHO_E010N0230_M175124932R.PYR.TIF`** — 915 066 604 bytes.
Reading its header and first IFD by range request shows it is far more convenient than the `.IMG`:

| Tag | Value |
|---|---|
| `ImageWidth` / `ImageLength` | 6387 x 105 056 |
| `BitsPerSample` / `SampleFormat` | 8 / unsigned |
| `Compression` | 1 (none) |
| `TileWidth` / `TileLength` | **256 x 256** |
| `TileOffsets` / `TileByteCounts` | 10 275 entries each (25 x 411 tiles) |
| `ModelPixelScale` | 0.230217629559956, 0.230187709412123 m |
| Next IFD | present (pyramid overviews) |

6387 x 0.2302 = 1470 m and 105 056 x 0.23019 = 24 183 m, identical to the `.IMG` extent, so the
browse is the same raster resampled to 0.2302 m/px, 8-bit. Because it is **tiled and
uncompressed**, any single 256x256 tile is a 65 536-byte contiguous byte range. The first IFD sits
at offset 673 382 408 (near the end, so this is not a cloud-optimized GeoTIFF and the index costs
one extra round trip), but after that each tile is one request.

Verified end to end: the site maps to pixel (3274, 57969), tile (12, 226), index 5662; its entry in
`TileOffsets` is 371 064 840 with byte count 65 536; fetching that range returns a 256x256 image in
which the lunar module and its shadow sit at in-tile pixel (202, 112), which is where the
arithmetic puts the Wagner coordinate.

**`DATA/BDR/NAC_ROI/APOLLO11HIB/NAC_ROI_APOLLO11HIB_E006N0235.IMG`** — "Apollo 11 Landing Site High
Sun controlled mosaic in Equirectangular projection centered at .6N, 23.5E. For more info, see
[KLEMETAL2014]." 73 407 x 6950, float32, 0.4 m/px, 2 040 742 400 bytes, built from `M175124932L`
and `M175124932R`. Its LROC product page gives bounds 23.45–23.54 E, 0.15–1.12 N and a bundle
adjustment sigma0 of 1.40. Companion products exist at 5 m and 20 m.

This is the **controlled** product: bundle-adjusted, so its georeferencing is the trustworthy one,
at the cost of being served on a 0.4 m grid rather than 0.2 m.

### 4.3 Other missions

- **Chandrayaan-2 OHRC.** The Orbiter High Resolution Camera has a 0.25 m ground sampling distance
  and a 3 km swath from 100 km. Nagori, Dagar and Rajasekhar, "Age estimation and boulder
  population analysis of the West crater at Apollo 11 landing site using Orbiter High Resolution
  Camera on board Chandrayaan-2 mission", *Planetary and Space Science* **240** (2024) 105828,
  doi:10.1016/j.pss.2023.105828, report an OHRC image of about 0.26 m covering the Apollo 11
  landing site and showing the lunar module, and map more than 8500 boulders around West crater
  with it. So OHRC imagery of this exact site exists at a resolution comparable with the best NAC
  frame. It lives in the ISRO Science Data Archive at `https://pradan.issdc.gov.in/ch2/`, which
  requires account registration; I could not retrieve the product identifier or the raster (see
  section 8). Given that it is at best comparable with `M175124932` and considerably harder to
  obtain, it is not the recommended route.
- **Kaguya Terrain Camera.** 10 m class, near-global, morning and evening illuminations. Archived at
  JAXA DARTS, and as Cloud Optimized GeoTIFFs in the AWS Open Data bucket `astrogeo-ard`
  (`moon/kaguya/terrain_camera/monoscopic/uncontrolled/`), CC0, cite
  `https://doi.org/10.5066/P9SH5YNV`. Relevant only to quadtree levels 6–7, where Trek under-serves
  Kaguya by 2.8x.
- **Chang'e-2.** The CCD stereo camera produced a 7 m global orthophoto map; the 1.5 m data covers
  the Chang'e-3 landing area in Sinus Iridum, not Mare Tranquillitatis. Not useful here.
- **Apollo panoramic and metric cameras.** The panoramic camera flew only on Apollo 15, 16 and 17,
  at 1–2 m below the spacecraft. Those missions did not overfly Tranquility Base, so there is no
  panoramic coverage of this site; secondary sources mention an off-nadir Apollo 16 metric camera
  frame that includes it, which at metric-camera scale and high obliquity is far coarser than NAC.
- **Lunar Orbiter V.** Photographed Apollo landing site 2 before the mission at roughly 1–2 m from
  film, and the LOIRP restorations improved the scans. Strictly worse than NAC and 1967 vintage.

## 5. The registration defect (found while measuring)

The site coordinate the pipeline is given, 0.67416 N 23.47314 E, is not arbitrary: it is the
`A11_LM` row of the LROC `ANTHROPOGENIC_OBJECTS` shapefile
(`EXTRAS/SHAPEFILE/ANTHROPOGENIC_OBJECTS/ANTHROPOGENIC_OBJECTS_360.ZIP`), latitude
0.674160000000000, longitude 23.473140000000001, radius 1 735 473.711 m, `UNCERTAIN` 0.3 m,
`COORD_SRC` "Averaged Mapping Orbit NACs + Laser Ranging", citing Wagner, R. V., Nelson, D. M.,
Plescia, J. B., Robinson, M. S., Speyerer, E. J., and Mazarico, E. (2017), "Coordinates of
anthropogenic features on the Moon", *Icarus* **283**, 92–103, doi:10.1016/j.icarus.2016.05.011.

Measuring where each product actually puts the lunar module, over a 256 m box centered on that
coordinate:

| Product | Offset of the LM / its shadow from the published coordinate |
|---|---|
| PDS `NAC_PHO_...M175124932R.IMG` (brightest pixel) | 0.6 m W, 3.6 m S |
| LunaServ `luna_apollo_11_high_resolution_nac_mosaic` | 4.0 m E, 2.2 m N |
| LunaServ `luna_apollo_11_high_sun_nac_mosaic` | 2.2 m W, 3.5 m N |
| LunaServ `luna_apollo_11_low_sun_nac_mosaic` | 15.5 m W, 2.0 m N |
| LunaServ `luna_apollo_11_moderate_sun_nac_mosaic` | 3.0 m W, 4.0 m N |
| **Trek `apollo11_26cm_mosaic_byte_geo_1_2_highContrast`** | 3.9 m E, **21.4 m S** |
| **Trek `LRO_NAC_Apollo11_Mosaic_p`** | 3.2 m E, **23.3 m S** |
| Trek `A11_60x60km.eq` | 4.5 m E, **12.3 m N** |

East–west scatter is illumination (these are shadow positions and the Sun azimuth differs between
mosaics); the north–south column is the signal. Phase correlation of the whole 256 m box confirms
it independently: the Trek imagery has to move **+25 m north and +8 m east** to register with the
LunaServ mosaics, at a correlation peak of 0.40.

Consequences for the viewer as it stands:

1. Quadtree levels 11 and 12 (from `a11_26cm` and `a11_nac`) are about 22 m south of where they
   should be, so the lander touches down about 22 m north of the imagery's depiction of its own
   landing site.
2. Levels 6–10 come from `A11_60x60km.eq`, which is about 12 m north. Levels 10 and 11 therefore
   disagree with each other by about 34 m, inside the region the camera looks at hardest.
3. The NAC DTM the pipeline already pulls (`NAC_DTM_APOLLO11`, 2 m/px) is a PDS product and is
   presumably registered with the PDS rasters, not with Trek, so the relief and the imagery may
   also be offset from each other.

This is worth fixing whether or not an extra level is added, and moving the fine levels to the PDS
raster fixes it as a side effect.

## 6. What resolution actually exists — measured, not claimed

Radially averaged power spectra, each normalized to its own first non-zero wavenumber (so compare
down a column, and across a row only within one table). "Wavelength" is the full period; the
smallest resolvable feature is about half of it.

**Same 256 m box, same normalization — Trek's output against the archive's, both resampled to
1024 px:**

| Wavelength | Trek zoom 15 (0.651 m/px) | archive source rendered at 0.25 m/px |
|---|---|---|
| 5.02 m | 1.84e-05 | 3.18e-05 |
| 2.00 m | 8.15e-08 | 6.22e-07 |
| 1.43 m | 1.64e-09 | 1.18e-06 |
| 1.11 m | 8.79e-10 | 1.64e-06 |
| 1.02 m | 7.08e-10 | 5.66e-07 |

**Separate run: a 128 m box, asking the archive for finer and finer sampling, to find where its own
information stops:**

| Wavelength | requested 0.25 m/px | 0.125 m/px | 0.0625 m/px |
|---|---|---|---|
| 4.0 m | 7.50e-04 | 7.10e-04 | 7.04e-04 |
| 2.0 m | 3.20e-05 | 2.69e-05 | 2.59e-05 |
| 1.0 m | beyond Nyquist | 1.14e-06 | 9.72e-07 |
| 0.7 m | beyond Nyquist | 2.51e-06 | 1.81e-06 |
| 0.5 m | beyond Nyquist | beyond Nyquist | 1.32e-07 |

Reading:

- **Trek's output really is band-limited at its grid.** Power falls by a factor of 50 between 2.0 m
  and 1.4 m and then sits at a flat 1e-09 floor, exactly as a 0.651 m/px sampling (Nyquist period
  1.30 m) requires. Nothing is being thrown away by the current pipeline's resampling; the cap is
  the whole loss.
- **The archive data carries real structure about three times finer.** At 1.11 m wavelength, in the
  same box and normalization, the archive source has 1.64e-06 against Trek's 8.79e-10 — three
  orders of magnitude more, and it is still rising rather than at a noise floor.
- **It stops at about 0.5–0.7 m.** In the sampling sweep, power holds around 1e-06 to 2e-06 at
  0.7 m and 1.0 m and then falls by an order of magnitude at 0.5 m. So the finest real feature in
  this data is roughly 0.25–0.35 m across, which is what a 0.24 m detector pixel should give.

The gain is anisotropic, as the 0.24 m by 0.56 m raw pixel predicts. One-dimensional cuts through
the two-dimensional spectrum of a 204.8 m patch read straight out of the PDS float raster, each
normalized to its own 4 m value:

| Wavelength | along samples (E–W) | along lines (N–S) | ratio |
|---|---|---|---|
| 2.00 m | 1.34e-01 | 2.08e-01 | 1.55 |
| 1.20 m | 9.00e-03 | 7.62e-03 | 0.85 |
| 0.80 m | 1.17e-02 | 3.91e-03 | 0.33 |
| 0.60 m | 3.46e-03 | 3.78e-05 | 0.01 |
| 0.45 m | 3.07e-03 | 1.11e-04 | 0.04 |

Below about 0.8 m the sample direction retains one to two orders of magnitude more power than the
line direction. The product is genuinely sharper across the frame than along it, by roughly the
2.3x the index record's `SCALED_PIXEL_WIDTH` / `SCALED_PIXEL_HEIGHT` pair implies. (An unexplained
spike at exactly 1.00 m in the sample direction is probably a resampling artifact; it does not
change the conclusion.)

A visual check over an 83 m box — the size of one level-13 node — confirms the numbers: at
0.325 m/px the lunar module, its shadow, the rim of the small crater west of it and the two
deployed experiment packages are all distinctly sharper than at 0.651 m/px.

Reproduce any of the above with:

```
# Trek, zoom 15, the tile containing the site
https://trek.nasa.gov/tiles/Moon/EQ/apollo11_26cm_mosaic_byte_geo_1_2_highContrast/1.0.0/default/default028mm/15/16261/37041.png

# LunaServ, any scale
https://wms.im-ldi.com/lunaserv.php?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap
  &LAYERS=luna_apollo_11_high_resolution_nac_mosaic&STYLES=&SRS=EPSG:4326
  &BBOX=23.468918,0.669938,23.477362,0.678382&WIDTH=1024&HEIGHT=1024&FORMAT=image/png

# PDS float raster, 256 rows around the site (7.5 MB, one request)
Range: bytes=1958308128-1965836575 on
https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/DATA/SDP/NAC_PHO/E010N0230/NAC_PHO_E010N0230_M175124932R.IMG
```

## 7. Feasibility, budget and the plan

### 7.1 How many levels are worth adding

Level L of the quadtree is Trek zoom 3 + L. The corridor coverage rule in `fetch_moon_site.py`
(`level_coverage`, `fine_lod_factor` 2.5) gives, for the existing geometry:

| Level | Trek zoom | m/px | Node side | Nodes | Ground extent of the level's node set |
|---|---|---|---|---|---|
| 11 | 14 | 1.3013 | 333.1 m | 56 | 3.33 x 2.00 km |
| 12 (current deepest) | 15 | 0.6507 | 166.6 m | 51 | 1.50 x 1.00 km |
| **13 (proposed)** | 16 | 0.3253 | 83.3 m | **46** | **0.67 x 0.50 km** |
| 14 | 17 | 0.1627 | 41.6 m | 40 | 0.29 x 0.25 km |
| 15 | 18 | 0.0813 | 20.8 m | 36 | 0.15 x 0.12 km |

Levels 13 and 14 are small enough that **every one of their nodes falls inside every candidate
product's footprint**, including the 1.47 km wide single-frame `NAC_PHO` strip. The narrow swath of
the low-altitude pass is not a constraint for the levels worth adding. (It would be for level 11,
where only 18 of 56 nodes fit inside the `NAC_PHO` strip — level 11 must keep a wider source.)

Level 13 at 0.325 m/px sits just inside the measured 0.5–0.6 m information limit. Level 14 at
0.163 m/px is past it and would only be interpolating.

### 7.2 Byte budget

Current build: 791 nodes, 3.414 MB of JPEG, which is 4.55 MB once base64-inlined into the page.
Measured cost of the new levels, by fetching real tiles at the exact node geometry and putting them
through the pipeline's own `detail_over_parent` normalization (`detail_std` 3.2, `detail_max_gain`
2.2) before encoding:

| Level | JPEG quality | Mean bytes/tile | Range | Nodes | Level total |
|---|---|---|---|---|---|
| 13 | 68 | 7 332 B | 6 872–7 813 | 46 | 0.337 MB |
| 14 | 70 | 6 194 B | 5 894–6 564 | 40 | 0.248 MB |

| Build | JPEG total | Inlined (x 4/3) |
|---|---|---|
| today | 3.414 MB | 4.552 MB |
| + level 13 | 3.751 MB | **4.999 MB** |
| + levels 13 and 14 | 3.999 MB | 5.332 MB |

**One extra level fits the 5 MB budget with about a kilobyte to spare. Two do not.** Adding level 14
would require taking roughly 0.25 MB back out of the coarse levels, and would buy interpolation
rather than detail. Level 13 only.

A useful secondary signal: the detail-band gain the pipeline computes at level 13 is 0.94, meaning
the source has slightly *more* detail than the target amplitude and is being attenuated — real data.
At level 14 the gain is 1.16, mild amplification, consistent with running out of signal.

### 7.3 Source data volume and access

For the level-13 and level-14 node sets (670 x 500 m and 290 x 250 m), at the browse product's
0.2302 m/px, the required source region is about 2911 x 2172 px, which is **108 GeoTIFF tiles of
65 536 bytes each, about 7.1 MB of HTTP range reads**. Alternatively the float `.IMG` can be read as
one contiguous row range: 2500 rows x 7352 samples x 4 bytes = 73.5 MB in a single request
(measured throughput: 7.5 MB in 1.1 s).

**Rate limiting is real and worth designing around.** Issuing roughly 500 small range requests in a
burst to `pds.lroc.im-ldi.com` produced `HTTP 429` from CloudFront with `retry-after: 3600` — a
one-hour block, which also applied to the `pds.mcp.nasa.gov` mirror. Fetch in a few large ranges,
not many small ones. 108 tile reads plus one index read is well inside what the service tolerates,
but the fetcher should cache them on disk exactly as it caches Trek tiles today.

### 7.4 Projection and reprojection

Both the PDS raster and the viewer's quadtree are simple lon/lat grids on a sphere of radius
1737.4 km, so there is no real reprojection, only a resample:

- Quadtree node (level, x, y) already maps to a lon/lat box in `node_pixel_origin` / `build_tiles`.
- The `NAC_PHO` product is Equirectangular with standard parallel 1.0 deg and central meridian
  180 deg, planetocentric, positive east, and its PDS4 label gives the grid directly as
  `west_bounding_coordinate` = 23.448279815425, `north_bounding_coordinate` = 1.1142051945617 and
  `pixel_resolution` = 6.595577243361499e-06 deg/px in both axes. Pixel from lon/lat is therefore
  a subtraction and a division, with no projection maths at all.
- The `.PYR.TIF` shares the extent but is on a 6387 x 105 056 grid, so index it from the same
  lon/lat bounds scaled to those dimensions rather than from its `ModelTiepoint`, which is in
  projected meters under a convention I did not verify (section 8).
- The resample is a decimation: 0.2302 m/px source to 0.3253 m/px node, a factor of 1.41, so an
  area-average or Lanczos downsample, unlike the current integer-pixel crop out of Trek tiles.
  This is the one place the new path is less exact than the old one, and it is unavoidable because
  the archive grid is not a power-of-two relative of Trek's.

### 7.5 Licensing and attribution

The LROC terms page (`https://lroc.im-ldi.com/about/terms`) draws a line that matters here:

> "Data products available through the LROC Data Node of the PDS or through NASA PDS interfaces are
> in the public domain."

> "Images available on the LROC website, including but not limited to Featured Images, Featured
> Sites, curated products, downloads and non-PDS interfaces listed above, are copyrighted."

> "All commercial use of LROC images not downloaded from the LROC PDS archive require prior
> permission."

> "For news media or educational purposes, credit LROC images with the text tag NASA/GSFC/Arizona
> State University. If space is limited, the shorter credit NASA/GSFC/ASU is also acceptable."

So the **PDS raster is public domain and the LunaServ WMS is not**. LunaServ is excellent for
prototyping and for the measurements in this note, but a shipped, redistributable page should take
its pixels from the PDS product. The `source` string in `tiles.json` and the viewer's attribution
should then carry:

- NASA/GSFC/Arizona State University.
- Robinson, M. S., et al. (2010), "Lunar Reconnaissance Orbiter Camera (LROC) Instrument Overview",
  *Space Science Reviews* **150**, 81–124 — the instrument reference the LROC terms page asks for.
- Klem, S. M., et al. (2014), "Controlled LROC Narrow Angle Camera High Resolution Mosaics",
  *Lunar and Planetary Science Conference* **45**, abstract #2885 — if a `NAC_ROI` controlled
  mosaic is used.
- Wagner, R. V., et al. (2017), *Icarus* **283**, 92–103, doi:10.1016/j.icarus.2016.05.011 — for the
  site coordinate, which the pipeline already depends on.

### 7.6 Implementation steps

1. **Fix the registration first**, independently of resolution. Either apply the measured
   +25 m north, +8 m east shift to the `a11_26cm` and `a11_nac` Trek layers (and the corresponding
   correction to `A11_60x60km.eq`), or, better, move levels 11 and 12 onto
   `NAC_ROI_APOLLO11HIB_E006N0235` (0.4 m/px, controlled, covers 23.45–23.54 E which spans the
   level-12 node set and most of level 11). Verify by re-measuring the LM position against the
   `A11_LM` shapefile row. This alone is a visible improvement.
2. **Add an archive source alongside `TrekTiles`** in `fetch_moon_site.py`: a small class that,
   given a lon/lat box and a pixel count, reads `NAC_PHO_E010N0230_M175124932R.PYR.TIF` by HTTP
   range. It needs the first IFD read once (offset 673 382 408, 22 entries) to recover
   `TileOffsets` (10 275 uint32 at 673 423 794) and `TileByteCounts` (at 673 382 694), then one
   65 536-byte read per tile. Cache tiles on disk under `pds_cache/` next to the existing
   `trek_cache/`, keyed by product and tile index. Do not hard-code the IFD offsets: read them from
   the TIFF header, and fail loudly if the layout changes.
3. **Extend `ZOOM_LAYERS`** with an entry for Trek zoom 16 that names the archive source rather than
   a Trek layer, and raise `--max-level` from 12 to 13. Keep zooms up to 15 exactly as they are so
   the existing levels do not change.
4. **Resample properly**: area-average the 0.2302 m/px source into the node's 256 px grid. The
   existing `detail_over_parent` machinery then applies unchanged; it already normalizes each
   level's own detail band, so the new level will blend with level 12 without a tone step provided
   step 1 has removed the geometric step.
5. **Do not add level 14.** It costs 0.25 MB, pushes the inlined page to 5.33 MB, and carries no
   information the measurement can find.
6. **Record the ceiling in the metadata**: add the source resolution and the measured information
   limit to `tiles.json` so the viewer can stop asking for more and switch to the treatment in
   section 8 below that threshold, rather than blurring an unlabelled texture.
7. **Update the attribution** in `tiles.json`'s `source` field and wherever the viewer surfaces it.

**Status of these steps.** 2, 3, 4, 5, 6 and 7 are done, as described in section 1: `tiles.json`
carries `archive_sources` (product, URL, credit, citation and the levels it feeds), `attribution`
and a `resolution` block, and the viewer's info panel now reads "finest 0.33 m/px sampling of
~0.5 m detail" rather than the sampling alone. Step 1 was done in the preceding round for the
Trek layers and this round extends the same `LAYER_REGISTRATION` mechanism to the archive raster
rather than moving levels 11 and 12 onto `NAC_ROI_APOLLO11HIB`: the measurement in section 1
shows the corrected Trek levels and the archive level agree to 0.05 m, so moving them would buy
nothing.

## 8. The honest ceiling, and what the page should do below it

At 4 m altitude a filled viewport asks for roughly **0.02 m/px**. The best map-projected orbital
imagery of Tranquility Base resolves features of about **0.5 m**, so the page is asking for about
25 times more detail than exists, and about 12 times more than the raw 0.24 m detector pixel could
ever have supplied. Adding level 13 closes a factor of two of that gap. Nothing closes the rest.

Where each thing stops:

| Scale | What exists |
|---|---|
| 0.65 m/px | what the viewer has today (Moon Trek cap) |
| 0.33 m/px | one more level, real, recommended |
| 0.24 m/px | the LROC NAC detector pixel from 24 km, across the frame only; 0.56 m along it |
| ~0.5 m feature size | the measured information limit of every orbital product over this site |
| 0.02 m/px | what the final seconds of the descent ask for — **no orbital data at any scale** |

Real imagery does exist at that scale, but it is not orbital and not map-projected:

- **Apollo 11 surface Hasselblad photography** (70 mm, magazines including AS11-40), digitized by
  the ASU Apollo Image Archive. It is sub-centimeter on the regolith in the foreground, and it
  looks at Tranquility Base from a meter and a half above it. But it is oblique, uncalibrated for
  photogrammetry outside the panoramas, covers a few tens of meters, and would need a full
  structure-from-motion reconstruction plus a terrain solution to become a draped texture.
- **The 16 mm sequence camera descent film**, which is the only imagery *of the descent itself*, but
  is low resolution, motion-blurred and also unprojected.

Neither turns into a quadtree tile without a research project of its own. The realistic answer for
the last two orders of magnitude is therefore **a blended procedural treatment**, and the page
should be explicit about the handover rather than upsampling silently:

1. **Down to 0.33 m/px**: real LROC NAC imagery, as above.
2. **From 0.33 m/px to roughly 0.05 m/px**: keep the real image as the low-frequency albedo and add
   a synthesized high-frequency band on top — a regolith micro-texture (fractal or noise-based
   normal and albedo detail, tuned so the crater and boulder statistics match what is measured in
   the level-13 tiles). This is the standard detail-texture approach and is honest as long as the
   viewer says so; the existing `detail_over_parent` code already thinks in terms of a per-level
   detail band, so the extension is natural.
3. **Below 0.05 m/px** (the last few meters of the descent): the high-frequency band dominates
   entirely and should be driven by the plume interaction rather than by imagery — the blast zone,
   the scoured area and the dust the descent engine raises are what a viewer at 4 m altitude is
   actually looking at, and those are being modeled elsewhere in this project.
4. **Label it.** A small annotation in the viewer when the camera goes below the data's resolution
   ("below 0.3 m/px: synthesized surface detail") costs nothing and keeps the page from implying
   that NASA photographed the regolith at 2 cm from orbit.

A further option, worth noting but not recommended now: the `NAC_DTM_APOLLO11` digital terrain
model at 2 m/px already in the pipeline could drive a physically-motivated shading model at fine
scales instead of a texture, but at 2 m it is coarser than the imagery and would not help below
0.33 m/px either.

## 9. What could not be verified

- **The Chandrayaan-2 OHRC image of this site.** The published result (doi:10.1016/j.pss.2023.105828)
  is solid and states about 0.26 m and coverage of the Apollo 11 landing site, but the paper is
  paywalled (ScienceDirect returned HTTP 403) and I could not read its methods section, so I have
  **no product identifier, no acquisition date and no footprint** from a primary source. PRADAN
  requires account registration, which this study did not do, so I did not confirm the image is
  actually downloadable, at what processing level, or under what terms. ISRO's stated policy is
  release without a lock-in period, but I did not verify that for OHRC specifically.
- **Which ground direction the 0.24 m and 0.56 m pixel dimensions correspond to.** The index record
  gives `SCALED_PIXEL_WIDTH` 0.24 m and `SCALED_PIXEL_HEIGHT` 0.56 m and a `NORTH_AZIMUTH` of
  268.22 deg, and the spectral measurement independently finds the map-projected product sharper
  east–west than north–south. I did not confirm the mapping from detector axes to ground axes from
  the instrument documentation; the two observations are consistent but the causal chain is
  inferred.
- **A Kaguya Terrain Camera COG covering this site** in the AWS `astrogeo-ard` bucket. The bucket's
  objects are named by latitude and longitude, and the first 6000 keys contained nothing near
  0.7 N 23.5 E. The bucket is described as uncontrolled per-observation products and may not be
  global. Trek remains the verified Kaguya route, at its 20.8 m/px cap.
- **Whether Moon Trek offers any bulk or subset download.** Four endpoint guesses returned 404.
  There may be an authenticated or UI-driven path I did not find; I only established that no
  obvious public REST route exists.
- ~~**The `.PYR.TIF` geotransform convention.**~~ **Resolved** while implementing (section 1):
  `x = R cos(1 deg) (lon - 180 deg)` and `y = R lat` reproduce `ModelTiepoint` from the label's
  bounding coordinates to 0.15 m and 0.17 m, so the roughly 50 m residual recorded here was an
  arithmetic error in this study, not a property of the product. The shipped reader takes the
  geometry from the file rather than from the label.
- **Apollo-era orbital coverage.** The statement that Apollo 15, 16 and 17 did not overfly
  Tranquility Base, and that an off-nadir Apollo 16 metric camera frame includes the site, comes
  from secondary sources (LPI and ASU outreach pages) rather than from an index of the frames. I did
  not search the Apollo photographic index for a frame covering 0.67 N 23.47 E. Since the panoramic
  camera is 1–2 m at best, this does not affect the recommendation.
- **The provenance of Trek's ~22 m latitude offset.** It is reproducible and large, and the
  corrected position agrees with a 0.3 m published coordinate, so the direction of the error is not
  in doubt. But I did not establish *why* Trek's copies are shifted — whether the mosaics were
  ingested uncontrolled, or a datum or tiepoint was mishandled in tiling.
- **Long-term stability of the hosts used here.** `lroc.asu.edu`, `pds.lroc.asu.edu`,
  `wms.lroc.asu.edu` and `ser.sese.asu.edu` all now redirect to `*.im-ldi.com` domains. The PDS
  content is mirrored at `pds.mcp.nasa.gov`, which is the more durable choice for a pipeline; I did
  not find a NASA-hosted equivalent for the LunaServ WMS, which is a further reason not to depend
  on it.
