# Viewer 3D models

Models the viewer can draw in place of a spacecraft's link boxes through
`export_visualization(prefix; models=Dict(id => path), model_scale=..., model_rotation_deg=...)`
or `spaceagora visualize --model=<id>=<file>`. STL, OBJ and glTF/GLB are
supported; a `.gltf` must embed its buffers, and glTF files must not require
Draco, meshopt or KTX2 (the exporter refuses them and says how to convert).

| File | Source | License | Notes |
|---|---|---|---|
| `iss_nasa_3d_resources_b.glb` | NASA 3D Resources, "International Space Station (ISS) (B)" (github.com/nasa/NASA-3D-Resources) | US government work, public domain (NASA media usage guidelines) | glTF 2.0, Y-up, 36 meshes, 19 materials, no textures, 3.5 MB. NASA's file is Draco-compressed, which the viewer cannot decode, so this copy was decompressed with `npx @gltf-transform/cli copy` (geometry unchanged). Extents 14 x 45 x 46 model units; `model_scale=2.4` puts the truss near 109 m and `model_rotation_deg=(-90, 0, -90)` lays the truss along body y with the modules along body x. |

## Mission display meshes

These six glTF 2.0 assets come from NASA 3D Resources. The source PR records
Draco decompression with glTF-Transform 4.5.0. They contain embedded geometry
and images and need no network or compression decoder. This extraction
adds missing core PNG fallback references alongside the existing WebP
references, preserving all geometry and encoded image bytes. Dimensions
below are decoded scene extents in model units, not a calibrated physical
scale. Select an appropriate `model_scale` and `model_rotation_deg` for the
mission. These visualization meshes are not certified engineering geometry.

| File | NASA source and credited contributor | Extents in model units |
|---|---|---|
| `apollo_lunar_module_nasa_3d_resources.glb` | [Apollo Lunar Module](https://science.nasa.gov/3d-resources/apollo-lunar-module/), NASA/Michael D. Carbajal | 6.427 x 5.013 x 6.427 |
| `cassini_huygens_nasa_3d_resources_a.glb` | [Cassini-Huygens (A)](https://science.nasa.gov/3d-resources/cassini-huygens-a/), NASA/Brian E. Kumanchik | 17.983 x 11.663 x 13.106 |
| `cassini_nasa_3d_resources_a_without_huygens.glb` | [Cassini-Huygens (A), without Huygens](https://science.nasa.gov/3d-resources/cassini-huygens-a/), NASA/Brian E. Kumanchik | 17.983 x 11.663 x 13.106 |
| `cygnss_nasa_3d_resources.glb` | [Cyclone Global Navigation Satellite System](https://science.nasa.gov/3d-resources/cyclone-global-navigation-satellite-system-cygnss/), NASA 3D Resources | 1.590 x 0.512 x 0.244 |
| `magellan_nasa_3d_resources.glb` | [Magellan](https://science.nasa.gov/3d-resources/magellan/), NASA/Brian E. Kumanchik; NASA/Christian A. Lopez | 8.899 x 4.525 x 3.696 |
| `mars_odyssey_nasa_3d_resources.glb` | [Mars Odyssey](https://science.nasa.gov/3d-resources/mars-odyssey/), NASA/JPL/Eyes on the Solar System | 5.669 x 3.360 x 8.221 |

NASA makes its 3D resources available for use subject to its
[media usage guidance](https://www.nasa.gov/nasa-brand-center/images-and-media/).
Acknowledge NASA and the credited contributors. NASA identifiers do not
imply endorsement; separately identified third-party material retains its
own rights. The provenance listing identifies these as NASA resources and
must be retained with redistributed copies.
