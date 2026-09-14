# Viewer 3D models

Models the viewer can draw in place of a spacecraft's link boxes through
`export_visualization(prefix; models=Dict(id => path), model_scale=..., model_rotation_deg=...)`
or `spaceagora visualize --model=<id>=<file>`. STL, OBJ and glTF/GLB are
supported; a `.gltf` must embed its buffers.

| File | Source | Licence | Notes |
|---|---|---|---|
| `iss_nasa_3d_resources_b.glb` | NASA 3D Resources, "International Space Station (ISS) (B)" (github.com/nasa/NASA-3D-Resources) | US government work, public domain (NASA media usage guidelines) | glTF 2.0, Y-up, 19 materials, no textures. Extents 14 x 45 x 46 model units; `model_scale=2.4` puts the truss near 109 m and `model_rotation_deg=(-90, 0, -90)` lays the truss along body y with the modules along body x. |
