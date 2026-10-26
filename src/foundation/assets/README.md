# Asset lookup

`asset_catalog.jl` (`AssetCatalog`) reads the repository's asset manifest,
`data/assets_manifest.toml`, resolves its paths, and checks whether each local
file or folder is ready: missing paths, empty folders, Git LFS pointer files,
zero-byte files and whether the expected bytes are present. It returns one
record per entry and never downloads anything. The command-line tool keeps its
own report types, argument parsing and printing, and calls these checks.

`odyssey_surrogate_assets.jl` (`OdysseySurrogateAssets`) finds and verifies the
assets of the Odyssey frozen-atmosphere example. `rpo_station_assets.jl`
(`RPOStationAssets`) finds the RPO station meshes and samples them into point
clouds. It decodes binary STL with the geometry reader in
`SimulationModel.Structure`; ASCII STL keeps its own parser and sampling.

A ready asset is not a validated one: these checks say a file is present, not
that a model built from it, or a native runtime, is scientifically sound.
Atmosphere preset schemas and their scientific validation stay with the
atmosphere models.
