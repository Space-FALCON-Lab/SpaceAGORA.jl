# Native-free grid adapter tests

These tests need the pinned GRAMSuite wrapper in `data/GRAMSuite.jl`, but not
its large LFS data. From the repository root, fetch the wrapper source if you
have not already, then prepare one environment for the wrapper and SpaceAGORA
and run the suite in a new Julia process:

```sh
GIT_LFS_SKIP_SMUDGE=1 git submodule update --init --recursive -- data/GRAMSuite.jl
julia --startup-file=no test/grid_atmosphere/setup.jl
julia --startup-file=no --threads=4 --project=test/grid_atmosphere test/grid_atmosphere/runtests.jl
```

Setup changes only this folder's project and its ignored manifest, and it does
not build or initialize native GRAM. The default tests generate small synthetic
grids. The optional Odyssey fine-grid regression needs its three input paths
and their hashes; [the test README](../README.md#grid-adapter-suite) lists the
variables.
