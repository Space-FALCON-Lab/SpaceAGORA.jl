# Native-free grid adapter tests

Prepare one dependency environment for the pinned local wrapper and SpaceAGORA,
then run the suite in a fresh Julia process:

```sh
julia --startup-file=no test/grid_atmosphere/setup.jl
julia --startup-file=no --threads=4 --project=test/grid_atmosphere test/grid_atmosphere/runtests.jl
```

Retrieve the pinned wrapper into `data/GRAMSuite.jl` first. The setup updates only
this test directory's project and ignored manifest. It does not build or initialize
native GRAM. The default tests generate small synthetic grids; the retained Odyssey
regression remains optional and requires the three explicit input paths and hashes.
