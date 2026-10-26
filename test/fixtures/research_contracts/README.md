# Research-derived regression fixtures

These small MIT-licensed fixtures come from SpaceAGORA `bfa32f6a44c5f242714b8cf182a451e863119a05`. Full drivers and record tests move to SpaceAGORAMissionReconstructions.

- `save_field_recipes.jl`: the body-rate getter and three output expressions from `benchmarks/studies/cygnss_reconstruction/cygnss_slew.jl:258-288`. Synthetic states and wheel models test the core output and actuator contracts.
- `cygnss_historical.tle`: the unchanged public catalogue fixture used by `scripts/dev/test_cygnss_environment.jl`; its source comments retain the catalogue provenance.
- `pass_exponential_profile.csv`: header and first two rows of `benchmarks/studies/telemetry_validation/manifests/vex_final/step2pp.csv`, retained only for the file-bearing atmosphere configuration regression. Full profile and scientific record assertions belong to the research project.

No private telemetry or mission assets are included.
