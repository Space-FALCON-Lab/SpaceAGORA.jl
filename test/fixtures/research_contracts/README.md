# Research-derived regression fixtures

Small MIT-licensed fixtures taken from SpaceAGORA commit
`bfa32f6a44c5f242714b8cf182a451e863119a05`. The full research drivers and their
record tests now live in the SpaceAGORAMissionReconstructions repository, so the
source paths below exist only at that commit.

- `save_field_recipes.jl`: the body-rate getter and three output expressions
  from `benchmarks/studies/cygnss_reconstruction/cygnss_slew.jl:258-288`. Its
  tests use synthetic states and wheel models to check the core output and
  actuator contracts.
- `cygnss_historical.tle`: the public catalogue fixture, unchanged, that
  `scripts/dev/test_cygnss_environment.jl` reads; its comments keep the
  catalogue provenance.
- `pass_exponential_profile.csv`: the header and first two rows of
  `benchmarks/studies/telemetry_validation/manifests/vex_final/step2pp.csv`,
  kept only for the regression test of a file-backed atmosphere configuration.
  The full profile and the scientific record checks belong to the research
  project.

No private telemetry or mission assets are included.
