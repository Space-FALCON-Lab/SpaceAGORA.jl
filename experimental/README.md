# Experimental scaffolds

Sketches of possible subsystem interfaces, kept outside `src/`. `using
SpaceAGORA` does not load them and they carry no stability promise.

- `vehicle/resources/`: battery, solar-array, load and power-bus resource
  models (module `Resources`).
- `mission/constellation/network/`: a constellation network topology and state
  (module `ConstellationNetwork`).
- `gnc/estimation/`: an estimator state and measurement packet (module
  `Estimation`).

Their functions only throw an error whose message starts with
`Not implemented:`. `test/gates/ci_scaffold_interface_gate.jl` loads the three
modules and checks that they still do. The supported power budget and onboard
data storage are in `src/models/vehicle/power/` and `src/models/vehicle/data/`,
not here.
