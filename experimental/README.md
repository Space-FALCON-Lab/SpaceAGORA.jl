## Experimental Subsystems

This tree holds typed subsystem scaffolds that are intentionally not part of the
package-owned `src/` surface.

Current experimental modules:

- `experimental/vehicle/resources/`
- `experimental/vehicle/laser_terminal/`
- `experimental/mission/constellation/`
- `experimental/gnc/estimation/`

These files may contain explicit `Not implemented:` failures while the API
shapes are still under evaluation. They are not loaded by `using SpaceAGORA`
and carry no stability guarantee.
