# SpaceAGORA Final-Clean Contract

This contract defines the canonical-clean end state after migration:

1. No `legacy_` identifiers in `src/`.
2. No GNC hook ownership under `src/gnc/*/effectors.jl`.
3. No GNC helper ownership under `src/gnc/control/bridge_helpers.jl` and `src/gnc/control/bridge_helpers.jl`.
4. Runtime IO ownership is canonical under `src/io/{config,serialization,logging,outputs}`.
5. `src/simulation/engine/*` orchestrates IO but does not own IO serialization/output implementations.

## Canonical GNC Hook Owners
1. `src/gnc/control/control_hooks.jl`
2. `src/gnc/guidance/guidance_hooks.jl`
3. `src/gnc/navigation/navigation_hooks.jl`

## Canonical Bridge Helper Owner
1. `src/gnc/internal/bridge_helpers.jl`

## Canonical Runtime IO Owners
1. `src/io/config/io_config.jl`
2. `src/io/serialization/io_serialization.jl`
3. `src/io/logging/io_logging.jl`
4. `src/io/outputs/io_outputs.jl`
