# Guidance, navigation and control packages

`SpaceAGORAEDG.jl` is the optional, experimental energy-depletion guidance
(EDG) package. It is pinned here as a submodule; its code, releases and tests
belong to its own repository. Install it in a separate environment with
`scripts/optional_packages.jl setup edg /path/to/edg-env` (see the
[package overview](../README.md) and the
[installation guide](../../docs/src/user/installation_environment.md)), and load
both `SpaceAGORA` and `SpaceAGORAEDG` before selecting an EDG model.

SpaceAGORA keeps the public EDG configuration and model types, the physical
dynamics, the run lifecycle and the panel actuation. The companion computes the
two decisions SpaceAGORA asks of it, `guidance_decision!` and
`control_decision!`, through the service module `SpaceAGORA.EDGServices`
(contract 3.0.0). It never sees the integrator or the vehicle's mutable
geometry. EDG has not been shown to reproduce the published results exactly and
is not approved for flight use; the
[EDG package boundary](../../docs/src/reference/edg_adapter_boundary.md) lists
what each side owns and the current limits.

HYPR also plans spacecraft rendezvous, but its checkout lives under
[robotics/](../robotics/README.md) and is not repeated here.
