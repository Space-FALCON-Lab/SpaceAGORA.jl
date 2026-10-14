# Guidance, navigation and control packages

`SpaceAGORAEDG.jl` is the optional, experimental energy-depletion guidance
implementation. The shared optional-package command installs its reviewed source
pair in a separate environment; the companion remains an independent repository.
Load both packages before selecting an EDG model. See the
[installation guide](../../docs/src/user/installation_environment.md) for the
`edg` target and the [package boundary](../../docs/src/maintainer/edg_adapter_boundary.md)
for ownership, compatibility and scientific limits.

SpaceAGORA retains public configuration/model types, physical dynamics, lifecycle
and actuation. The companion owns the eight numerical guidance/control kernels.
Their service boundary does not expose the physical integrator or mutable vehicle
geometry, and it does not claim journal reproduction or flight qualification.

HYPR's source checkout is under [robotics/](../robotics/README.md), although it
also supplies spacecraft rendezvous planning. That remains its supported entry
point; domain grouping does not duplicate its implementation under GNC.
