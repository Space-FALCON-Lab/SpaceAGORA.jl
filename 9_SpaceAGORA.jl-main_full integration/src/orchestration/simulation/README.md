# Simulation responsibilities

`config/` defines the physical scenario and its declared run settings.
`engine/config/` assembles execution settings and reuses the scenario's solver
configuration type. `engine/` prepares and integrates a run; `callbacks/` applies
scheduled or event-driven updates; `campaigns/` coordinates repeated runs.

`composition/` assembles the existing `SimulationModel` module. `state/` owns
concrete run buffers and caches. Their defining module identities remain stable
for callers and Julia serialization.

`callback_services.jl` declares the internal state and sampling operations used
by callbacks. The engine implements those operations by delegating to its
existing methods, preserving the state representations, frames, policies and
mutation rules. Callbacks can use the declared contract without looking up the
engine module dynamically.

Mission policies live in `../mission/`; parallel routing, cost and workers live
in `../parallel/`. Result formats and checkpoint storage belong to Interfaces.
