# Constellation Type Design

## Link representation

A link is `(emitter_id, receiver_id)::Tuple{Int, Int}`, where IDs match `SpacecraftModel.id`. If you later need per-link metadata (range, power, schedule, etc.), graduate to a named tuple or small struct without changing the `constellation_struct` interface.

---

## `constellation_struct`

```julia
mutable struct constellation_struct
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}}  # all pairs that geometry could ever allow
    active_links::Vector{Tuple{Int, Int}}    # subset currently firing; mutated at runtime
end

`possible_links` is fixed at construction (every pair the caller's spacecraft could ever link). `active_links` starts as a copy and is the only field the runtime scheduler mutates.

---

## Building one: `build_constellation(spacecraft_model_list)`

Takes any list/array/tuple of pre-built `SpacecraftModel`s — each already carrying its own orbit (`initial_condition`), mass, and bus properties. It reassigns `id` 1..N by position (whatever `id` you passed in is overwritten), computes `possible_links` as every directed pair between distinct satellites, and returns a `constellation_struct` with `active_links` starting as a full copy of `possible_links`.

---

## Mutating links at runtime

- `activate_link!(c, emitter_id, receiver_id)` — adds the pair to `active_links` if it's in `possible_links`; errors otherwise.
- `deactivate_link!(c, emitter_id, receiver_id)` — removes the pair; no-op if already absent.
- `reset_active_links!(c)` — restores `active_links` to a full copy of `possible_links`.

These are called from a `DiffEq` scheduler callback as satellite geometry changes (e.g. a helper crossing into/out of line-of-sight).

---

## Key Design Decisions

| Decision | Rationale |
|---|---|
| Tuple element = spacecraft `id`, not array index | IDs survive reordering and deletion; indices don't |
| `spacecraft` is built entirely by the caller | Keeps `constellation_struct` independent of bus geometry, which differs between ORACLE, CYGNSS, etc. |
| `build_constellation` (re)assigns `id` by position | Callers don't need to hand-number satellites |
| `possible_links` is immutable by convention | Only `active_links` mutates at runtime |

---

## Usage Example (ORACLE-style: 1 target + N helpers)

```julia
make_sc(ic; mass_kg=227.0) = begin
    bus = Link(root=true, m=mass_kg)
    SpacecraftModel([bus], bus, ic=ic, id=0)  # id is reassigned by build_constellation
end

target_ic  = InitialCondition(ra=500e3 + 6_371_000.0, rp=500e3 + 6_371_000.0, i=0.0, ω=0.0, Ω=0.0, ν=0.0)
helper_ics = [InitialCondition(ra=600e3 + 6_371_000.0, rp=600e3 + 6_371_000.0, i=0.0, ω=0.0, Ω=0.0, ν=360.0*(k-1)/10) for k in 1:10]

c = build_constellation([make_sc(target_ic); make_sc.(helper_ics)])

deactivate_link!(c, 3, 1)
activate_link!(c, 5, 1)
```

---

## End-to-End Flow

1. **Build spacecraft** — construct one `SpacecraftModel` per satellite yourself (orbit, mass, bus).
2. **Build the constellation** — `build_constellation(spacecraft_model_list)` assigns IDs and computes `possible_links`/`active_links`.
3. **Wire into simulation** — `DynamicsModel(constellation.spacecraft, effectors...)` feeds directly into `SimulationConfiguration`; the constellation and integrator share the same spacecraft objects.
4. **Propagate** — `run_constellation_ensemble` integrates each spacecraft independently; the RHS reads `constellation.active_links` each step to know which laser links are firing.
5. **Schedule links** — a `DiffEq` callback calls `activate_link!`/`deactivate_link!`/`reset_active_links!` as geometry changes.
6. **Save results** — SpaceAGORA's IO layer writes telemetry via `save_csv`/`generate_plots`. The `constellation_struct` itself is never serialized — only per-spacecraft time-series matter.
