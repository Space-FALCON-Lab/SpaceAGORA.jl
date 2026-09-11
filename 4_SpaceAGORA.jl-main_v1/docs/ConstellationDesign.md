# Constellation Type Design

## The Link Tuple

Each link needs to know which two satellites it connects. The simplest representation is:

```julia
# (emitter_id, receiver_id)
# IDs match SpacecraftModel.id
const LinkPair = Tuple{Int, Int}
```

If you later need to store metadata per link (range, power, schedule, etc.), you can graduate to a named tuple or a small struct without changing the `Constellation` interface.

---

## The `Constellation` Struct

```julia
# abstract type lets multiple dispatch drive the right constructor
abstract type ConstellationPattern end

mutable struct Constellation
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}}  # all pairs that geometry could ever allow
    active_links::Vector{Tuple{Int, Int}}    # subset currently firing; mutated at runtime
    pattern::ConstellationPattern            # retains the geometry parameters for introspection
end
```

`possible_links` is populated at construction from the geometry (e.g. every pair within range, or every adjacent pair in a Walker plane). `active_links` starts as a copy of `possible_links` and is mutated by the scheduler during the simulation.

---

## Pattern Types and Constructors

```julia
Base.@kwdef struct WalkerPattern <: ConstellationPattern
    altitude_m::Float64
    inclination_deg::Float64
    T::Int              # total satellites
    P::Int              # orbital planes
    F::Int              # phasing parameter
    mass_kg::Float64 = 100.0
    ecc::Float64 = 0.0
end

Base.@kwdef struct OraclePattern <: ConstellationPattern
    target_altitude_m::Float64
    helper_altitude_m::Float64
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    n_helpers::Int = 10
    mass_kg::Float64 = 227.0
end

Base.@kwdef struct FlowerPattern <: ConstellationPattern
    altitude_m::Float64
    inclination_deg::Float64
    petals::Int
    loops::Int
    mass_kg::Float64 = 100.0
end
```

```julia
# Builds a Walker T/P/F constellation: evenly distributes T satellites across P planes with phasing F.
# Inputs: WalkerPattern (geometry params), make_spacecraft (factory returning SpacecraftModel).
# Output: Constellation with spacecraft, possible_links (adjacent in-plane + cross-plane pairs), active_links.
function Constellation(p::WalkerPattern, make_spacecraft::Function)
    sats_per_plane = p.T ÷ p.P
    raan_spacing   = 360.0 / p.P
    nu_spacing     = 360.0 / sats_per_plane
    phase_offset   = p.F * (360.0 / p.T)

    spacecraft = SpacecraftModel[]
    id = 1
    for plane in 0:(p.P - 1)
        raan = plane * raan_spacing
        for slot in 0:(sats_per_plane - 1)
            nu = mod(slot * nu_spacing + plane * phase_offset, 360.0)
            ic = InitialCondition(
                ra = p.altitude_m,
                rp = p.altitude_m,
                i  = p.inclination_deg,
                ω  = 0.0,
                Ω  = raan,
                ν  = nu
            )
            push!(spacecraft, make_spacecraft(id, ic, p))
            id += 1
        end
    end

    possible = _walker_adjacent_pairs(spacecraft, p)
    return Constellation(spacecraft, possible, copy(possible), p)
end

# Builds an ORACLE constellation: one debris-target satellite plus N evenly spaced helper satellites at a different altitude.
# Inputs: OraclePattern (target/helper altitudes, inclinations, helper count), make_spacecraft factory.
# Output: Constellation where possible_links are all helper→target pairs; helpers orbit the same plane.
function Constellation(p::OraclePattern, make_spacecraft::Function)
    spacecraft = SpacecraftModel[]
    ic_target = InitialCondition(
        ra = p.target_altitude_m, rp = p.target_altitude_m,
        i  = p.target_inclination_deg, ω = 0.0, Ω = 0.0, ν = 0.0
    )
    push!(spacecraft, make_spacecraft(1, ic_target, p))

    for k in 1:p.n_helpers
        nu = 360.0 * (k - 1) / p.n_helpers
        ic = InitialCondition(
            ra = p.helper_altitude_m, rp = p.helper_altitude_m,
            i  = p.helper_inclination_deg, ω = 0.0, Ω = 0.0, ν = nu
        )
        push!(spacecraft, make_spacecraft(k + 1, ic, p))
    end

    # all helpers can potentially link to the target
    possible = [(spacecraft[k+1].id, spacecraft[1].id) for k in 1:p.n_helpers]
    return Constellation(spacecraft, possible, copy(possible), p)
end

# Flower: places `p.petals * p.loops` satellites on a repeating ground-track orbit.
# Petals set the number of ground-track lobes; loops distribute satellites along each lobe.
# Adjacent satellites within the same petal are linked as possible_links.
function Constellation(p::FlowerPattern, make_spacecraft::Function)
    spacecraft = SpacecraftModel[]
    id = 1
    n_sats = p.petals * p.loops
    for petal in 0:(p.petals - 1)
        raan = petal * (360.0 / p.petals)
        for loop in 0:(p.loops - 1)
            nu = loop * (360.0 / p.loops)
            ic = InitialCondition(
                ra = p.altitude_m,
                rp = p.altitude_m,
                i  = p.inclination_deg,
                ω  = 0.0,
                Ω  = raan,
                ν  = nu
            )
            push!(spacecraft, make_spacecraft(id, ic, p))
            id += 1
        end
    end

    possible = Tuple{Int,Int}[]
    for petal in 0:(p.petals - 1)
        base = petal * p.loops
        for loop in 0:(p.loops - 1)
            i = spacecraft[base + loop + 1].id
            j = spacecraft[base + mod(loop + 1, p.loops) + 1].id
            push!(possible, (i, j))
        end
    end
    return Constellation(spacecraft, possible, copy(possible), p)
end
```

---

## Link Update Helpers

```julia
# Adds (emitter_id, receiver_id) to active_links if it exists in possible_links; errors otherwise.
# Inputs: Constellation (mutated), emitter/receiver satellite IDs (Int).
# Output: nothing; mutates c.active_links in place.
function activate_link!(c::Constellation, emitter_id::Int, receiver_id::Int)
    pair = (emitter_id, receiver_id)
    pair in c.possible_links || throw(ArgumentError("($emitter_id, $receiver_id) is not a possible link."))
    pair in c.active_links   || push!(c.active_links, pair)
end

# Removes (emitter_id, receiver_id) from active_links; silently a no-op if the pair is already absent.
# Inputs: Constellation (mutated), emitter/receiver satellite IDs (Int).
# Output: nothing; mutates c.active_links in place.
function deactivate_link!(c::Constellation, emitter_id::Int, receiver_id::Int)
    filter!(p -> p != (emitter_id, receiver_id), c.active_links)
end

# Restores active_links to the full set of possible_links, undoing any runtime activations/deactivations.
# Input: Constellation (mutated).
# Output: nothing; active_links is a fresh copy of possible_links after the call.
function reset_active_links!(c::Constellation)
    resize!(c.active_links, length(c.possible_links))
    copyto!(c.active_links, c.possible_links)
end
```

---

## Private Helper: Walker Adjacency

```julia
# Computes within-plane and cross-plane adjacent satellite ID pairs for a Walker constellation.
# Inputs: ordered spacecraft vector, WalkerPattern (for T and P counts).
# Output: Vector{Tuple{Int,Int}} of (emitter_id, receiver_id) pairs representing all adjacent links.
function _walker_adjacent_pairs(spacecraft::Vector{SpacecraftModel}, p::WalkerPattern)
    sats_per_plane = p.T ÷ p.P
    pairs = Tuple{Int, Int}[]
    for plane in 0:(p.P - 1)
        base = plane * sats_per_plane
        for slot in 0:(sats_per_plane - 1)
            # within-plane: link each satellite to the next in its orbital plane
            i = spacecraft[base + slot + 1].id
            j = spacecraft[base + mod(slot + 1, sats_per_plane) + 1].id
            push!(pairs, (i, j))
            # cross-plane: link to the co-slot satellite in the adjacent plane
            next_plane = mod(plane + 1, p.P)
            k = spacecraft[next_plane * sats_per_plane + slot + 1].id
            push!(pairs, (i, k))
        end
    end
    return pairs
end
```

---

## Usage Example

```julia
# provide your own spacecraft factory matching your bus geometry
make_sc = (id, ic, p) -> begin
    bus = Link(root=true, m=p.mass_kg)
    SpacecraftModel([bus], bus, ic=ic, id=id)
end

# Walker 24/3/1 at 700 km, 53° inclination
c = Constellation(
    WalkerPattern(altitude_m=700e3 + 6_371_000.0, inclination_deg=53.0, T=24, P=3, F=1),
    make_sc
)

# mutate active links at runtime (e.g. inside a DiffEq callback)
deactivate_link!(c, 3, 1)
activate_link!(c, 5, 1)
```

---
## Key Design Decisions

| Decision | Rationale |
|---|---|
| Tuple element = spacecraft `id`, not array index | IDs survive reordering and deletion; indices don't |
| `make_spacecraft` is caller-supplied | Keeps `Constellation` independent of bus geometry, which differs between Walker, ORACLE, CYGNSS, etc. |
| `possible_links` is immutable by convention | Only `active_links` is mutated at runtime. If range-gating should also shrink `possible_links`, that can be added explicitly |
---

## How a Simulation Actually Runs — A Plain-Language Walkthrough

This section traces one complete ORACLE-style run from Julia `main` to final CSV output, naming which function does what at each step.

### Step 1 — Describe your constellation geometry (before any physics)

You pick a pattern type and fill in the numbers. For ORACLE that is one debris target at 500 km and ten helpers at 600 km:

```julia
pattern = OraclePattern(
    target_altitude_m  = 500e3 + 6_371_000.0,
    helper_altitude_m  = 600e3 + 6_371_000.0,
    n_helpers          = 10,
    mass_kg            = 227.0
)
```

Nothing moves yet. This is just a data record describing the shape of the constellation.

---

### Step 2 — Build the spacecraft objects and the link table

`Constellation(p::OraclePattern, make_spacecraft)` is called. Internally it:

1. Calls your `make_spacecraft` factory once for the target (id = 1) and ten times for the helpers (ids 2–11), each with a different `InitialCondition` (true anomaly spaced 36° apart).
2. Builds `possible_links = [(2,1), (3,1), …, (11,1)]` — every helper can potentially illuminate the target.
3. Sets `active_links = copy(possible_links)` so all links start open.

At this point you have a `Constellation` object: a vector of `SpacecraftModel` plus two link lists. No orbit has been propagated.

---

### Step 3 — Wire the constellation into a `SimulationConfiguration`

SpaceAGORA's main configuration struct wraps your constellation's spacecraft list together with the physics models you want:

```julia
dyn = DynamicsModel(constellation.spacecraft, (J2Model(), SRPModel()))
args = SimulationConfiguration(dynamics_model=dyn, …)
```

The `DynamicsModel` constructor accepts the spacecraft vector directly from `Constellation.spacecraft`, so the constellation and the integrator share the same objects.

---

### Step 4 — Run the integrator (one satellite at a time, or all together)

For uncoupled satellites SpaceAGORA calls `run_constellation_ensemble`, which propagates each spacecraft in `constellation.spacecraft` as an independent ODE from `t = 0` to `t = T_sim`. At each integration step the dynamics RHS evaluates forces (gravity, drag, SRP, laser) on that satellite and advances its Keplerian + attitude state.

For ORACLE the laser force terms need to know which links are currently active. The integrator reads `constellation.active_links` at every callback timestep.

---

### Step 5 — The scheduler mutates active links at runtime (inside a callback)

Inside a `DiffEq` discrete callback — which fires every scheduler period (e.g. every 60 s) — you call the link-management helpers:

```julia
# Geometry check says helper 5 is now behind the limb → turn it off
deactivate_link!(constellation, 5, 1)

# Helper 7 just crossed into line-of-sight → switch it on
activate_link!(constellation, 7, 1)
```

`deactivate_link!` just `filter!`s the pair out of `active_links`. `activate_link!` pushes it back in after confirming it is in `possible_links`. These are O(n_links) vector operations; they are fast enough for a 10-satellite case and cheap relative to the ODE step cost.

If a new orbit pass starts and you want to reset to "all links open" without listing each pair manually, call:

```julia
reset_active_links!(constellation)
```

---

### Step 6 — Post-processing and output

After the integrator finishes, SpaceAGORA collects telemetry (orbital elements, ΔV, mass, laser power delivered) and writes results via `save_csv` / `generate_plots`. The `Constellation` object itself is not serialised — only the per-spacecraft time-series matter for analysis.

---

### How the functions map to these steps

| Step | Function(s) used |
|---|---|
| 1 — Describe geometry | `WalkerPattern`, `OraclePattern`, `FlowerPattern` structs |
| 2 — Build spacecraft + links | `Constellation(p::OraclePattern, …)`, `Constellation(p::WalkerPattern, …)`, `_walker_adjacent_pairs` |
| 3 — Wire into simulation | `DynamicsModel(constellation.spacecraft, …)` (SpaceAGORA core) |
| 4 — Propagate orbits | `run_constellation_ensemble` / `SimulationEngine` (SpaceAGORA core) |
| 5 — Schedule links at runtime | `activate_link!`, `deactivate_link!`, `reset_active_links!` |
| 6 — Save results | SpaceAGORA IO layer (not part of this file) |
