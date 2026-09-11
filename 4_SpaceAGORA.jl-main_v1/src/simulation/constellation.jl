module Constellations

using ..SpacecraftModels: SpacecraftModel, InitialCondition

export ConstellationPattern, WalkerPattern, OraclePattern, FlowerPattern
export Constellation
export activate_link!, deactivate_link!, reset_active_links!

# ── Pattern types ──────────────────────────────────────────────────────────────

abstract type ConstellationPattern end

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

# ── Constellation struct ───────────────────────────────────────────────────────

mutable struct Constellation
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}}  # all pairs that geometry could ever allow
    active_links::Vector{Tuple{Int, Int}}    # subset currently firing; mutated at runtime
    pattern::ConstellationPattern            # retains the geometry parameters for introspection
end

# ── Constructors ───────────────────────────────────────────────────────────────

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

# ── Link update helpers ────────────────────────────────────────────────────────

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

# ── Private helpers ────────────────────────────────────────────────────────────

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

end # module Constellations
