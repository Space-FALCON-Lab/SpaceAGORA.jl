# The third-body sample's remaining per-call allocation (WS11f item 2).
#
# After WS10c's `map`-over-`body_names` fix (docs/architecture/third_body_cost.md,
# "Deliverable 1"), `sample_third_body_ephemerides` still allocated about
# 4.3 KB per call. Root cause: two of the `do`-block closure's captured
# variables were not hoisted to locals the way `et`/`primary_body_name`/
# `spice_rhs_memo`/etc. already were.
#
#   1. `p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls`
#      was read from inside the closure, which captures the whole
#      `p::ODEParams{...}` -- an inlinable immutable holding the entire
#      `SimulationConfiguration` -- to reach one `Threads.Atomic{Int64}`
#      three levels down. That forced the compiler to materialise a full
#      boxed copy of `p` once per spacecraft per derivative evaluation.
#   2. `cache_entry::Union{Nothing, NBodyEphemerisCache}` was captured
#      directly, so the closure's own type was a `Union` of two closures (one
#      per branch), and `map` over a `Union`-typed closure cannot infer a
#      concrete result type: `positions_ii` came back `Tuple{Any, Any}`,
#      boxing both `SVector{3, Float64}` positions on every call.
#
# The fix (src/simulation/engine/effector_sampling.jl): hoist the counter to
# a local before the closure, exactly like every other shared-buffer read
# already did, and take the `cache_entry isa NBodyEphemerisCache` branch
# OUTSIDE the `map`, once, so each of the two closures below it captures a
# concretely-typed `cache_entry` (or none at all) instead of the `Union`.
#
# Run: `julia --project=. --threads=1 benchmarks/studies/small_wins/third_body_sample_allocation.jl`
#
# Measured on space-falcon-1, 1 thread, back to back in one process (this
# script), 200000-call loops, warmed up first:
#
#   before (WS10c tip, adb343566)  4272.0 B/call
#   after (this change)              192.0 B/call   (22.3x fewer bytes)
#
# The "before" row above is a historical point (this script measures the
# current tree only); reproduce it by checking out
# src/simulation/engine/effector_sampling.jl at adb343566 and re-running.
# Values (positions, order, cache-then-backend fallback) are unchanged --
# see the SRP + third-body dump/cmp evidence in
# docs/architecture/small_wins_20260923.md.

using Test
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PE = SM.DynamicEffectors.PerturbationEffectors

function tbsa_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=planet.Rp_e + 550e3 + 1e3 * id,
        rp=planet.Rp_e + 500e3 + 1e3 * id,
        i=53.0, ω=24.0, Ω=10.0, ν=(360.0 * id) / 32,
    )
    return SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

function tbsa_configuration(planet, effectors::Tuple, n_sats::Int)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=10.0, orientation_sim=false, num_steps_to_save=10, data_rate=1.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel([tbsa_spacecraft(planet, i) for i in 1:n_sats], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0,
        ),
    )
end

function tbsa_nbody_cache(body_query_names::Vector{String})
    ets = collect(0.0:600.0:6000.0)
    positions = Matrix{SVector{3, Float64}}(undef, length(ets), length(body_query_names))
    for (i, et) in pairs(ets), k in eachindex(body_query_names)
        angle = 1e-6 * et + 0.7 * k
        radius = k == 1 ? 1.49e11 : 3.84e8
        positions[i, k] = SVector{3, Float64}(
            radius * cos(angle), radius * sin(angle), 0.1 * radius * sin(0.5 * angle)
        )
    end
    return SM.NBodyEphemerisCache(
        PE._spice_query_name("Earth"),
        body_query_names,
        Dict(name => k for (k, name) in pairs(body_query_names)),
        ets,
        positions,
    )
end

function build_scenario()
    planet = Earth()
    nbody = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)
    n_sats = 4
    args = tbsa_configuration(planet, (nbody,), n_sats)
    p = SM.ODEParams(n_sats=n_sats, args=args)
    p.shared_buffers.nbody_ephemeris_cache[] = tbsa_nbody_cache(
        [PE._spice_query_name("Sun"), PE._spice_query_name("Moon")]
    )
    p.shared_buffers.et_start[] = 0.0
    u = SE.build_initial_conditions(args)
    return nbody, p, u.sc[1]
end

function loop_sample(nbody, p, sc1, n::Int)
    last = SE.sample_third_body_ephemerides(nbody, sc1, p, 1, 2750.0)
    for _ in 2:n
        last = SE.sample_third_body_ephemerides(nbody, sc1, p, 1, 2750.0)
    end
    return last
end

nbody, p, sc1 = build_scenario()
n = 200_000

loop_sample(nbody, p, sc1, 10)   # warm up compilation
bytes_total = @allocated loop_sample(nbody, p, sc1, n)
bytes_per_call = bytes_total / n
println("sample_third_body_ephemerides: ", bytes_per_call, " bytes/call over ", n, " calls")

@testset "sample_third_body_ephemerides allocates a small, bounded amount per call" begin
    # 4272 B/call was the measured pre-fix (WS10c tip) cost; the fix must stay
    # well under it. 512 B/call leaves headroom for compiler/version drift
    # while still failing hard if the closure regresses back toward capturing
    # `p` or a `Union`-typed `cache_entry`.
    @test bytes_per_call < 512.0
end
