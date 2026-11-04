module AeroBatchParityTests

# Guards for the aerodynamic and density batching measured in
# docs/architecture/aero_batch.md.
#
# The claim under test is not "close enough". It is that a constellation's
# aerodynamic force and its atmosphere sample are the same Float64 values
# whether they are produced one spacecraft at a time on the per-satellite RHS
# route or in a batch pre-pass on the flat constellation route, and that the
# number of workers the pre-pass is split across cannot move a bit. So every
# comparison below is `===` on the raw Float64 (or on its bit pattern), never
# `≈` and never a tolerance.
#
# THREADS. At an inner thread budget of one the flat constellation queue is
# admitted only when every effector is served by a pre-pass; the aerodynamic
# pre-pass counts, so the (gravity, aero) stack below takes the flat route in
# any process and the route comparison always runs. The worker-count sweep
# needs more than one thread and is skipped, visibly, on a single-threaded
# process. Run this file with
# `julia --project=. --threads=4 test/unit/dynamics/aero_batch_parity_tests.jl`
# to exercise every testset.

using Test
using StaticArrays
using ComponentArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM_ABP = SpaceAGORA.SimulationModel
const SE_ABP = SpaceAGORA.SimulationEngine
const EM_ABP = SM_ABP.EnvironmentModels

# ── Fixtures ─────────────────────────────────────────────────────────────────

# Perigee altitudes are spread low enough that the exponential model returns a
# nonzero density for most spacecraft and high enough for the last few that it
# underflows -- both the "there is drag" and the "there is no drag" branch of
# `_aero_pure_wrench` are exercised in one constellation.
function abp_spacecraft(planet, id::Int, n::Int; with_panel::Bool=false)
    root = Link(root=true, m=500.0, ref_area=12.0)
    links = Link[root]
    if with_panel
        for panel_idx in 1:2
            theta = pi * (panel_idx - 1)
            push!(links, Link(
                root=false, m=8.0, ref_area=3.0,
                r=MVector{3, Float64}(1.8 * cos(theta), 1.8 * sin(theta), 0.4),
            ))
        end
    end
    rp_alt = 180e3 + (700e3 - 180e3) * (id - 1) / max(1, n - 1)
    ic = InitialCondition(
        ra=planet.Rp_e + rp_alt + 40e3,
        rp=planet.Rp_e + rp_alt,
        i=53.0, ω=24.0, Ω=10.0, ν=(360.0 * (id - 1)) / n,
    )
    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(Joint[], links, root, true, dry_mass, 0.0, root.inertia, 0, 0, ic, id)
end

function abp_configuration(
    planet,
    n_sats::Int;
    incidence::Symbol=:max_drag,
    with_panel::Bool=false,
    density_model=EM_ABP.ExponentialAtmosphereModel(planet),
    per_link_atmosphere::Bool=false,
)
    effectors = (
        InverseSquaredJ2GravityModel(),
        AerodynamicCoefficientfM(
            fixed_attitude_incidence=incidence, per_link_atmosphere=per_link_atmosphere,
        ),
    )
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=10.0, orientation_sim=false, num_steps_to_save=10, data_rate=1.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            # No SPICE: a fresh checkout has no kernels, and neither the aero
            # wrench nor the exponential atmosphere needs an ephemeris.
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel(
            [abp_spacecraft(planet, i, n_sats; with_panel=with_panel) for i in 1:n_sats],
            effectors,
        ),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0
        ),
    )
end

# The env a solve would present, with everything that could make two calls
# differ for a reason other than the route pinned off.
function abp_env(mode::String, budget::Int)
    return [
        "SPACEAGORA_RHS_EXECUTION_MODE" => mode,
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
        "SPACEAGORA_PARALLEL_PROFILE" => nothing,
        "SPACEAGORA_RHS_PLAN_STEP_CACHE" => nothing,
        "SPACEAGORA_EFFECTOR_FLAT_MIN_SATS" => "8",
        "SPACEAGORA_EFFECTOR_FLAT_MIN_EFFECTORS" => "1",
        "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => nothing,
        "SPACEAGORA_VACUUM_GRAM_CACHE" => nothing,
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
    ]
end

# One derivative evaluation of `args` at time `t`, on the route `mode` names.
# Returns the raw derivative data and the route actually taken, so a testset
# can assert it got the route it asked for instead of silently comparing a
# route with itself.
function abp_derivative(args, mode::String, budget::Int, t::Float64)
    n_sats = length(args.dynamics_model.spacecraft)
    u = SE_ABP.build_initial_conditions(args)
    du = zero(u)
    p = ODEParams(n_sats=n_sats, args=args)
    SE_ABP._initialize_save_cache_buffers!(p)
    SE_ABP._initialize_heat_rate_buffers!(p)
    taken = withenv(abp_env(mode, budget)...) do
        plan = SE_ABP._rhs_execution_plan(args, p, args.dynamics_model.dynamic_effectors, n_sats)
        SE_ABP.spacecraft_dynamics!(du, u, p, t)
        plan.mode
    end
    return (
        data = copy(ComponentArrays.getdata(du)),
        mode = taken,
        drag = copy(p.save_cache.drag_cache),
        lift = copy(p.save_cache.lift_cache),
        cross = copy(p.save_cache.cross_cache),
    )
end

# Bit-for-bit, including the sign of a zero and the payload of a NaN: `===` on
# Float64 already distinguishes -0.0 from 0.0, and reinterpreting makes a NaN
# comparable too, which `===` on NaN would not be.
abp_bits(x::Float64)::UInt64 = reinterpret(UInt64, x)

function abp_count_mismatches(a::AbstractVector{Float64}, b::AbstractVector{Float64})::Int
    length(a) == length(b) || return typemax(Int)
    n = 0
    @inbounds for i in eachindex(a)
        abp_bits(a[i]) === abp_bits(b[i]) || (n += 1)
    end
    return n
end

function abp_count_vec_mismatches(
    a::AbstractVector{SVector{3, Float64}},
    b::AbstractVector{SVector{3, Float64}},
)::Int
    length(a) == length(b) || return typemax(Int)
    n = 0
    @inbounds for i in eachindex(a), k in 1:3
        abp_bits(a[i][k]) === abp_bits(b[i][k]) || (n += 1)
    end
    return n
end

const ABP_MULTITHREADED = Threads.nthreads() > 1

# ── Testsets ─────────────────────────────────────────────────────────────────

@testset "aerodynamic and density batching" begin

@testset "the flat route reproduces the per-satellite route bit for bit" begin
    planet = Earth()
    # Budget 1 is the serial admission of the flat route, which only the auto
    # route makes (an explicit flat request at budget 1 still falls back to the
    # per-satellite batch); the wider budget is the threaded one, forced. Both
    # must reproduce the per-satellite route.
    budgets = unique((1, min(4, Threads.nthreads())))
    for incidence in (:max_drag, :attitude, :tumbling_average),
        with_panel in (false, true),
        n_sats in (9, 32),
        budget in budgets

        args = abp_configuration(planet, n_sats; incidence=incidence, with_panel=with_panel)
        serial = abp_derivative(args, "satellite", 1, 0.0)
        flat = abp_derivative(args, budget == 1 ? "auto" : "flat", budget, 0.0)
        @test serial.mode === :satellite_batch
        @test flat.mode === :flat_constellation_effector_queue
        @test abp_count_mismatches(serial.data, flat.data) == 0
        @test abp_count_vec_mismatches(serial.drag, flat.drag) == 0
        @test abp_count_vec_mismatches(serial.lift, flat.lift) == 0
        @test abp_count_vec_mismatches(serial.cross, flat.cross) == 0
    end

    # The auto route at budget 1 admits the flat queue for this stack, and a
    # per-link-atmosphere model, which stays on the queue, keeps it out.
    auto = abp_derivative(abp_configuration(planet, 32), "auto", 1, 0.0)
    @test auto.mode === :flat_constellation_effector_queue
    per_link = abp_derivative(abp_configuration(planet, 32; per_link_atmosphere=true), "auto", 1, 0.0)
    @test per_link.mode === :satellite_batch
end

@testset "the flat route's answer does not depend on the worker count" begin
    planet = Earth()
    if !ABP_MULTITHREADED
        @test Threads.nthreads() == 1
        @info "aero_batch_parity: single-threaded process, worker-count sweep not run"
    else
        # Satellite counts that do not divide evenly into worker slices, so a
        # short last slice is exercised as well as a full one.
        for n_sats in (9, 17, 32, 33)
            args = abp_configuration(planet, n_sats)
            reference = abp_derivative(args, "flat", min(2, Threads.nthreads()), 0.0)
            @test reference.mode === :flat_constellation_effector_queue
            for budget in (2, 3, 4, 8)
                budget <= Threads.nthreads() || continue
                other = abp_derivative(args, "flat", budget, 0.0)
                @test other.mode === :flat_constellation_effector_queue
                @test abp_count_mismatches(reference.data, other.data) == 0
                @test abp_count_vec_mismatches(reference.drag, other.drag) == 0
            end
        end
    end
end

@testset "a spacecraft the atmosphere does not reach contributes exactly zero" begin
    planet = Earth()
    # NoAtmosphereModel returns rho = 0, so `_aero_pure_wrench`'s vacuum guard
    # fires for every spacecraft and the aero effector must contribute nothing
    # at all -- not a denormal, not a signed zero of the wrong sign.
    args = abp_configuration(planet, 32; density_model=NoAtmosphereModel())
    r = abp_derivative(args, "satellite", 1, 0.0)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    @test all(v === zero3 for v in r.drag[1:32])
    @test all(v === zero3 for v in r.lift[1:32])
    @test all(v === zero3 for v in r.cross[1:32])

    if ABP_MULTITHREADED
        flat = abp_derivative(args, "flat", min(4, Threads.nthreads()), 0.0)
        @test flat.mode === :flat_constellation_effector_queue
        @test abp_count_mismatches(r.data, flat.data) == 0
        @test all(v === zero3 for v in flat.drag[1:32])
    end
end

@testset "getDensityBatch! reproduces getDensity bit for bit" begin
    planet = Earth()
    args = abp_configuration(planet, 4)
    p = ODEParams(n_sats=4, args=args)

    # The uniform-light route of the RHS atmosphere pre-sample answers N
    # queries with one getDensityBatch! call instead of N getDensity calls.
    # That substitution is only legitimate if the two are the same arithmetic,
    # which is asserted here rather than read off the source.
    exponential = EM_ABP.ExponentialAtmosphereModel(planet)
    piecewise = EM_ABP.PiecewiseExponentialAtmosphereModel(
        [0.0, 100e3, 300e3, 1000e3],
        [1.2, 5.6e-7, 2.4e-11],
        [8.5e3, 22e3, 60e3],
    )
    # Altitudes that straddle the piecewise model's breakpoints exactly, plus a
    # sweep across and beyond its advisory band.
    alts = Float64[
        -5e3, 0.0, 1.0, 99_999.0, 100e3, 100_001.0, 180e3, 299_999.0, 300e3,
        300_001.0, 550e3, 999_999.0, 1000e3, 1_000_001.0, 2000e3, 1e7,
    ]
    lats = [0.3 * sin(0.7 * i) for i in eachindex(alts)]
    lons = [0.9 * cos(0.4 * i) for i in eachindex(alts)]

    for model in (exponential, piecewise, NoAtmosphereModel())
        n = length(alts)
        rhos = zeros(Float64, n)
        Ts = zeros(Float64, n)
        winds = fill(SVector{3, Float64}(0.0, 0.0, 0.0), n)
        getDensityBatch!(rhos, Ts, winds, model, alts, lats, lons, 0.0, true, p)
        mismatches = 0
        for i in 1:n
            rho_i, T_i, wind_i = getDensity(model, alts[i], lats[i], lons[i], 0.0, true, p)
            abp_bits(rhos[i]) === abp_bits(rho_i) || (mismatches += 1)
            abp_bits(Ts[i]) === abp_bits(T_i) || (mismatches += 1)
            for k in 1:3
                abp_bits(winds[i][k]) === abp_bits(wind_i[k]) || (mismatches += 1)
            end
        end
        @test mismatches == 0
    end
end

@testset "the flat route's atmosphere buffers match the per-satellite samples" begin
    planet = Earth()
    n_sats = 32
    args = abp_configuration(planet, n_sats)
    p = ODEParams(n_sats=n_sats, args=args)
    SE_ABP._initialize_save_cache_buffers!(p)
    SE_ABP._initialize_heat_rate_buffers!(p)
    u = SE_ABP.build_initial_conditions(args)

    # What the RHS pre-sample writes into shared_buffers, against the
    # per-satellite sampler called directly on the same state. These are the
    # values the aero effector then reads, so a difference here would be a
    # difference in the force even if the wrench code were untouched.
    withenv(abp_env("flat", max(1, min(4, Threads.nthreads())))...) do
        SE_ABP._prefill_atmosphere_samples!(p, 0.0, u.sc)
    end
    batched_rho = copy(p.shared_buffers.densities[1:n_sats])
    batched_T = copy(p.shared_buffers.temperatures[1:n_sats])
    batched_wind = copy(p.shared_buffers.winds[1:n_sats])

    mismatches = 0
    for i in 1:n_sats
        sample = withenv(abp_env("satellite", 1)...) do
            SE_ABP.sample_atmosphere(u.sc[i], p, i, 0.0; write_buffers=false)
        end
        abp_bits(batched_rho[i]) === abp_bits(sample.rho_kg_m3) || (mismatches += 1)
        abp_bits(batched_T[i]) === abp_bits(sample.temperature_k) || (mismatches += 1)
        for k in 1:3
            abp_bits(batched_wind[i][k]) === abp_bits(sample.wind_pp[k]) || (mismatches += 1)
        end
    end
    @test mismatches == 0
end

end # testset "aerodynamic and density batching"

end # module
