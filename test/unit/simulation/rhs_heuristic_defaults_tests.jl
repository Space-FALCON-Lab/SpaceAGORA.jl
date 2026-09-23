module RhsHeuristicDefaultsTests

# The DEFAULT RHS execution plan for a single-harmonics constellation.
#
# What this pins is the routing width floor measured in
# `docs/architecture/rhs_heuristic_defaults.md`: the heuristic used to size the
# flat route's worker count from `fld(active_sats, 4)` and took it whenever that
# was two or more, where the 4 is the
# harmonics pre-pass's SIMD slice floor and has nothing to do with what a worker
# costs to wake. On the archived TRX50 P1 ladder that default ran 10.6x, 7.7x
# and 2.5x slower than the plan the calibration sweep measured at 64, 256 and
# 1024 spacecraft, and every user who does not enable an adaptive profile paid
# it.
#
# The fix gates the route with a satellites-per-worker floor at the full
# budget; below it the default returns the satellite batch.
#
# The tests below are about the ROUTING DECISION, not about timing: they ask
# which plan `_rhs_execution_plan` returns for a given (active satellites,
# thread budget), which is deterministic and cheap to check.

using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SE_RHD = SpaceAGORA.SimulationEngine

const RHD_HARMONICS_FILE = joinpath(
    dirname(dirname(dirname(@__DIR__))), "data", "Gravity_harmonics_data", "EarthGGM05C.csv"
)

function rhd_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra = planet.Rp_e + 540e3 + 2e3 * (id - 1),
        rp = planet.Rp_e + 500e3 + 1e3 * (id - 1),
        i = 35.0, ω = 40.0, Ω = 10.0, ν = 120.0 + 12.0 * (id - 1),
    )
    return SpacecraftModel(Joint[], Link[root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

# A single gravitational-harmonics effector is what makes this the
# `single_harmonics_flat` shape. Degree 8 rather than 50: the route decision
# does not read the degree, and the file parse is the expensive part.
function rhd_config(n_sats::Int; aero::Bool=false)
    planet = Earth()
    gravity = isfile(RHD_HARMONICS_FILE) ?
        GravitationalHarmonicsModel(8, 8, RHD_HARMONICS_FILE, planet) :
        InverseSquaredJ2GravityModel()
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=600.0, orientation_sim=false, num_steps_to_save=10, data_rate=10.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0,
            density_model=aero ? ExponentialAtmosphereModel(planet) : NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false
        ),
        dynamics_model=DynamicsModel([rhd_spacecraft(planet, i) for i in 1:n_sats],
            aero ? (gravity, AerodynamicCoefficientfM()) : (gravity,)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0
        ),
    )
end

# The route decision as a solve would take it, with everything the harness or a
# leftover profile could otherwise inject cleared.
rhd_env(budget::Int, extra...) = vcat(Pair{String, Union{Nothing, String}}[
    "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
    "SPACEAGORA_RHS_EXECUTION_MODE" => nothing,
    "SPACEAGORA_PARALLEL_PROFILE" => nothing,
    "SPACEAGORA_RHS_PLAN_STEP_CACHE" => nothing,
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => nothing,
    "SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER" => nothing,
    "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => nothing,
    "SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER" => nothing,
    "SPACEAGORA_EFFECTOR_FLAT_MIN_SATS" => nothing,
    "SPACEAGORA_RHS_BATCH_PARALLEL" => nothing,
], collect(Pair{String, Union{Nothing, String}}, extra))

function rhd_plan(n_sats::Int, budget::Int, extra...; aero::Bool=false)
    args = rhd_config(n_sats; aero=aero)
    p = SE_RHD.ODEParams(n_sats=n_sats, args=args)
    SE_RHD._initialize_heat_rate_buffers!(p)
    SE_RHD._initialize_harmonics_workspace_buffers!(p)
    SE_RHD._initialize_save_cache_buffers!(p)
    return withenv(rhd_env(budget, extra...)...) do
        SE_RHD._rhs_execution_plan(args, p, args.dynamics_model.dynamic_effectors, n_sats)
    end
end

@testset "single-harmonics default routing width floor" begin
    floor_per_worker = withenv("SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER" => nothing) do
        SE_RHD._rhs_harmonics_flat_min_sats_per_worker()
    end
    # The measured default, not merely "some positive number": a regression that
    # put it back at the pre-pass slice floor of 4 is exactly what this file
    # exists to catch.
    @test floor_per_worker == 64

    # SPACEAGORA_INNER_THREAD_BUDGET is capped at the pool, so the budgets
    # below need that many threads to mean what they say.
    if Threads.nthreads() < 2
        @test_skip "needs julia --threads>=2 to reach a multi-worker route"
        return
    end
    wide = min(8, Threads.nthreads())

    # Below two workers' worth of satellites the default must not open a flat
    # worker team; the satellite batch is the route that costs nothing to
    # dispatch.
    plan_small = rhd_plan(64, wide)
    @test plan_small.mode == :satellite_batch
    @test plan_small.allotment == 1
    @test plan_small.effector_decision.use_threads == false

    # At the floor the flat route is admitted again, at the full width it always
    # had: the floor gates the route, not the width.
    # A budget of 2 keeps the constellation small: the route opens at
    # floor * budget satellites, one below that it does not.
    plan_at = rhd_plan(2 * floor_per_worker, 2)
    @test plan_at.mode == :flat_constellation_effector_queue
    @test plan_at.allotment == 2
    plan_below = rhd_plan(2 * floor_per_worker - 1, 2)
    @test plan_below.mode == :satellite_batch

    # The floor is an override, and a floor of 1 restores the pre-change
    # decision exactly (the width is still the pre-pass slice floor's).
    plan_legacy = rhd_plan(64, wide, "SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER" => "1")
    @test plan_legacy.mode == :flat_constellation_effector_queue
    @test plan_legacy.allotment == min(wide, fld(64, 4))
end

@testset "the width floor leaves the one-thread flat admission alone" begin
    # At a budget of one the flat route spawns no tasks at all -- the pre-pass
    # runs its slice inline -- so it is admitted for the batched coefficient
    # sweep, not for threading, and the width floor has no business gating it.
    # This is the admission `third_body_route_parity_tests.jl` pins for
    # multi-effector pre-pass stacks; the single-harmonics branch has its own.
    plan = rhd_plan(64, 1)
    @test plan.mode == :flat_constellation_effector_queue
    @test plan.allotment == 1
end

@testset "an enclosing outer split still collapses the width to one" begin
    plan = rhd_plan(64, max(2, min(8, Threads.nthreads())), "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1")
    @test plan.mode == :flat_constellation_effector_queue
    @test plan.allotment == 1
end

@testset "harmonics plus aerodynamics takes the same floor" begin
    # The generic flat branch, reached by a two-effector stack with no batched
    # pre-pass kernel. Measured 2.89x behind the satellite batch at 64
    # spacecraft on 8 threads before the floor.
    if Threads.nthreads() < 4
        # Below 4 the flat queue's own thread-budget floor routes to the batch
        # before the branch under test is reached.
        @test_skip "needs julia --threads>=4"
    else
        wide = min(8, Threads.nthreads())
        plan = rhd_plan(64, wide; aero=true)
        @test plan.mode == :satellite_batch
        plan_legacy = rhd_plan(64, wide, "SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER" => "1"; aero=true)
        @test plan_legacy.mode == :flat_constellation_effector_queue
    end
end

end # module
