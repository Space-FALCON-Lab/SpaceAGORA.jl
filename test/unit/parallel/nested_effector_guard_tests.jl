using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays

# The `:per_satellite_effector_reduce` route is the only auto route that keeps a
# threaded effector_decision. When the satellite axis is already threaded -- the
# RHS dispatch runs a Polyester `@batch` over the spacecraft whenever
# `_rhs_batch_parallel_enabled` and the batch is wider than one -- that decision
# nests a second split of the same pool inside every satellite of the first, on
# every RHS call, with a fresh per-effector slot buffer and a task spawn for
# each.
#
# `_satellite_batch_saturates_pool` does not catch it: it asks whether the
# satellite count reaches the thread budget, and a constellation with no outer
# split advertised sees the WHOLE pool as its budget precisely because nothing
# is splitting it. Eight spacecraft against a budget of twelve or thirty-two
# fall straight through to the nested route.
#
# Measured on 8 spacecraft, three effectors, 12 threads, same process,
# alternating: 341904 B per RHS call nested against 230464 B with the inner
# axis serial.

const SE = SpaceAGORA.SimulationEngine
const PPol = SpaceAGORA.SimulationModel.ParallelPolicy

const HARMONICS_FILE = joinpath(
    dirname(dirname(dirname(@__DIR__))), "data", "Gravity_harmonics_data", "EarthGGM05C.csv"
)

function guard_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra = planet.Rp_e + 540e3 + 2e3 * (id - 1),
        rp = planet.Rp_e + 500e3 + 1e3 * (id - 1),
        i = 35.0, ω = 40.0, Ω = 10.0, ν = 120.0 + 12.0 * (id - 1),
    )
    return SpacecraftModel(Joint[], Link[root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

function guard_config(n_sats::Int)
    planet = Earth()
    gravity = isfile(HARMONICS_FILE) ?
        GravitationalHarmonicsModel(8, 8, HARMONICS_FILE, planet) :
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
            density_model=ExponentialAtmosphereModel(planet),
            # No SPICE: the baseline no-GRAM path is what a fresh checkout has,
            # and the routing this file pins does not depend on the ephemeris.
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false
        ),
        dynamics_model=DynamicsModel(
            [guard_spacecraft(planet, i) for i in 1:n_sats],
            # Two effectors, both thread-safe and neither needing an ephemeris:
            # SolarRadiationPressureModel reaches SPICE for the Sun whatever the
            # configured ephemerides model is, and a fresh checkout has no
            # kernels.
            (gravity, AerodynamicCoefficientfM())
        ),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0
        ),
    )
end

# The effector policy refuses to thread light work (`heavy_only`, on by default)
# and the satellite batch only engages from `SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD`
# spacecraft up. Both are lifted here so the nested route is reachable on a
# small constellation -- the guard is what this file is about, not the gates
# that happen to keep the shipped defaults off that route today.
const NESTING_ENV = [
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
    "SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
    "SPACEAGORA_EFFECTOR_PARALLEL" => "on",
    "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "0",
    "SPACEAGORA_RHS_BATCH_PARALLEL" => "on",
    "SPACEAGORA_RHS_EXECUTION_MODE" => "auto",
]

const N_SATS = 8

function guard_plan(pairs)
    args = guard_config(N_SATS)
    p = SE.ODEParams(n_sats=N_SATS, args=args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_save_cache_buffers!(p)
    return withenv(pairs...) do
        plan = SE._rhs_execution_plan(args, p, args.dynamics_model.dynamic_effectors, N_SATS)
        (
            plan = plan,
            batch_parallel = SE._rhs_batch_parallel_enabled(p, N_SATS),
            batch_workers = SE._rhs_batch_workers(p),
        )
    end
end

@testset "a threaded satellite batch does not nest a threaded effector team" begin
    if Threads.nthreads() < 2
        @test_skip "needs julia --threads>=2 to exercise the nested route"
    else
        probe = guard_plan(NESTING_ENV)
        # Preconditions: this is the route the guard is about, and the satellite
        # axis really is threaded. If either stops holding the assertion below
        # would pass for the wrong reason.
        @test probe.plan.mode == :per_satellite_effector_reduce
        @test probe.batch_parallel
        @test probe.batch_workers > 1
        # The guard itself.
        @test probe.plan.effector_decision.use_threads == false
        @test probe.plan.effector_decision.allotment == 1
        @test probe.plan.allotment == 1
    end
end

@testset "the effector team survives where the satellite axis is serial" begin
    if Threads.nthreads() < 2
        @test_skip "needs julia --threads>=2 to exercise the effector team"
    else
        # Same shape, satellite batch off: the effector team is then the only
        # parallelism in the RHS and the guard must leave it alone. This is the
        # half that keeps the guard from degenerating into "never thread
        # effectors when no outer split is advertised".
        pairs = copy(NESTING_ENV)
        pairs[5] = "SPACEAGORA_RHS_BATCH_PARALLEL" => "off"
        probe = guard_plan(pairs)
        @test probe.plan.mode == :per_satellite_effector_reduce
        @test probe.batch_parallel == false
        @test probe.plan.effector_decision.use_threads == true
        @test probe.plan.effector_decision.allotment > 1
    end
end

@testset "per-RHS-call allocation does not grow with the thread budget" begin
    if Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4 to separate the budgets"
    else
        args = guard_config(N_SATS)
        u = SE.build_initial_conditions(args)
        p = SE.ODEParams(n_sats=N_SATS, args=args)
        SE._initialize_heat_rate_buffers!(p)
        SE._initialize_harmonics_workspace_buffers!(p)
        SE._initialize_save_cache_buffers!(p)
        du = copy(u)
        du .= 0.0

        function bytes_per_call(pairs; calls::Int = 40)
            return withenv(pairs...) do
                SE.spacecraft_dynamics!(du, u, p, 0.0)   # warm the route
                GC.gc()
                best = typemax(Int)
                for _ in 1:3
                    allocated = @allocated begin
                        for _ in 1:calls
                            SE.spacecraft_dynamics!(du, u, p, 0.0)
                        end
                    end
                    best = min(best, allocated ÷ calls)
                end
                return best
            end
        end

        narrow = bytes_per_call(vcat(NESTING_ENV, ["SPACEAGORA_INNER_THREAD_BUDGET" => "2"]))
        wide = bytes_per_call(NESTING_ENV)
        # Unguarded, the wide budget buys a wider effector team inside every
        # satellite of the batch and the per-call allocation climbs with it.
        # Guarded, both budgets take the same serial inner axis, so the two
        # differ only by whatever the rest of the step does.
        @test wide <= 1.25 * narrow
    end
end

@testset "the guard does not change results" begin
    args = guard_config(N_SATS)
    u = SE.build_initial_conditions(args)

    function rhs_under(pairs)
        p = SE.ODEParams(n_sats=N_SATS, args=args)
        SE._initialize_heat_rate_buffers!(p)
        SE._initialize_harmonics_workspace_buffers!(p)
        SE._initialize_save_cache_buffers!(p)
        du = copy(u)
        du .= 0.0
        withenv(pairs...) do
            SE.spacecraft_dynamics!(du, u, p, 0.0)
        end
        return du
    end

    serial = rhs_under([
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1",
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
    ])
    guarded = rhs_under(NESTING_ENV)
    @test guarded == serial
end
