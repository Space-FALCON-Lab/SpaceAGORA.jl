using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays

# The density DiscreteCallback used to thread its per-satellite loop on an item
# count alone. It asked `thread_policy_decision` "are there at least
# SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD satellites and is there a thread
# budget", and never passed the light-work guard (`heavy_only`) that the
# effector policy has always used -- so a constellation on a closed-form
# atmosphere dispatched eight workers, once per accepted step, for a loop body
# that is a 3x3 rotation, a lat/lon conversion and one `exp`.
#
# Worse, on the batch route -- which is the route a constellation sharing one
# density model takes -- the threaded region does not evaluate density at all.
# It stages altitude, latitude and longitude into
# shared_buffers.density_batch_*, and `getDensityBatch!` then runs on one
# thread afterwards. There is nothing there for a worker pool to win.
#
# Measured on 8 spacecraft, degree-50 harmonics, exponential atmosphere, a
# one-hour mission at 24 threads with no outer split: 3.64 s with
# SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto against 0.64 s with it off, and
# bit-identical final states. After the guard, auto measures 0.59 s -- the
# `off` number.

const SE = SpaceAGORA.SimulationEngine
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks

const N_SATS = 8

function light_work_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra = planet.Rp_e + 540e3 + 2e3 * (id - 1),
        rp = planet.Rp_e + 500e3 + 1e3 * (id - 1),
        i = 35.0, ω = 40.0, Ω = 10.0, ν = 120.0 + 12.0 * (id - 1),
    )
    return SpacecraftModel(Joint[], Link[root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

# No SPICE and no harmonics file: the routing this file pins does not depend on
# either, and the baseline no-GRAM path is what a fresh checkout has.
function light_work_config(n_sats::Int)
    planet = Earth()
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
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false
        ),
        dynamics_model=DynamicsModel(
            [light_work_spacecraft(planet, i) for i in 1:n_sats],
            (InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM())
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

# A wide inner budget with the item-count gate and the auto budget floor lifted,
# so the *only* thing left that can keep the callback serial is the light-work
# guard this file is about.
const WIDE_ENV = [
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
    "SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
    "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
    "SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2",
    "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET" => "2",
]

const OFF_ENV = vcat(WIDE_ENV, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off"])

function light_work_params(args)
    p = SE.ODEParams(n_sats=N_SATS, args=args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_save_cache_buffers!(p)
    p.shared_buffers.et_start[] =
        SimulationModel.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    return p
end

# The callback's affect! only reads `p`, `u` and `t` off its integrator, so a
# NamedTuple standing in for one keeps this a callback test rather than a solve.
light_work_integrator(p, u, t::Float64) = (p=p, u=u, t=t)

@testset "a closed-form atmosphere is light work for the density callback" begin
    planet = Earth()
    @test CB.density_model_work_is_heavy(ExponentialAtmosphereModel(planet)) == false
    @test CB.density_model_work_is_heavy(NoAtmosphereModel()) == false
    @test CB.density_model_work_is_heavy(
        PiecewiseExponentialAtmosphereModel([0.0, 1.0e5], [1.2], [8.5e3])
    ) == false
end

@testset "the density callback refuses a threaded dispatch for light work" begin
    args = light_work_config(N_SATS)
    if Threads.nthreads() < 2
        @test_skip "needs julia --threads>=2 to exercise the threaded dispatch"
    else
        withenv(WIDE_ENV...) do
            light = CB._density_callback_thread_decision(args, N_SATS; heavy_work=false)
            heavy = CB._density_callback_thread_decision(args, N_SATS; heavy_work=true)
            # The guard.
            @test light.use_threads == false
            @test light.allotment == 1
            # ... and it must not swallow the case the threaded path exists for.
            # If this half stops holding, the assertion above would pass for the
            # wrong reason.
            @test heavy.use_threads == true
            @test heavy.allotment > 1
        end
        # An explicit `on` is an override, not a request the guard may decline.
        withenv(vcat(WIDE_ENV, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on"])...) do
            forced = CB._density_callback_thread_decision(args, N_SATS; heavy_work=false)
            @test forced.use_threads == true
        end
    end
end

@testset "the staged density callback does not dispatch on a wide budget" begin
    args = light_work_config(N_SATS)
    u = SE.build_initial_conditions(args)
    if Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4 to separate the budgets"
    else
        function bytes_per_affect(pairs; calls::Int = 40)
            p = light_work_params(args)
            affect! = CB.get_density_callback(N_SATS, args.dynamics_model.dynamic_effectors, args).affect!
            integrator = light_work_integrator(p, u, 0.0)
            return withenv(pairs...) do
                affect!(integrator)   # warm the route
                GC.gc()
                best = typemax(Int)
                for _ in 1:3
                    allocated = @allocated begin
                        for _ in 1:calls
                            affect!(integrator)
                        end
                    end
                    best = min(best, allocated ÷ calls)
                end
                return best
            end
        end

        # Unguarded, the wide budget sends eight workers through the persistent
        # pool on every call and the per-call allocation climbs with it.
        # Guarded, both budgets run the same serial loop.
        narrow = bytes_per_affect(vcat(WIDE_ENV, ["SPACEAGORA_INNER_THREAD_BUDGET" => "1"]))
        wide = bytes_per_affect(WIDE_ENV)
        @test wide <= 1.25 * max(narrow, 1)
    end
end

@testset "the guard does not change the staged density samples" begin
    args = light_work_config(N_SATS)
    u = SE.build_initial_conditions(args)

    function staged_under(pairs)
        p = light_work_params(args)
        affect! = CB.get_density_callback(N_SATS, args.dynamics_model.dynamic_effectors, args).affect!
        withenv(pairs...) do
            affect!(light_work_integrator(p, u, 0.0))
        end
        return (
            copy(p.shared_buffers.densities),
            copy(p.shared_buffers.temperatures),
            copy(p.shared_buffers.winds),
        )
    end

    # `==` and not `isapprox`: both routes evaluate the same closed-form
    # atmosphere at the same states in the same order, so any difference at all
    # would be a race on the staging buffers rather than a rounding difference.
    @test staged_under(WIDE_ENV) == staged_under(OFF_ENV)
    @test all(isfinite, staged_under(WIDE_ENV)[1])
end
