using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays

# The thermal DiscreteCallback threads its per-satellite loop purely on
# satellite count and thread budget, the same shape the density callback had
# before it gained a light-work guard (see
# test/unit/parallel/density_callback_light_work_tests.jl). Its per-satellite
# body, `_compute_stage_heat_rates!` called with `use_buffered_density=true`,
# does one `sample_planet_frame` plus one `getHeatRate` per thermal link --
# light for a single-link spacecraft, and only becomes heavy enough to be
# worth a threaded dispatch once a spacecraft has enough links. Measured
# driving the affect body directly on an 8-spacecraft shape at thread widths
# 1/2/4/8: below THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD links the extra work
# does not reliably amortize the persistent-pool dispatch; at and above it,
# widths 4 and 8 are consistently faster than serial. See
# thermal_callbacks.jl's `_thermal_callback_work_is_heavy` for the constant
# and its provenance.

const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const SE = SpaceAGORA.SimulationEngine

const N_SATS = 8

function thermal_work_spacecraft(planet, id::Int, n_links::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    links = Link[root]
    for _ in 2:n_links
        push!(links, Link(root=false, m=10.0, ref_area=1.0))
    end
    ic = InitialCondition(
        ra = planet.Rp_e + 540e3 + 2e3 * (id - 1),
        rp = planet.Rp_e + 500e3 + 1e3 * (id - 1),
        i = 35.0, ω = 40.0, Ω = 10.0, ν = 120.0 + 12.0 * (id - 1),
    )
    return SpacecraftModel(Joint[], links, root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

# No SPICE and no harmonics file, matching the density light-work test: the
# baseline no-GRAM path is what a fresh checkout has.
function thermal_work_config(n_sats::Int, n_links::Int)
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
            [thermal_work_spacecraft(planet, i, n_links) for i in 1:n_sats],
            (InverseSquaredJ2GravityModel(),)
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

# A wide inner budget with the item-count gate and the auto budget floor
# lifted, so the *only* thing left that can keep the callback serial is the
# light-work guard this file is about.
const WIDE_ENV = [
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
    "SPACEAGORA_INNER_THREAD_BUDGET" => nothing,
    "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "auto",
    "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "1",
    "SPACEAGORA_THERMAL_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2",
]

const OFF_ENV = vcat(WIDE_ENV, ["SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "off"])

function thermal_work_params(args)
    p = SE.ODEParams(n_sats=N_SATS, args=args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_save_cache_buffers!(p)
    p.shared_buffers.et_start[] =
        SimulationModel.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    # Buffered-atmosphere path (matches update_thermal_sat!'s
    # use_buffered_density=true), pre-populated so the callback under test
    # never has to evaluate a density model.
    for i in 1:N_SATS
        p.shared_buffers.densities[i] = 1e-11
        p.shared_buffers.temperatures[i] = 250.0
        p.shared_buffers.winds[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
        p.shared_buffers.density_sample_t[i] = 0.0
    end
    return p
end

@testset "a single-link spacecraft is light work for the thermal callback" begin
    args_light = thermal_work_config(N_SATS, 1)
    p_light = thermal_work_params(args_light)
    @test CB._thermal_callback_work_is_heavy(p_light, N_SATS) == false

    args_heavy = thermal_work_config(N_SATS, CB.THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD)
    p_heavy = thermal_work_params(args_heavy)
    @test CB._thermal_callback_work_is_heavy(p_heavy, N_SATS) == true

    # nothing stands in for the no-run-state overload used only by direct
    # unit tests of _thermal_callback_thread_decision(num_sats); it must not
    # change that call path's behavior.
    @test CB._thermal_callback_work_is_heavy(nothing, N_SATS) == true
end

@testset "the thermal callback refuses a threaded dispatch for light work" begin
    if Threads.nthreads() < 2
        @test_skip "needs julia --threads>=2 to exercise the threaded dispatch"
    else
        light = CB._thermal_callback_thread_decision(N_SATS; heavy_work=false)
        # The guard: a single-link vehicle (pre-fix, this always threaded).
        @test light.use_threads == false
        @test light.allotment == 1
        # ... and it must not swallow the case the threaded path exists for.
        # If this half stops holding, the assertion above would pass for the
        # wrong reason. Computed inside WIDE_ENV: the default thread-count
        # gate and auto-budget floor would otherwise pin allotment at 1
        # regardless of heavy_work, and the assertion below would fail for a
        # reason unrelated to the guard.
        withenv(WIDE_ENV...) do
            heavy = CB._thermal_callback_thread_decision(N_SATS; heavy_work=true)
            @test heavy.use_threads == true
            @test heavy.allotment > 1
        end
    end
end

@testset "the thermal callback does not dispatch a single-link vehicle on a wide budget" begin
    args = thermal_work_config(N_SATS, 1)
    u = SE.build_initial_conditions(args)
    if Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4 to separate the budgets"
    else
        function bytes_per_affect(pairs; calls::Int=40)
            p = thermal_work_params(args)
            affect! = CB.get_thermal_callback(N_SATS, args).affect!
            integrator = (p=p, u=u, t=0.0)
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

        # Pre-fix, the wide budget sent eight workers through the persistent
        # pool on every call for a single-link vehicle, and the per-call
        # allocation climbed with it. Post-fix, both budgets run the same
        # serial loop.
        narrow = bytes_per_affect(vcat(WIDE_ENV, ["SPACEAGORA_INNER_THREAD_BUDGET" => "1"]))
        wide = bytes_per_affect(WIDE_ENV)
        @test wide <= 1.25 * max(narrow, 1)
    end
end

@testset "the thermal callback still threads a many-link vehicle on a wide budget" begin
    args = thermal_work_config(N_SATS, CB.THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD)
    u = SE.build_initial_conditions(args)
    if Threads.nthreads() < 4
        @test_skip "needs julia --threads>=4 to exercise the threaded dispatch"
    else
        p = thermal_work_params(args)
        affect! = CB.get_thermal_callback(N_SATS, args).affect!
        integrator = (p=p, u=u, t=0.0)
        withenv(WIDE_ENV...) do
            decision = CB._thermal_callback_thread_decision(
                p, N_SATS;
                heavy_work=CB._thermal_callback_work_is_heavy(p, N_SATS)
            )
            @test decision.use_threads == true
            affect!(integrator)  # must not throw, and must exercise the threaded path
        end
        @test all(isfinite, p.shared_buffers.heat_rates[1])
    end
end

@testset "the guard does not change the staged heat rates" begin
    args = thermal_work_config(N_SATS, CB.THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD)
    u = SE.build_initial_conditions(args)

    function staged_under(pairs)
        p = thermal_work_params(args)
        affect! = CB.get_thermal_callback(N_SATS, args).affect!
        withenv(pairs...) do
            affect!((p=p, u=u, t=0.0))
        end
        return [copy(hr) for hr in p.shared_buffers.heat_rates]
    end

    # `==` and not `isapprox`: both routes evaluate the same closed-form
    # thermal model at the same states in the same order, so any difference
    # at all would be a race on the staging buffer rather than a rounding
    # difference.
    @test staged_under(WIDE_ENV) == staged_under(OFF_ENV)
    @test all(hr -> all(isfinite, hr), staged_under(WIDE_ENV))
end
