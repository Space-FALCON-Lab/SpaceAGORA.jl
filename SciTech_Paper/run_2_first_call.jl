using Dates
using LinearAlgebra
using TOML

isempty(ARGS) || error("Usage: run_2_first_call.jl (fixed 50-helper, 100-orbit SMA case)")
const PAPER_DIR = @__DIR__
const SOURCE_ROOT = normpath(joinpath(PAPER_DIR, "..", "2_SpaceAGORA.jl"))
const LOAD_WALL_S = @elapsed include(joinpath(SOURCE_ROOT, "examples", "common.jl"))
using DiffEqBase

normpath(pathof(SpaceAGORA)) == joinpath(SOURCE_ROOT, "src", "SpaceAGORA.jl") ||
    error("Loaded the wrong SpaceAGORA snapshot: $(pathof(SpaceAGORA))")

const CASE_DIR = "50_helpers_100_orbits_gve_sma_target_inc_1.0_deg"
const BENCH_ROOT = joinpath(PAPER_DIR, "runtime_benchmarks")

function run_case()
    planet = make_no_gram_planet(:earth)
    planet.μ == 3.98600436233e14 || error("Gravitational parameter differs from the comparison.")
    target_radius = 6_378_136.6 + 1050e3
    helper_radius = 6_378_136.6 + 1000e3
    target_period = 2pi * sqrt(target_radius^3 / planet.μ)
    duration = 100 * target_period
    spacecraft = map(1:51) do index
        bus = Link(root=true, m=227.0, ref_area=1.0)
        ic = InitialCondition(index == 1 ? target_radius : helper_radius,
            0.0, index == 1 ? 1.0 : 0.0, 0.0, 0.0,
            index == 1 ? 0.0 : 360.0 * (index - 2) / 50)
        SpacecraftModel(Joint[], [bus], bus, true, 227.0, 0.0,
            bus.inertia, 0, 0, ic, index)
    end
    laser = OpenCavityLaserLinkModel(target_idx=1, helper_indices=collect(2:51),
        range_m=200e3, power_w=10_000.0, magnification=100.0,
        beta=1.0, eta=1.0, schedule=:gve_sma)
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false,
            verbose=false, results_directory=joinpath(BENCH_ROOT, "first_call", "runs", CASE_DIR, "2_spaceagora")),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(), laser)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9,
            abstol_orbit=1e-11, dt_max_orbit=10.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))

    tracker = LaserImpulseTracker()
    impulse_callback = laser_impulse_callback(laser, tracker, 227.0)
    scheduler_callback = laser_link_scheduler_callback(laser)
    schedule = Tuple{Float64, Int}[]
    function record!(integrator)
        push!(schedule, (Float64(integrator.t), laser.active_helper_idx))
        return nothing
    end
    recorder = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, record!;
        initialize=(callback, u, t, integrator) -> record!(integrator),
        save_positions=(false, false))
    # Preserve the native example's kick-before-rescheduling order and model identity.
    elapsed = @elapsed solution = run_simulation(args; isolate_state=false,
        return_solution=true, extra_callbacks=(impulse_callback, scheduler_callback, recorder))
    string(solution.retcode) == "Success" || error("Simulation failed: $(solution.retcode)")
    solution.t[end] == duration || error("Simulation did not reach the requested duration.")
    solution.prob.p.args.dynamics_model.dynamic_effectors[2] === laser ||
        error("Callbacks did not share the simulated laser model.")
    first(schedule)[1] == 0.0 && last(schedule)[1] == duration ||
        error("Recorded schedule does not span the simulation.")
    issorted(first.(schedule)) || error("Schedule times are not monotonic.")
    all(entry -> entry[2] == 0 || 2 <= entry[2] <= 51, schedule) ||
        error("Schedule contains an invalid helper index.")
    count(entry -> entry[2] > 0, schedule) == laser.active_link_step_count ||
        error("Recorded schedule does not match the native scheduler count.")
    tracker.t_prev == duration || error("Impulse callback did not reach the final time.")
    inclination(state) = let h = cross(state.pos, state.vel)
        rad2deg(atan(hypot(h[1], h[2]), h[3]))
    end
    initial_inc = inclination(first(solution.u).sc[1])
    isapprox(initial_inc, 1.0; atol=1e-10) || error("Incorrect initial inclination.")
    state = last(solution.u).sc[1]
    final_sma = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / planet.μ)
    final_inc = inclination(state)
    final_ecc = norm(cross(state.vel, cross(state.pos, state.vel)) / planet.μ -
        state.pos / norm(state.pos))
    all(isfinite, (elapsed, final_sma, final_inc, final_ecc)) && elapsed > 0 ||
        error("Nonfinite simulation output or invalid timing.")
    result = (
        code="2_SpaceAGORA.jl", helper_count=50, schedule="gve_sma",
        requested_target_orbits=100, target_period_s=target_period, duration_s=duration,
        initial_target_altitude_m=1050e3, initial_helper_altitude_m=1000e3,
        initial_helper_inclination_deg=0.0, initial_inclination_deg=initial_inc,
        final_inclination_deg=final_inc, delta_inclination_deg=final_inc - initial_inc,
        initial_sma_m=target_radius, final_sma_m=final_sma,
        delta_sma_m=final_sma - target_radius, final_eccentricity=final_ecc,
        scheduler_updates=length(schedule), link_activations=laser.link_activation_count,
        active_updates=laser.active_link_step_count, final_helper=laser.active_helper_idx,
        runtime_s=elapsed, retcode=string(solution.retcode))
    result_dir = joinpath(BENCH_ROOT, "first_call", "results", CASE_DIR, "2_spaceagora")
    mkpath(result_dir)
    open(joinpath(result_dir, "summary.toml"), "w") do io
        TOML.print(io, Dict(string(k) => v for (k, v) in pairs(result)); sorted=true)
    end
    open(joinpath(result_dir, "schedule.csv"), "w") do io
        println(io, "time_s,target_orbits,active_helper")
        for (time, helper) in schedule
            println(io, time, ",", time / target_period, ",", helper)
        end
    end
    println(result)
    return result
end

function benchmark_first_call()
    GC.gc()
    wrapper_wall = @elapsed result = run_case()
    wrapper_wall >= result.runtime_s || error("Invalid nested timer measurements.")
    summary = Dict(
        "code" => "2_spaceagora", "helper_count" => 50, "orbits" => 100,
        "schedule" => "gve_sma", "target_inclination_deg" => 1.0,
        "trial_count" => 1, "first_call_only" => true,
        "completed_at" => string(now()), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "cpu_model" => first(Sys.cpu_info()).model,
        "runner_load_wall_s" => LOAD_WALL_S,
        "first_simulation_call_wall_s" => result.runtime_s,
        "first_wrapper_wall_s" => wrapper_wall,
        "trials" => [Dict(string(k) => v for (k, v) in pairs(result))])
    mkpath(joinpath(BENCH_ROOT, CASE_DIR))
    open(joinpath(BENCH_ROOT, CASE_DIR, "2_spaceagora.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    println("FIRST CALL: ", result.runtime_s, " s; wrapper: ", wrapper_wall, " s")
end

abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_first_call()
