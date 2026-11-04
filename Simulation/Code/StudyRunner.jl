module OracleStudy

include("StudyConfig.jl")
using .StudyConfig
include(joinpath(LIBRARY_ROOT, "examples", "common.jl"))
using Arrow, DataFrames, Dates, DiffEqBase, LinearAlgebra, SHA, StaticArrays, TOML
include("EncounterTimingRecorder.jl")
using .EncounterTimingRecorder
include("TrajectoryWriter.jl")
include("EncounterProducts.jl")
using .EncounterProducts

export StudyCase, case_directory, run_case, study_configuration, StudyConfig

function study_configuration(c::StudyCase)
    planet = make_no_gram_planet(:earth)
    (planet.μ, planet.Rp_e, planet.J2) == (MU, EARTH_RADIUS_M, J2) ||
        error("Library Earth constants differ from the recorded study constants")
    conditions = [InitialCondition(a=planet.Rp_e + c.helper_altitude_km * 1000, e=0.0, i=0.0,
        Ω=0.0, ω=0.0, ν=360.0 * (k - 1) / c.helpers) for k in 1:c.helpers]
    push!(conditions, InitialCondition(a=planet.Rp_e + c.altitude_km * 1000,
        e=0.0, i=c.inclination_deg, Ω=0.0, ω=0.0, ν=1.0))
    spacecraft = map(eachindex(conditions)) do i
        bus = Link(root=true, m=MASS_KG)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=conditions[i], id=i, n_terminal=1)
    end
    links = ConstellationInterLinkModel()
    for helper in 1:c.helpers
        register_allowed_link!(links, spacecraft, (helper, 1), (c.helpers + 1, 1);
            parameters=InterLinkParameters(P=c.power_W, B=B, range=RANGE_M))
    end
    SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=mission_duration(c), orientation_sim=false,
            data_rate=mission_duration(c) / (TRAJECTORY_ROWS - 1)),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-12,
            abstol_orbit=1e-12, dt_max_orbit=STEP_S),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=links, interlink_diagnostics=true,
        scheduling_policy_model=SchedulingPolicyModel(c.policy; target_idx=c.helpers + 1))
end

function check_initial_geometry(state)
    nearest = Inf
    for a in 1:length(state.sc)-1, b in a+1:length(state.sc)
        distance = norm(state.sc[a].pos - state.sc[b].pos)
        nearest = min(nearest, distance)
    end
    nearest > 1e-3 || throw(ArgumentError("Initial spacecraft overlap or separation below 1 mm"))
    nearest
end

function completed_step_pairs(integrator, previous_time, power_W)
    history = integrator.p.args.interlink_model.history
    index = last(history).time > previous_time ? length(history) - 1 : length(history)
    index > 0 && history[index].time == previous_time ||
        error("Missing scheduler state for the completed integration step")
    power_W == 0 && return Tuple{Int,Int}[]
    pairs = [minmax(edge[1][1], edge[2][1]) for edge in history[index].active]
    length(pairs) <= 1 || error("Single-terminal target has multiple simultaneous links")
    pairs
end

function pipeline_callback(initial, c, directory, trajectory, spool; record_events=true, observer=nothing)
    duration = mission_duration(c)
    recorder = Recorder(snapshot(initial), RANGE_M)
    previous = snapshot(initial)
    previous_time = 0.0
    times = StudyConfig.trajectory_times(c)
    save_index = 2
    next_dense = ENCOUNTER_SAMPLE_S
    last_progress = time()
    start_wall = last_progress
    steps = 0
    minimum_radius = Inf
    maximum_step = 0.0

    function progress(simulation_time)
        atomic_toml(joinpath(directory, "progress.toml"), Dict(
            "pid" => getpid(), "simulation_time_s" => simulation_time,
            "duration_s" => duration, "fraction" => simulation_time / duration,
            "wall_elapsed_s" => time() - start_wall, "accepted_steps" => steps,
            "saved_rows" => trajectory.rows, "peak_rss_bytes" => peak_rss_bytes(),
            "updated_utc" => string(now(UTC))))
    end

    function affect!(integrator)
        stop = Float64(integrator.t)
        dt = stop - previous_time
        0 < dt <= STEP_S + 2eps(duration) || error("Invalid accepted step: $dt seconds")
        current = snapshot(integrator.u)
        pairs = completed_step_pairs(integrator, previous_time, c.power_W)
        all(pair -> pair[2] == c.helpers + 1, pairs) || error("Non-target firing link")
        dense(t) = snapshot(integrator(t))
        if record_events
            observe_firing!(recorder, pairs, previous, previous_time; initial=steps == 0)
            observe_geometry!(recorder, previous, current, previous_time, stop, dense)
        end
        observer === nothing || observer(integrator, previous_time, pairs)
        while save_index <= length(times) && times[save_index] <= stop
            saved_time = times[save_index]
            save_sample!(trajectory, saved_time, saved_time == stop ? integrator.u : integrator(saved_time))
            save_index += 1
        end
        if spool !== nothing
            while next_dense <= stop && next_dense < duration
                write_dense_sample(spool, next_dense, dense(next_dense))
                next_dense += ENCOUNTER_SAMPLE_S
            end
            stop == duration && write_dense_sample(spool, stop, current)
        end
        for sc in integrator.u.sc
            minimum_radius = min(minimum_radius, norm(sc.pos))
        end
        minimum_radius > EARTH_RADIUS_M + RANGE_M ||
            error("Trajectory left the verified nonocculted encounter domain")
        maximum_step = max(maximum_step, dt)
        steps += 1
        previous, previous_time = current, stop
        if time() - last_progress >= 30 || stop == duration
            require_disk_space(directory)
            progress(stop)
            println("progress=$(round(100stop / duration; digits=2))% t_s=$stop steps=$steps")
            flush(stdout)
            last_progress = time()
        end
        nothing
    end

    initialize = function (callback, state, time, integrator)
        time == 0.0 || error("This runner resumes completed cases, not partial trajectories")
        save_sample!(trajectory, 0.0, state)
        spool === nothing || write_dense_sample(spool, 0.0, snapshot(state))
        progress(0.0)
    end
    callback = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, affect!;
        initialize, save_positions=(false, false))
    metrics() = Dict("accepted_steps_observed" => steps, "maximum_step_s" => maximum_step,
        "minimum_radius_m" => minimum_radius,
        "maximum_encounter_boundary_residual_m" => recorder.max_boundary_residual_m)
    (; callback, recorder, metrics)
end

function write_table(path, table, run_id)
    Arrow.write(path * ".partial", table; file=true, compress=:zstd,
        metadata=Dict("run_id" => run_id, "schema_version" => string(SCHEMA_VERSION)))
    mv(path * ".partial", path)
end

function schedule_fingerprint(history)
    context = SHA.SHA2_256_CTX()
    for sample in history
        SHA.update!(context, codeunits(repr((sample.time, sample.active)) * "\n"))
    end
    bytes2hex(SHA.digest!(context))
end

function check_firing_history(laser, history, power)
    expected = Tuple{Int,Float64,Float64}[]
    for i in 1:length(history)-1
        sample, stop = history[i], history[i + 1].time
        power == 0 && continue
        isempty(sample.active) && continue
        pair = only(sample.active)
        helper = min(pair[1][1], pair[2][1])
        if !isempty(expected) && expected[end][1] == helper && expected[end][3] == sample.time
            expected[end] = (helper, expected[end][2], stop)
        else
            push!(expected, (helper, sample.time, stop))
        end
    end
    actual = [(row.sc_a, row.start_time_s, row.end_time_s) for row in eachrow(laser)]
    actual == expected || error("Actual firing intervals do not match the accepted-step force history")
    nothing
end

function run_case(c::StudyCase; output_root=DEFAULT_OUTPUT, record_events=true,
    observer=nothing, save_schedule=false)
    directory = case_directory(output_root, c)
    ispath(directory) && error("Refusing to overwrite $directory; resume through run_grid.jl")
    mkpath(dirname(directory))
    mkdir(directory)
    hashes = source_hashes()
    info = Dict{String,Any}("status" => "running", "run_id" => case_signature(c),
        "settings" => case_settings(c), "source_hashes" => hashes,
        "source_signature" => source_signature(hashes), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "record_events" => record_events,
        "started_utc" => string(now(UTC)), "pid" => getpid(),
        "trajectory_chunk_rows" => 512, "rhs_execution_mode" => "serial")
    info_path = joinpath(directory, "run_info.toml")
    atomic_toml(info_path, info)
    trajectory = nothing
    spool = nothing
    spool_path = joinpath(directory, "encounter_states.bin.partial")
    start_wall = time()
    try
        require_disk_space(directory)
        args = study_configuration(c)
        initial = SimulationEngine.build_initial_conditions(args)
        info["minimum_initial_separation_m"] = check_initial_geometry(initial)
        info["initial_states"] = [Dict("id" => sc.id, "index" => i,
            "position_m" => collect(initial.sc[i].pos),
            "velocity_mps" => collect(initial.sc[i].vel), "mass_kg" => initial.sc[i].mass)
            for (i, sc) in enumerate(args.dynamics_model.spacecraft)]
        atomic_toml(info_path, info)
        trajectory = TrajectoryWriter(joinpath(directory, "trajectory.feather.partial"),
            c.helpers + 1, info["run_id"])
        spool = record_events ? open(spool_path, "w") : nothing
        pipeline = pipeline_callback(initial, c, directory, trajectory, spool; record_events, observer)
        println("Starting $(c.policy): helper$(c.helper_altitude_km)km target$(c.altitude_km)km i=$(c.inclination_deg) N=$(c.helpers), $(c.orbits) orbits")
        println("output=$directory")
        flush(stdout)
        solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false",
            "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
            "SPACEAGORA_PARALLEL_PROFILE" => "R0", "SPACEAGORA_VISUALIZATION" => "0") do
            run_simulation(args; return_solution=true, save_fields=SM.SaveField[],
                extra_callbacks=(pipeline.callback,))
        end
        string(solution.retcode) == "Success" || error("Solver returned $(solution.retcode)")
        last(solution.t) == mission_duration(c) || error("Simulation ended before the requested duration")
        trajectory.last_time == mission_duration(c) || error("Missing final trajectory sample")
        trajectory.rows == TRAJECTORY_ROWS || error("Incorrect trajectory row count")
        close_trajectory!(trajectory)
        spool === nothing || close(spool)
        mv(joinpath(directory, "trajectory.feather.partial"), joinpath(directory, "trajectory.feather"))
        record_events && finish!(pipeline.recorder, snapshot(solution.u[end]), mission_duration(c))
        geometry, laser = pipeline.recorder.geometry, pipeline.recorder.laser
        geometry.run_id = fill(info["run_id"], nrow(geometry))
        laser.run_id = fill(info["run_id"], nrow(laser))
        laser.power_W = fill(c.power_W, nrow(laser))
        laser.force_N = fill(B * c.power_W / LIGHT_SPEED, nrow(laser))
        history = solution.prob.p.args.interlink_model.history
        record_events && check_firing_history(laser, history, c.power_W)
        tracks = build_encounter_products!(geometry, spool_path, c.helpers + 1,
            info["run_id"], mission_duration(c); cadence=ENCOUNTER_SAMPLE_S)
        segments = range_segments(geometry, laser; run_id=info["run_id"])
        files = ["trajectory.feather", "encounters.feather", "laser_firing.feather",
            "firing_range_audit.feather", "longest_encounter.feather"]
        write_table(joinpath(directory, files[2]), geometry, info["run_id"])
        write_table(joinpath(directory, files[3]), laser, info["run_id"])
        write_table(joinpath(directory, files[4]), segments, info["run_id"])
        write_table(joinpath(directory, files[5]), tracks, info["run_id"])
        info["temporary_dense_state_bytes"] = record_events ? filesize(spool_path) : 0
        record_events && rm(spool_path)
        info["schedule_sha256"] = schedule_fingerprint(history)
        if save_schedule
            table = DataFrame(time=[sample.time for sample in history],
                active_helper=[isempty(sample.active) ? 0 : only(sample.active)[1][1] for sample in history])
            write_table(joinpath(directory, "schedule.feather"), table, info["run_id"])
            push!(files, "schedule.feather")
        end
        merge!(info, pipeline.metrics())
        info["solver_accepted_steps"] = solution.stats.naccept
        info["solver_rejected_steps"] = solution.stats.nreject
        info["trajectory_rows"] = trajectory.rows
        info["encounters"] = nrow(geometry)
        info["firing_intervals"] = nrow(laser)
        info["actual_firing_s"] = sum(laser.duration_s)
        info["outside_range_firing_s"] = sum(segments.duration_s[.!segments.within_range])
        info["checksums"] = Dict(name => filehash(joinpath(directory, name)) for name in files)
        info["output_bytes"] = sum(filesize(joinpath(directory, name)) for name in files)
        info["wall_elapsed_s"] = time() - start_wall
        info["peak_rss_bytes"] = peak_rss_bytes()
        info["finished_utc"] = string(now(UTC))
        info["status"] = "complete"
        atomic_toml(info_path, info)
        println("COMPLETE output=$directory bytes=$(info["output_bytes"]) peak_rss_bytes=$(info["peak_rss_bytes"])")
        flush(stdout)
        return directory
    catch err
        info["status"] = "failed"
        info["error"] = sprint(showerror, err)
        info["failed_utc"] = string(now(UTC))
        atomic_toml(info_path, info)
        @error "Study case failed; partial output was retained" directory exception=(err, catch_backtrace())
        rethrow()
    finally
        trajectory === nothing || close_trajectory!(trajectory)
        spool === nothing || !isopen(spool) || close(spool)
    end
end

end
