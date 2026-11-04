include(joinpath(@__DIR__, "validate_encounter_recording.jl"))

const PREFLIGHT_ROOT = joinpath(@__DIR__, "results", "study_preflight")

function checked_initials(altitude, inclination, helpers; phase=1.0)
    helpers > 0 || throw(ArgumentError("Helper count must be positive"))
    planet = make_no_gram_planet(:earth)
    conditions = [InitialCondition(a=planet.Rp_e + 1000e3, e=0.0, i=0.0,
        Ω=0.0, ω=0.0, ν=360.0 * (i - 1) / helpers) for i in 1:helpers]
    push!(conditions, InitialCondition(a=planet.Rp_e + altitude * 1000,
        e=0.0, i=inclination, Ω=0.0, ω=0.0, ν=phase))
    # Check the requested circular-orbit positions before starting the solver.
    positions = [SVector(c.a * cos(c.ν), c.a * sin(c.ν) * cos(c.i),
        c.a * sin(c.ν) * sin(c.i)) for c in conditions]
    minimum_separation = minimum(norm(positions[a] - positions[b])
        for a in 1:length(positions)-1 for b in a+1:length(positions))
    minimum_separation > 1e-3 ||
        throw(ArgumentError("Initial spacecraft overlap or separation below 1 mm"))
    return conditions, minimum_separation
end

function preflight_config(; altitude=1050.0, inclination=1.0, helpers=20,
    policy=:gve_sma, dtmax=10.0, save_interval=10.0, orbits=100.0, duration=nothing,
    power=10_000.0, j2=true, links_enabled=true, output="", save=true)
    conditions, _ = checked_initials(altitude, inclination, helpers)
    planet = make_no_gram_planet(:earth)
    spacecraft = map(eachindex(conditions)) do i
        bus = Link(root=true, m=MASS_KG)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=conditions[i], id=i, n_terminal=1)
    end
    links = ConstellationInterLinkModel()
    for helper in 1:helpers
        register_allowed_link!(links, spacecraft, (helper, 1), (helpers + 1, 1);
            parameters=InterLinkParameters(P=power, B=100.0, range=RANGE_M))
    end
    end_time = duration === nothing ? orbits * 2pi *
        sqrt(conditions[end].a^3 / planet.μ) : duration
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=save, generate_plots=false,
            save_csv=false, results_directory=output),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=end_time, orientation_sim=false, data_rate=save_interval),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft,
            (j2 ? InverseSquaredJ2GravityModel() : InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-12,
            abstol_orbit=1e-12, dt_max_orbit=dtmax),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=links_enabled ? links : nothing,
        scheduling_policy_model=SchedulingPolicyModel(policy; target_idx=helpers + 1))
end

function endpoint_run(args)
    return withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(args; return_solution=true)
    end
end

function energy_and_momentum(state, planet, j2)
    energy = 0.0
    momentum = zero(SVector{3, Float64})
    for i in axes(state, 2)
        r, v = SVector{3,Float64}(state[1:3, i]), SVector{3,Float64}(state[4:6, i])
        radius = norm(r)
        potential = -planet.μ / radius
        if j2
            potential += planet.μ * planet.J2 * planet.Rp_e^2 / (2radius^3) *
                (3(r[3] / radius)^2 - 1)
        end
        energy += MASS_KG * (dot(v, v) / 2 + potential)
        momentum += MASS_KG * cross(r, v)
    end
    return energy, momentum
end

function check_segments(geometry, laser, segments)
    @test all(segments.end_time_s .> segments.start_time_s)
    @test all(segments.duration_s .== segments.end_time_s .- segments.start_time_s)
    @test all(segments.within_range .== .!ismissing.(segments.encounter_id))
    @test all(segments.laser_id .∈ Ref(Set(laser.laser_id)))
    for firing in eachrow(laser)
        parts = sort(filter(s -> s.laser_id == firing.laser_id, segments), :start_time_s)
        @test !isempty(parts)
        @test first(parts.start_time_s) == firing.start_time_s
        @test last(parts.end_time_s) == firing.end_time_s
        @test parts.end_time_s[1:end-1] == parts.start_time_s[2:end]
        @test all(parts.sc_a .== firing.sc_a) && all(parts.sc_b .== firing.sc_b)
        @test isapprox(sum(parts.duration_s), firing.duration_s; atol=1e-8)
        for part in eachrow(parts)
            if part.within_range
                e = only(eachrow(filter(e -> e.encounter_id == part.encounter_id, geometry)))
                @test (e.sc_a, e.sc_b) == (part.sc_a, part.sc_b)
                @test e.start_time_s <= part.start_time_s < part.end_time_s <= e.end_time_s
            else
                @test !any(e -> (e.sc_a, e.sc_b) == (part.sc_a, part.sc_b) &&
                    e.start_time_s < part.end_time_s && e.end_time_s > part.start_time_s,
                    eachrow(geometry))
            end
        end
    end
    encounters, audit = analyze_intervals(geometry, laser)
    @test isapprox(sum(segments.duration_s[segments.within_range]),
        sum(encounters.laser_on_s); atol=1e-8)
    @test isapprox(sum(segments.duration_s[.!segments.within_range]),
        sum(audit.outside_range_s); atol=1e-8)
end

function fixture_checks()
    mkpath(PREFLIGHT_ROOT)
    synthetic_tests()
    results = Dict{String,Any}()
    @testset "Phase grid and LOS" begin
        separations = Float64[]
        for h in (850.0, 950.0, 1000.0, 1050.0, 1150.0), i in (0.0, 0.5, 1.0),
            n in (1, 50, 100, 150, 200, 250, 300)
            conditions, minimum_distance = checked_initials(h, i, n)
            @test isapprox(conditions[end].ν, deg2rad(1.0))
            @test all(c.a == 6_378_136.6 + 1000e3 for c in conditions[1:end-1])
            push!(separations, minimum_distance)
        end
        @test_throws ArgumentError checked_initials(1000.0, 0.0, 50; phase=0.0)
        # For any point on a <=200 km segment, the reverse triangle inequality
        # gives radius >= endpoint radius - 200 km.
        clearance = 850e3 - RANGE_M
        @test clearance == 650e3
        results["minimum_initial_pair_separation_m"] = minimum(separations)
        results["initial_domain_LOS_clearance_lower_bound_m"] = clearance
    end
    @testset "Force size and pair balances" begin
        args = preflight_config(helpers=1, save=false, duration=120.0)
        initial = SimulationEngine.build_initial_conditions(args)
        args.interlink_model.graph.active .= true
        f1 = SM.laser_force_on_spacecraft(args.interlink_model, initial, 1)
        f2 = SM.laser_force_on_spacecraft(args.interlink_model, initial, 2)
        @test isapprox(norm(f1), FORCE_N; rtol=1e-14)
        @test isapprox(norm(f1), 0.0033356409519815205; rtol=1e-14)
        @test f1 == -f2
        torque = cross(initial.sc[1].pos, f1) + cross(initial.sc[2].pos, f2)
        @test norm(torque) < 1e-8
        initial.sc[2].mass = 2MASS_KG
        derivative = copy(initial)
        fill!(derivative, 0.0)
        SimulationEngine._apply_interlink_rhs!(derivative, initial,
            (args=args, is_active=trues(2)))
        @test isapprox(norm(derivative.sc[1].vel), FORCE_N / MASS_KG; rtol=1e-14)
        @test isapprox(derivative.sc[1].vel, -2derivative.sc[2].vel; rtol=1e-14)
        @test derivative.sc[1].laser_dv == derivative.sc[1].vel
        @test derivative.sc[2].laser_dv == derivative.sc[2].vel
        results["force_N"] = norm(f1)
        results["pair_torque_residual_Nm"] = norm(torque)
    end
    @testset "Laser on/off and zero power" begin
        for policy in (:gve_sma, :gve_inc)
            base = (; helpers=1, duration=120.0, save=false, policy, dtmax=1.0)
            on = endpoint_run(preflight_config(; base...))
            off = endpoint_run(preflight_config(; base..., power=0.0))
            gravity = endpoint_run(preflight_config(; base..., links_enabled=false))
            @test all(string(s.retcode) == "Success" for s in (on, off, gravity))
            @test norm(on.u[end].sc[2].pos - off.u[end].sc[2].pos) > 0.01
            @test norm(on.u[end].sc[2].laser_dv) > 1e-5
            @test all(iszero(sc.laser_dv) for sc in off.u[end].sc)
            @test maximum(abs, snapshot(off.u[end]) - snapshot(gravity.u[end])) < 1e-5
            results["$(policy)_on_off_position_difference_m"] =
                norm(on.u[end].sc[2].pos - off.u[end].sc[2].pos)
        end
        equatorial = endpoint_run(preflight_config(helpers=1, inclination=0.0,
            policy=:gve_inc, duration=120.0, save=false))
        @test all(iszero(sc.laser_dv) for sc in equatorial.u[end].sc)
    end
    @testset "Simultaneous contacts and range-audit partitions" begin
        state(t) = [0.0 -2.0+4t 2.0-4t; zeros(2,3); 0.0 4.0 -4.0; zeros(2,3)]
        r = Recorder(state(0.0), 1.0)
        observe_geometry!(r, state(0.0), state(1.0), 0.0, 1.0, state)
        finish!(r, state(1.0), 1.0)
        target_contacts = filter(e -> e.sc_a == 1, r.geometry)
        @test nrow(target_contacts) == 2
        @test all(isapprox.(target_contacts.start_time_s, 0.25; atol=1e-8))
        @test all(isapprox.(target_contacts.end_time_s, 0.75; atol=1e-8))

        # A cubic crosses both sides of the sphere repeatedly during one step.
        multi(t) = [0.0 40(t-0.2)*(t-0.5)*(t-0.8); zeros(2,2);
            0.0 40(3t^2-3t+0.66); zeros(2,2)]
        r = Recorder(multi(0.0), 0.2)
        observe_firing!(r, [(1, 2)], multi(0.0), 0.0; initial=true)
        observe_geometry!(r, multi(0.0), multi(1.0), 0.0, 1.0, multi)
        finish!(r, multi(1.0), 1.0)
        @test nrow(r.geometry) == 3
        segments = range_segments(r.geometry, r.laser; run_id="multiple-crossings")
        @test nrow(segments) == 7
        check_segments(r.geometry, r.laser, segments)
        Arrow.write(joinpath(PREFLIGHT_ROOT, "multiple_crossing_segments.feather"), segments)
        @test_throws ArgumentError range_segments(vcat(r.geometry, r.geometry), r.laser; run_id="bad")

        outside = [0.0 2.0; zeros(5, 2)]
        r = Recorder(outside, 1.0)
        observe_firing!(r, [(1, 2)], outside, 0.0; initial=true)
        finish!(r, outside, 1.0)
        segments = range_segments(r.geometry, r.laser; run_id="outside")
        @test nrow(segments) == 1
        @test !only(segments.within_range)
        check_segments(r.geometry, r.laser, segments)
        r = Recorder(outside, 1.0)
        finish!(r, outside, 1.0)
        segments = range_segments(r.geometry, r.laser; run_id="no-firing")
        @test isempty(segments)
        check_segments(r.geometry, r.laser, segments)
        Arrow.write(joinpath(PREFLIGHT_ROOT, "no_firing_segments.feather"), segments)
    end
    open(joinpath(PREFLIGHT_ROOT, "fixtures.toml"), "w") do io
        TOML.print(io, results; sorted=true)
    end
    TOML.print(stdout, results; sorted=true)
end

function diagnostics_observer(args; record_events=true)
    initial = snapshot(SimulationEngine.build_initial_conditions(args))
    count = size(initial, 2)
    target = count
    planet = args.environment_model.planet
    j2 = only(args.dynamics_model.dynamic_effectors) isa InverseSquaredJ2GravityModel
    duration = args.mission_configuration.mission_time
    interval = args.mission_configuration.data_rate
    times = collect(0.0:interval:duration)
    last(times) < duration && push!(times, duration)
    next_save = 2
    previous, previous_time = initial, 0.0
    totals = zeros(5 + count)
    rows = DataFrame(time_s=Float64[], target_dv_R_mps=Float64[], target_dv_T_mps=Float64[],
        target_dv_N_mps=Float64[], target_dv_T_positive_mps=Float64[],
        target_dv_T_negative_mps=Float64[])
    for i in 1:count
        rows[!, "sc$(i)_work_J"] = Float64[]
    end
    push!(rows, (0.0, totals...))
    recorder = Recorder(initial, RANGE_M)
    E0, H0 = energy_and_momentum(initial, planet, j2)
    max_energy_error, max_momentum_error = Ref(0.0), Ref(0.0)
    min_radius = Ref(minimum(norm(initial[1:3, i]) for i in 1:count))
    max_quadrature_error = Ref(0.0)
    function rates(state, pairs)
        values = zeros(5 + count)
        isempty(pairs) && return values
        length(pairs) == 1 || error("Multiple selected target links")
        a, b = only(pairs)
        b == target || error("Unexpected selected target")
        r = SVector{3,Float64}(state[1:3, target])
        v = SVector{3,Float64}(state[4:6, target])
        force = FORCE_N * normalize(r - SVector{3,Float64}(state[1:3, a]))
        R = normalize(r)
        N = normalize(cross(r, v))
        T = cross(N, R)
        values[1:3] = [dot(force, R), dot(force, T), dot(force, N)] / MASS_KG
        values[4] = max(values[2], 0)
        values[5] = min(values[2], 0)
        values[5 + target] = dot(force, v)
        values[5 + a] = -dot(force, state[4:6, a])
        return values
    end
    function affect!(integrator)
        time = Float64(integrator.t)
        current = snapshot(integrator.u)
        history = integrator.p.args.interlink_model.history
        index = history[end].time > previous_time ? length(history) - 1 : length(history)
        sample = history[index]
        sample.time == previous_time || error("Scheduler sample does not cover this step")
        pairs = sample_pairs(sample)
        dense(t) = t == previous_time ? previous : t == time ? current : snapshot(integrator(t))
        if record_events
            observe_firing!(recorder, pairs, previous, previous_time; initial=previous_time == 0)
            observe_geometry!(recorder, previous, current, previous_time, time, dense)
        end
        function integrate_to(stop)
            h = stop - previous_time
            h == 0 && return zeros(length(totals))
            q0 = rates(previous, pairs)
            q1 = rates(dense(previous_time + h / 4), pairs)
            q2 = rates(dense(previous_time + h / 2), pairs)
            q3 = rates(dense(previous_time + 3h / 4), pairs)
            q4 = rates(dense(stop), pairs)
            fine = h / 12 * (q0 + 4q1 + 2q2 + 4q3 + q4)
            coarse = h / 6 * (q0 + 4q2 + q4)
            max_quadrature_error[] = max(max_quadrature_error[], norm(fine[1:3] - coarse[1:3]))
            return fine
        end
        while next_save <= length(times) && times[next_save] <= time
            saved_time = times[next_save]
            push!(rows, (saved_time, (totals + integrate_to(saved_time))...))
            next_save += 1
        end
        totals .+= integrate_to(time)
        energy, momentum = energy_and_momentum(current, planet, j2)
        max_energy_error[] = max(max_energy_error[], abs(energy - E0 - sum(totals[6:end])))
        defect = j2 ? abs(momentum[3] - H0[3]) : norm(momentum - H0)
        max_momentum_error[] = max(max_momentum_error[], defect / norm(H0))
        min_radius[] = min(min_radius[], minimum(norm(current[1:3, i]) for i in 1:count))
        previous, previous_time = current, time
    end
    callback = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, affect!;
        save_positions=(false, false))
    return (; callback, recorder, rows, totals, max_energy_error, max_momentum_error,
        min_radius, max_quadrature_error)
end

function refinement_run(case, dtmax, save_interval; control=false, j2=true, short=false)
    case in ("sma850", "inc1050") || error("Unknown case $case")
    policy = case == "sma850" ? :gve_sma : :gve_inc
    altitude = case == "sma850" ? 850.0 : 1050.0
    label = "$(case)_dt$(dtmax)_save$(save_interval)" *
        (control ? "_control" : "") * (j2 ? "" : "_noJ2") * (short ? "_short" : "")
    output = joinpath(PREFLIGHT_ROOT, label)
    isdir(output) && !isempty(readdir(output)) && error("Refusing to overwrite $output")
    mkpath(output)
    args = preflight_config(; altitude, policy, dtmax, save_interval, output, j2,
        orbits=short ? 2.0 : 100.0)
    observer = diagnostics_observer(args; record_events=!control)
    fields = filter(f -> f.name in (:position, :velocity, :laser_dv, :laser_delta_sma),
        SM.default_save_fields(args))
    println("Running $label for $(args.mission_configuration.mission_time) s")
    flush(stdout)
    solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(args; return_solution=true, save_fields=fields,
            extra_callbacks=(observer.callback,))
    end
    duration = args.mission_configuration.mission_time
    @test string(solution.retcode) == "Success"
    @test solution.t[end] == duration
    @test all(isfinite, solution.u[end])
    @test maximum(abs, Vector(observer.rows[end, 2:end]) .- observer.totals) < 1e-8
    @test all(isapprox.(observer.rows.target_dv_T_mps,
        observer.rows.target_dv_T_positive_mps + observer.rows.target_dv_T_negative_mps; atol=1e-12))
    @test observer.min_radius[] - RANGE_M > args.environment_model.planet.Rp_e
    @test observer.max_momentum_error[] < 1e-9
    work = sum(observer.totals[6:end])
    energy_tolerance = 20.0 + 1e-4 * abs(work)
    @test observer.max_energy_error[] < energy_tolerance
    @test norm(sum(MASS_KG * SVector{3,Float64}(s.laser_dv) for s in solution.u[end].sc)) < 1e-8
    history = solution.prob.p.args.interlink_model.history
    schedule = DataFrame(time_s=[s.time for s in history],
        active_helper=[isempty(s.active) ? 0 : only(s.active)[1][1] for s in history])
    @test length(schedule.time_s) == solution.stats.naccept + 1
    @test maximum(diff(schedule.time_s)) <= dtmax + 2eps(duration)
    Arrow.write(joinpath(output, "diagnostics.feather"), observer.rows)
    Arrow.write(joinpath(output, "schedule.feather"), schedule)
    planet = args.environment_model.planet
    function elements(state)
        r = SVector{3,Float64}(state.sc[end].pos)
        v = SVector{3,Float64}(state.sc[end].vel)
        h = cross(r, v)
        return inv(2 / norm(r) - dot(v, v) / planet.μ), rad2deg(atan(hypot(h[1], h[2]), h[3]))
    end
    a0, i0 = elements(solution.u[1])
    af, inf = elements(solution.u[end])
    summary = Dict{String,Any}(
        "case" => case, "label" => label, "retcode" => string(solution.retcode),
        "target_altitude_km" => altitude, "helper_altitude_km" => 1000.0,
        "inclination_deg" => 1.0, "target_nu_deg" => 1.0, "helper_count" => 20,
        "target_idx" => 21, "duration_s" => duration, "orbits" => short ? 2 : 100,
        "dtmax_s" => dtmax, "save_interval_s" => save_interval, "J2_enabled" => j2,
        "scheduler" => string(policy), "range_m" => RANGE_M, "mass_kg" => MASS_KG,
        "power_W" => 10_000.0, "B" => 100.0, "reltol" => 1e-12, "abstol" => 1e-12,
        "mu" => planet.μ, "earth_radius_m" => planet.Rp_e, "J2" => planet.J2,
        "delta_sma_m" => af - a0, "delta_inc_deg" => inf - i0,
        "dv_RTN_mps" => observer.totals[1:3], "optical_work_total_J" => work,
        "max_energy_work_error_J" => observer.max_energy_error[],
        "energy_work_tolerance_J" => energy_tolerance,
        "max_relative_angular_momentum_error" => observer.max_momentum_error[],
        "angular_momentum_quantity" => j2 ? "axial component" : "full vector",
        "min_radius_m" => observer.min_radius[],
        "max_step_RTN_quadrature_difference_mps" => observer.max_quadrature_error[],
        "diagnostic_method" => "Composite Simpson quadrature on solver interpolation; NOT RTN solver states",
        "accepted_steps" => solution.stats.naccept, "rejected_steps" => solution.stats.nreject,
        "firing_s" => sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)),
        "test_sha256" => filehash(@__FILE__), "recorder_sha256" => filehash(joinpath(@__DIR__, "EncounterTimingRecorder.jl")))
    if !control
        recorder = finish!(observer.recorder, snapshot(solution.u[end]), duration)
        segments = range_segments(recorder.geometry, recorder.laser; run_id=label)
        check_segments(recorder.geometry, recorder.laser, segments)
        @test observer.recorder.max_boundary_residual_m < 1e-3
        @test isapprox(sum(recorder.laser.duration_s), summary["firing_s"]; atol=1e-8)
        for (name, data) in (("geometry", recorder.geometry), ("laser", recorder.laser),
            ("range_segments", segments))
            Arrow.write(joinpath(output, "$name.feather"), data)
        end
        summary["encounter_count"] = nrow(recorder.geometry)
        summary["firing_count"] = nrow(recorder.laser)
        summary["geometry_s"] = sum(recorder.geometry.duration_s)
        summary["outside_firing_s"] = sum(segments.duration_s[.!segments.within_range])
    end
    open(joinpath(output, "summary.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    TOML.print(stdout, summary; sorted=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Study preflight checks" begin
        if ARGS == ["fixtures"]
            fixture_checks()
        elseif length(ARGS) >= 4 && ARGS[1] == "run"
            refinement_run(ARGS[2], parse(Float64, ARGS[3]), parse(Float64, ARGS[4]);
                control="control" in ARGS, j2=!("noJ2" in ARGS), short="short" in ARGS)
        else
            error("Usage: preflight_study_checks.jl fixtures | run CASE DTMAX SAVE_INTERVAL [control] [noJ2] [short]")
        end
    end
end
