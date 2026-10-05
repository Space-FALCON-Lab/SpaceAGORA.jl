using SpaceAGORA
using SpaceAGORA.SimulationModel
using LinearAlgebra
using DiffEqBase
using DiffEqCallbacks
using StaticArrays
using CSV
using DataFrames
using Printf
using Test

const CASE1_MODEL = SpaceAGORA.SimulationModel
const CASE1_ENGINE = SpaceAGORA.SimulationEngine
const CASE1_OUTPUT = joinpath(@__DIR__, "..", "III_output", "pdf_case1_comparison")

case1_sma(state, mu) = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)

function case1_configuration(mode::Symbol; orbits::Float64=10.0, dtmax::Float64=10.0,
    laser_enabled::Bool=true)
    # Step 1: Validate the implementation mode and define the Earth/orbit geometry.
    is_v5 = mode === :v5
    mode in (:v5, :legacy) || throw(ArgumentError("Mode must be :v5 or :legacy."))
    @assert is_v5 == isdefined(CASE1_MODEL, :InterLinkModel)
    planet = make_no_gram_planet(:earth)
    target_radius = planet.Rp_e + 1050e3
    helper_radius = planet.Rp_e + 1000e3
    period = 2pi * sqrt(target_radius^3 / planet.μ)

    # Step 2: Create one target and twenty evenly spaced helper spacecraft.
    spacecraft = map(1:21) do index
        bus = Link(root=true, m=227.0)
        radius = index == 1 ? target_radius : helper_radius
        phase = index == 1 ? 0.0 : 360.0 * (index - 2) / 20
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, id=index,
            initial_condition=InitialCondition(radius, 0.0, 0.0, 0.0, 0.0, phase))
    end

    # Step 3: Attach either the v5 scheduler or the legacy laser model when enabled.
    extra_config = NamedTuple()
    dynamic_effectors = (InverseSquaredJ2GravityModel(),)
    laser = nothing
    if is_v5 && laser_enabled
        laser = InterLinkModel()
        for helper in 2:21
            register_candidate!(laser, spacecraft, (1, 1), (helper, 1);
                parameters=InterLinkParameters(P=10_000.0, B=100.0, range=200e3))
        end
        extra_config = (interlink_model=laser, scheduling_policy_model=SchedulingPolicyModel(:gve_sma; target_idx=1))
    elseif laser_enabled
        laser = OpenCavityLaserLinkModel(target_idx=1, helper_indices=collect(2:21),
            range_m=200e3, power_w=10_000.0, magnification=100.0, beta=1.0,
            eta=1.0, schedule=:gve_sma)
        dynamic_effectors = (dynamic_effectors..., laser)
    end

    # Step 4: Assemble the mission, J2 environment, dynamics, and solver settings.
    args = SimulationConfiguration(;
        simulation_settings=SimulationSettings(results=false, generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=orbits * period, orientation_sim=false, data_rate=10.0),
        environment_model=EnvironmentModel(planet=planet, EI=120.0,
            density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-12,
            abstol_orbit=1e-12, dt_max_orbit=dtmax),
        solver_config=CASE1_MODEL.SolverConfig(solver_mode=:tsit5), extra_config...)
    return args, laser, period
end

function run_case1(mode::Symbol; orbits::Float64=10.0, dtmax::Float64=10.0,
    laser_enabled::Bool=true)
    # Step 1: Build the selected case and prepare initial semimajor axes and sample times.
    args, laser, period = case1_configuration(mode; orbits, dtmax, laser_enabled)
    duration = args.mission_configuration.mission_time
    initial = CASE1_ENGINE.build_initial_conditions(args)
    mu = args.environment_model.planet.μ
    initial_sma = [case1_sma(current, mu) for current in initial.sc]
    sample_times = collect(0.0:10.0:duration)
    sample_times[end] < duration && push!(sample_times, duration)
    values = SavedValues(Float64, Vector{Float64})

    # Step 2: Save semimajor-axis samples and, for legacy, record helper switches.
    save_callback = SavingCallback(
        (state, time, integrator) -> [case1_sma(current, mu) for current in state.sc],
        values; saveat=sample_times, save_everystep=false, save_end=false)
    legacy_history = NamedTuple{(:time_s, :active_helper), Tuple{Float64, Int}}[]
    extra_callbacks = (save_callback,)
    if mode === :legacy && laser_enabled
        tracker = LaserImpulseTracker()
        function record_schedule!(integrator)
            push!(legacy_history, (time_s=Float64(integrator.t), active_helper=laser.active_helper_idx))
            DiffEqBase.u_modified!(integrator, false)
        end
        recorder = DiffEqBase.DiscreteCallback((state, time, integrator) -> true,
            record_schedule!; initialize=(callback, state, time, integrator) -> record_schedule!(integrator),
            save_positions=(false, false))
        extra_callbacks = (laser_impulse_callback(laser, tracker, 227.0),
            laser_link_scheduler_callback(laser), recorder, save_callback)
    end

    # Step 3: Run the simulation and verify its endpoint and saved samples.
    println("Running ", mode, ", laser=", laser_enabled, ", ", orbits,
        " initial target periods, dtmax=", dtmax, " s; package=", pathof(SpaceAGORA))
    elapsed = @elapsed solution = withenv("SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false", "SPACEAGORA_SOLVER_DENSE" => "false") do
        run_simulation(args; isolate_state=false, return_solution=true,
            extra_callbacks=extra_callbacks, save_fields=SaveField[])
    end
    @test solution.t[end] == duration
    @test values.t == sample_times
    @test all(isfinite, solution.u[end])
    @test all(all(isfinite, sample) for sample in values.saveval)
    @test all(isapprox(values.saveval[end][index], case1_sma(solution.u[end].sc[index], mu); atol=1e-5, rtol=0) for index in 1:21)

    # Step 4: Build and save per-spacecraft SMA histories and link schedules.
    table = DataFrame(time_s=values.t)
    for satellite in 1:21
        table[!, Symbol("satellite_$(satellite)_sma_m")] = [sample[satellite] for sample in values.saveval]
        table[!, Symbol("satellite_$(satellite)_delta_sma_m")] = table[!, Symbol("satellite_$(satellite)_sma_m")] .- initial_sma[satellite]
    end
    suffix = "$(mode)_$(laser_enabled ? "laser" : "no_laser")_$(orbits)orbits_dt$(dtmax)s"
    mkpath(CASE1_OUTPUT)
    CSV.write(joinpath(CASE1_OUTPUT, "$(suffix)_semimajor_axis.csv"), table)
    if laser_enabled
        history = mode === :v5 ? DataFrame(time_s=[sample.time for sample in laser.history],
            active_helper=[isempty(sample.active) ? 0 : only(sample.active)[2][1] for sample in laser.history]) : DataFrame(legacy_history)
        CSV.write(joinpath(CASE1_OUTPUT, "$(suffix)_schedule.csv"), history)
    end

    # Step 5: Save a compact run summary and report the target's total SMA change.
    summary = (implementation=String(mode), laser_enabled=laser_enabled,
        helpers=20, target_altitude_km=1050.0, helper_altitude_km=1000.0,
        mass_kg=227.0, power_w=10_000.0, magnification=100.0, range_m=200e3,
        force_per_endpoint_N=laser_enabled ? 1e6 / 299_792_458.0 : 0.0,
        objective="target_sma_rate",
        gravity="J2", mu_m3_s2=mu, earth_equatorial_radius_m=args.environment_model.planet.Rp_e,
        target_period_s=period, orbits=orbits, duration_s=duration, dtmax_s=dtmax,
        reltol_orbit=1e-12, abstol_orbit=1e-12, accepted_steps=solution.destats.naccept,
        rejected_steps=solution.destats.nreject, elapsed_wall_s=elapsed,
        initial_target_sma_m=initial_sma[1], final_target_sma_m=case1_sma(solution.u[end].sc[1], mu),
        target_delta_sma_m=case1_sma(solution.u[end].sc[1], mu) - initial_sma[1],
        target_integrated_laser_delta_sma_m=mode === :v5 && laser_enabled ? solution.u[end].sc[1].laser_delta_sma : missing,
        retcode=string(solution.retcode), source_file=pathof(SpaceAGORA))
    CSV.write(joinpath(CASE1_OUTPUT, "$(suffix)_summary.csv"), DataFrame([summary]))
    @printf("%s laser=%s: target delta SMA = %.9f m; %d accepted / %d rejected steps; wall %.1f s\n",
        String(mode), string(laser_enabled), summary.target_delta_sma_m,
        summary.accepted_steps, summary.rejected_steps, elapsed)
    return table, summary
end

function compare_case1(mode::Symbol; orbits::Float64=10.0, dtmax::Float64=10.0)
    # Step 1: Run matched laser-on and laser-off cases for the selected implementation.
    @testset "PDF Case 1 $(mode)" begin
        laser_table, _ = run_case1(mode; orbits, dtmax)
        baseline_table, _ = run_case1(mode; orbits, dtmax, laser_enabled=false)

        # Step 2: Confirm matching sample times and export the target comparison.
        @test laser_table.time_s == baseline_table.time_s
        target = DataFrame(time_s=laser_table.time_s,
            target_sma_m=laser_table.satellite_1_sma_m,
            target_delta_sma_m=laser_table.satellite_1_delta_sma_m,
            no_laser_target_sma_m=baseline_table.satellite_1_sma_m,
            target_sma_gain_vs_no_laser_m=laser_table.satellite_1_sma_m - baseline_table.satellite_1_sma_m)
        CSV.write(joinpath(CASE1_OUTPUT, "$(mode)_$(orbits)orbits_dt$(dtmax)s_target.csv"), target)
        @printf("%s target SMA difference from laser-off reference: %.9f m\n",
            String(mode), target.target_sma_gain_vs_no_laser_m[end])
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Step 1: Read the implementation and run-size options from the command line.
    mode = Symbol(isempty(ARGS) ? "v5" : ARGS[1])
    orbits = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 10.0
    dtmax = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 10.0

    # Step 2: Compare laser-on and laser-off simulations and write the CSV outputs.
    compare_case1(mode; orbits, dtmax)
end