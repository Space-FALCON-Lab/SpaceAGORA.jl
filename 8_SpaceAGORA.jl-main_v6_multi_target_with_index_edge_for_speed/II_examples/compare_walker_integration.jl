using SpaceAGORA
using SpaceAGORA.SimulationModel
using LinearAlgebra
using DiffEqBase
using DiffEqCallbacks
using CSV
using DataFrames
using Printf
using Test

const WALKER_MODEL = SpaceAGORA.SimulationModel
const WALKER_ENGINE = SpaceAGORA.SimulationEngine
const WALKER_FIXTURE = normpath(joinpath(@__DIR__, "..", "III_output", "walker_interlinks", "seed_20261003"))
const WALKER_COMPARISON = joinpath(WALKER_FIXTURE, "force_coupled_20.0s", "version_comparison")

walker_comparison_sma(state, mu) = inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)

function walker_comparison_configuration(mode::Symbol)
    mode in (:legacy, :v5, :v6) || throw(ArgumentError("Mode must be legacy, v5, or v6."))
    modern = mode !== :legacy
    @assert modern == isdefined(WALKER_MODEL, :InterLinkModel)
    project = mode === :legacy ? "2_SpaceAGORA.jl" : "8_SpaceAGORA.jl-main_$(mode)"
    @assert dirname(dirname(pathof(SpaceAGORA))) == normpath(joinpath(@__DIR__, "..", "..", project))
    metadata = only(eachrow(CSV.read(joinpath(WALKER_FIXTURE, "constellation.csv"), DataFrame)))
    @assert metadata.satellites == 4000 && metadata.seed == 20261003
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:metadata.satellites) do index
        plane, slot = divrem(index - 1, metadata.per_plane)
        bus = Link(root=true, m=227.0)
        initial_condition = InitialCondition(a=planet.Rp_e + metadata.altitude, i=metadata.inclination,
            Ω=mod(metadata.raan_offset + 360plane / metadata.planes, 360),
            ν=mod(metadata.anomaly_offset + 360slot / metadata.per_plane +
                360metadata.phasing * plane / metadata.satellites, 360))
        terminal_config = modern ? (n_terminal=1,) : NamedTuple()
        SpacecraftModel(; root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial_condition, id=index, terminal_config...)
    end
    effectors = (InverseSquaredGravityModel(),)
    extra_config = NamedTuple()
    if modern
        laser = InterLinkModel()
        for helper in 2:metadata.satellites
            register_candidate!(laser, spacecraft, (1, 1), (helper, 1);
                parameters=InterLinkParameters(P=10_000.0, B=100.0, range=metadata.range_m))
        end
        extra_config = (interlink_model=laser, scheduling_policy_model=SchedulingPolicyModel(:gve_sma; target_idx=1))
    else
        laser = OpenCavityLaserLinkModel(target_idx=1, helper_indices=collect(2:metadata.satellites),
            range_m=metadata.range_m, power_w=10_000.0, magnification=100.0, beta=1.0, eta=1.0, schedule=:gve_sma)
        effectors = (effectors..., laser)
    end
    args = SimulationConfiguration(;
        simulation_settings=SimulationSettings(results=false, generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=20.0, orientation_sim=false, data_rate=20.0),
        environment_model=make_no_gram_environment(; planet),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-11, dt_max_orbit=5.0),
        solver_config=WALKER_MODEL.SolverConfig(solver_mode=:tsit5), extra_config...)
    return args, laser
end

function run_walker_comparison(mode::Symbol, phase::String, reference_initial)
    args, laser = walker_comparison_configuration(mode)
    initial = WALKER_ENGINE.build_initial_conditions(args)
    mu = args.environment_model.planet.μ
    sample_times = collect(0.0:1.0:20.0)
    samples = SavedValues(Float64, Vector{Float64})
    sampler = SavingCallback((state, time, integrator) ->
        [walker_comparison_sma(current, mu) for current in state.sc], samples;
        saveat=sample_times, save_everystep=false, save_end=false)
    history = DataFrame(time_s=Float64[], active_helper=Int[])
    extra_callbacks = (sampler,)
    if mode === :legacy
        function record_schedule!(integrator)
            push!(history, (Float64(integrator.t), laser.active_helper_idx))
            DiffEqBase.u_modified!(integrator, false)
        end
        recorder = DiffEqBase.DiscreteCallback((state, time, integrator) -> true, record_schedule!;
            initialize=(callback, state, time, integrator) -> record_schedule!(integrator),
            save_positions=(false, false))
        extra_callbacks = (laser_impulse_callback(laser, LaserImpulseTracker(), 227.0),
            laser_link_scheduler_callback(laser), recorder, sampler)
    end
    println("Running ", mode, " ", phase, ": 4000 spacecraft, 20 s; package=", pathof(SpaceAGORA))
    flush(stdout)
    measured = @timed withenv("SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false", "SPACEAGORA_SOLVER_DENSE" => "false") do
        run_simulation(args; isolate_state=false, return_solution=true, extra_callbacks, save_fields=SaveField[])
    end
    solution = measured.value
    if mode !== :legacy
        history = DataFrame(time_s=[sample.time for sample in laser.history],
            active_helper=[isempty(sample.active) ? 0 : only(sample.active)[2][1] for sample in laser.history])
    end
    states = DataFrame(time_s=Float64[], satellite=Int[], x_m=Float64[], y_m=Float64[], z_m=Float64[],
        vx_m_s=Float64[], vy_m_s=Float64[], vz_m_s=Float64[], mass_kg=Float64[], sma_m=Float64[])
    for (time, state) in ((0.0, initial), (20.0, solution.u[end]))
        for (satellite, current) in enumerate(state.sc)
            push!(states, (time, satellite, current.pos..., current.vel..., current.mass, walker_comparison_sma(current, mu)))
        end
    end
    initial_sma = states.sma_m[1:4000]
    final_sma = states.sma_m[4001:8000]
    @testset "Walker $(mode) $(phase)" begin
        @test string(solution.retcode) == "Success"
        @test solution.t[end] == 20.0
        @test length(solution.u[end].sc) == 4000
        @test all(isfinite, solution.u[end])
        @test states[1:4000, names(reference_initial)] == reference_initial
        @test samples.t == sample_times
        @test all(all(isfinite, sample) for sample in samples.saveval)
        @test all(isapprox.(samples.saveval[end], final_sma; atol=1e-5, rtol=0))
        @test nrow(history) == solution.destats.naccept + 1
        @test history.time_s[1] == 0.0 && history.time_s[end] == 20.0
        @test all(diff(history.time_s) .> 0)
        @test all(1 .< history.active_helper .<= 4000)
        @test final_sma[1] > initial_sma[1]
    end
    summary = (implementation=String(mode), phase=phase, satellites=4000, candidates=3999,
        seed=20261003, duration_s=20.0, dtmax_s=5.0, reltol=1e-9, abstol=1e-11,
        gravity="two_body", solver="Tsit5", elapsed_s=measured.time, allocated_bytes=measured.bytes,
        gc_s=measured.gctime, accepted_steps=solution.destats.naccept, rejected_steps=solution.destats.nreject,
        schedule_samples=nrow(history), initial_target_sma_m=initial_sma[1], final_target_sma_m=final_sma[1],
        target_delta_sma_m=final_sma[1] - initial_sma[1],
        target_integrated_laser_delta_sma_m=mode === :legacy ? missing : solution.u[end].sc[1].laser_delta_sma,
        retcode=string(solution.retcode), julia_threads=Threads.nthreads(), source_file=pathof(SpaceAGORA))
    mkpath(WALKER_COMPARISON)
    CSV.write(joinpath(WALKER_COMPARISON, "$(mode)_initial_final_states.csv"), states)
    CSV.write(joinpath(WALKER_COMPARISON, "$(mode)_schedule.csv"), history)
    CSV.write(joinpath(WALKER_COMPARISON, "$(mode)_sma.csv"),
        DataFrame(time_s=repeat(samples.t; inner=4000), satellite=repeat(1:4000; outer=length(samples.t)),
            sma_m=reduce(vcat, samples.saveval)))
    @printf("%s %s: %.6f s, %d accepted / %d rejected, target SMA change %.9f m\n",
        String(mode), phase, measured.time, solution.destats.naccept, solution.destats.nreject, summary.target_delta_sma_m)
    return summary
end

function benchmark_walker_comparison(mode::Symbol; warm_runs::Int=5)
    reference = CSV.read(joinpath(WALKER_FIXTURE, "force_coupled_20.0s", "trajectory.csv"), DataFrame)
    reference_initial = reference[reference.time_s .== 0.0,
        [:time_s, :satellite, :x_m, :y_m, :z_m, :vx_m_s, :vy_m_s, :vz_m_s, :mass_kg, :sma_m]]
    timings = DataFrame()
    for trial in 0:warm_runs
        phase = iszero(trial) ? "cold" : "warm_$(trial)"
        push!(timings, run_walker_comparison(mode, phase, reference_initial))
        CSV.write(joinpath(WALKER_COMPARISON, "$(mode)_timings.csv"), timings)
    end
    return timings
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = Symbol(isempty(ARGS) ? "v6" : ARGS[1])
    warm_runs = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
    benchmark_walker_comparison(mode; warm_runs)
end