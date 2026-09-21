#!/usr/bin/env julia
#
# Runs the Nihal MPC firing plan (Input/Schedule_v1/initial_conditions.csv +
# mpc_schedule.csv) through the SpaceAGORA OpenCavityLaserLinkModel, using the
# schedule's exact per-interval helper on/off table (supports simultaneous
# multi-helper firing) instead of one of the built-in heuristic schedules.

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
println("Loading SpaceAGORA (precompilation may take a few minutes on first run)...")
include(joinpath(REPO_ROOT, "examples", "common.jl"))

using Arrow
using CSV
using DataFrames
using LinearAlgebra
using Printf
using DiffEqBase
using StaticArrays
using Plots
using .SimulationModel

const FIRING_PLAN_DIR = joinpath(@__DIR__, "Input", "Schedule_v1")
const RESULTS_DIR     = joinpath(@__DIR__, "Output", "Schedule_v1")
const MASS_KG          = 300.0  # matches Nihal's actual src/ run (verified against his ground-truth history), not the 227 kg paper default
const LASER_POWER_W    = 10_000.0
const MAGNIFICATION    = 100.0
const BETA             = 1.0
const ETA              = 1.0

# Minimal stand-in for ORACLE/run_case2_laser_links.jl's option struct — only
# used so the ported plotting pipeline (functions/9_Runners.jl) type-checks;
# only the fields it actually reads (helpers, mass_kg, eta, laser_range_km) matter.
Base.@kwdef struct OracleCase2Options
    helpers::Int = 10
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    target_nu_deg::Float64 = 0.0
    target_ecc::Float64 = 0.0
    orbits::Float64 = 10.0
    schedule::Symbol = :external_schedule
    laser_range_km::Float64 = 200.0
    laser_power_w::Float64 = 10_000.0
    magnification::Float64 = 100.0
    beta::Float64 = 1.0
    eta::Float64 = 1.0
    mass_kg::Float64 = 227.0
    dt_max_s::Float64 = 10.0
    paper_grid::Bool = false
    feather_only::Bool = false
    output_dir::String = "output"
    timeseries_points::Int = 1001
    animate::Bool = false
end

include(joinpath(REPO_ROOT, "ORACLE", "functions", "0_Spacecraft.jl"))
include(joinpath(REPO_ROOT, "ORACLE", "functions", "4_Diagnostics.jl"))
include(joinpath(REPO_ROOT, "ORACLE", "functions", "5_OE_Converters.jl"))
include(joinpath(REPO_ROOT, "ORACLE", "functions", "6_OE_and_dv_in_RTN.jl"))
include(joinpath(REPO_ROOT, "ORACLE", "functions", "7_Plots.jl"))
include(joinpath(REPO_ROOT, "ORACLE", "functions", "9_Runners.jl"))
const DT_MAX_S         = 10.0

# --------- 1. Parse initial_conditions.csv + mpc_schedule.csv ---------
function _load_firing_plan(dir::String)
    ic_df    = CSV.read(joinpath(dir, "initial_conditions.csv"), DataFrame)
    sched_df = CSV.read(joinpath(dir, "mpc_schedule.csv"), DataFrame)

    target_rows = filter(:is_target => identity, ic_df)
    nrow(target_rows) == 1 || error("Expected exactly one target satellite in initial_conditions.csv, found $(nrow(target_rows)).")
    target_csv_id  = target_rows.satellite[1]
    helper_csv_ids = sort(ic_df.satellite[ic_df.satellite .!= target_csv_id])

    # spacecraft build order: position 1 = target, positions 2..N+1 = helpers (sorted by csv id)
    ordered_csv_ids = vcat([target_csv_id], helper_csv_ids)
    helper_slot = Dict(cid => k for (k, cid) in enumerate(helper_csv_ids))  # csv id -> 1-based helper slot

    n_intervals   = length(unique(sched_df.interval))
    interval_pos  = Dict(iv => k for (k, iv) in enumerate(sort(unique(sched_df.interval))))
    schedule_starts = Vector{Float64}(undef, n_intervals)
    schedule_active = falses(n_intervals, length(helper_csv_ids))
    for row in eachrow(sched_df)
        row.target == target_csv_id || error("mpc_schedule.csv edge target $(row.target) does not match initial_conditions.csv target $(target_csv_id).")
        k = interval_pos[row.interval]
        schedule_starts[k] = row.start_time_s
        schedule_active[k, helper_slot[row.source]] = (row.u_full_power == 1)
    end
    mission_time_s = maximum(sched_df.end_time_s)

    return (
        ordered_csv_ids=ordered_csv_ids, ic_df=ic_df,
        schedule_starts=schedule_starts, schedule_active=BitMatrix(schedule_active),
        mission_time_s=mission_time_s, n_helpers=length(helper_csv_ids),
    )
end

# --------- 2. Build SimulationConfiguration from the parsed plan ---------
function _build_firing_plan_config(plan)
    planet = make_no_gram_planet(:earth)

    spacecraft = SpacecraftModel[]
    for (pos, cid) in enumerate(plan.ordered_csv_ids)
        row = only(findall(==(cid), plan.ic_df.satellite))
        r   = plan.ic_df[row, :]
        radius_m = planet.Rp_e + r.altitude_km * 1e3
        push!(spacecraft, _spacecraft(pos, MASS_KG, InitialCondition(
            radius_m, 0.0, r.inclination_deg, r.raan_deg, 0.0, r.true_anomaly_deg,
        )))
    end

    laser_model = OpenCavityLaserLinkModel(
        target_idx=1,
        helper_indices=collect(2:(plan.n_helpers + 1)),
        range_m=Inf,   # trust the MPC schedule fully — no extra range cutoff
        power_w=LASER_POWER_W,
        magnification=MAGNIFICATION,
        beta=BETA,
        eta=ETA,
        schedule=:external_schedule,
        schedule_starts=plan.schedule_starts,
        schedule_active=plan.schedule_active,
    )

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=RESULTS_DIR,
            generate_plots=false, save_csv=true, normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=plan.mission_time_s, orientation_sim=false,
            num_steps_to_save=1000, data_rate=max(10.0, plan.mission_time_s / 1000.0),
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(), laser_model)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12, abstol_orbit=1e-12, dt_max_orbit=DT_MAX_S,
        ),
        solver_config=SolverConfig(solver_mode=:tsit5),
    )
    return args, laser_model
end

# --------- 2b. Build a no-laser reference SimulationConfiguration (same ICs/gravity) ---------
function _build_reference_config(plan)
    planet = make_no_gram_planet(:earth)

    spacecraft = SpacecraftModel[]
    for (pos, cid) in enumerate(plan.ordered_csv_ids)
        row = only(findall(==(cid), plan.ic_df.satellite))
        r   = plan.ic_df[row, :]
        radius_m = planet.Rp_e + r.altitude_km * 1e3
        push!(spacecraft, _spacecraft(pos, MASS_KG, InitialCondition(
            radius_m, 0.0, r.inclination_deg, r.raan_deg, 0.0, r.true_anomaly_deg,
        )))
    end

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=joinpath(RESULTS_DIR, "reference"),
            generate_plots=false, save_csv=false, normalize=false,
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=plan.mission_time_s, orientation_sim=false,
            num_steps_to_save=1000, data_rate=max(10.0, plan.mission_time_s / 1000.0),
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
        ),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),  # no laser effector
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12, abstol_orbit=1e-12, dt_max_orbit=DT_MAX_S,
        ),
        solver_config=SolverConfig(solver_mode=:tsit5),
    )
    return args
end

# --------- 3. Run ---------
function main()
    plan = _load_firing_plan(FIRING_PLAN_DIR)
    println(@sprintf("Loaded %d satellites (1 target + %d helpers), %d schedule intervals, mission time %.1f s (%.2f hours).",
        plan.n_helpers + 1, plan.n_helpers, length(plan.schedule_starts), plan.mission_time_s, plan.mission_time_s / 3600.0))

    args, laser_model = _build_firing_plan_config(plan)

    impulse_tracker = LaserImpulseTracker()
    impulse_cb   = laser_impulse_callback(laser_model, impulse_tracker, MASS_KG)
    scheduler_cb = laser_link_scheduler_callback(laser_model)

    laser_fields = SaveField[
        SaveField(:dv_r_accumulated, (u, t, integrator) -> impulse_tracker.dv_R; per_satellite=false, column_prefix="dv_r_accumulated"),
        SaveField(:dv_t_accumulated, (u, t, integrator) -> impulse_tracker.dv_T; per_satellite=false, column_prefix="dv_t_accumulated"),
        SaveField(:dv_n_accumulated, (u, t, integrator) -> impulse_tracker.dv_N; per_satellite=false, column_prefix="dv_n_accumulated"),
        SaveField(:laser_active_helper_count, (u, t, integrator) -> Float64(length(laser_model.active_helper_indices)); per_satellite=false, column_prefix="laser_active_helper_count"),
    ]
    all_save_fields = vcat(SimulationModel.default_save_fields(args), laser_fields)

    elapsed = @elapsed begin
        result = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
            run_simulation(
                args; isolate_state=false, return_solution=true, return_solver_metadata=true,
                extra_callbacks=(impulse_cb, scheduler_cb), save_fields=all_save_fields,
            )
        end
    end
    sol = result.solution
    @printf("\nSimulation runtime: %.2f s (%.2f min)  retcode=%s\n", elapsed, elapsed / 60.0, sol.retcode)

    # --------- 4. Post-analysis: Δa/Δe/Δi/ΔRAAN per satellite ---------
    feather_df = DataFrame(Arrow.Table(joinpath(RESULTS_DIR, "simulation_results.feather")))
    N   = plan.n_helpers + 1
    mu  = args.environment_model.planet.μ

    function _elements_at(row)
        [_rv_to_elements(
            SVector{3,Float64}(feather_df[row, "sc$(i)_pos_1"], feather_df[row, "sc$(i)_pos_2"], feather_df[row, "sc$(i)_pos_3"]),
            SVector{3,Float64}(feather_df[row, "sc$(i)_vel_1"], feather_df[row, "sc$(i)_vel_2"], feather_df[row, "sc$(i)_vel_3"]),
            mu,
        ) for i in 1:N]
    end
    oe0 = _elements_at(1)
    oef = _elements_at(nrow(feather_df))

    summary_rows = DataFrame(
        csv_satellite_id=plan.ordered_csv_ids,
        internal_index=1:N,
        is_target=[i == 1 for i in 1:N],
        da_m=[oef[i].a - oe0[i].a for i in 1:N],
        de=[oef[i].e - oe0[i].e for i in 1:N],
        di_deg=[rad2deg(oef[i].i - oe0[i].i) for i in 1:N],
        draan_deg=[rad2deg(oef[i].raan - oe0[i].raan) for i in 1:N],
    )

    println("\n=============== Δ Orbital Elements (start \u2192 end) ===============")
    for r in eachrow(summary_rows)
        tag = r.is_target ? "TARGET" : "helper"
        @printf("  sat csv#%-2d [%-6s]  Δa=% .4f m   Δe=% .4e   Δi=% .4e deg   ΔRAAN=% .4e deg\n",
            r.csv_satellite_id, tag, r.da_m, r.de, r.di_deg, r.draan_deg)
    end

    println("\n=============== Laser link stats ===============")
    println("  link_activation_count  : ", laser_model.link_activation_count)
    println("  active_link_step_count : ", laser_model.active_link_step_count)
    @printf("  dV_RTN tracker [m/s]    : R=% .6e  T=% .6e  N=% .6e\n",
        impulse_tracker.dv_R, impulse_tracker.dv_T, impulse_tracker.dv_N)

    mkpath(RESULTS_DIR)
    CSV.write(joinpath(RESULTS_DIR, "summary.csv"), summary_rows)

    # --------- 5. Plots (same layout as the prototype's output/images/) ---------
    flat_sol   = _make_flat_sol_from_feather(feather_df, N)
    plot_result = (flat_sol=flat_sol, impulse_tracker=impulse_tracker, mu=mu)
    plot_opts   = OracleCase2Options(
        helpers=plan.n_helpers, mass_kg=MASS_KG, eta=ETA, beta=BETA,
        magnification=MAGNIFICATION, laser_power_w=LASER_POWER_W, laser_range_km=Inf,
    )
    img_dir = joinpath(RESULTS_DIR, "images")
    plot_open_cavity_results(plot_result, plot_opts; IMG_DIR=img_dir, target_only=false)

    # --------- 6. Target Δa vs. a no-laser reference trajectory ---------
    # Reference: same ICs + gravity model, propagated with zero laser commands.
    ref_args = _build_reference_config(plan)
    ref_elapsed = @elapsed begin
        ref_sol = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
            run_simulation(ref_args; isolate_state=false, return_solution=true)
        end
    end
    @printf("Reference (no-laser) simulation runtime: %.2f s  retcode=%s\n",
        ref_elapsed, ref_sol.retcode)
    ref_feather_df = DataFrame(Arrow.Table(joinpath(RESULTS_DIR, "reference", "simulation_results.feather")))

    target_elements_series(df) = [
        _rv_to_elements(
            SVector{3,Float64}(df[k, "sc1_pos_1"], df[k, "sc1_pos_2"], df[k, "sc1_pos_3"]),
            SVector{3,Float64}(df[k, "sc1_vel_1"], df[k, "sc1_vel_2"], df[k, "sc1_vel_3"]),
            mu,
        ) for k in 1:nrow(df)
    ]
    oe_controlled = target_elements_series(feather_df)
    oe_reference  = target_elements_series(ref_feather_df)
    length(oe_controlled) == length(oe_reference) || error(
        "Controlled ($(length(oe_controlled))) and reference ($(length(oe_reference))) time grids differ in length.")

    a_controlled = getfield.(oe_controlled, :a)
    a_reference  = getfield.(oe_reference, :a)
    e_controlled = getfield.(oe_controlled, :e)
    e_reference  = getfield.(oe_reference, :e)

    Δa_km  = (a_controlled .- a_reference) ./ 1e3
    Δe     = e_controlled .- e_reference
    T_ref  = 2π * sqrt(a_reference[1]^3 / mu)
    x_orbits = feather_df.time ./ T_ref

    plt_da_ref = plot(x_orbits, Δa_km, legend=false,
        title="Target Δa vs. no-laser reference",
        xlabel="t / T_ref", ylabel="Δa, km")
    savefig(plt_da_ref, joinpath(img_dir, "target_da_vs_no_laser_reference.png"))
    @printf("Final target Δa vs. no-laser reference: %+.4f km\n", Δa_km[end])

    plt_de_ref = plot(x_orbits, Δe, legend=false,
        title="Target Δe vs. no-laser reference",
        xlabel="t / T_ref", ylabel="Δe")
    savefig(plt_de_ref, joinpath(img_dir, "target_de_vs_no_laser_reference.png"))
    @printf("Final target Δe vs. no-laser reference: %+.6e\n", Δe[end])

    println("\nOutput directory: $RESULTS_DIR")
    return (result=result, summary=summary_rows, laser_model=laser_model, impulse_tracker=impulse_tracker)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
