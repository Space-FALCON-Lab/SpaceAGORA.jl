#!/usr/bin/env julia
#
# Runs each MPC firing schedule in Input/Schedule_v2/
# (schedule_*.csv) against the new initial_conditions.csv in that same folder,
# one scheduler at a time. Reuses run_nihal_firing_plan_v1.jl's setup (package
# loading, spacecraft/plotting helpers) — including it here does not run its
# main() because PROGRAM_FILE points at this file instead.
#
# New initial_conditions.csv columns: satellite, shell_index, plane_index,
# slot_index, a_km, e, i_deg, raan_deg, omega_deg, M_deg (full osculating
# elements, incl. nonzero eccentricity — unlike the old altitude/true-anomaly
# format). Schedules are wide-format: one "<helper>_to_<target>" column per
# helper, exactly one helper active per row, target identified by the shared
# "_to_<target>" suffix.

include(joinpath(@__DIR__, "run_nihal_firing_plan_v1.jl"))

const NEW_FIRING_PLAN_DIR = joinpath(@__DIR__, "Input", "Schedule_v2")
const NEW_RESULTS_DIR     = joinpath(@__DIR__, "Output", "Schedule_v2")
const NEW_IC_PATH         = joinpath(NEW_FIRING_PLAN_DIR, "initial_conditions.csv")

# --------- 1. Parse a "<helper>_to_<target>" wide-format schedule ---------
function _load_new_firing_plan(ic_df::DataFrame, sched_path::String)
    sched_df = CSV.read(sched_path, DataFrame)

    helper_cols = names(sched_df)[2:end]
    col_matches = match.(r"^(\d+)_to_(\d+)$", helper_cols)
    any(isnothing, col_matches) && error("Schedule columns in $sched_path must match '<helper>_to_<target>'.")
    col_helper_ids = [parse(Int, m.captures[1]) for m in col_matches]
    target_ids     = unique(parse(Int, m.captures[2]) for m in col_matches)
    length(target_ids) == 1 || error("Schedule $sched_path references multiple targets: $target_ids.")
    target_csv_id  = target_ids[1]

    helper_csv_ids = sort(unique(col_helper_ids))
    helper_slot    = Dict(cid => k for (k, cid) in enumerate(helper_csv_ids))
    ordered_csv_ids = vcat([target_csv_id], helper_csv_ids)

    n_intervals     = nrow(sched_df)
    schedule_starts = Float64.(sched_df.start_time_s)
    schedule_active = falses(n_intervals, length(helper_csv_ids))
    for (ci, cid) in enumerate(col_helper_ids)
        slot = helper_slot[cid]
        col  = sched_df[!, helper_cols[ci]]
        schedule_active[:, slot] .= (col .== 1.0)
    end
    # Schedule only gives interval start times; step is uniform, so the mission
    # ends one more step past the last recorded interval start.
    dt = n_intervals >= 2 ? (schedule_starts[end] - schedule_starts[1]) / (n_intervals - 1) : 0.0
    mission_time_s = schedule_starts[end] + dt

    return (
        ordered_csv_ids=ordered_csv_ids, ic_df=ic_df,
        schedule_starts=schedule_starts, schedule_active=BitMatrix(schedule_active),
        mission_time_s=mission_time_s, n_helpers=length(helper_csv_ids),
    )
end

# --------- 2. Build SimulationConfigurations (controlled + no-laser reference) ---------
function _new_plan_spacecraft(plan)
    spacecraft = SpacecraftModel[]
    for (pos, cid) in enumerate(plan.ordered_csv_ids)
        row = only(findall(==(cid), plan.ic_df.satellite))
        r   = plan.ic_df[row, :]
        ν_deg = _mean_to_true_anomaly_deg(r.M_deg, r.e)
        push!(spacecraft, _spacecraft(pos, MASS_KG, InitialCondition(
            r.a_km * 1e3, r.e, r.i_deg, r.omega_deg, r.raan_deg, ν_deg,
        )))
    end
    return spacecraft
end

function _build_new_firing_plan_config(plan, results_dir::String)
    planet = make_no_gram_planet(:earth)
    spacecraft = _new_plan_spacecraft(plan)

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
            results=true, verbose=false, results_directory=results_dir,
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

function _build_new_reference_config(plan, results_dir::String)
    planet = make_no_gram_planet(:earth)
    spacecraft = _new_plan_spacecraft(plan)

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=joinpath(results_dir, "reference"),
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

# --------- 3. Run each scheduler in Input/Schedule_v2/, one at a time ---------
function main_new()
    ic_df = CSV.read(NEW_IC_PATH, DataFrame)
    ic_df.satellite = Int.(round.(ic_df.satellite))

    schedule_paths = sort(filter(
        f -> occursin(r"^schedule_.+\.csv$", basename(f)),
        readdir(NEW_FIRING_PLAN_DIR; join=true),
    ))
    isempty(schedule_paths) && error("No schedule_*.csv files found in $NEW_FIRING_PLAN_DIR.")

    for sched_path in schedule_paths
        name = match(r"^schedule_(.+)\.csv$", basename(sched_path)).captures[1]
        println("\n=============== Scheduler: $name ===============")

        plan = _load_new_firing_plan(ic_df, sched_path)
        println(@sprintf("Loaded %d satellites (1 target + %d helpers), %d schedule intervals, mission time %.1f s (%.2f hours).",
            plan.n_helpers + 1, plan.n_helpers, length(plan.schedule_starts), plan.mission_time_s, plan.mission_time_s / 3600.0))

        results_dir = joinpath(NEW_RESULTS_DIR, name)
        args, laser_model = _build_new_firing_plan_config(plan, results_dir)

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
        @printf("[%s] controlled runtime: %.2f s  retcode=%s\n", name, elapsed, sol.retcode)

        feather_df = DataFrame(Arrow.Table(joinpath(results_dir, "simulation_results.feather")))
        N  = plan.n_helpers + 1
        mu = args.environment_model.planet.μ

        # --------- no-laser reference trajectory (same ICs + gravity model) ---------
        ref_args = _build_new_reference_config(plan, results_dir)
        ref_elapsed = @elapsed begin
            ref_sol = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
                run_simulation(ref_args; isolate_state=false, return_solution=true)
            end
        end
        @printf("[%s] reference runtime:  %.2f s  retcode=%s\n", name, ref_elapsed, ref_sol.retcode)
        ref_feather_df = DataFrame(Arrow.Table(joinpath(results_dir, "reference", "simulation_results.feather")))

        # --------- plots (same layout as run_nihal_firing_plan_v1.jl) ---------
        flat_sol    = _make_flat_sol_from_feather(feather_df, N)
        plot_result = (flat_sol=flat_sol, impulse_tracker=impulse_tracker, mu=mu)
        plot_opts   = OracleCase2Options(
            helpers=plan.n_helpers, mass_kg=MASS_KG, eta=ETA, beta=BETA,
            magnification=MAGNIFICATION, laser_power_w=LASER_POWER_W, laser_range_km=Inf,
        )
        img_dir = joinpath(results_dir, "images")
        plot_open_cavity_results(plot_result, plot_opts; IMG_DIR=img_dir, target_only=false)

        # --------- target orbital elements: controlled vs. no-laser reference ---------
        target_elements_series(df) = [
            rv2coe(
                SVector{3,Float64}(df[k, "sc1_pos_1"], df[k, "sc1_pos_2"], df[k, "sc1_pos_3"]),
                SVector{3,Float64}(df[k, "sc1_vel_1"], df[k, "sc1_vel_2"], df[k, "sc1_vel_3"]),
                mu,
            ) for k in 1:nrow(df)
        ]
        oe_laser    = target_elements_series(feather_df)
        oe_no_laser = target_elements_series(ref_feather_df)
        length(oe_laser) == length(oe_no_laser) || error(
            "Controlled ($(length(oe_laser))) and reference ($(length(oe_no_laser))) time grids differ in length.")

        a_laser,     a_no_laser     = getfield.(oe_laser, :a), getfield.(oe_no_laser, :a)
        e_laser,     e_no_laser     = getfield.(oe_laser, :e), getfield.(oe_no_laser, :e)
        i_laser,     i_no_laser     = rad2deg.(getfield.(oe_laser, :i)), rad2deg.(getfield.(oe_no_laser, :i))
        Omega_laser, Omega_no_laser = rad2deg.(getfield.(oe_laser, :Ω)), rad2deg.(getfield.(oe_no_laser, :Ω))
        omega_laser, omega_no_laser = rad2deg.(getfield.(oe_laser, :ω)), rad2deg.(getfield.(oe_no_laser, :ω))

        Δa_km    = (a_laser .- a_no_laser) ./ 1e3
        Δe       = e_laser .- e_no_laser
        T_ref    = 2π * sqrt(a_no_laser[1]^3 / mu)
        x_orbits = feather_df.time ./ T_ref

        plt_da_ref = plot(x_orbits, Δa_km, legend=false,
            title="Target Δa vs. no-laser reference ($name)", xlabel="t / T_ref", ylabel="Δa, km")
        savefig(plt_da_ref, joinpath(img_dir, "target_da_vs_no_laser_reference_$name.png"))

        plt_de_ref = plot(x_orbits, Δe, legend=false,
            title="Target Δe vs. no-laser reference ($name)", xlabel="t / T_ref", ylabel="Δe")
        savefig(plt_de_ref, joinpath(img_dir, "target_de_vs_no_laser_reference_$name.png"))

        @printf("[%s] final target Δa vs. no-laser reference: %+.4f km   Δe: %+.6e\n", name, Δa_km[end], Δe[end])

        oe_csv = DataFrame(
            time=feather_df.time,
            a_laser=a_laser, e_laser=e_laser, i_laser=i_laser, Omega_laser=Omega_laser, omega_laser=omega_laser,
            a_no_laser=a_no_laser, e_no_laser=e_no_laser, i_no_laser=i_no_laser, Omega_no_laser=Omega_no_laser, omega_no_laser=omega_no_laser,
        )
        mkpath(results_dir)
        CSV.write(joinpath(results_dir, "oe_timeseries_$name.csv"), oe_csv)

        println("[$name] output directory: $results_dir")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_new()
end
