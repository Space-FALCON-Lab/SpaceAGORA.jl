#!/usr/bin/env julia
#
# Runs each of the 8 MPC firing schedules in
# Input/Schedule_v3/{corrected,non_corrected}/*_schedule.csv
# against the same initial_conditions.csv used by run_nihal_firing_plan_v2.jl
# (identical wide "<helper>_to_<target>" schedule format and satellite
# numbering), one scheduler at a time. Reuses run_nihal_firing_plan_v2.jl's
# setup (package loading, plan parsing, config builders) — including it here
# does not run its main_new() because PROGRAM_FILE points at this file instead.

include(joinpath(@__DIR__, "run_nihal_firing_plan_v2.jl"))

const NEWNEW_FIRING_PLAN_DIR = joinpath(@__DIR__, "Input", "Schedule_v3")
const NEWNEW_RESULTS_DIR     = joinpath(@__DIR__, "Output", "Schedule_v3")
const NEWNEW_VARIANTS        = ("corrected", "non_corrected")

# Same as _build_new_reference_config, but with save_csv=true so the reference
# run also produces csv+feather+toml (not just feather+toml).
function _build_newnew_reference_config(plan, results_dir::String)
    planet = make_no_gram_planet(:earth)
    spacecraft = _new_plan_spacecraft(plan)

    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=joinpath(results_dir, "reference"),
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

# --------- Run each of the 8 schedulers (corrected/ and non_corrected/), one at a time ---------
function main_newnew()
    ic_df = CSV.read(NEW_IC_PATH, DataFrame)
    ic_df.satellite = Int.(round.(ic_df.satellite))

    for variant in NEWNEW_VARIANTS
        variant_dir = joinpath(NEWNEW_FIRING_PLAN_DIR, variant)
        schedule_paths = sort(filter(
            f -> occursin(r"^.+_schedule\.csv$", basename(f)),
            readdir(variant_dir; join=true),
        ))
        isempty(schedule_paths) && error("No *_schedule.csv files found in $variant_dir.")

        for sched_path in schedule_paths
            name = match(r"^(.+)_schedule\.csv$", basename(sched_path)).captures[1]
            println("\n=============== [$variant] Scheduler: $name ===============")

            plan = _load_new_firing_plan(ic_df, sched_path)
            println(@sprintf("Loaded %d satellites (1 target + %d helpers), %d schedule intervals, mission time %.1f s (%.2f hours).",
                plan.n_helpers + 1, plan.n_helpers, length(plan.schedule_starts), plan.mission_time_s, plan.mission_time_s / 3600.0))

            results_dir = joinpath(NEWNEW_RESULTS_DIR, variant, name)
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
            @printf("[%s/%s] controlled runtime: %.2f s  retcode=%s\n", variant, name, elapsed, sol.retcode)

            feather_df = DataFrame(Arrow.Table(joinpath(results_dir, "simulation_results.feather")))
            N  = plan.n_helpers + 1
            mu = args.environment_model.planet.μ

            # --------- no-laser reference trajectory (same ICs + gravity model) ---------
            ref_args = _build_newnew_reference_config(plan, results_dir)
            ref_elapsed = @elapsed begin
                ref_sol = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
                    run_simulation(ref_args; isolate_state=false, return_solution=true)
                end
            end
            @printf("[%s/%s] reference runtime:  %.2f s  retcode=%s\n", variant, name, ref_elapsed, ref_sol.retcode)
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
                title="Target Δa vs. no-laser reference ($variant/$name)", xlabel="t / T_ref", ylabel="Δa, km")
            savefig(plt_da_ref, joinpath(img_dir, "target_da_vs_no_laser_reference.png"))

            plt_de_ref = plot(x_orbits, Δe, legend=false,
                title="Target Δe vs. no-laser reference ($variant/$name)", xlabel="t / T_ref", ylabel="Δe")
            savefig(plt_de_ref, joinpath(img_dir, "target_de_vs_no_laser_reference.png"))

            @printf("[%s/%s] final target Δa vs. no-laser reference: %+.4f km   Δe: %+.6e\n", variant, name, Δa_km[end], Δe[end])

            oe_csv = DataFrame(
                time=feather_df.time,
                a_laser=a_laser, e_laser=e_laser, i_laser=i_laser, Omega_laser=Omega_laser, omega_laser=omega_laser,
                a_no_laser=a_no_laser, e_no_laser=e_no_laser, i_no_laser=i_no_laser, Omega_no_laser=Omega_no_laser, omega_no_laser=omega_no_laser,
            )
            mkpath(results_dir)
            CSV.write(joinpath(results_dir, "oe_timeseries.csv"), oe_csv)

            println("[$variant/$name] output directory: $results_dir")
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_newnew()
end
