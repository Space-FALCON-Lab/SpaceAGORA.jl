"""
Odyssey P20 panel-control exercise with one named frozen atmosphere.

Run this file for the 90-degree/30-degree comparison. Include it to use
`OdysseySurrogateExample.run_case` or `compare_panel_caps` from a notebook.
"""
module OdysseySurrogateExample

using SpaceAGORA, GRAMSuite, StaticArrays, LinearAlgebra, TOML, CSV, Printf, Libdl
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CB = SM.SimulationCallbacks
const SP = SpaceAGORA.TelemetryVerification.SPICE
const AE = SM.DynamicEffectors.AerodynamicEffectors
const PRESET = "odyssey_p20_frozen_v1"
const VERSION = "1.0.0"
const ENTRY_UTC = "2001-11-07T11:45:45.690115Z"
const COMMON_TIME_S = 600.0

function validate_panel_cap(panel_cap_deg)
    cap = Float64(panel_cap_deg)
    isfinite(cap) && 0 < cap <= 90 || throw(ArgumentError("panel_cap_deg must be finite and in (0, 90]."))
    return cap
end

native_gram_loaded() = any(path -> occursin("libgram", lowercase(path)), Libdl.dllist())

function initial_time()
    # InitialTime stores Float32 seconds. Resolve the engine clock first, then
    # query the kernel at that same instant so the clock and state agree.
    return SM.InitialTime(year=2001, month=11, day=7, hour=11, minute=45,
        second=45.690115)
end

function furnish_assets!(assets)
    length(assets.kernels) == 4 || throw(ArgumentError("Odyssey P20 requires the four identified scenario kernels."))
    all(isfile, assets.kernels) || throw(ArgumentError("A required Odyssey kernel is missing."))
    isfile(assets.gravity) || throw(ArgumentError("The identified Mars gravity coefficients are missing."))
    SP.kclear()
    foreach(SP.furnsh, assets.kernels)
    return nothing
end

function check_kernel_pool(assets)
    SP.ktotal("ALL") == 4 || error("The Odyssey four-kernel pool changed during simulation.")
    actual = Set(realpath(SP.kdata(i, "ALL")[1]) for i in 1:4)
    actual == Set(realpath.(assets.kernels)) || error("Unexpected kernel in the Odyssey scenario.")
    native_gram_loaded() && error("This example must run without native GRAM.")
    return nothing
end

function configure_case(model, assets, panel_cap_deg)
    cap = validate_panel_cap(panel_cap_deg)
    planet = SM.Mars()
    epoch = initial_time()
    start_et = SM.ephemerides_time_seconds(epoch, SM.SpiceEphemeridesModel())
    state = SP.spkezr("-53", start_et, "J2000", "NONE", "499")[1] .* 1000.0
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7), panel_dims=(0.01, 1.945, 1.7),
        bus_mass=391.0, panel_mass_each=10.0, panel_offset_y=2.2725,
        prop_mass=50.0, id=1001,
        ic=SM.CartesianInitialCondition(state[1:3], state[4:6]))
    effectors = (
        SM.GravitationalHarmonicsModel(20, 20, assets.gravity, planet;
            coefficients_normalized=true, j2_source=:file_c20),
        SM.NBodyGravityModel(body_names=("sun", "earth", "moon", "jupiter"), primary_body_name="Mars", planet=planet),
        SM.SolarRadiationPressureModel(1.3, 3.74),
        SM.AerodynamicCoefficientfM(fixed_attitude_incidence=:attitude),
    )
    config = SM.AerobrakingEnergyDepletionConfig(
        guidance_modes=(:max_energy_depletion,), max_energy_submodes=(:heat_rate,),
        controlled_panel_links=(2, 3), max_alpha_rad=deg2rad(cap),
        min_alpha_rad=min(1e-4, deg2rad(cap)), heat_rate_limit_w_cm2=Inf)
    control_state = SM.AerobrakingEnergyDepletionState(num_sats=1)
    guidance = SM.AerobrakingEnergyDepletionGuidanceModel(config, control_state)
    control = SM.AerobrakingEnergyDepletionControlModel(config, control_state)
    base = make_example_config(planet=planet, spacecraft=spacecraft,
        mission_time=900.0, initial_time=epoch, dynamic_effectors=effectors,
        density_model=model, orientation_sim=false, keplerian=false, EI_km=250.0,
        verbose=false, results=false,
        solver_config=SM.SolverConfig(solver_mode=:tsit5, maxiters=20_000_000))
    environment = SM.EnvironmentModel(planet=planet, EI=250.0, density_model=model,
        ephemerides_model=SM.SpiceEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, thermal_contact=false, planet=planet),
        topography=false, wind=true)
    args = SM.SimConfig._with_configuration(base;
        environment_model=environment,
        guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[3.0]),
        control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.1]),
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false,
            normalize=false, save_csv=false, checkpoint_enabled=false),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-7, abstol_orbit=1e-9,
            reltol_atmosphere=1e-7, abstol_atmosphere=1e-9,
            dt_max=0.2, dt_max_orbit=0.2, dt_max_atmosphere=0.2))
    return (; args, start_et, initial_state=state, panel_cap_deg=cap)
end

function engine_config()
    SE.SimulationEngineConfig(
        parallel=SE.ParallelConfig(effector_parallel_mode="serial", rhs_batch_parallel_mode="serial",
            density_callback_parallel_mode="serial", thermal_callback_parallel_mode="serial", control_callback_parallel_mode="serial"),
        runtime_policy=SE.RuntimePolicyConfig(srp_ephemeris_cache=false, nbody_ephemeris_cache=false,
            planet_frame_cache=false, spice_rhs_memo=true),
        env_overrides=Dict("SPACEAGORA_DENSITY_FREEZE_PER_STEP"=>"0", "SPACEAGORA_VACUUM_GRAM_CACHE"=>"0",
            "SPACEAGORA_ENTRY_TARGET_COUNT"=>"0", "SPACEAGORA_RHS_CALIBRATE"=>"off",
            "SPACEAGORA_GRAM_TRACK_CACHE"=>"off", "SPACEAGORA_VISUALIZATION"=>"0",
            "SPACEAGORA_SOLVER_SAVE_EVERYSTEP"=>"1", "SPACEAGORA_SOLVER_SAVE_ON"=>"1",
            "SPACEAGORA_SOLVER_SAVE_START"=>"1", "SPACEAGORA_SOLVER_SAVE_END"=>"1"))
end

function sample_record(u, p, t)
    t = Float64(t)
    sc = u.sc[1]
    spacecraft = p.args.dynamics_model.spacecraft[1]
    frame = SE.sample_planet_frame(sc, p, 1, t)
    atmosphere = SE.sample_atmosphere(sc, p, 1, t; write_buffers=false)
    state = SE.build_state_sample(sc, spacecraft, false)
    env = SpaceAGORA.EnvironmentSample(p.args.environment_model.planet;
        planet_frame=frame, atmosphere=atmosphere)
    _, _, drag, _, _ = AE._aero_pure_wrench(:fm, state, env, nothing, :attitude)
    rates = copy(CB._compute_stage_heat_rates!(p, sc, 1, t))
    control = p.args.control_model.control_effectors[1]
    left_incidence = AE._aero_link_angles(spacecraft, spacecraft.links[2], 1,
        false, SVector(0.0,0.0,0.0), 0.0, :attitude)[1]
    right_incidence = AE._aero_link_angles(spacecraft, spacecraft.links[3], 1,
        false, SVector(0.0,0.0,0.0), 0.0, :attitude)[1]
    return (elapsed_s=t, height_km=frame.alt_m/1000, latitude_deg=rad2deg(frame.lat_rad),
        longitude_deg=rad2deg(frame.lon_rad), density_kg_m3=atmosphere.rho_kg_m3,
        temperature_K=atmosphere.temperature_k, drag_N=norm(drag),
        panel_left_deg=rad2deg(spacecraft.links[2].α), panel_right_deg=rad2deg(spacecraft.links[3].α),
        panel_left_incidence_deg=rad2deg(left_incidence), panel_right_incidence_deg=rad2deg(right_incidence),
        mode=string(control.state.selected_mode[1]),
        heat_left_W_cm2=rates[2], heat_right_W_cm2=rates[3],
        heat_left_J_cm2=Float64(sc.heat_loads[2]), heat_right_J_cm2=Float64(sc.heat_loads[3]),
        x_m=Float64(sc.pos[1]), y_m=Float64(sc.pos[2]), z_m=Float64(sc.pos[3]),
        vx_m_s=Float64(sc.vel[1]), vy_m_s=Float64(sc.vel[2]), vz_m_s=Float64(sc.vel[3]))
end

"""
    run_case(; panel_cap_deg=90, grid_file="", offline=false, output_dir, assets=nothing,
               allow_unreleased=false)

Run one bounded P20 passage. `panel_cap_deg` changes the active panel controller.
The thermal threshold is infinite: this exercise tests commanded panel incidence,
not closed-loop thermal protection or mission accuracy. `assets` is normally
resolved by `odyssey_surrogate_assets`; explicit values support controlled tests.
"""
function run_case(; panel_cap_deg=90.0, grid_file="", offline=false,
        output_dir=joinpath(pwd(), "odyssey_surrogate_$(panel_cap_deg)"),
        assets=nothing, allow_unreleased=false)
    validate_panel_cap(panel_cap_deg)
    ispath(output_dir) && throw(ArgumentError("Choose a new output directory: $output_dir"))
    model = SpaceAGORA.surrogate_preset_model(PRESET; version=VERSION, planet="Mars",
        file=grid_file, offline=offline, allow_unreleased=allow_unreleased)
    assets === nothing && (assets = SpaceAGORA.odyssey_surrogate_assets(; offline=offline))
    furnish_assets!(assets)
    setup = configure_case(model, assets, panel_cap_deg)
    check_kernel_pool(assets)
    mkpath(output_dir)
    rows = NamedTuple[]
    function capture!(integrator)
        t = Float64(integrator.t)
        if isempty(rows) || t - rows[end].elapsed_s >= 1.0 - 1e-9
            push!(rows, sample_record(integrator.u, integrator.p, t))
        end
        # Preserve invalidation from guidance/control callbacks: articulating
        # mutable link geometry can change the RHS without changing the ODE u.
    end
    recorder = CB.DiscreteCallback((u,t,i)->true, capture!;
        initialize=(c,u,t,i)->capture!(i), save_positions=(false,false))
    exit_found = Ref(false)
    function exit_passage!(integrator)
        integrator.t > 1 || return nothing
        exit_found[] = true
        CB.terminate!(integrator)
    end
    exit_condition(u,t,i) = SE.sample_planet_frame(u.sc[1], i.p, 1, Float64(t)).alt_m - 250_000.0
    exit_callback = CB.ContinuousCallback(exit_condition, exit_passage!, nothing;
        rootfind=CB.DiffEqBase.SciMLBase.RightRootFind, abstol=1e-6, reltol=0.0, save_positions=(true,true))
    @printf("Running the %g-degree panel-cap passage with the frozen Odyssey atmosphere.\n", panel_cap_deg)
    flush(stdout)
    result = nothing
    wall_s = @elapsed result = SE.run_simulation(engine_config(), setup.args;
        return_solution=true, return_solver_metadata=true, save_fields=[],
        extra_callbacks=(recorder, exit_callback), visualization=false)
    sol = result.solution
    check_kernel_pool(assets)
    exit_found[] && SE.SciMLBase.successful_retcode(sol.retcode) || error("Passage did not reach its outbound 250 km event.")
    sol.t[end] > COMMON_TIME_S || error("Passage ended before the comparison time.")
    endpoint = sample_record(sol.u[end], sol.prob.p, sol.t[end])
    common = sample_record(sol(COMMON_TIME_S), sol.prob.p, COMMON_TIME_S)
    push!(rows, endpoint)
    active = filter(r->r.elapsed_s >= 1.0, rows)
    all(r->r.mode == "max_energy_depletion", active) || error("Guidance was inactive.")
    all(r->isapprox(r.panel_left_deg, panel_cap_deg; atol=1e-10, rtol=0) &&
        isapprox(r.panel_right_deg, panel_cap_deg; atol=1e-10, rtol=0), active) || error("Panel cap was not applied.")
    all(r->isapprox(r.panel_left_incidence_deg, panel_cap_deg; atol=1e-10, rtol=0) &&
        isapprox(r.panel_right_incidence_deg, panel_cap_deg; atol=1e-10, rtol=0), active) || error("Applied panel geometry does not match the command.")
    all(r->100 <= r.height_km <= 260 && 40 <= r.latitude_deg <= 90, rows) || error("Recorded trajectory left the preset domain.")
    CSV.write(joinpath(output_dir, "trajectory.csv"), rows)
    summary = Dict{String,Any}("status"=>"completed_native_free_control_exercise", "panel_cap_deg"=>Float64(panel_cap_deg),
        "start_et_tdb_s"=>setup.start_et, "requested_start_utc"=>ENTRY_UTC, "start_utc"=>SP.et2utc(setup.start_et,"ISOC",9),
        "wall_s"=>wall_s, "solver_retcode"=>string(sol.retcode), "solver"=>"Tsit5", "exit_elapsed_s"=>Float64(sol.t[end]),
        "initial_state_m_m_s"=>setup.initial_state, "native_gram_loaded"=>false,
        "atmosphere"=>SpaceAGORA.atmosphere_provenance(model), "scenario_assets"=>assets.provenance,
        "panel_cap_active"=>true, "thermal_feedback_tested"=>false,
        "comparison_time_s"=>COMMON_TIME_S,
        "common_state_m_m_s"=>[common.x_m, common.y_m, common.z_m, common.vx_m_s, common.vy_m_s, common.vz_m_s],
        "common_panel_heat_load_J_cm2"=>[common.heat_left_J_cm2, common.heat_right_J_cm2],
        "maximum_recorded_drag_N"=>maximum(r.drag_N for r in rows),
        "minimum_recorded_height_km"=>minimum(r.height_km for r in rows),
        "maximum_recorded_latitude_deg"=>maximum(r.latitude_deg for r in rows),
        "scope"=>"Algorithm usability exercise with a frozen atmosphere; numerical release limits remain unset.")
    open(joinpath(output_dir, "summary.toml"), "w") do io; TOML.print(io, summary); end
    @printf("Panel cap %g degrees: passage %.3f s, maximum sampled drag %.6f N, run %.2f s.\n",
        panel_cap_deg, sol.t[end], summary["maximum_recorded_drag_N"], wall_s)
    flush(stdout)
    return (; summary, rows, common, endpoint)
end

function compare_panel_caps(; output_dir=joinpath(pwd(), "odyssey_surrogate_results"), kwargs...)
    ispath(output_dir) && throw(ArgumentError("Choose a new output directory: $output_dir"))
    first = run_case(; panel_cap_deg=90.0, output_dir=joinpath(output_dir,"cap_90_deg"), kwargs...)
    second = run_case(; panel_cap_deg=30.0, output_dir=joinpath(output_dir,"cap_30_deg"), kwargs...)
    a = first.summary["common_state_m_m_s"]; b = second.summary["common_state_m_m_s"]
    position_difference = norm(a[1:3]-b[1:3])
    velocity_difference = norm(a[4:6]-b[4:6])
    heat_difference = norm(first.summary["common_panel_heat_load_J_cm2"] - second.summary["common_panel_heat_load_J_cm2"])
    position_difference > 1e-3 && velocity_difference > 1e-6 && heat_difference > 1e-6 ||
        error("The active panel-cap change did not produce a measurable trajectory and heating effect.")
    comparison = Dict("comparison_time_s"=>COMMON_TIME_S, "position_difference_m"=>position_difference,
        "velocity_difference_m_s"=>velocity_difference, "panel_heat_load_difference_norm_J_cm2"=>heat_difference,
        "effect_assertions_passed"=>true, "native_gram_loaded"=>native_gram_loaded(),
        "interpretation"=>"Different active panel commands under the same atmosphere; this is not an accuracy tolerance.")
    open(joinpath(output_dir,"comparison.toml"),"w") do io; TOML.print(io, comparison); end
    @printf("At %.0f s: position difference %.3f m, velocity difference %.6f m/s.\n", COMMON_TIME_S, position_difference, velocity_difference)
    return (; first, second, comparison)
end

function main(argv=ARGS)
    argv == ["--help"] && return println("julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl [--offline] [--output=NEW_DIRECTORY]\nRuns the active 90-degree/30-degree panel-cap comparison with the named frozen Odyssey preset.")
    offline = false; output = joinpath(pwd(), "odyssey_surrogate_results")
    seen = Set{String}()
    for arg in argv
        key = first(split(arg,'=';limit=2))
        key in seen && throw(ArgumentError("Repeated option $key")); push!(seen,key)
        if arg == "--offline"
            offline = true
        elseif startswith(arg,"--output=")
            output = last(split(arg,'=';limit=2)); isempty(output) && throw(ArgumentError("--output needs a directory."))
        else
            throw(ArgumentError("Unknown option $arg; use --help."))
        end
    end
    compare_panel_caps(; output_dir=output, offline=offline)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    OdysseySurrogateExample.main()
end
