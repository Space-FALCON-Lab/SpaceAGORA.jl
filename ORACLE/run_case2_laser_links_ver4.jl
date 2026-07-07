#!/usr/bin/env julia

# 1. load common.jl
const REPO_ROOT = normpath(joinpath(@__DIR__, "..")) # find path to the repository root
include(joinpath(REPO_ROOT, "examples", "common.jl")) # load common.jl for utility functions and types

# 2. load dependencies 
using CSV
using DataFrames
using LinearAlgebra
using Printf
using DiffEqBase
using StaticArrays
using .SimulationModel

# 3. -animate flag & import "10_Animation_ver2.jl"
const _HAS_GLMAKIE = "--animate" in ARGS && (try; @eval using GLMakie; true; catch; false; end)
_HAS_GLMAKIE && include(joinpath(@__DIR__, "10_Animation_ver2.jl"))

# 4. define output path
const DEFAULT_SUMMARY_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_summary.csv")
const DEFAULT_TIMESERIES_CSV = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_timeseries.csv")

# 5. define paper grid parameters (overridden by --paper-grid flag)
const PAPER_HELPER_ALTITUDES_KM = (1150.0, 1050.0, 1000.0, 950.0, 850.0)
const PAPER_HELPER_INCLINATION_DELTAS_DEG = (0.0, 0.5, 1.0)

# 6. settings container & default values
Base.@kwdef struct OracleCase2Options
    helpers::Int = 200
    helper_altitude_km::Float64 = 1050.0
    target_altitude_km::Float64 = 1000.0
    inclination_deg::Float64 = 0.0
    helper_inclination_delta_deg::Float64 = 0.0
    orbits::Float64 = 80.0
    schedule::Symbol = :naive_next_entering
    laser_range_km::Float64 = 200.0
    laser_power_w::Float64 = 10_000.0
    magnification::Float64 = 100.0
    beta::Float64 = 1.0
    eta::Float64 = 1.0
    mass_kg::Float64 = 227.0
    dt_max_s::Float64 = 10.0
    paper_grid::Bool = false
    output_csv::String = DEFAULT_SUMMARY_CSV
    timeseries_csv::String = DEFAULT_TIMESERIES_CSV
    append_output::Bool = false
    timeseries_points::Int = 1001
    animate::Bool = false
end


# 7. help text printed when you run the script with --help.
function _usage()
    return """
    Usage:
      julia --project=. ORACLE/run_case2_laser_links.jl [options]

    Options:
      --helpers N
      --helper-altitude-km KM
      --target-altitude-km KM
      --inclination-deg DEG
      --helper-inclination-delta-deg DEG
      --orbits N
      --schedule naive_next_entering|positive_along_track
      --laser-range-km KM
      --laser-power-w W
      --magnification B
      --beta VALUE
      --eta VALUE
      --mass-kg KG
      --dt-max-s SEC
      --paper-grid
      --output-csv PATH
      --timeseries-csv PATH
      --append-output
      --timeseries-points N
      --animate             Show 3D animation after the simulation (requires GLMakie)
    """
end

# 8. group the command-line option names by types
const _FLAG_OPTS   = (:paper_grid, :append_output, :animate)
const _INT_OPTS    = (:helpers, :timeseries_points)
const _SYMBOL_OPTS = (:schedule,)
const _PATH_OPTS   = (:output_csv, :timeseries_csv)
const _FLOAT_OPTS  = (
    :helper_altitude_km, :target_altitude_km, :inclination_deg,
    :helper_inclination_delta_deg, :orbits,
    :laser_range_km, :laser_power_w, :magnification, :beta, :eta,
    :mass_kg, :dt_max_s,
)

# 9. Reads what the user typed on the command line and turns it into an OracleCase2Options struct.
function _parse_options(argv)::OracleCase2Options
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(argv) # example of argv: ["--helpers", "200", "--laser-range-km", "150.0", "--schedule", "positive_along_track"]
        arg = argv[i] # Get the current argument string, e.g. "--helpers" or "--laser-range-km"
        arg in ("--help", "-h") && (println(_usage()); exit(0)) # If the user typed --help or -h, print the usage text and quit
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'.\n$(_usage())")) # Every argument must start with "--", otherwise crash with an error
        key = Symbol(replace(arg[3:end], '-' => '_')) # Strip the leading "--" and convert dashes to underscores
        if key in _FLAG_OPTS
            opts[key] = true; i += 1; continue # For flags like --animate where there is no value after them, just mark them as true, step forward by 1, and go to the next argument
        end
        i < length(argv) || throw(ArgumentError("Missing value for $arg.")) # Make sure there IS a next input to read as the value for this arg
        val = argv[i + 1] # Grab the next word as the value, e.g. "200" for --helpers 200
        if     key in _INT_OPTS    opts[key] = parse(Int, val) # if the key is in the integer options, parse the value as an integer
        elseif key in _SYMBOL_OPTS opts[key] = Symbol(val) # if the key is in the symbol options, convert the value to a symbol
        elseif key in _PATH_OPTS   opts[key] = abspath(val) # if the key is in the path options, convert the value to an absolute path
        elseif key in _FLOAT_OPTS  opts[key] = parse(Float64, val) # if the key is in the float options, parse the value as a float
        else   throw(ArgumentError("Unknown option $arg.\n$(_usage())"))
        end
        i += 2 # Step forward by 2 to move past the key and its value
    end
    return OracleCase2Options(; opts...)
end

# 10. Validate that the inputs are within reasonable bounds, and throw an error if not.
function _validate_options(opts::OracleCase2Options)
    opts.helpers >= 1 || throw(ArgumentError("--helpers must be >= 1."))
    opts.helper_altitude_km > 0.0 || throw(ArgumentError("--helper-altitude-km must be positive."))
    opts.target_altitude_km > 0.0 || throw(ArgumentError("--target-altitude-km must be positive."))
    opts.orbits > 0.0 || throw(ArgumentError("--orbits must be positive."))
    opts.schedule in (:naive_next_entering, :positive_along_track) ||
        throw(ArgumentError("--schedule must be naive_next_entering or positive_along_track."))
    opts.laser_range_km >= 0.0 || throw(ArgumentError("--laser-range-km must be nonnegative."))
    opts.laser_power_w >= 0.0 || throw(ArgumentError("--laser-power-w must be nonnegative."))
    opts.magnification >= 0.0 || throw(ArgumentError("--magnification must be nonnegative."))
    opts.beta >= 0.0 || throw(ArgumentError("--beta must be nonnegative."))
    opts.eta >= 0.0 || throw(ArgumentError("--eta must be nonnegative."))
    opts.mass_kg > 0.0 || throw(ArgumentError("--mass-kg must be positive."))
    opts.dt_max_s > 0.0 || throw(ArgumentError("--dt-max-s must be positive."))
    opts.timeseries_points >= 2 || throw(ArgumentError("--timeseries-points must be >= 2."))
    return nothing
end

################################################################ DEFINE SPACECRAFT HERE ################################################################
# 11. This function creates a SpacecraftModel with the given id, mass, and initial condition.
function _spacecraft(id::Int, mass_kg::Float64, ic::InitialCondition)
# Takes 3 inputs: a numeric ID, the mass, and an initial orbital condition (position/velocity)

    bus = Link(root=true, m=mass_kg, ref_area=1.0)
    # Creates the spacecraft's single structural body — the "bus" (main body).
    # root=true means it's the top of the body tree (no parent).
    # ref_area=1.0 is the reference area (m²) used for drag/pressure calculations.

    return SpacecraftModel(
        Joint[],       # No joints — this is a rigid single-body spacecraft, no articulated parts
        [bus],         # List of all Links (bodies); just the one bus
        bus,           # The root link (the bus itself)
        true,          # keplerian=true — use Keplerian orbital dynamics
        mass_kg,       # Total mass
        0.0,           # Initial drag coefficient or similar scalar — set to zero here
        bus.inertia,   # Moment of inertia tensor, taken from the bus body
        0,             # Number of reaction wheels = 0
        0,             # Number of thrusters = 0
        ic,            # Initial orbital condition (altitude, inclination, true anomaly, etc.)
        id,            # Unique integer ID for this spacecraft (1 = target, 2..N = helpers)
    )
end
########################################################################### END SPACECRAFT DEFINITION ###############################################################

# RTN frame unit vectors
@inline function _rtn_basis(r::SVector{3, Float64}, v::SVector{3, Float64})
    rhat = r / norm(r)
    nhat = cross(r, v)
    nhat = nhat / norm(nhat)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

# rv to orbit elements
function _rv_to_elements(r::SVector{3, Float64}, v::SVector{3, Float64}, mu::Float64)
    rmag = norm(r)
    v2 = dot(v, v)
    h = cross(r, v)
    hmag = norm(h)
    n = cross(SVector{3, Float64}(0.0, 0.0, 1.0), h)
    nmag = norm(n)
    evec = cross(v, h) / mu - r / rmag
    e = norm(evec)
    a = -mu / (2.0 * (0.5 * v2 - mu / rmag))
    inc = acos(clamp(h[3] / hmag, -1.0, 1.0))
    raan = nmag <= 1e-12 ? 0.0 : mod(atan(n[2], n[1]), 2pi)
    return (a=a, e=e, i=inc, raan=raan)
end

# Laser impulse tracker: integrates F/m·dt at every accepted ODE step in RTN frame
Base.@kwdef mutable struct _LaserImpulseTracker
    # ---------------------------------------------------------------------------
    # Laser impulse tracker: integrates F/m·dt at every accepted ODE step so that
    # the reported dV_RTN is the true accumulated laser impulse, not an ECI
    # velocity difference that is contaminated by orbital-phase drift.
    # ---------------------------------------------------------------------------
    t_prev::Float64            = 0.0
    dv_R::Float64              = 0.0
    dv_T::Float64              = 0.0
    dv_N::Float64              = 0.0
    t_hist::Vector{Float64}    = Float64[]
    dv_R_hist::Vector{Float64} = Float64[]
    dv_T_hist::Vector{Float64} = Float64[]
    dv_N_hist::Vector{Float64} = Float64[]
end

# Creates a DiscreteCallback that accumulates the laser impulse in RTN frame at every accepted ODE step.
function _make_laser_impulse_callback(
    model::OpenCavityLaserLinkModel, # the live laser model (isolate_state=false, so no deep copy)
    tracker::_LaserImpulseTracker,   # mutable struct that accumulates dV; lives outside the integrator
    mass_kg::Float64,                # target satellite mass, needed to convert force → acceleration
)
    function affect!(integrator)            # called by DiffEq at every accepted ODE step
        dt = integrator.t - tracker.t_prev  # time elapsed since the last accepted step
        if dt > 0.0                         # skip the very first call where t_prev == t (dt = 0)
            helper_idx = model.active_helper_idx  # which helper is currently firing (0 = none)
            if helper_idx > 0                     # only accumulate when a link is active
                sc = integrator.u.sc              # array of all spacecraft states at this step
                tgt_pos = SVector{3, Float64}(sc[model.target_idx].pos)  # target ECI position [m]
                tgt_vel = SVector{3, Float64}(sc[model.target_idx].vel)  # target ECI velocity [m/s]
                hlp_pos = SVector{3, Float64}(sc[helper_idx].pos)        # active helper ECI position [m]

                # Recompute the laser force inline (laser_link_pair_force is not exported)
                rel = tgt_pos - hlp_pos   # vector from helper → target [m]
                rho = norm(rel)            # separation distance [m]
                if rho > 0.0 && rho <= model.range_m  # link is physically valid
                    # Scalar force magnitude: F = η·β·M·P / c  [N]
                    F_mag = model.eta * model.beta *
                            model.magnification * model.power_w /
                            299_792_458.0
                    force = F_mag * rel / rho              # force vector along line-of-sight [N]
                    rhat, that, nhat = _rtn_basis(tgt_pos, tgt_vel)  # RTN unit vectors at target
                    accel = force / mass_kg                # acceleration vector [m/s²]
                    # Trapezoidal-style: a·dt ≈ ΔV for this step (left-endpoint rectangle rule)
                    tracker.dv_R += dot(accel, rhat) * dt  # radial component [m/s]
                    tracker.dv_T += dot(accel, that) * dt  # along-track component [m/s]
                    tracker.dv_N += dot(accel, nhat) * dt  # cross-track component [m/s]
                end
            end
        end
        tracker.t_prev = integrator.t  # advance the "last step" timestamp for the next call

        # Append current time and cumulative dV to the history vectors for post-run lookup
        push!(tracker.t_hist,    integrator.t)
        push!(tracker.dv_R_hist, tracker.dv_R)
        push!(tracker.dv_T_hist, tracker.dv_T)
        push!(tracker.dv_N_hist, tracker.dv_N)
    end
    return DiscreteCallback(
        (u, t, integrator) -> true,  # condition: fire at every accepted step
        affect!;
        save_positions=(false, false), # don't save an extra solution snapshot at trigger/affect time
    )
end

# Look up the accumulated (dv_R, dv_T, dv_N) at a given time from the history.
function _tracked_dv_at(tracker::_LaserImpulseTracker, t::Float64)
    isempty(tracker.t_hist) && return (0.0, 0.0, 0.0)
    idx = clamp(searchsortedlast(tracker.t_hist, t), 1, length(tracker.t_hist))
    return (tracker.dv_R_hist[idx], tracker.dv_T_hist[idx], tracker.dv_N_hist[idx])
end

# Compute cumulative orbit count for the TARGET satellite using the instantaneous
# semi-major axis a(t) at each saved time step. Correctly accounts for orbital period
# changes caused by the laser (a grows → T grows → fewer orbits per second).
# orbit_count(t) = ∫₀ᵗ dt'/T(t')  ≈  cumsum(Δt / T_mid(tₖ))
function _orbit_count_from_sol(sol, mu::Float64)
    T_series = [begin
        sc = u.sc[1]
        r  = SVector{3, Float64}(sc.pos)
        v  = SVector{3, Float64}(sc.vel)
        a  = _rv_to_elements(r, v, mu).a
        2pi * sqrt(a^3 / mu)
    end for u in sol.u]
    dt    = diff(Float64.(sol.t))
    T_mid = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
    return cumsum([0.0; dt ./ T_mid])  # length == length(sol.t)
end

# Builds the configuration for the case 2 of the Oracle simulation.
function build_case_config(opts::OracleCase2Options)
    # 1. Sanity-checks every number the user passed in (altitudes > 0, mass > 0, etc.) 
    _validate_options(opts)
    
    # 2. Planet + radii
    planet = make_no_gram_planet(:earth)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3

    # 3. Build the spacecraft array
    spacecraft = SpacecraftModel[]
    push!(spacecraft, _spacecraft(1, opts.mass_kg, InitialCondition(
        target_radius_m, 0.0, opts.inclination_deg, 0.0, 0.0, 0.0
    ))) # 3.1. Creates one target satellite (ID 1)

    helper_inclination_deg = opts.inclination_deg + opts.helper_inclination_delta_deg
    for helper in 1:opts.helpers
        nu_deg = 360.0 * (helper - 1) / opts.helpers # evenly space them around the orbit
        push!(spacecraft, _spacecraft(helper + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, helper_inclination_deg, 0.0, 0.0, nu_deg
        )))
    end # 3.2. Create N helper satellites (IDs 2…N+1)

    # 4. Build the laser model
    laser_model = OpenCavityLaserLinkModel(
        1,  # target is spacecraft #1
        collect(2:(opts.helpers + 1)); # helpers are #2..N+1
        range_m=opts.laser_range_km * 1e3,
        power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        schedule=opts.schedule,
    )

    # 5. Mission time duration in seconds
    target_period_s = 2pi * sqrt(target_radius_m^3 / planet.μ)  # initial period [s] (for reference)
    mission_time_s = opts.orbits * target_period_s

    # 6. The big constructor — wires everything together into one config object
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
        ), # 6.1. Turns off results files, verbose logging, and plots — this script handles output itself

        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000,
            data_rate=max(10.0, mission_time_s / 1000.0),
        ), # 6.2. Sets the stop condition (MissionTime), how many solution snapshots to save (~1000), and the data output rate

        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
        ), # 6.3. Earth gravity + no atmosphere + simple Sun/Moon ephemerides + thermal model

        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(), laser_model)), # 6.4. forces acting on all spacecraft: J2 gravity + the laser link model

        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]), 
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]), # 6.5, no attitude control, no navigation

        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0), # 6.6. Calendar start date (Jan 1 2026) used for ephemerides
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-12,
            abstol_orbit=1e-12,
            dt_max_orbit=opts.dt_max_s,
        ), # 6.7. Very tight ODE tolerances (1e-12) and a max step size so the integrator doesn't skip over fast events

        solver_config=SolverConfig(solver_mode=:tsit5), # 6.8. 	Use the Tsitouras 5th-order Runge-Kutta solver (:tsit5)
    )
    return args, laser_model, target_period_s
end

# 7. get rv at any time t from the solution object (sol) for the target satellite (sc[1])
@inline function _target_rv_at(sol, t::Float64)
    u = if isapprox(t, Float64(sol.t[1]); atol=0.0, rtol=0.0)
        sol.u[1]
    elseif isapprox(t, Float64(sol.t[end]); atol=0.0, rtol=0.0)
        sol.u[end]
    else
        sol(t)
    end
    sc = u.sc[1]
    return SVector{3, Float64}(sc.pos), SVector{3, Float64}(sc.vel)
end 

# 8. Generate a unique case ID string based on the options (used for labeling output)
function _case_id(opts::OracleCase2Options)::String
    return @sprintf(
        "Nh%d_h%.0fkm_dinc%.1fdeg",
        opts.helpers,
        opts.helper_altitude_km,
        opts.helper_inclination_delta_deg,
    )
end

# 9. Build time-series DataFrame of cumulative ΔP and ΔE for each (helper, target) pair over the simulation.
function _build_timeseries_dataframe(
    opts::OracleCase2Options,
    sol,
    oe0,
    orbit_counts::Vector{Float64},   # cumulative orbit count at each sol.t point (from _orbit_count_from_sol)
    mu::Float64,
    tracker::_LaserImpulseTracker,
)::DataFrame
    tf = Float64(sol.t[end])
    times = collect(range(0.0, tf; length=opts.timeseries_points))
    sol_times = Float64.(sol.t)      # saved time points corresponding to orbit_counts
    rows = NamedTuple[]
    case_id = _case_id(opts)
    helper_inclination_deg = opts.inclination_deg + opts.helper_inclination_delta_deg
    for t in times
        r, v = _target_rv_at(sol, Float64(t))
        dv_r_ts, dv_t_ts, dv_n_ts = _tracked_dv_at(tracker, Float64(t))
        oe = _rv_to_elements(r, v, mu)
        oc_idx = clamp(searchsortedlast(sol_times, Float64(t)), 1, length(sol_times))
        push!(rows, (
            case_id=case_id,
            time_s=Float64(t),
            orbit=orbit_counts[oc_idx],
            helpers=opts.helpers,
            helper_altitude_km=opts.helper_altitude_km,
            target_altitude_km=opts.target_altitude_km,
            target_inclination_deg=opts.inclination_deg,
            helper_inclination_delta_deg=opts.helper_inclination_delta_deg,
            helper_inclination_deg=helper_inclination_deg,
            schedule=opts.schedule,
            laser_range_km=opts.laser_range_km,
            laser_power_w=opts.laser_power_w,
            magnification=opts.magnification,
            beta=opts.beta,
            eta=opts.eta,
            mass_kg=opts.mass_kg,
            dv_r_mps=dv_r_ts,
            dv_t_mps=dv_t_ts,
            dv_n_mps=dv_n_ts,
            da_m=oe.a - oe0.a,
            de=oe.e - oe0.e,
            di_deg=rad2deg(oe.i - oe0.i),
            draan_deg=rad2deg(oe.raan - oe0.raan),
        ))
    end
    return DataFrame(rows)
end

# push!() is a Julia standard library function that appends one or more items to the end of a collection (array, vector, etc.) in-place (the ! convention means it mutates its first argument).
# v = [1, 2, 3]
# push!(v, 4)   
# v is now [1, 2, 3, 4]

# 10. Run the simulation for the given options, returning a summary, time-series DataFrame, and the full solution object.
function run_open_cavity_case(opts::OracleCase2Options)
    # --- Block 1: build config ---
    # Build the SimulationConfiguration, laser model, and initial orbital period from opts
    args, laser_model, target_period_s = build_case_config(opts)

    # --- Block 2: initial orbital elements ---
    u0 = SimulationEngine.build_initial_conditions(args)  # convert args → initial state vector
    r0 = SVector{3, Float64}(u0.sc[1].pos)                # target initial position [m]
    v0 = SVector{3, Float64}(u0.sc[1].vel)                # target initial velocity [m/s]
    oe0 = _rv_to_elements(r0, v0, args.environment_model.planet.μ)  # initial orbital elements (a,e,i,raan)

    # --- Block 3: callbacks ---
    # impulse_cb must come FIRST so it reads the active_helper_idx that was set
    # for the step just completed, before the scheduler updates it for the next.
    impulse_tracker = _LaserImpulseTracker()               # accumulates dV in RTN at every ODE step
    impulse_cb  = _make_laser_impulse_callback(laser_model, impulse_tracker, opts.mass_kg)
    scheduler_cb = laser_link_scheduler_callback(laser_model)  # updates which helper fires next

    # --- Block 4: run the ODE solver ---
    result = run_simulation(
        args;
        isolate_state=false,       # no deep copy — laser_model IS the live model, so the callback can use it directly
        return_solution=true,      # keep the full solution object (dense output + saved snapshots)
        return_solver_metadata=true,  # keep solver name and step stats for the summary
        extra_callbacks=(impulse_cb, scheduler_cb),
    )

    # --- Block 5: post-simulation quantities ---
    sol = result.solution
    orbit_counts   = _orbit_count_from_sol(sol, args.environment_model.planet.μ)  # varying-period orbit counter
    orbits_elapsed = orbit_counts[end]   # total orbits completed (accounts for changing period)
    final_state = sol.u[end].sc[1]       # target state at the last saved time step
    rf = SVector{3, Float64}(final_state.pos)  # final position [m]
    vf = SVector{3, Float64}(final_state.vel)  # final velocity [m/s]
    oef = _rv_to_elements(rf, vf, args.environment_model.planet.μ)  # final orbital elements
    # laser_model is the live model (isolate_state=false), so its final scheduler state is already up to date

    # --- Block 6: pack results into a summary named tuple ---
    # Each field becomes one column when written to the summary CSV
    summary = (
        case_id=_case_id(opts),                          # unique string label for this run
        helpers=opts.helpers,                            # number of helper satellites
        helper_altitude_km=opts.helper_altitude_km,      # helper orbital altitude
        schedule=opts.schedule,                          # which helper-selection algorithm was used
        target_period_s=target_period_s,                 # initial orbital period of the target [s]
        active_helper=laser_model.active_helper_idx,     # which helper was firing at end of sim (0 = none)
        activations=laser_model.link_activation_count,   # total number of helper switches
        active_steps=laser_model.active_link_step_count, # total ODE steps where the laser was on
        target_altitude_km=opts.target_altitude_km,
        inclination_deg=opts.inclination_deg,
        helper_inclination_delta_deg=opts.helper_inclination_delta_deg,
        helper_inclination_deg=opts.inclination_deg + opts.helper_inclination_delta_deg,
        orbits=opts.orbits,                              # requested number of orbits (input)
        orbits_elapsed=orbits_elapsed,                   # actual orbits completed (varying period)
        laser_range_km=opts.laser_range_km,
        laser_power_w=opts.laser_power_w,
        magnification=opts.magnification,
        beta=opts.beta,
        eta=opts.eta,
        mass_kg=opts.mass_kg,
        dt_max_s=opts.dt_max_s,
        dv_r_mps=impulse_tracker.dv_R,   # total radial ΔV from laser [m/s]
        dv_t_mps=impulse_tracker.dv_T,   # total along-track ΔV from laser [m/s]
        dv_n_mps=impulse_tracker.dv_N,   # total cross-track ΔV from laser [m/s]
        da_m=oef.a - oe0.a,              # change in semi-major axis [m]
        de=oef.e - oe0.e,                # change in eccentricity
        di_deg=rad2deg(oef.i - oe0.i),   # change in inclination [deg]
        draan_deg=rad2deg(oef.raan - oe0.raan),  # change in RAAN [deg]
        retcode=sol.retcode,             # ODE solver exit status (Success, MaxIters, etc.)
        solver=result.solver_trace[end].solver,  # name of the solver used
    )

    # --- Block 7: build the time-series DataFrame and return ---
    timeseries = _build_timeseries_dataframe(
        opts,
        sol,
        oe0,
        orbit_counts,
        args.environment_model.planet.μ,
        impulse_tracker,
    )
    return (summary=summary, timeseries=timeseries, sol=sol, helper_num=opts.helpers)
end

# 11. Write a DataFrame to CSV, creating the directory if needed, and optionally appending to an existing file.
function _write_csv!(path::String, data; append::Bool=false)
    isempty(path) && return ""
    mkpath(dirname(path))
    do_append = append && isfile(path)
    df = data isa DataFrame ? data : DataFrame(data)
    CSV.write(path, df; append=do_append, header=!do_append)
    return path
end

# 12. Print a summary of the simulation results to the console.
function _print_summary(s)
    println("ORACLE Open Cavity Case laser-link run")
    println("  helpers              : $(s.helpers)")
    println("  helper altitude       : $(s.helper_altitude_km) km")
    println("  schedule              : $(s.schedule)")
    println("  solver                : $(s.solver)")
    println("  retcode               : $(s.retcode)")
    @printf("  mission time          : %.2f s (%.4f orbits elapsed)\n", s.orbits * s.target_period_s, s.orbits_elapsed)
    println("  link activations      : $(s.activations)")
    println("  accepted active steps : $(s.active_steps)")
    @printf("  dV_RTN [m/s]          : R=% .6e  T=% .6e  N=% .6e\n", s.dv_r_mps, s.dv_t_mps, s.dv_n_mps)
    @printf("  da [m]                : % .6e\n", s.da_m)
    @printf("  de                    : % .6e\n", s.de)
    @printf("  di [deg]              : % .6e\n", s.di_deg)
    @printf("  dRAAN [deg]           : % .6e\n", s.draan_deg)
end

# 13. Copy opts, assist publication grid search
_with(opts::OracleCase2Options; kwargs...) =
    OracleCase2Options(; (f => getfield(opts, f) for f in fieldnames(OracleCase2Options))..., kwargs...)

function main(argv=ARGS)
    opts = _parse_options(argv)  # parse command-line arguments into an OracleCase2Options struct

    if opts.paper_grid
        # --- Paper-grid mode: sweep over all (helper_altitude, inclination_delta) combinations ---

        # Delete old output files unless the user asked to append
        opts.append_output || (isfile(opts.output_csv)      && rm(opts.output_csv))
        opts.append_output || (isfile(opts.timeseries_csv)  && rm(opts.timeseries_csv))

        for helper_alt in PAPER_HELPER_ALTITUDES_KM               # outer loop: 5 altitudes
            for helper_inclination_delta in PAPER_HELPER_INCLINATION_DELTAS_DEG  # inner loop: 3 inclination deltas
                # Run one simulation with this (altitude, inclination_delta) combination,
                # keeping all other opts unchanged
                result = run_open_cavity_case(_with(opts;
                    helper_altitude_km=helper_alt,
                    helper_inclination_delta_deg=helper_inclination_delta,
                    append_output=true,   # always append inside the grid loop
                    animate=false,        # no animation during batch runs
                ))
                s = result.summary

                # Append this run's one-row summary and full timeseries to the CSVs
                _write_csv!(opts.output_csv,     [s];              append=true)
                _write_csv!(opts.timeseries_csv, result.timeseries; append=true)

                # Print a one-line progress update to the terminal
                @printf(
                    "helpers=%d helper_alt_km=%.1f helper_dinc_deg=%.1f dv_R=%.6e dv_T=%.6e dv_N=%.6e activations=%d\n",
                    s.helpers, s.helper_altitude_km, s.helper_inclination_delta_deg,
                    s.dv_r_mps, s.dv_t_mps, s.dv_n_mps, s.activations
                )
            end
        end
        println("Summary CSV: $(opts.output_csv)")
        println("Time-series CSV: $(opts.timeseries_csv)")

    else
        # --- Single-run mode: run once with exactly what the user typed ---
        result = run_open_cavity_case(opts)
        s = result.summary

        _print_summary(s)  # print human-readable results to the terminal

        # Write summary and timeseries CSVs (returns "" if path is empty, so we skip the print)
        out    = _write_csv!(opts.output_csv,     [s];              append=opts.append_output)
        ts_out = _write_csv!(opts.timeseries_csv, result.timeseries; append=opts.append_output)
        isempty(out)    || println("Summary CSV: $out")
        isempty(ts_out) || println("Time-series CSV: $ts_out")

        # --- Optional 3D animation (only if --animate was passed and GLMakie loaded) ---
        if opts.animate
            if _HAS_GLMAKIE
                p_anim = Dict{Symbol, Any}(:N => opts.helpers + 1)  # total spacecraft count
                _, anim_controls = animate_all_satellites_3d_smooth_helper_target(
                    result.sol, p_anim, opts.helpers)
                println("Animation window opened. Close the window to exit.")
                wait(anim_controls.screen)  # block until the user closes the window
            else
                @warn "GLMakie is not installed — animation skipped. " *
                      "Install it with: using Pkg; Pkg.add(\"GLMakie\")"
            end
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
