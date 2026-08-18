#!/usr/bin/env julia
# =============================================================================
# Section IV.F — Conservation-Law Numerical Validation
# =============================================================================
# Runs two 2-spacecraft OCL simulations (J2 inactive and J2 active) and writes
# a Feather file for each run.  All OCL parameters match Table 2 of the paper:
#
#   Laser range:       200 km      (L_max)
#   Laser power:    10 000 W       (P_ij)
#   OCL magnification: 100          (B)
#   Geometric loss:      1.0        (β_ij)
#
# The target is initialized at 1100 km altitude; the helper at 1000 km.  Both
# satellites are in equatorial circular orbits (i = 0°).  The 100 km altitude
# difference keeps the helper in and out of the 200 km laser range, producing
# the intermittent-link scenario described in Sec. IV.F.
#
# Usage:
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/run_verification.jl
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/run_verification.jl --orbits 2.0
#   julia --project=. ORACLE/Paper_Numerical_Verification_Code/run_verification.jl --no-j2
# =============================================================================

# ── 1. Bootstrap SpaceAGORA ──────────────────────────────────────────────────
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))

# ── 2. Dependencies ──────────────────────────────────────────────────────────
using Arrow
using DataFrames
using LinearAlgebra
using Printf
using StaticArrays
using DiffEqBase
using .SimulationModel

# ── 3. Shared orbital-mechanics helpers (only OE converters + RTN basis needed) ──
# Note: 3_Dynamics.jl is NOT included here because it defines build_case_config
# with an OracleCase2Options type annotation that is only in scope when loaded via
# run_case2_laser_links.jl.  The _spacecraft function from 0_Spacecraft.jl and
# OE converters are the only shared dependencies.
include(joinpath(@__DIR__, "..", "functions", "0_Spacecraft.jl"))
include(joinpath(@__DIR__, "..", "functions", "5_OE_Converters.jl"))  # _rtn_basis, rv2coe, ...

# ── 4. Default output directory ───────────────────────────────────────────────
const DEFAULT_VER_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "Paper_Numerical_Verification_Code")

# ─────────────────────────────────────────────────────────────────────────────
# Options container
# ─────────────────────────────────────────────────────────────────────────────
Base.@kwdef struct VerificationOptions
    helpers::Int                = 1             # number of helper spacecraft
    helper_altitude_km::Float64 = 1000.0
    target_altitude_km::Float64 = 1100.0
    inclination_deg::Float64    = 0.0
    orbits::Float64             = 50
    use_j2::Bool                = true        # InverseSquaredJ2GravityModel vs InverseSquaredGravityModel
    laser_range_km::Float64     = 200.0
    laser_power_w::Float64      = 10_000.0
    magnification::Float64      = 100.0
    beta::Float64               = 1.0
    eta::Float64                = 1.0
    mass_kg::Float64            = 227.0
    # dt_max_s controls conservation quality of H_z in the saved output.
    # Tsit5 (non-symplectic, 5th-order) accumulates H_z integration error at O(h^4).
    # dt_max=10s → |ΔH_z|/|H_z| ≈ 10⁻¹¹; dt_max=1s → ≈ 10⁻¹⁵ (near machine epsilon).
    dt_max_s::Float64           = 1.0
    output_dir::String          = DEFAULT_VER_OUTPUT_DIR
    timeseries_points::Int      = 1001
end

# ─────────────────────────────────────────────────────────────────────────────
# Verification tracker — extends the RTN impulse tracker with:
#   • W_ocl   : accumulated OCL work  ∫ (F_t·v_t + F_h·v_h) dt  [J]
#   • eps_H   : angular-momentum exchange residual
#               ∫ (r_t × F_{t←h} + r_h × F_{h←t}) dt            [kg·m²/s]
# Both are accumulated at every accepted ODE step via DiscreteCallback.
# ─────────────────────────────────────────────────────────────────────────────
Base.@kwdef mutable struct _VerificationTracker
    # ── RTN impulse (mirrors _LaserImpulseTracker) ──
    t_prev::Float64             = 0.0
    dv_R::Float64               = 0.0
    dv_T::Float64               = 0.0
    dv_N::Float64               = 0.0
    t_hist::Vector{Float64}     = Float64[]
    dv_R_hist::Vector{Float64}  = Float64[]
    dv_T_hist::Vector{Float64}  = Float64[]
    dv_N_hist::Vector{Float64}  = Float64[]

    # ── OCL energy and angular-momentum accumulators ──
    W_ocl::Float64              = 0.0          # accumulated OCL work [J]
    eps_H_x::Float64            = 0.0          # angular-momentum exchange residual x [kg·m²/s]
    eps_H_y::Float64            = 0.0
    eps_H_z::Float64            = 0.0
end

# Build the DiscreteCallback that accumulates all verification quantities.
function _make_verification_callback(
    model::OpenCavityLaserLinkModel,
    tracker::_VerificationTracker,
    mass_kg::Float64,
)
    c_light = 299_792_458.0   # speed of light [m/s]

    function affect!(integrator)
        dt = integrator.t - tracker.t_prev
        if dt > 0.0
            helper_idx = model.active_helper_idx
            if helper_idx > 0
                sc = integrator.u.sc

                # Target state
                r_t = SVector{3,Float64}(sc[model.target_idx].pos)
                v_t = SVector{3,Float64}(sc[model.target_idx].vel)

                # Helper state
                r_h = SVector{3,Float64}(sc[helper_idx].pos)
                v_h = SVector{3,Float64}(sc[helper_idx].vel)

                # OCL force on target (from helper → target direction)
                rel = r_t - r_h         # displacement from helper to target
                rho = norm(rel)

                if rho > 0.0 && rho <= model.range_m
                    F_mag = model.eta * model.beta * model.magnification *
                            model.power_w / c_light      # [N]
                    F_t = F_mag * rel / rho              # force on target
                    F_h = -F_t                           # Newton's 3rd law

                    # ── RTN dV accumulation ──────────────────────────────────
                    rhat, that, nhat = _rtn_basis(r_t, v_t)
                    accel_t = F_t / mass_kg
                    tracker.dv_R += dot(accel_t, rhat) * dt
                    tracker.dv_T += dot(accel_t, that) * dt
                    tracker.dv_N += dot(accel_t, nhat) * dt

                    # ── OCL mechanical power W_OCL ───────────────────────────
                    # W_dot = F_t·v_t + F_h·v_h
                    tracker.W_ocl += (dot(F_t, v_t) + dot(F_h, v_h)) * dt

                    # ── Angular-momentum exchange residual ε_H ───────────────
                    # ε_H = ∫ (r_t × F_{t←h} + r_h × F_{h←t}) dt
                    # Since F_{h←t} = -F_{t←h}:
                    #   r_t × F_t + r_h × (-F_t) = (r_t - r_h) × F_t
                    # This should be ≈ 0 (forces are collinear with displacement).
                    tau = cross(r_t - r_h, F_t)
                    tracker.eps_H_x += tau[1] * dt
                    tracker.eps_H_y += tau[2] * dt
                    tracker.eps_H_z += tau[3] * dt
                end
            end
        end
        tracker.t_prev = integrator.t

        push!(tracker.t_hist,    integrator.t)
        push!(tracker.dv_R_hist, tracker.dv_R)
        push!(tracker.dv_T_hist, tracker.dv_T)
        push!(tracker.dv_N_hist, tracker.dv_N)
    end

    return DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions=(false, false),
    )
end

# SaveField entries for extra verification columns.
# `laser_model` is captured in closures so its live state is read at each output point.
function _build_verification_save_fields(
    tracker::_VerificationTracker,
    laser_model::OpenCavityLaserLinkModel,
)
    return SaveField[
        # Accumulated OCL work [J] — from tracker, integrated over all ODE steps
        SaveField(
            :W_ocl_accumulated,
            (u, t, integrator) -> tracker.W_ocl;
            per_satellite=false, column_prefix="W_ocl_accumulated",
        ),
        # Angular-momentum exchange residual components [kg·m²/s] — from tracker
        SaveField(
            :eps_H_x,
            (u, t, integrator) -> tracker.eps_H_x;
            per_satellite=false, column_prefix="eps_H_x",
        ),
        SaveField(
            :eps_H_y,
            (u, t, integrator) -> tracker.eps_H_y;
            per_satellite=false, column_prefix="eps_H_y",
        ),
        SaveField(
            :eps_H_z,
            (u, t, integrator) -> tracker.eps_H_z;
            per_satellite=false, column_prefix="eps_H_z",
        ),
        # Active-link indicator: 0 = no link active, 2 = helper spacecraft 2 is linked
        SaveField(
            :laser_active_helper,
            (u, t, integrator) -> Float64(laser_model.active_helper_idx);
            per_satellite=false, column_prefix="laser_active_helper",
        ),
    ]
end

# ─────────────────────────────────────────────────────────────────────────────
# SimulationConfiguration builder
# ─────────────────────────────────────────────────────────────────────────────
function _build_verification_config(opts::VerificationOptions, results_dir::String)
    planet          = make_no_gram_planet(:earth)
    target_radius_m = planet.Rp_e + opts.target_altitude_km * 1e3
    helper_radius_m = planet.Rp_e + opts.helper_altitude_km * 1e3
    target_period_s = 2π * sqrt(target_radius_m^3 / planet.μ)
    mission_time_s  = opts.orbits * target_period_s

    # N+1 spacecraft: target (id=1) + opts.helpers helpers (ids 2..N+1)
    # Helpers are uniformly phased in true anomaly over [0°, 360°)
    spacecraft = SpacecraftModel[
        _spacecraft(1, opts.mass_kg, InitialCondition(
            target_radius_m, 0.0, opts.inclination_deg, 0.0, 0.0, 0.0,
        ))
    ]
    for k in 1:opts.helpers
        nu_deg = 360.0 * (k - 1) / opts.helpers   # uniform spacing
        push!(spacecraft, _spacecraft(k + 1, opts.mass_kg, InitialCondition(
            helper_radius_m, 0.0, opts.inclination_deg, 0.0, 0.0, nu_deg,
        )))
    end

    laser_model = OpenCavityLaserLinkModel(
        1, collect(2:(opts.helpers + 1));
        range_m       = opts.laser_range_km * 1e3,
        power_w       = opts.laser_power_w,
        magnification = opts.magnification,
        beta          = opts.beta,
        eta           = opts.eta,
        schedule      = :naive_next_entering,
    )

    # J2 toggle — use InverseSquaredGravityModel for two-body only
    gravity_model = opts.use_j2 ? InverseSquaredJ2GravityModel() : InverseSquaredGravityModel()

    args = SimulationConfiguration(
        simulation_settings = SimulationSettings(
            results           = true,
            verbose           = false,
            results_directory = results_dir,
            generate_plots    = false,
            save_csv          = false,
            normalize         = false,
        ),
        mission_configuration = MissionConfiguration(
            mission_type       = MissionTime,
            keplerian          = true,
            number_of_orbits   = 1,
            mission_time       = mission_time_s,
            orientation_sim    = false,
            num_steps_to_save  = opts.timeseries_points - 1,
            data_rate          = max(1.0, mission_time_s / (opts.timeseries_points - 1)),
        ),
        environment_model = EnvironmentModel(
            planet             = planet,
            EI                 = 120.0,
            density_model      = NoAtmosphereModel(),
            ephemerides_model  = SimpleEphemeridesModel(),
            thermal_model      = MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography         = false,
            wind               = false,
        ),
        dynamics_model = DynamicsModel(spacecraft, (gravity_model, laser_model)),
        guidance_model    = GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model  = NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model     = ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time      = InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances = IntegrationTolerances(
            reltol_orbit = 1e-12,
            abstol_orbit = 1e-12,
            dt_max_orbit = opts.dt_max_s,
        ),
        solver_config = SolverConfig(solver_mode=:tsit5),
    )

    return args, laser_model, target_period_s
end

# ─────────────────────────────────────────────────────────────────────────────
# Main runner
# ─────────────────────────────────────────────────────────────────────────────
"""
    run_verification_case(opts; label="") → results_dir::String

Run a single 2-spacecraft validation simulation and write a Feather file.
`label` is appended to the output directory stem to distinguish J2-on / J2-off runs.
"""
function run_verification_case(opts::VerificationOptions; label::String="")
    j2_tag      = opts.use_j2 ? "j2_on" : "j2_off"
    stem        = label == "" ? j2_tag : "$(label)_$(j2_tag)"
    results_dir = joinpath(opts.output_dir, stem)
    mkpath(results_dir)

    # Build config
    args, laser_model, target_period_s = _build_verification_config(opts, results_dir)

    # Build callbacks
    ver_tracker  = _VerificationTracker()
    ver_cb       = _make_verification_callback(laser_model, ver_tracker, opts.mass_kg)
    scheduler_cb = laser_link_scheduler_callback(laser_model)

    # Build save fields (SpaceAGORA defaults + verification extras)
    extra_fields = _build_verification_save_fields(ver_tracker, laser_model)
    all_fields   = vcat(SimulationModel.default_save_fields(args), extra_fields)

    # Run
    println("  Running: $(stem)  (J2=$(opts.use_j2), $(opts.orbits) orbits, " *
            "h_target=$(opts.target_altitude_km) km, h_helper=$(opts.helper_altitude_km) km, " *
            "N_helpers=$(opts.helpers))")
    elapsed = @elapsed begin
        result = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
            run_simulation(
                args;
                isolate_state       = false,
                return_solution     = true,
                return_solver_metadata = true,
                extra_callbacks     = (ver_cb, scheduler_cb),
                save_fields         = all_fields,
            )
        end
    end

    feather_path = joinpath(results_dir, "simulation_results.feather")
    println("  → Feather saved: $feather_path  ($(round(elapsed; digits=1)) s)")
    println("  → OCL activations: $(laser_model.link_activation_count), " *
            "active steps: $(laser_model.active_link_step_count)")
    println("  → W_OCL final: $(round(ver_tracker.W_ocl; sigdigits=5)) J")
    println("  → |ε_H| final: $(round(norm([ver_tracker.eps_H_x,
                                            ver_tracker.eps_H_y,
                                            ver_tracker.eps_H_z]); sigdigits=5)) kg·m²/s")

    # Write a small metadata TOML alongside the feather for the post-processor
    meta_path = joinpath(results_dir, "verification_meta.toml")
    open(meta_path, "w") do io
        println(io, "# Verification run metadata")
        println(io, "use_j2              = $(opts.use_j2)")
        println(io, "helper_altitude_km  = $(opts.helper_altitude_km)")
        println(io, "target_altitude_km  = $(opts.target_altitude_km)")
        println(io, "inclination_deg     = $(opts.inclination_deg)")
        println(io, "orbits              = $(opts.orbits)")
        println(io, "laser_range_km      = $(opts.laser_range_km)")
        println(io, "laser_power_w       = $(opts.laser_power_w)")
        println(io, "magnification       = $(opts.magnification)")
        println(io, "beta                = $(opts.beta)")
        println(io, "eta                 = $(opts.eta)")
        println(io, "mass_kg             = $(opts.mass_kg)")
        println(io, "dt_max_s            = $(opts.dt_max_s)")
        println(io, "target_period_s     = $(target_period_s)")
        println(io, "W_ocl_final_J       = $(ver_tracker.W_ocl)")
        println(io, "eps_H_x_final       = $(ver_tracker.eps_H_x)")
        println(io, "eps_H_y_final       = $(ver_tracker.eps_H_y)")
        println(io, "eps_H_z_final       = $(ver_tracker.eps_H_z)")
        println(io, "activations         = $(laser_model.link_activation_count)")
        println(io, "active_link_steps   = $(laser_model.active_link_step_count)")
        println(io, "retcode             = \"$(result.solution.retcode)\"")
    end

    return results_dir
end

# ─────────────────────────────────────────────────────────────────────────────
# Usage / argument parsing
# ─────────────────────────────────────────────────────────────────────────────
function _ver_usage()
    return """
    Usage:
      julia --project=. ORACLE/Paper_Numerical_Verification_Code/run_verification.jl [options]

    Options:
      --helpers N               Number of helper spacecraft [default: 10]
      --helper-altitude-km KM   Helper altitude [default: 1000.0]
      --target-altitude-km KM   Target altitude [default: 1100.0]
      --inclination-deg DEG     Common inclination [default: 0.0]
      --orbits N                Number of target orbits to simulate [default: 1.5]
      --laser-range-km KM       Maximum OCL range [default: 200.0]
      --laser-power-w W         OCL laser power [default: 10000.0]
      --magnification B         OCL thrust magnification [default: 100.0]
      --beta VALUE              Geometric loss factor [default: 1.0]
      --eta VALUE               Additional efficiency factor [default: 1.0]
      --mass-kg KG              Spacecraft mass [default: 227.0]
      --dt-max-s SEC            Maximum ODE step size [default: 10.0]
      --output-dir PATH         Output directory [default: output/Paper_Numerical_Verification_Code]
      --no-j2                   Run only the J2-inactive case (default: run both)
      --j2-only                 Run only the J2-active case
    """
end

function _ver_parse_options(argv)
    opts = Dict{Symbol,Any}(:run_no_j2 => true, :run_with_j2 => true)
    i = 1
    while i <= length(argv)
        arg = argv[i]
        arg in ("--help", "-h") && (println(_ver_usage()); exit(0))
        if arg == "--no-j2"
            opts[:run_with_j2] = false; i += 1; continue
        elseif arg == "--j2-only"
            opts[:run_no_j2] = false; i += 1; continue
        end
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'."))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        val = argv[i + 1]
        float_keys = (:helper_altitude_km, :target_altitude_km, :inclination_deg,
                      :orbits, :laser_range_km, :laser_power_w, :magnification,
                      :beta, :eta, :mass_kg, :dt_max_s)
        int_keys   = (:timeseries_points, :helpers)
        if key in float_keys
            opts[key] = parse(Float64, val)
        elseif key in int_keys
            opts[key] = parse(Int, val)
        elseif key == :output_dir
            opts[key] = abspath(val)
        else
            throw(ArgumentError("Unknown option: $arg"))
        end
        i += 2
    end
    return opts
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
function main(argv=ARGS)
    parsed = _ver_parse_options(argv)
    run_no_j2   = pop!(parsed, :run_no_j2)
    run_with_j2 = pop!(parsed, :run_with_j2)

    # Build base options from parsed args
    base_pairs = [(k, v) for (k, v) in parsed
                  if k in fieldnames(VerificationOptions)]
    base_opts  = VerificationOptions(; base_pairs...)

    output_dirs = String[]

    println("\n=== Section IV.F Conservation-Law Numerical Validation ===\n")
    println("Parameters (Table 2):")
    @printf("  Laser range:      %.0f km\n", base_opts.laser_range_km)
    @printf("  Laser power:  %6.0f W\n",   base_opts.laser_power_w)
    @printf("  Magnification:    %.0f\n",     base_opts.magnification)
    @printf("  β (geometric):    %.1f\n",     base_opts.beta)
    @printf("  Spacecraft mass:  %.0f kg\n",  base_opts.mass_kg)
    @printf("  Helper altitude:  %.0f km\n",  base_opts.helper_altitude_km)
    @printf("  Target altitude:  %.0f km\n",  base_opts.target_altitude_km)
    @printf("  Number of helpers:%d\n",       base_opts.helpers)
    @printf("  Orbits:           %.2f\n\n",   base_opts.orbits)

    if run_no_j2
        opts_no_j2 = VerificationOptions(;
            (f => getfield(base_opts, f) for f in fieldnames(VerificationOptions))...,
            use_j2 = false,
        )
        push!(output_dirs, run_verification_case(opts_no_j2))
    end

    if run_with_j2
        opts_j2 = VerificationOptions(;
            (f => getfield(base_opts, f) for f in fieldnames(VerificationOptions))...,
            use_j2 = true,
        )
        push!(output_dirs, run_verification_case(opts_j2))
    end

    println("\nDone.  Feather files written to:")
    for d in output_dirs
        println("  $d")
    end
    println("\nRun postprocess_verification.jl to compute residuals and produce plots.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
