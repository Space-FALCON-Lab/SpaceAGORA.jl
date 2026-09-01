# Shared case construction for the contention probes.
#
# Extracted so the thread-ladder probe measures the SAME workload the matrix
# does. A second copy of a case definition drifts, and this repo has already
# recorded a probe that reported a regression which had been fixed an hour
# earlier because it carried its own copy of the profile table.

const _SUPPORT_REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const _SUPPORT_HARM_FILE =
    joinpath(_SUPPORT_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
const _SUPPORT_N = parse(Int, get(ENV, "PCM_N", "256"))
const _SUPPORT_L = parse(Int, get(ENV, "PCM_L", "20"))
# One hour, matching heavy_1024sat_fullstack_1hr. The span matters: the
# ephemeris look-ahead caches are built over the mission, so a probe that
# steps the epoch past mission_time misses every one of them and measures live
# SPICE traffic no solve would pay.
const _SUPPORT_MISSION_S = parse(Float64, get(ENV, "PCM_MISSION_S", "3600.0"))

function _support_make_sc(planet; ra_alt_m, rp_alt_m, ν_deg)
    root = SM.Link{0}(root = true, m = 500.0, ref_area = 12.0)
    panel = SM.Link{0}(root = false, m = 30.0, ref_area = 6.0,
                       r = MVector{3, Float64}(0.0, 1.2, 0.0))
    ic = SM.InitialCondition(ra = planet.Rp_e + ra_alt_m, rp = planet.Rp_e + rp_alt_m,
                             i = 35.0, ω = 40.0, Ω = 10.0, ν = ν_deg)
    return SM.SpacecraftModel(SM.Joint[], [root, panel], root, true,
                              root.m + panel.m, 0.0, root.inertia, 0, 0, ic, 1)
end

function build_probe_case(planet, effectors, density_model, ephemerides_model;
                          n_sats = _SUPPORT_N, ra_alt_m = 200e3, rp_alt_m = 160e3)
    scs = [_support_make_sc(planet; ra_alt_m = ra_alt_m + 10e3 * (i % 17),
                            rp_alt_m = rp_alt_m + 5e3 * (i % 13),
                            ν_deg = (120.0 + 3.0 * i) % 360.0) for i in 1:n_sats]
    env = SM.EnvironmentModel(
        planet = planet, EI = 120.0,
        density_model = density_model,
        ephemerides_model = ephemerides_model,
        thermal_model = SM.MaxwellianHeat(thermal_accomodation_factor = 1.0, planet = planet),
        topography = false, wind = false)
    args = SM.SimulationConfiguration(
        simulation_settings = SM.SimulationSettings(results = false, verbose = false,
                                                    generate_plots = false, normalize = false),
        mission_configuration = SM.MissionConfiguration(
            mission_type = SM.MissionTime, keplerian = true, number_of_orbits = 1,
            mission_time = _SUPPORT_MISSION_S, orientation_sim = false, num_steps_to_save = 10),
        environment_model = env,
        dynamics_model = SM.DynamicsModel(scs, effectors),
        guidance_model = SM.GuidanceModel(guidance_effectors = (), guidance_rates = Float64[]),
        navigation_model = SM.NavigationModel(navigation_effectors = (), navigation_rates = Float64[]),
        control_model = SM.ControlModel(control_effectors = (), control_rates = Float64[]),
        initial_time = SM.InitialTime(year = 2020, month = 1, day = 1, hour = 0, minute = 0, second = 0.0),
        integration_tolerances = SM.IntegrationTolerances())
    return probe_params_from_args(args)
end

"""
    probe_params_from_args(args) -> (args, p, u0)

Build the ODE params and initial state for `args`, running the same
initialisation `execution.jl` does before its first RHS evaluation.

One place, because omitting any of it does not under-warm the probe -- it makes
the probe measure a configuration no solve runs. Two omissions have already been
caught this way: the save-cache buffers (whose absence races the flat queue) and
the ephemeris look-ahead caches (whose absence sent every frame transform live
to SPICE and reported the resulting lock traffic as a property of the workload).
Anything driving the RHS directly, including a config built by the benchmark
harness, goes through here.
"""
function probe_params_from_args(args)
    n_sats = length(args.dynamics_model.spacecraft)
    p = SM.ODEParams(n_sats = n_sats, args = args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_density_model_instances!(p)
    SE._initialize_density_cache_buffers!(p)
    SE._initialize_gram_isolated_pool_buffers!(p)
    SE._initialize_save_cache_buffers!(p)

    SE._initialize_aero_workspace_buffers!(p)
    SE._initialize_nbody_workspace_buffers!(p)
    SE._initialize_nbody_ephemeris_cache_buffer!(p)
    SE._initialize_srp_sun_cache_buffer!(p)
    SE._initialize_planet_frame_cache_buffer!(p)
    SE._initialize_spice_rhs_memo_mode!(p)
    SE._reset_spice_rhs_memo!(p)
    SE._initialize_runtime_env_config!(p)

    # The ephemeris look-ahead caches, populated exactly as execution.jl does.
    #
    # Omitting these does not merely under-warm the probe -- it measures a
    # configuration no solve runs. A real solve pre-samples the planet frame,
    # the Sun vector and third-body states onto a grid and interpolates, so its
    # RHS touches SPICE rarely; a probe without them calls pxform/spkpos per
    # satellite per pass and reports lock occupancy that belongs to the probe
    # rather than to the workload.
    ephemerides_model = args.environment_model.ephemerides_model
    et_start = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    p.shared_buffers.et_start[] = et_start
    args.environment_model.planet.L_PI .=
        SM.planet_frame_lpi(args.environment_model.planet, et_start, ephemerides_model)
    mission_end = Float64(args.mission_configuration.mission_time)
    SE._initialize_nbody_ephemeris_cache!(p, et_start, mission_end)
    SE._initialize_srp_sun_ephemeris_cache!(p, et_start, mission_end)
    SE._initialize_planet_frame_ephemeris_cache!(p, et_start, mission_end)

    return (args = args, p = p, u0 = SE.build_initial_conditions(args))
end


# The fullstack shape, named once so the matrix and the thread ladder cannot
# drift apart. This is the configuration `heavy_1024sat_fullstack_1hr` uses:
# harmonics + SRP + aero over an analytic atmosphere on a SPICE-backed planet.
function build_fullstack_case(; n_sats = _SUPPORT_N)
    planet = SM.Earth(_SUPPORT_HARM_FILE)
    return build_probe_case(
        planet,
        (SM.GravitationalHarmonicsModel(_SUPPORT_L, _SUPPORT_L, _SUPPORT_HARM_FILE, planet),
         SM.SolarRadiationPressureModel(1.2, 12.0),
         SM.AerodynamicCoefficientfM()),
        SM.ExponentialAtmosphereModel(planet),
        SM.SpiceEphemeridesModel();
        n_sats = n_sats)
end


# The two shapes the scaling question is actually about: one that plateaus and
# one that inverts, differing by 48x in allocation per RHS pass (16 KiB against
# 790 KiB at N=256) and by almost nothing else -- same planet, same harmonics
# degree, same satellite count, same working set relative to cache.
function build_gravity_case(; n_sats = _SUPPORT_N)
    planet = SM.Earth(_SUPPORT_HARM_FILE)
    return build_probe_case(
        planet,
        (SM.GravitationalHarmonicsModel(_SUPPORT_L, _SUPPORT_L, _SUPPORT_HARM_FILE, planet),),
        SM.NoAtmosphereModel(), SM.SpiceEphemeridesModel();
        n_sats = n_sats, ra_alt_m = 400e3, rp_alt_m = 380e3)
end

const PROBE_CASES = Dict(
    "fullstack" => build_fullstack_case,
    "gravity"   => build_gravity_case,
)


# The benchmark harness's own cases, driven through the same probe.
#
# Reconstructing a case by hand is how a comparison against a published number
# turns into a comparison against an approximation of it: a hand-built aero
# shape came out 1.9x cheaper per RHS call than the figure it was being checked
# against, which is enough to make a non-reproduction meaningless. This builds
# the harness's configuration with the harness's own code.
function build_harness_case(case_name::String; profile::String = "full")
    root = normpath(joinpath(@__DIR__, "..", "benchmarks", "studies",
                             "parallelization_performance"))
    isdefined(Main, :ppc_single_config) || begin
        Base.include(Main, joinpath(root, "cli.jl"))
        Base.include(Main, joinpath(root, "cases.jl"))
    end
    cfg = Base.invokelatest(getfield(Main, :PPCConfig); profile = profile)
    args = Base.invokelatest(getfield(Main, :ppc_single_config), case_name, cfg)
    return probe_params_from_args(args)
end


# GRAMSuite is a weak dependency living in a vendored submodule, so it needs its
# path pushed onto LOAD_PATH before `import` resolves. Callers must invoke this
# at TOP LEVEL, before any function that constructs a GRAM model: the ext's
# keyword constructors are added to a newer world than a call already in flight,
# and `@eval import` inside the same top-level statement leaves them
# unreachable. Two probe scripts hit that separately.
function ensure_gramsuite_for_probe!()
    isdefined(SpaceAGORA, :GRAMSuite) && return true
    vendored = normpath(joinpath(@__DIR__, "..", "data", "GRAMSuite.jl"))
    try
        if Base.find_package("GRAMSuite") === nothing && isdir(vendored)
            pushfirst!(LOAD_PATH, vendored)
        end
        @eval Main import GRAMSuite
        return true
    catch err
        @info "GRAMSuite not loadable." exception = err
        return false
    end
end

# A live-GRAM density shape: the one workload measured to hold the shared native
# lock for most of an RHS pass (rho = 0.89), and therefore the one where an
# Amdahl width cap should bind.
function build_gram_case(; n_sats = _SUPPORT_N)
    planet = SM.Earth(_SUPPORT_HARM_FILE)
    return build_probe_case(
        planet,
        (SM.GravitationalHarmonicsModel(_SUPPORT_L, _SUPPORT_L, _SUPPORT_HARM_FILE, planet),
         SM.AerodynamicCoefficientfM()),
        Base.invokelatest(SM.GRAMAtmosphereModel; planet_name = "earth"),
        SM.SpiceEphemeridesModel();
        n_sats = n_sats, ra_alt_m = 200e3, rp_alt_m = 160e3)
end
