#!/usr/bin/env julia
#
# Cold-cache trajectory of a routing profile's pre-solve RHS calibration.
#
#   julia --project=. --threads=8 scripts/r6_cold_cache_trajectory.jl \
#       --cases=gravity_4096_l50,interact_256_aero --modes=policy_v2,full_smart,outer_threads --runs=6
#
# The paper harness runs its modes in a fixed order against one shared
# calibration store, so by the time policy_v2 (R6) runs a case, the R4 and R5
# sweeps before it have already cast the consecutive "heuristic won" votes that
# R6's cached-verdict rule needs. A user running R6 on its own gets none of that
# help. This script models that user: every (case, mode) pair starts from an
# EMPTY calibration store and is solved `--runs` times in succession, with the
# in-process calibration memo dropped between solves so each one re-enters
# calibration exactly as a fresh `julia` invocation would. What is printed per
# run is the wall time and where the installed plan came from (sweep or cache,
# and which plan), so the trajectory -- sweep, sweep, sweep, then cache -- is
# visible run by run rather than averaged away.
#
# Each pair is warmed once beforehand against a throwaway store, so run 1 pays
# the sweep but not JIT. Thread count is fixed at process start; run once per
# budget.

using SpaceAGORA
using StaticArrays
using Printf

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const EARTH = SM.Earth()
const HARM = joinpath(@__DIR__, "..", "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

const _COMMON = Dict(
    "SPACEAGORA_PERF_PARALLEL_BACKEND" => "none",
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
    "SPACEAGORA_SAVE_BUNDLE" => "0",
    "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
)

const _MODE_PROFILE = Dict(
    "serial" => "R0",
    "outer_threads" => "R1_a",
    "inner_only" => "R2",
    "outer_inner_static" => "R3",
    "outer_inner_adaptive" => "R4",
    "full_smart" => "R5",
    "policy_v2" => "R6",
)

# Same derivation as scripts/paired_profile_probe.jl: the shipped profile
# definition, not a local copy of its env table.
function mode_env(name::String)::Dict{String, String}
    haskey(_MODE_PROFILE, name) || error("unknown mode '$name'")
    profile = _MODE_PROFILE[name]
    pairs = SpaceAGORA.ParallelProfiles.profile_env_pairs(profile; preserve_existing=false)
    envd = Dict{String, String}(String(k) => String(v) for (k, v) in pairs)
    adaptive = profile in ("R4", "R5", "R6")
    envd["SPACEAGORA_RHS_CALIBRATE"] = adaptive ? "auto" : "off"
    return merge(_COMMON, envd)
end

function make_sc(i, planet)
    root = SM.Link(root=true, m=500.0, ref_area=12.0)
    panel = SM.Link(root=false, m=30.0, ref_area=6.0, r=MVector{3,Float64}(0.0, 1.2, 0.0))
    ic = SM.InitialCondition(ra=planet.Rp_e + 500e3 + 10e3 * (i % 17),
                             rp=planet.Rp_e + 480e3 + 5e3 * (i % 13),
                             i=35.0, ω=40.0, Ω=10.0, ν=(120.0 + 3.0 * i) % 360.0)
    return SM.SpacecraftModel(SM.Joint[], [root, panel], root, true,
                              root.m + panel.m, 0.0, root.inertia, 0, 0, ic, 1)
end

function build(case::String)
    n, L, M, mission, effs_kind = if case == "gravity_4096_l50"
        4096, 50, 50, 3600.0, :harm
    elseif case == "heavy_1024_l50_6hr"
        1024, 50, 50, 21600.0, :harm
    elseif case == "interact_256_aero"
        256, 50, 50, 3600.0, :harm_aero
    elseif case == "light_256_aero_j2"
        256, 50, 50, 3600.0, :harm_aero_j2
    else
        error("unknown case '$case'")
    end
    harm = SM.GravitationalHarmonicsModel(L, M, HARM, EARTH)
    effs, dens = if effs_kind === :harm_aero
        ((harm, SM.AerodynamicCoefficientfM()), SM.ExponentialAtmosphereModel(EARTH))
    elseif effs_kind === :harm_aero_j2
        ((harm, SM.AerodynamicCoefficientfM(), SM.InverseSquaredJ2GravityModel()),
         SM.ExponentialAtmosphereModel(EARTH))
    else
        ((harm,), SM.NoAtmosphereModel())
    end
    env = SM.EnvironmentModel(planet=EARTH, EI=120.0, density_model=dens,
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
                                                  generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=mission, orientation_sim=false, num_steps_to_save=10),
        environment_model=env,
        dynamics_model=SM.DynamicsModel([make_sc(i, EARTH) for i in 1:n], effs),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances())
end

# Drop the in-process calibration memo so the next solve re-enters calibration
# as a fresh process would; the on-disk store is left alone, which is what
# carries the vote count between runs.
function _drop_calibration_memo!()
    empty!(SE._rhs_calib_cache)
    SE._rhs_calib_loaded[] = false
    return nothing
end

# One solve under `envd`. Returns wall time and the installed plan's origin,
# read from the solve's own metadata (the harness reads it the same way; a
# snapshot taken from this task afterwards would see the global context, which
# the engine's scoped context never wrote to).
function solve_once(args, envd::Dict{String, String})
    _drop_calibration_memo!()
    return withenv((k => v for (k, v) in envd)...) do
        t0 = time_ns()
        res = SpaceAGORA.run_simulation(args; return_solution=true, return_solver_metadata=true)
        el = (time_ns() - t0) / 1e9
        pp = res.parallel_policy
        src = pp === nothing ? "?" : string(pp.rhs_plan_source)
        mode = pp === nothing ? "?" : string(pp.rhs_plan_mode)
        allot = pp === nothing ? 0 : Int(pp.rhs_plan_allotment)
        return (el, src, mode, allot)
    end
end

cases = ["gravity_4096_l50", "interact_256_aero"]
modes = ["policy_v2", "full_smart", "outer_threads"]
runs = 6
for arg in ARGS
    startswith(arg, "--cases=") && (global cases = String.(split(split(arg, "=", limit=2)[2], ",")))
    startswith(arg, "--modes=") && (global modes = String.(split(split(arg, "=", limit=2)[2], ",")))
    startswith(arg, "--runs=") && (global runs = parse(Int, split(arg, "=", limit=2)[2]))
end
runs >= 4 || error("--runs must be at least 4 so the post-vote regime has at least one run")

const STATE_DIR = mktempdir(; prefix="spaceagora_cold_cache_")
# Isolated inner-policy and outer-route stores too, so this run neither reads
# nor writes the machine's persistent state.
const _ISOLATION = Dict(
    "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath(STATE_DIR, "inner_policy_state.toml"),
    "SPACEAGORA_OUTER_ROUTE_STATE_PATH"     => joinpath(STATE_DIR, "outer_route_state.toml"),
)

@printf("cold-cache trajectory: %d threads, %d runs per (case, mode), state under %s\n\n",
        Threads.nthreads(), runs, STATE_DIR)

for case in cases
    args = build(case)
    for mode in modes
        base = merge(mode_env(mode), _ISOLATION)
        # Warm-up against a throwaway store: JIT and one sweep, discarded.
        warm = merge(base, Dict("SPACEAGORA_RHS_CALIBRATION_PATH" =>
                                joinpath(STATE_DIR, "warm_$(case)_$(mode).toml")))
        solve_once(args, warm)
        fresh = merge(base, Dict("SPACEAGORA_RHS_CALIBRATION_PATH" =>
                                 joinpath(STATE_DIR, "trajectory_$(case)_$(mode).toml")))
        @printf("%-20s %-22s\n", case, mode)
        times = Float64[]
        for r in 1:runs
            el, src, pm, allot = solve_once(args, fresh)
            push!(times, el)
            @printf("    run %2d  %7.3fs  %-6s %-32s allot=%d\n", r, el, src, pm, allot)
        end
        early = sum(times[1:3]) / 3
        late  = sum(times[4:end]) / (runs - 3)
        @printf("    -> runs 1-3 mean %.3fs, runs 4-%d mean %.3fs (%+.1f%%)\n\n",
                early, runs, late, 100 * (late / early - 1))
    end
end
