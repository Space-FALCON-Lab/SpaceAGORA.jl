#!/usr/bin/env julia
#
# Paired in-process A/B between two routing profiles.
#
#   julia --project=. --threads=12 scripts/paired_profile_probe.jl \
#       --cases=gravity_4096_l50,heavy_1024_l50_6hr --a=full_smart --b=inner_only --pairs=15
#
# Exists because the benchmark harness cannot resolve the differences that are
# left. It spawns a fresh Julia process per (case, mode, thread) point -- ~40 s
# of startup against a few seconds of solving -- and runs each mode as a block,
# so drift between blocks is free to masquerade as effect. On this machine the
# per-solve spread is 15-45% while the residuals under investigation are 3-10%,
# and medians of block-ordered samples cannot separate those.
#
# Two changes fix that without needing more samples:
#
#   INTERLEAVING. Profiles alternate A,B,A,B and the reduction is over per-pair
#   ratios, so anything that drifts slowly -- thermal state, frequency, a
#   background process -- scales both halves of a pair together and cancels. The
#   order within each pair alternates too, so a systematic first-position or
#   second-position bias cancels across pairs rather than accumulating. This is
#   not a refinement: the same thermal-callback comparison measured 11.1% block
#   ordered and 1.7% paired, and the 11.1% was entirely ordering.
#
#   A SIGNIFICANCE VERDICT. The reduction is a two-sided sign test on pair wins,
#   not a difference of medians. It needs no distributional assumption, is
#   immune to a single wild pair, and -- the point -- can return "too close to
#   call", which a median cannot. A point estimate of +6% invites a fix; "not
#   distinguishable from zero at n=15" correctly does not.
#
# One process, so setup is paid once. Thread count is fixed at process start, so
# run once per budget rather than passing a ladder.

using SpaceAGORA
using StaticArrays
using Printf

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PC = SM.ParallelCost
const EARTH = SM.Earth()
const HARM = joinpath(@__DIR__, "..", "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

# Profile environments come from the shipped definition, not a local copy.
#
# The first version of this file duplicated the env table from
# benchmarks/studies/parallelization_performance/modes.jl, on the reasoning that
# importing that file drags in the harness's PPCConfig machinery and the point
# here is a process that starts fast. That duplication drifted within the hour:
# R5's scheduler was changed in profile_definitions.jl and this probe kept
# measuring the old value, reporting a regression that had already been fixed.
#
# ParallelProfiles.profile_env_pairs is what the shipped profile actually
# resolves to, so deriving from it cannot drift. preserve_existing=false so the
# profile's own values win rather than inheriting whatever the caller's shell
# happens to have set.
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
    "outer_process" => "R1_b",
    "inner_only" => "R2",
    "outer_inner_static" => "R3",
    "outer_inner_adaptive" => "R4",
    "full_smart" => "R5",
    "policy_v2" => "R6",
)

function mode_env(name::String)::Dict{String, String}
    haskey(_MODE_PROFILE, name) || error("unknown mode '$name'")
    profile = _MODE_PROFILE[name]
    pairs = SpaceAGORA.ParallelProfiles.profile_env_pairs(profile; preserve_existing=false)
    envd = Dict{String, String}(String(k) => String(v) for (k, v) in pairs)
    # Pre-solve plan calibration follows the profile's adaptive flag, matching
    # modes.jl: the static profiles measure one fixed route, the adaptive ones
    # are meant to exercise everything they ship with.
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

# The three original cases all land on calibration signatures that HAVE a stored
# plan, so none of them re-sweeps and none of them exercises the no-regret
# floor's retain-the-heuristic outcome. The light_* cases below sit in the
# signature buckets that carry no stored entry -- small constellations and
# three-effector shapes -- which is where a per-run sweep is both largest
# relative to the solve and least likely to buy anything.
function build(case::String)
    n, L, M, mission, effs_kind = if case == "gravity_4096_l50"
        4096, 50, 50, 3600.0, :harm
    elseif case == "heavy_1024_l50_6hr"
        1024, 50, 50, 21600.0, :harm
    elseif case == "interact_256_aero"
        256, 50, 50, 3600.0, :harm_aero
    elseif case == "light_64_aero"          # 33_64 | effs=2
        64, 50, 50, 3600.0, :harm_aero
    elseif case == "light_16_harm"          # 9_16 | effs=1 | harm=1
        16, 50, 50, 3600.0, :harm
    elseif case == "light_256_aero_j2"      # 129_256 | effs=3
        256, 50, 50, 3600.0, :harm_aero_j2
    # Collides with interact_256_aero under the v3 signature: same satellite
    # bucket, same effector COUNT, same harm=0 flag, entirely different work.
    elseif case == "invsq_256_aero"         # 129_256 | effs=2, collides
        256, 50, 50, 3600.0, :invsq_aero
    else
        error("unknown case '$case'")
    end
    harm = SM.GravitationalHarmonicsModel(L, M, HARM, EARTH)
    effs, dens = if effs_kind === :harm_aero
        ((harm, SM.AerodynamicCoefficientfM()), SM.ExponentialAtmosphereModel(EARTH))
    elseif effs_kind === :invsq_aero
        ((SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientfM()),
         SM.ExponentialAtmosphereModel(EARTH))
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

# Per-arm calibration stores, wiped at start.
#
# The calibration TOML is shared mutable state that BOTH arms read and write,
# and the sweep's verdict is not stable run to run (near-ties resolve
# differently under noise), so a single shared file lets one arm's stored plan
# silently become the other arm's cache hit. Pointing each arm at its own file
# under a scratch dir -- and deleting both before the first sample -- makes the
# two arms independent and makes every run start from the same state, which is
# what a machine that has never seen this shape does.
function _arm_calibration_paths(tag::String)
    dir = mktempdir(; prefix = "spaceagora_probe_calib_")
    return (
        joinpath(dir, "calib_a_$(tag).toml"),
        joinpath(dir, "calib_b_$(tag).toml"),
    )
end

# Drop the in-process RHS-calibration memo so the next solve re-enters
# _calibrate_rhs_plan_if_needed! exactly as a fresh process would.
#
# The probe runs every sample in ONE process, and _rhs_calib_cache is a module
# global. So the pre-solve sweep runs on the first sample and every later sample
# takes a cache hit -- the probe structurally cannot see per-run calibration
# cost, it amortises it across all the pairs. The benchmark harness spawns a
# fresh process per point and so does every real user, so the regime the probe
# measures is the one regime in which that cost does not exist.
#
# This models the fresh-process regime honestly: the in-process memo is dropped
# but the on-disk TOML is left alone, so a signature that HAS a stored plan
# reloads it cheaply (one file read, shared across all cases) while a signature
# that has none re-sweeps. That is exactly what a second `julia` invocation on a
# previously-calibrated machine does.
function _drop_calibration_memo!()
    SE = SpaceAGORA.SimulationEngine
    empty!(SE._rhs_calib_cache)
    SE._rhs_calib_loaded[] = false
    return nothing
end

function solve_once(args, envd::Dict{String, String}; cold_calib::Bool=false)::Float64
    cold_calib && _drop_calibration_memo!()
    return withenv((k => v for (k, v) in envd)...) do
        t0 = time_ns()
        SpaceAGORA.run_simulation(args)
        return Float64(time_ns() - t0)
    end
end

# ── CLI ───────────────────────────────────────────────────────────────────────
cases = ["gravity_4096_l50"]
a_name, b_name, pairs = "full_smart", "inner_only", 15
# Applied on top of profile B, so a single mechanism can be isolated by running
# the same profile on both sides with one setting reverted on B. If B wins
# significantly, that setting is what the profile is paying for.
b_override = Dict{String, String}()
cold_calib = false
isolate_calib = false
for arg in ARGS
    arg == "--cold-calib" && (global cold_calib = true)
    arg == "--isolate-calib" && (global isolate_calib = true)
    startswith(arg, "--cases=") && (global cases = String.(split(split(arg, "=")[2], ",")))
    startswith(arg, "--a=")     && (global a_name = String(split(arg, "=")[2]))
    startswith(arg, "--b=")     && (global b_name = String(split(arg, "=")[2]))
    startswith(arg, "--pairs=") && (global pairs = parse(Int, split(arg, "=")[2]))
    if startswith(arg, "--b-override=")
        for kv in split(split(arg, "=", limit=2)[2], ";")
            isempty(kv) && continue
            k, v = split(kv, "=", limit=2)
            global b_override[String(k)] = String(v)
        end
    end
end

# A two-sided sign test on n pairs cannot return below 2/2^n, so fewer than six
# pairs can never reach 0.05 however lopsided the result: five wins out of five
# gives p = 0.0625. Refusing is better than reporting "not distinguishable" for
# a comparison the test was never able to resolve.
if pairs < 6
    error("--pairs must be at least 6; a sign test on $(pairs) pairs bottoms out at " *
          "p = $(round(2.0^(1 - pairs), digits=4)) and can never clear 0.05.")
end

env_a = mode_env(a_name)
env_b = merge(mode_env(b_name), b_override)
if isolate_calib
    calib_a, calib_b = _arm_calibration_paths(string(getpid()))
    env_a["SPACEAGORA_RHS_CALIBRATION_PATH"] = calib_a
    env_b["SPACEAGORA_RHS_CALIBRATION_PATH"] = calib_b
    @printf("isolated calibration stores: A=%s B=%s\n", calib_a, calib_b)
end
isempty(b_override) || @printf("B overrides: %s\n", join(["$k=$v" for (k,v) in b_override], " "))
@printf("paired probe: %s (A) vs %s (B), %d pairs, %d threads%s\n\n",
        a_name, b_name, pairs, Threads.nthreads(),
        cold_calib ? ", cold calibration memo" : "")

for case in cases
    args = build(String(case))
    solve_once(args, env_a); solve_once(args, env_b)      # warm both paths once

    ratios = Float64[]; wins_a = 0; wins_b = 0
    for i in 1:pairs
        ta, tb = if isodd(i)
            x = solve_once(args, env_a; cold_calib=cold_calib)
            y = solve_once(args, env_b; cold_calib=cold_calib)
            (x, y)
        else
            y = solve_once(args, env_b; cold_calib=cold_calib)
            x = solve_once(args, env_a; cold_calib=cold_calib)
            (x, y)
        end
        push!(ratios, ta / tb)
        ta < tb ? (wins_a += 1) : (tb < ta ? (wins_b += 1) : nothing)
    end

    sort!(ratios)
    med = ratios[cld(length(ratios), 2)]
    decided = wins_a + wins_b
    p = PC._sign_test_two_sided(max(wins_a, wins_b), decided)
    verdict = p <= 0.05 ? (wins_a > wins_b ? "A FASTER" : "B FASTER") : "NOT DISTINGUISHABLE"
    @printf("%-22s  A wins %2d / B wins %2d   median A/B = %.3f (%+.1f%%)  p = %.4f  -> %s\n",
            case, wins_a, wins_b, med, 100 * (med - 1), p, verdict)
end
