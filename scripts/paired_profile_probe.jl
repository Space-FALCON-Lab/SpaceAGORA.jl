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

# Profile environments, mirroring benchmarks/studies/parallelization_performance/modes.jl.
# Kept local rather than imported because that file needs the harness's PPCConfig
# machinery, and the point here is a process that starts fast.
const _COMMON = Dict(
    "SPACEAGORA_PERF_PARALLEL_BACKEND" => "none",
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
    "SPACEAGORA_RHS_BATCH_PARALLEL" => "auto",
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
    "SPACEAGORA_SAVE_BUNDLE" => "0",
    "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
)

function mode_env(name::String)::Dict{String, String}
    inner(v) = Dict(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => v,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => v,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => v,
        "SPACEAGORA_MULTIBODY_PARALLEL" => v,
        "SPACEAGORA_EFFECTOR_PARALLEL" => v,
    )
    base = merge(_COMMON, Dict("SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
                               "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
                               "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
                               "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => "0",
                               "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "static",
                               "SPACEAGORA_RHS_CALIBRATE" => "off"))
    if name == "serial"
        return merge(base, inner("off"), Dict("SPACEAGORA_PARALLEL_PROFILE" => "R0"))
    elseif name == "inner_only"
        return merge(base, inner("auto"), Dict("SPACEAGORA_PARALLEL_PROFILE" => "R2"))
    elseif name == "outer_inner_static"
        return merge(base, inner("auto"), Dict("SPACEAGORA_PARALLEL_PROFILE" => "R3"))
    elseif name == "outer_inner_adaptive"
        return merge(base, inner("auto"),
                     Dict("SPACEAGORA_PARALLEL_PROFILE" => "R4",
                          "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
                          "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
                          "SPACEAGORA_RHS_CALIBRATE" => "auto"))
    elseif name == "full_smart"
        return merge(base, inner("auto"),
                     Dict("SPACEAGORA_PARALLEL_PROFILE" => "R5",
                          "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
                          "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "on",
                          "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "1",
                          "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "1",
                          "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => "1",
                          "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "1",
                          "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => "dynamic",
                          "SPACEAGORA_RHS_CALIBRATE" => "auto"))
    end
    error("unknown mode '$name'")
end

function make_sc(i, planet)
    root = SM.Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = SM.Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3,Float64}(0.0, 1.2, 0.0))
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
    else
        error("unknown case '$case'")
    end
    harm = SM.GravitationalHarmonicsModel(L, M, HARM, EARTH)
    effs, dens = effs_kind === :harm_aero ?
        ((harm, SM.AerodynamicCoefficientfM()), SM.ExponentialAtmosphereModel(EARTH)) :
        ((harm,), SM.NoAtmosphereModel())
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

function solve_once(args, envd::Dict{String, String})::Float64
    return withenv((k => v for (k, v) in envd)...) do
        t0 = time_ns()
        SpaceAGORA.run_simulation(args)
        return Float64(time_ns() - t0)
    end
end

# ── CLI ───────────────────────────────────────────────────────────────────────
cases = ["gravity_4096_l50"]
a_name, b_name, pairs = "full_smart", "inner_only", 15
for arg in ARGS
    startswith(arg, "--cases=") && (global cases = String.(split(split(arg, "=")[2], ",")))
    startswith(arg, "--a=")     && (global a_name = String(split(arg, "=")[2]))
    startswith(arg, "--b=")     && (global b_name = String(split(arg, "=")[2]))
    startswith(arg, "--pairs=") && (global pairs = parse(Int, split(arg, "=")[2]))
end

# A two-sided sign test on n pairs cannot return below 2/2^n, so fewer than six
# pairs can never reach 0.05 however lopsided the result: five wins out of five
# gives p = 0.0625. Refusing is better than reporting "not distinguishable" for
# a comparison the test was never able to resolve.
if pairs < 6
    error("--pairs must be at least 6; a sign test on $(pairs) pairs bottoms out at " *
          "p = $(round(2.0^(1 - pairs), digits=4)) and can never clear 0.05.")
end

env_a, env_b = mode_env(a_name), mode_env(b_name)
@printf("paired probe: %s (A) vs %s (B), %d pairs, %d threads\n\n", a_name, b_name, pairs, Threads.nthreads())

for case in cases
    args = build(String(case))
    solve_once(args, env_a); solve_once(args, env_b)      # warm both paths once

    ratios = Float64[]; wins_a = 0; wins_b = 0
    for i in 1:pairs
        ta, tb = if isodd(i)
            x = solve_once(args, env_a); y = solve_once(args, env_b); (x, y)
        else
            y = solve_once(args, env_b); x = solve_once(args, env_a); (x, y)
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
