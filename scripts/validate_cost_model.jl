#!/usr/bin/env julia
#
# Cold validation of the analytic cost model against the REAL harmonics kernel.
#
#   julia --project=. --threads=12 scripts/validate_cost_model.jl
#
# The machine constants are measured entirely from synthetic kernels
# (scripts/calibrate_machine.jl). Nothing in this grid is fitted, which is the
# point: a good prediction here is evidence the synthetic constants transfer to
# real code, and a bad one is a localised residual that says where the model is
# wrong rather than being absorbed into a coefficient.
#
# The headline number is NOT prediction error. It is decision accuracy: how
# often argmin over predicted cost matches argmin over measured cost. A model
# with 30% error on absolute time but the right ranking everywhere is a complete
# success, because ranking is all the router consumes.

using SpaceAGORA
using StaticArrays
using Printf

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PC = SM.ParallelCost

const EARTH = SM.Earth()
const HARM_FILE = joinpath(@__DIR__, "..", "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

function make_sc(; ra_alt_m, rp_alt_m, ν_deg)
    root = SM.Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = SM.Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))
    ic = SM.InitialCondition(ra=EARTH.Rp_e + ra_alt_m, rp=EARTH.Rp_e + rp_alt_m,
                             i=35.0, ω=40.0, Ω=10.0, ν=ν_deg)
    return SM.SpacecraftModel(SM.Joint[], [root, panel], root, true,
                              root.m + panel.m, 0.0, root.inertia, 0, 0, ic, 1)
end

function build_case(n_sats::Int, L::Int, M::Int)
    scs = [make_sc(ra_alt_m=500e3 + 10e3 * (i % 17), rp_alt_m=480e3 + 5e3 * (i % 13),
                   ν_deg=(120.0 + 3.0 * i) % 360.0) for i in 1:n_sats]
    harm = SM.GravitationalHarmonicsModel(L, M, HARM_FILE, EARTH)
    env = SM.EnvironmentModel(
        planet=EARTH, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false, wind=false)
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
                                                  generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=30.0, orientation_sim=false, num_steps_to_save=10),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(scs, (harm,)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances())
    p = SM.ODEParams(n_sats=n_sats, args=args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_density_model_instances!(p)
    SE._initialize_density_cache_buffers!(p)
    SE._initialize_gram_isolated_pool_buffers!(p)
    SE._initialize_aero_workspace_buffers!(p)
    SE._initialize_nbody_workspace_buffers!(p)
    SE._initialize_runtime_env_config!(p)
    u0 = SE.build_initial_conditions(args)
    return (args=args, p=p, u0=u0, harm=harm)
end

# Time one candidate by pinning it and evaluating the real RHS.
function measure_candidate(case, cand::PC.PlanCandidate)
    plan = cand.mode === :satellite_batch ?
        SE._make_calib_satellite_batch_plan() :
        SE._make_calib_flat_plan(cand.allotment, cand.scheduler)
    du = zero(case.u0)
    case.p.shared_buffers.rhs_plan_override[] = plan
    try
        ns = PC.timed_min(k = 9, target_window_ns = 200_000) do
            SE.spacecraft_dynamics!(du, case.u0, case.p, 0.0)
            @inbounds du.sc[1][1]
        end
        return ns
    finally
        case.p.shared_buffers.rhs_plan_override[] = nothing
    end
end

mc = PC.load_machine_constants()
if mc === nothing
    println("No cached constants; calibrating now.")
    mc = PC.calibrate_machine()
    PC.save_machine_constants(mc)
end
canary = PC.constants_are_current(mc)
@printf("constants fingerprint=%s  canary drift=%+.1f%%  ok=%s\n\n",
        mc.fingerprint, 100 * canary.drift, canary.ok)

budget = Threads.nthreads()
grid = [(N, L, M) for N in (8, 64, 256, 1024)
                  for (L, M) in ((4, 4), (20, 0), (20, 20), (50, 0), (50, 50))]

@printf("%-6s %-4s %-4s | %-22s %-22s | %-8s %s\n",
        "N", "L", "M", "predicted best", "measured best", "ratio", "agree")
println(repeat("-", 100))

agree = 0; total = 0; ratios = Float64[]
for (N, L, M) in grid
    case = build_case(N, L, M)
    counts = PC.constellation_work_counts(case.args, N)
    cands = PC.plan_candidates(counts; budget=budget, n_active_sats=N)

    preds = [(c, PC.predict_plan_ns(counts, mc, c; n_active_sats=N, budget=budget).ns) for c in cands]
    meas  = [(c, measure_candidate(case, c)) for c in cands]

    pbest = preds[argmin(map(last, preds))][1]
    mbest = meas[argmin(map(last, meas))][1]
    mbest_ns = minimum(map(last, meas))
    # Cost of following the model's pick, relative to the true best.
    pbest_measured = meas[findfirst(x -> x[1] === pbest, meas)][2]
    regret = pbest_measured / mbest_ns

    ok = pbest == mbest
    global agree += ok ? 1 : 0
    global total += 1
    push!(ratios, regret)
    lbl(c) = c.mode === :satellite_batch ? "satellite_batch" :
        "flat(a=$(c.allotment),$(c.scheduler))"
    @printf("%-6d %-4d %-4d | %-22s %-22s | %-8.3f %s\n",
            N, L, M, lbl(pbest), lbl(mbest), regret, ok ? "yes" : "NO")
end

println()
@printf("decision accuracy      : %d/%d (%.0f%%)\n", agree, total, 100*agree/total)
@printf("median regret vs best  : %+.1f%%\n", 100*(sort(ratios)[cld(length(ratios),2)] - 1))
@printf("worst regret vs best   : %+.1f%%\n", 100*(maximum(ratios) - 1))
