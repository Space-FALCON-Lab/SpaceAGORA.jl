# Runs one ppc case once and prints a sha256 of the full final state plus the
# RHS plan the run used. Usage: julia --project=. rhs_seed_bitident_probe.jl <case> <outdir>
# Arm is selected by env (SPACEAGORA_RHS_CALIBRATE, ..._SEED_TOPK, ..._SEED_ACCEPT_WINNER).
using SpaceAGORA, SHA
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PPC = joinpath(ROOT, "benchmarks", "studies", "parallelization_performance")
for f in ("cli.jl", "cases.jl", "modes.jl", "execution.jl"); include(joinpath(PPC, f)); end
case, outdir = ARGS[1], ARGS[2]
mkpath(outdir)
cfg = PPCConfig(profile = "test")
args = ppc_single_config(case, cfg)
t0 = time()
sol = ppc_solve_once(args, cfg).solution
wall = time() - t0
state = Float64[x for x in sol.u[end]]
tvec = Float64.(sol.t)
write(joinpath(outdir, "state.bin"), state)
snap = SpaceAGORA.SimulationModel.ParallelPolicy.policy_telemetry_snapshot()
println("RESULT case=", case, " retcode=", sol.retcode, " wall_s=", round(wall; digits = 2),
        " nsteps=", length(tvec), " nstate=", length(state), " nsat=", length(sol.u[end].sc))
println("RESULT sha_state=", bytes2hex(sha256(reinterpret(UInt8, state))))
println("RESULT sha_t=", bytes2hex(sha256(reinterpret(UInt8, tvec))))
println("RESULT plan source=", snap.rhs_plan_source, " mode=", snap.rhs_plan_mode,
        " allotment=", snap.rhs_plan_allotment, " scheduler=", snap.rhs_plan_scheduler,
        " sweep_s=", round(snap.rhs_plan_sweep_ns / 1e9; digits = 3), " sweeps=", snap.rhs_plan_sweeps_total)
