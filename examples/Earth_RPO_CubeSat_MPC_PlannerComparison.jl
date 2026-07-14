include(joinpath(@__DIR__, "Earth_RPO_CubeSat_MPC_Batch.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    run_rpo_cubesat_mpc_planner_comparison_batch_cli(copy(ARGS))
end
