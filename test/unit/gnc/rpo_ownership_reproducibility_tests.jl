using Test, SpaceAGORA, LinearAlgebra, Random, StaticArrays, JSON
const RC = SpaceAGORA.SimulationModel.ControlHooks
const RG = SpaceAGORA.SimulationModel.GuidanceHooks
const RM = SpaceAGORA.SimulationModel

fresh_controller() = RC.init_rpo_lqmpc(0.001, 1.0, Matrix{Float64}(I, 6, 6),
    10.0 .* Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 6, 6), 4)

@testset "RPO controller owns its copied solver" begin
    original = fresh_controller()
    original.U_prev .= range(0.0, 0.01; length=length(original.U_prev))
    graph = deepcopy((a=original, b=original, warm=original.U_prev))
    copied = graph.a
    @test copied === graph.b
    @test copied !== original
    @test copied.U_prev === graph.warm
    @test copied.U_prev == original.U_prev
    @test copied.U_prev !== original.U_prev
    @test copied.qp_model !== original.qp_model
    @test copied.qp_model.workspace != original.qp_model.workspace
    @test copied.qp_results !== original.qp_results
    @test copied.H == original.H && copied.H !== original.H
    @test copied.W == original.W && copied.W !== original.W

    # Explicit finalization models the reported lifetime failure without
    # depending on when a particular host's garbage collector happens to run.
    finalize(original.qp_model)
    GC.gc(true)
    x = [0.3, -0.2, 0.1, 0.01, -0.02, 0.03]
    ref = zeros(6, copied.horizon + 1)
    analytical = -(Matrix(copied.H) \ (copied.E * x))
    command = RC.rpo_lqmpc_control(copied, x, ref)
    @test all(isfinite, command)
    @test command ≈ analytical[1:3] atol=1e-4 rtol=1e-4
    @test copied.qp_results.info.status in (:Solved, :Solved_inaccurate)
    second = deepcopy(copied)
    @test second.qp_model.workspace != copied.qp_model.workspace
    saved_warm = copy(copied.U_prev)
    RC.rpo_lqmpc_control(second, -x, ref)
    @test copied.U_prev == saved_warm
    finalize(copied.qp_model)
    GC.gc(true)
    @test all(isfinite, RC.rpo_lqmpc_control(second, x, ref))
end

@testset "Particle random streams are independent of thread scheduling" begin
    geometry = RM.RPOReferenceGeometry(RM.RPOStationGeometry(
        [0.0 0.0 0.0; -0.5 0.0 0.5; 0.0 0.0 0.0]; keepout_radius_m=0.25);
        chaser=RM.RPOCubeSatGeometry(dims_m=(0.2, 0.2, 0.3)))
    config = RM.RPOPSOConfig(n_waypoints=2, n_particles=12, n_iters=5,
        adaptive_enable=false, sample_ds_m=0.25, safe_distance_m=0.5,
        search_margin_m=5.0, cost_ref_distance_m=10.0,
        iteration_runtime_limit_s=Inf)
    snapshots = []
    for seed in (741, 742, 743)
        run() = RM.rpo_pso_plan_path(SVector(-5.0, -1.0, 0.0),
            SVector(5.0, 1.0, 0.0), geometry, config; rng=MersenneTwister(seed))
        a, b = run(), run()
        @test a.path == b.path
        @test a.cost == b.cost
        @test a.cost_history == b.cost_history
        @test all(isfinite, a.path)
        @test !a.iteration_timed_out
        # Hexadecimal bits make cross-process/thread-count checks exact.
        bits(x) = string(reinterpret(UInt64, Float64(x)); base=16, pad=16)
        push!(snapshots, (seed=seed, shape=size(a.path), path=bits.(vec(a.path)),
            cost=bits(a.cost), history=bits.(a.cost_history)))
    end
    output = get(ENV, "SPACEAGORA_RPO_REPRO_SNAPSHOT", "")
    isempty(output) || write(output, JSON.json(snapshots))
end
