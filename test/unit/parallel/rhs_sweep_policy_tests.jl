using Test
using SpaceAGORA
const SEng = SpaceAGORA.SimulationEngine

@testset "skip-by-prediction decision is pure and off by default" begin
    skip = SEng._rhs_calib_skip_by_prediction
    @test skip(0.0, 1e9, 1e8) === nothing           # ratio off
    @test skip(5.0, 0.0, 1e8) === nothing           # no solve prediction
    @test skip(5.0, 1e9, 0.0) === nothing           # no sweep prediction
    @test skip(5.0, 4e8, 1e8) === true              # 4e8 < 5e8
    @test skip(5.0, 5e8, 1e8) === false             # boundary sweeps
    @test skip(5.0, 1e9, 1e8) === false
end

@testset "predicted solve time uses solve_ns, then evals x recorded default call" begin
    pred = SEng._rhs_calib_predicted_solve_ns
    cands = Any[Dict("default" => false, "ns" => 9.0), Dict("default" => true, "ns" => 2.0)]
    @test pred(Dict("solve_ns" => 3e9), 100.0) == 3e9
    @test pred(Dict("candidates" => cands), 100.0) == 200.0
    @test pred(Dict("candidates" => cands), 0.0) == 0.0
    @test pred(Dict{String, Any}(), 100.0) == 0.0
end

@testset "skip ratio env gates the solve threshold" begin
    entry = Dict{String, Any}("mode" => "heuristic", "solve_ns" => 4e8, "sweep_ns" => 1e8)
    mktempdir() do dir
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(dir, "s.toml"),
                "SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S" => "1.0") do
            lock(SEng._rhs_calib_lock) do
                SEng._rhs_calib_loaded[] = true; SEng._rhs_calib_loaded_path[] = SEng._rhs_calib_path()
                empty!(SEng._rhs_calib_cache); SEng._rhs_calib_cache["sig"] = entry
            end
            try
                # 0.4 s solve < MIN_SOLVE_S: short either way with the ratio off.
                @test !SEng._rhs_calib_solve_exceeds_threshold("sig")
                # Ratio 2: 4e8 >= 2 * 1e8 -> long, which the fixed rule would call short.
                withenv("SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO" => "2") do
                    @test SEng._rhs_calib_solve_exceeds_threshold("sig")
                end
                withenv("SPACEAGORA_RHS_CALIBRATE_SKIP_RATIO" => "5") do
                    @test !SEng._rhs_calib_solve_exceeds_threshold("sig")
                end
            finally
                lock(SEng._rhs_calib_lock) do
                    empty!(SEng._rhs_calib_cache); SEng._rhs_calib_loaded[] = false; SEng._rhs_calib_loaded_path[] = ""
                end
            end
        end
    end
end
