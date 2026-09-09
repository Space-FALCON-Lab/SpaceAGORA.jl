using Test
using SpaceAGORA
using TOML
const SEng = SpaceAGORA.SimulationEngine

# The calibration cache mirrors ONE store file, and _rhs_calib_save! rewrites
# that file from the cache wholesale. Before the loaded-path guard, a test that
# redirected SPACEAGORA_RHS_CALIBRATION_PATH to a temp file, populated the cache
# from it and let the path revert would have the next real-path save overwrite
# the real store with the temp file's contents -- which truncated the real
# store from 82 verdicts to a handful twice on 2026-09-08.

_reset!() = lock(SEng._rhs_calib_lock) do
    empty!(SEng._rhs_calib_cache)
    SEng._rhs_calib_loaded[] = false
    SEng._rhs_calib_loaded_path[] = ""
end

function _write_store(path, sigs)
    rows = [Dict("signature" => s, "mode" => "satellite_batch", "allotment" => 1,
                 "scheduler" => "auto", "elapsed_mean_ns" => 1.0) for s in sigs]
    open(path, "w") do io
        TOML.print(io, Dict("schema_version" => 1, "calibrations" => rows))
    end
end
_sigs(path) = Set(String(r["signature"]) for r in TOML.parsefile(path)["calibrations"])

@testset "a store path change reloads the cache instead of carrying the old file along" begin
    mktempdir() do dir
        real = joinpath(dir, "real.toml"); tmp = joinpath(dir, "tmp.toml")
        _write_store(real, ["real_a", "real_b", "real_c"])
        _write_store(tmp, ["tmp_x"])
        _reset!()
        try
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => real) do
                @test SEng._rhs_calib_lookup("real_a") !== nothing
            end
            # Redirect: the cache follows the path, so the real entries are not
            # visible here and the temp entry is.
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => tmp) do
                @test SEng._rhs_calib_lookup("real_a") === nothing
                @test SEng._rhs_calib_lookup("tmp_x") !== nothing
                SEng._rhs_calib_store!("tmp_y", SEng._make_calib_satellite_batch_plan(), 2.0)
                SEng._rhs_calib_save!()
                @test _sigs(tmp) == Set(["tmp_x", "tmp_y"])
            end
            # Back on the real path: a store + save adds to the real file rather
            # than replacing it with the temp cache.
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => real) do
                @test SEng._rhs_calib_lookup("tmp_x") === nothing
                SEng._rhs_calib_store!("real_d", SEng._make_calib_satellite_batch_plan(), 3.0)
                SEng._rhs_calib_save!()
                @test _sigs(real) == Set(["real_a", "real_b", "real_c", "real_d"])
            end
            @test _sigs(tmp) == Set(["tmp_x", "tmp_y"])
        finally
            _reset!()
        end
    end
end

@testset "a save never writes a cache loaded from another file" begin
    mktempdir() do dir
        real = joinpath(dir, "real.toml"); tmp = joinpath(dir, "tmp.toml")
        _write_store(real, ["real_a", "real_b"])
        _write_store(tmp, ["tmp_x"])
        _reset!()
        try
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => tmp) do
                @test SEng._rhs_calib_lookup("tmp_x") !== nothing
            end
            # The cache mirrors tmp; a save with the path pointing at real is refused.
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => real) do
                SEng._rhs_calib_save!()
                @test _sigs(real) == Set(["real_a", "real_b"])
            end
        finally
            _reset!()
        end
    end
end

@testset "a cache marked loaded by hand belongs to the current path" begin
    mktempdir() do dir
        path = joinpath(dir, "injected.toml")
        _reset!()
        try
            withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => path) do
                lock(SEng._rhs_calib_lock) do
                    SEng._rhs_calib_loaded[] = true
                    SEng._rhs_calib_cache["inj"] = Dict{String, Any}("mode" => "satellite_batch",
                        "allotment" => 1, "scheduler" => "auto", "elapsed_mean_ns" => 1.0)
                end
                @test SEng._rhs_calib_lookup("inj") !== nothing      # not wiped by the load
                SEng._rhs_calib_save!()
                @test _sigs(path) == Set(["inj"])
            end
        finally
            _reset!()
        end
    end
end
