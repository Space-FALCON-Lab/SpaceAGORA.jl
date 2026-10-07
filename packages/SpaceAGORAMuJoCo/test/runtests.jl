# Run with the SpaceAGORAMuJoCo project active (see the package README). Skips cleanly when the pinned
# libmujoco cannot be loaded (unsupported platform or no network for the lazy artifact).
using Test
using SpaceAGORAMuJoCo

const _mujoco_ok = try
    SpaceAGORAMuJoCo.Binding.version() == 3_011_000
catch err
    # CI sets SPACEAGORA_MUJOCO_REQUIRED=1 so a failed artifact download fails the job instead of skipping.
    get(ENV, "SPACEAGORA_MUJOCO_REQUIRED", "") == "1" && rethrow()
    @warn "MuJoCo 3.11.0 is not loadable; skipping SpaceAGORAMuJoCo tests" exception = (err, catch_backtrace())
    false
end

if _mujoco_ok
    @testset "SpaceAGORAMuJoCo" begin
        include("common.jl")
        include("binding_tests.jl")
        include("scene_tests.jl")
        include("engine_tests.jl")
    end
else
    @testset "SpaceAGORAMuJoCo (skipped: MuJoCo unavailable)" begin
        @test_skip false
    end
end
