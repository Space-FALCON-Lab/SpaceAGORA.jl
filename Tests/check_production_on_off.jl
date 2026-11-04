include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyConfig.jl"))
using .StudyConfig
import Pkg
Pkg.activate(LIBRARY_ROOT; io=devnull)
using Arrow, LinearAlgebra, Test, TOML

root = joinpath(@__DIR__, "results", "production_pipeline_compact")
outputs = TOML.parsefile(joinpath(root, "short_pilot_outputs.toml"))
on = Arrow.Table(joinpath(outputs["final_on"], "trajectory.feather"))
off = Arrow.Table(joinpath(outputs["final_zero"], "trajectory.feather"))
position(table) = [table[Symbol("sc2_pos_$k")][end] for k in 1:3]
difference = norm(position(on) - position(off))
@testset "Persisted production on/off response" begin
    @test on.time == off.time
    @test difference > 0.01
    @test any(!iszero, on.sc2_laser_dv_1)
    for name in propertynames(off)
        occursin("laser", string(name)) || startswith(string(name), "target_dv") || continue
        @test all(iszero, off[name])
    end
end
atomic_toml(joinpath(root, "on_off_verification.toml"), Dict(
    "passed" => true, "target_position_difference_m" => difference))
