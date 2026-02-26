const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const EXAMPLES_DIR = joinpath(REPO_ROOT, "src", "examples")
const PROJECT_PATH = joinpath(REPO_ROOT, ".AGORA")

using CSV
using DataFrames
using LinearAlgebra
using StaticArrays
using Test

Base.@kwdef struct RegressionCase
    file::String
    mission_time::Float64 = 120.0
    expected_rows::Int
    rows_atol::Int = 0
    expected_t_final::Float64
    t_atol::Float64 = 1e-9
    expected_pos_norm::Float64
    pos_rtol::Float64 = 1e-6
    expected_vel_norm::Float64
    vel_rtol::Float64 = 1e-6
    expected_mass::Float64
    mass_rtol::Float64 = 1e-8
    expected_q_norm::Union{Nothing, Float64} = nothing
    q_atol::Float64 = 1e-8
end

const CASES = [
    RegressionCase(
        file="Earth_RW_Test.jl",
        mission_time=120.0,
        expected_rows=1447,
        rows_atol=8,
        expected_t_final=120.0,
        t_atol=1e-8,
        expected_pos_norm=7.377853166890683e6,
        pos_rtol=1e-6,
        expected_vel_norm=7273.739232146093,
        vel_rtol=1e-6,
        expected_mass=328.0,
        mass_rtol=1e-10,
        expected_q_norm=1.0,
        q_atol=1e-8
    ),
    RegressionCase(
        file="AGORA_Odyssey.jl",
        mission_time=120.0,
        expected_rows=67,
        rows_atol=4,
        expected_t_final=120.0,
        t_atol=1e-8,
        expected_pos_norm=2.819489491576329e7,
        pos_rtol=1e-6,
        expected_vel_norm=604.0083429788058,
        vel_rtol=1e-6,
        expected_mass=461.0,
        mass_rtol=1e-10
    )
]

function run_case(case::RegressionCase)
    output = IOBuffer()
    example_path = joinpath(EXAMPLES_DIR, case.file)
    isfile(example_path) || error("Missing example: $(case.file)")

    mktempdir() do tmp
        cmd = `$(Base.julia_cmd()) --startup-file=no --depwarn=error --project=$(PROJECT_PATH) $(example_path)`
        cmd = Cmd(cmd; dir=tmp)
        cmd = addenv(
            cmd,
            "SPACEAGORA_EXAMPLE_SMOKE" => "1",
            "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME" => string(case.mission_time),
            "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => "1",
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )

        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        if !success(proc)
            println(text)
            error("Example execution failed for $(case.file)")
        end
        if occursin("init_NaN", text) || occursin("First function call produced NaNs", text)
            println(text)
            error("NaN initialization warning found for $(case.file)")
        end

        csv_path = joinpath(tmp, "simulation_results.csv")
        isfile(csv_path) || error("Missing simulation_results.csv for $(case.file)")

        df = CSV.read(csv_path, DataFrame)
        n = nrow(df)
        t_final = Float64(df.time[end])
        pos_norm = norm(SVector{3, Float64}(Float64(df.sc1_pos_1[end]), Float64(df.sc1_pos_2[end]), Float64(df.sc1_pos_3[end])))
        vel_norm = norm(SVector{3, Float64}(Float64(df.sc1_vel_1[end]), Float64(df.sc1_vel_2[end]), Float64(df.sc1_vel_3[end])))
        mass = Float64(df.sc1_mass[end])

        @test abs(n - case.expected_rows) <= case.rows_atol
        @test isapprox(t_final, case.expected_t_final; atol=case.t_atol, rtol=0.0)
        @test isapprox(pos_norm, case.expected_pos_norm; rtol=case.pos_rtol, atol=0.0)
        @test isapprox(vel_norm, case.expected_vel_norm; rtol=case.vel_rtol, atol=0.0)
        @test isapprox(mass, case.expected_mass; rtol=case.mass_rtol, atol=0.0)

        if case.expected_q_norm !== nothing
            @test all(col -> col in names(df), ("sc1_q_1", "sc1_q_2", "sc1_q_3", "sc1_q_4"))
            q_norm = norm(SVector{4, Float64}(Float64(df.sc1_q_1[end]), Float64(df.sc1_q_2[end]), Float64(df.sc1_q_3[end]), Float64(df.sc1_q_4[end])))
            @test isapprox(q_norm, case.expected_q_norm; atol=case.q_atol, rtol=0.0)
        end

        println("regression_ok $(case.file): rows=$(n), tf=$(t_final), |r|=$(pos_norm), |v|=$(vel_norm), m=$(mass)")
    end
end

@testset "Example Regression Cases" begin
    for case in CASES
        run_case(case)
    end
end

println("examples_regression_ok")
