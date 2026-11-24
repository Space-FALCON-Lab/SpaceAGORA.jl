using Test, TOML

const CODE = normpath(joinpath(@__DIR__, "..", "..",
    "Output_8_Simulation_absolute_RN_range_based", "Code"))
include(joinpath(CODE, "run_grid.jl"))
include(joinpath(CODE, "VerifyOutput.jl"))
using .StudyVerification: verify_output

@testset "Range-based campaign grid and launcher" begin
    cases = grid_cases()
    expected = Set((h, i, n, :range_based) for h in (900, 950, 1000, 1050, 1100),
        i in (0, 0.5, 1), n in (1, 100, 200, 300, 400, 500, 600, 700, 800))
    @test Set((c.altitude_km, c.inclination_deg, c.helpers, c.policy) for c in cases) == expected
    @test length(unique(case_directory.(Ref(DEFAULT_OUTPUT), cases))) == length(cases) == 135
    @test all(c -> c.orbits == 100 && c.power_W == 10_000, cases)

    mktempdir(@__DIR__) do root
        launcher = joinpath(CODE, "launch_study.sh")
        run(`bash $launcher stop-after-case $root`)
        @test isfile(joinpath(root, "_batch_range_based", "STOP_AFTER_CASE"))
        run(`bash $launcher run $root`)
        status = TOML.parse(read(`bash $launcher status $root`, String))
        @test status["status"] == "stopped_after_case"
        @test status["campaign"] == "range_based" && status["total_cases"] == 135
        @test !isdir(joinpath(root, "paper_plot_mode"))
    end
end

@testset "Short range-based case writes verifiable results" begin
    mktempdir(@__DIR__) do root
        c = StudyCase(1000, 0, 1, :range_based; orbits=0.002)
        run(simulation_command(root, c))
        directory = case_directory(root, c)
        report = verify_output(directory)
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        @test report["passed"] && report["rows"] == 1001
        @test info["status"] == "verified" && info["settings"]["scheduler"] == "range_based"
        # The only helper stays in range; even the coplanar case must fire throughout.
        @test report["firing_intervals"] == 1 && info["actual_firing_s"] ≈ mission_duration(c)
    end
end
