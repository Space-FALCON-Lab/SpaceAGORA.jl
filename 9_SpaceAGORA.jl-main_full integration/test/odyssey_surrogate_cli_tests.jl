using Test

@testset "Odyssey CLI preserves downstream diagnostics" begin
    script = abspath(joinpath(@__DIR__, "..", "examples", "odyssey_surrogate.jl"))
    project = dirname(Base.active_project())
    julia = `$(Base.julia_cmd()) --startup-file=no --project=$project`
    mktempdir() do directory
        function capture(command)
            stdout_file = joinpath(directory, "stdout.txt")
            stderr_file = joinpath(directory, "stderr.txt")
            child = run(pipeline(ignorestatus(command); stdout=stdout_file, stderr=stderr_file))
            return child.exitcode, read(stdout_file, String), read(stderr_file, String)
        end
        for (argv, message) in ((["--cap=0"], "must be finite"),
                (["--caps=60"], "Unknown option"),
                (["--cap=90"], "caps must differ"),
                (["--output=" * directory], "Choose a new output directory"))
            code, out, err = capture(`$julia $script $argv`)
            @test code == 1
            @test occursin(message, err)
            # Package startup may emit compilation progress on a fresh cache.
            # The CLI diagnostic itself must remain one error line.
            @test count(line -> startswith(line, "ERROR: "), split(chomp(err), '\n')) == 1
            @test !occursin("Stacktrace:", err)
            @test isempty(out)
        end
        code, out, err = capture(`$julia $script --help`)
        @test code == 0
        @test occursin("Usage:", out)
        @test isempty(err)

        # An unavailable version reaches the real preset resolver after parsing.
        # Confine the changed version constant to this subprocess; no downloads
        # or simulation should start, and its ArgumentError must retain a trace.
        output = joinpath(directory, "unused")
        probe = """
            include($(repr(script)))
            @eval OdysseySurrogateExample const VERSION = "0.0.0-cli-regression"
            exit(OdysseySurrogateExample.cli_main(["--offline", $(repr("--output=" * output))]))
            """
        code, _, err = capture(`$julia -e $probe`)
        @test code != 0
        @test occursin("ArgumentError", err)
        @test occursin("0.0.0-cli-regression", err)
        @test occursin("Stacktrace:", err)
        @test occursin("surrogate_presets.jl", err)
        @test !ispath(output)
    end
end
