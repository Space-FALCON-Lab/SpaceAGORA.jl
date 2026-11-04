using Test

@testset "Development launcher executes the requested script" begin
    launcher = joinpath(@__DIR__, "run.jl")
    mktempdir() do tmp
        output = joinpath(tmp, "result.txt")
        guarded = joinpath(tmp, "guarded script.jl")
        write(guarded, """
        if abspath(PROGRAM_FILE) == @__FILE__
            length(ARGS) == 2 || error("arguments not forwarded")
            write(ARGS[1], ARGS[2])
            println("GUARDED_MAIN_EXECUTED")
        end
        """)
        cmd = `$(Base.julia_cmd()) --startup-file=no $launcher $guarded $output "argument with spaces"`
        @test occursin("GUARDED_MAIN_EXECUTED", read(cmd, String))
        @test isfile(output)
        if isfile(output)
            @test read(output, String) == "argument with spaces"
        end
        unguarded = joinpath(tmp, "unguarded.jl")
        write(unguarded, "println(\"UNGUARDED_EXECUTED\")")
        @test strip(read(`$(Base.julia_cmd()) --startup-file=no $launcher $unguarded`, String)) == "UNGUARDED_EXECUTED"
        failure = joinpath(tmp, "failure.jl")
        write(failure, "if abspath(PROGRAM_FILE) == @__FILE__; exit(7); end")
        failed = run(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $launcher $failure`))
        @test failed.exitcode == 7
    end
end
