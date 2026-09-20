using Test, SpaceAGORA

@testset "Preset CLI validates user arguments before retrieval" begin
    listing = sprint(io -> @test SpaceAGORA.run_cli(["assets", "list"]; io) == 0)
    @test occursin("odyssey_p20_frozen_v1@1.0.0", listing)
    @test occursin("mars, available", listing)
    @test_throws ArgumentError SpaceAGORA.run_cli(["assets", "list", "extra"])

    for command in ("fetch", "check")
        for args in (
            ["--unknown=value"],
            ["--preset=x", "--preset=y"],
            ["--preset"],
            ["--preset", "--version=1.0.0"],
            ["--preset="],
            ["unexpected"],
            ["--preset=x"],
            ["--version=1.0.0"],
        )
            @test_throws ArgumentError SpaceAGORA.run_cli(["assets", command, args...])
        end
        # Both accepted option spellings must reach resolution, not fail in the
        # parser. An unknown preset guarantees that no network access can occur.
        for args in (
            ["--preset=missing_preset", "--version=1.0.0", "--offline"],
            ["--preset", "missing_preset", "--version", "1.0.0", "--file", "absent"],
        )
            error = try
                SpaceAGORA.run_cli(["assets", command, args...])
                nothing
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("Unknown surrogate preset 'missing_preset@1.0.0'", sprint(showerror, error))
        end
    end
end
