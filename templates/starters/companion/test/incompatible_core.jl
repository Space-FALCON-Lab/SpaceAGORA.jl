# Run in a fresh process without loading the companion first: its source, with
# the contract version changed to one the core does not provide, must refuse
# before it defines any method of the facade's functions.
using Test, SpaceAGORA

@testset "an incompatible contract version is refused before any method is defined" begin
    @test !isdefined(Main, :StarterName)
    path = joinpath(@__DIR__, "..", "src", "StarterName.jl")
    source = read(path, String)
    pattern = r"const CORE_CONTRACT_VERSION = v\"[^\"]*\""
    @test length(collect(eachmatch(pattern, source))) == 1
    incompatible = replace(source, pattern => "const CORE_CONTRACT_VERSION = v\"999.0.0\"")
    failure = try
        Base.include_string(Main, incompatible, path)
        nothing
    catch error
        while error isa LoadError
            error = error.error
        end
        error
    end
    services = SpaceAGORA.ExternalPropagationServices
    @test failure isa services.CompatibilityError
    @test isdefined(Main, :StarterName)
    @test !isdefined(Main.StarterName, :CircularOrbitPropagator)
    for name in services.IMPLEMENTED_FUNCTIONS
        @test !any(method -> method.module === Main.StarterName, methods(getproperty(services, name)))
    end
end
