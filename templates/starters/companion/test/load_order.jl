# Run in a fresh process, with core-first or companion-first: either load order
# gives the same facade, contract version and methods.
using Test
length(ARGS) == 1 || error("Pass core-first or companion-first.")
if only(ARGS) == "core-first"
    using SpaceAGORA
    @test !any(package -> package.name == "StarterName", keys(Base.loaded_modules))
    using StarterName
elseif only(ARGS) == "companion-first"
    using StarterName
    using SpaceAGORA
else
    error("Pass core-first or companion-first.")
end

@testset "load order $(only(ARGS))" begin
    services = SpaceAGORA.ExternalPropagationServices
    @test StarterName.EP === services
    @test services.CONTRACT_VERSION == StarterName.CORE_CONTRACT_VERSION
    @test StarterName.CircularOrbitPropagator <: services.AbstractExternalPropagator
    counts = map(name -> length(methods(getproperty(services, name))), services.IMPLEMENTED_FUNCTIONS)
    @test StarterName.__init__() === nothing
    @test map(name -> length(methods(getproperty(services, name))), services.IMPLEMENTED_FUNCTIONS) == counts
    for name in services.IMPLEMENTED_FUNCTIONS
        @test any(method -> method.module === StarterName, methods(getproperty(services, name)))
    end
end
