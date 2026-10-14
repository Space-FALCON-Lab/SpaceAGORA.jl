module AdapterContractTests
using Test, SpaceAGORA

binding_path(name::Symbol) = (name,)
binding_path(path::Tuple{Vararg{Symbol}}) = path
resolve_binding(facade, path) = foldl(getproperty, binding_path(path); init=facade)

@testset "Optional adapter facade conformance" begin
    for facade in (SpaceAGORA.HYPRServices, SpaceAGORA.ExternalPropagationServices, SpaceAGORA.EDGServices)
        @testset "$(nameof(facade))" begin
            @test facade.CONTRACT_VERSION == v"1.0.0"
            @test facade.CompatibilityError <: Exception
            @test facade.require_version(facade.CONTRACT_VERSION) === nothing
            for incompatible in (v"0.0.0", v"1.0.1", v"2.0.0")
                @test_throws facade.CompatibilityError facade.require_version(incompatible)
            end
            # Reject duplicate/overlapping declarations, missing bindings and a
            # consumed value accidentally advertised as an extension function.
            implemented = binding_path.(facade.IMPLEMENTED_FUNCTIONS)
            consumed = binding_path.(facade.CONSUMED_SERVICES)
            @test length(unique(implemented)) == length(implemented)
            @test length(unique(consumed)) == length(consumed)
            @test isempty(intersect(implemented, consumed))
            for path in implemented
                @test !isempty(path)
                @test resolve_binding(facade, path) isa Function
            end
            for path in consumed
                @test !isempty(path)
                @test resolve_binding(facade, path) !== nothing
            end
        end
    end
    for path in SpaceAGORA.HYPRServices.CLOTH_SERVICES
        @test resolve_binding(SpaceAGORA.HYPRServices.Cloth, path) isa Function
    end
end
end
