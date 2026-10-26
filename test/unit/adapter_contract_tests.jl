module AdapterContractTests
using Test, SpaceAGORA
using SpaceAGORA.Conformance: test_service_facade

binding_path(name::Symbol) = (name,)
binding_path(path::Tuple{Vararg{Symbol}}) = path
resolve_binding(facade, path) = foldl(getproperty, binding_path(path); init=facade)

# Each facade's reviewed contract version; a bump changes its entry here.
const CONTRACT_VERSIONS = (HYPRServices=v"3.0.0", ExternalPropagationServices=v"2.0.0", EDGServices=v"3.0.0")

@testset "Optional adapter facade conformance" begin
    # Version, exact require_version, unique and disjoint inventories that resolve:
    # SpaceAGORA.Conformance.check_service_facade, which companions call on their facade.
    for facade in (SpaceAGORA.HYPRServices, SpaceAGORA.ExternalPropagationServices, SpaceAGORA.EDGServices)
        test_service_facade(facade; version=CONTRACT_VERSIONS[nameof(facade)])
    end
    for path in SpaceAGORA.HYPRServices.CLOTH_SERVICES
        @test resolve_binding(SpaceAGORA.HYPRServices.Cloth, path) isa Function
    end
end
end
