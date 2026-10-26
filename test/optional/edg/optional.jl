using Test, Pkg
length(ARGS) == 1 || error("Pass core-first or edg-first")
if only(ARGS) == "core-first"
    using SpaceAGORA
    @test !SpaceAGORA.EDGServices.available()
    const configuration_before = SpaceAGORA.AerobrakingEnergyDepletionConfig
    @test_throws SpaceAGORA.EDGServices.EDGUnavailableError SpaceAGORA.EDGServices.require_provider()
    using SpaceAGORAEDG
    @test configuration_before === SpaceAGORA.AerobrakingEnergyDepletionConfig
elseif only(ARGS) == "edg-first"
    using SpaceAGORAEDG
    using SpaceAGORA
else
    error("Pass core-first or edg-first")
end

@testset "Public free-molecular aerodynamic model" begin
    @test :AerodynamicCoefficientfM in names(SpaceAGORA)
    @test SpaceAGORA.AerodynamicCoefficientfM === SpaceAGORA.SimulationModel.AerodynamicCoefficientfM
    model = SpaceAGORA.AerodynamicCoefficientfM()
    @test model.per_link_atmosphere === false
    @test model.fixed_attitude_incidence === :max_drag
    @test SpaceAGORA.AerodynamicCoefficientfM(fixed_attitude_incidence=:attitude).fixed_attitude_incidence === :attitude
    documentation = Base.Docs.doc(SpaceAGORA.AerodynamicCoefficientfM)
    @test documentation !== nothing
    @test occursin("fixed_attitude_incidence", sprint(show, MIME"text/plain"(), documentation))
end
@testset "EDG integration in $(only(ARGS)) order" begin
    ES = SpaceAGORA.EDGServices
    @test Base.PkgId(SpaceAGORAEDG).uuid == Base.UUID("a782e8ee-96d7-4dd6-8231-bb64649e6ad4")
    @test Pkg.dependencies()[Base.PkgId(SpaceAGORAEDG).uuid].version == v"0.5.1"
    @test ES.CONTRACT_VERSION == SpaceAGORAEDG.REQUIRED_CONTRACT_VERSION == v"3.0.0"
    @test ES.available()
    @test ES.require_provider() === Base.get_extension(SpaceAGORA, :SpaceAGORAEDGExt)
    @test_throws ES.CompatibilityError SpaceAGORAEDG.require_version(v"1.0.0")
    @test_throws ES.CompatibilityError SpaceAGORAEDG.require_version(v"2.0.0")
    @test SpaceAGORAEDG.require_version(v"3.0.0") === nothing
    @test SpaceAGORAEDG.AerobrakingEnergyDepletionConfig === SpaceAGORA.AerobrakingEnergyDepletionConfig
    @test SpaceAGORAEDG.AerobrakingEnergyDepletionState === SpaceAGORA.AerobrakingEnergyDepletionState
end
