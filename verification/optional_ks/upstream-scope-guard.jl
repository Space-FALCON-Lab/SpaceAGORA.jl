include("baseline-fixtures.jl")
using KSPropagator, Test
@testset "KS gravity preflight uses the effective spacecraft scope" begin
    cfg=baseline_config(:two_body,:ks;duration=2.)
    sc=only(cfg.dynamics_model.spacecraft)
    for mask in (trues(1),falses(1))
        positions=Union{Nothing,Vector{Int}}[findall(mask)]
        masks=Union{Nothing,BitVector}[mask]
        table=S.SpacecraftModels.EffectorScopes(positions,masks,Any[SpaceAGORA.AllSpacecraft()],1,nothing)
        scoped_cfg=S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc],cfg.dynamics_model.dynamic_effectors,table))
        resolved=SE._resolve_effector_scopes(scoped_cfg)
        err=try SE._validate_ks_config!(resolved);nothing catch e;e end
        println("SCOPED_CENTRAL_GUARD active=",only(mask)," error=",err);flush(stdout)
        if only(mask)
            @test isnothing(err)
        else
            @test err isa ArgumentError
            @test err isa ArgumentError && occursin("exactly one configured",sprint(showerror,err))
        end
    end
    public_cfg=S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc],(SpaceAGORA.scoped(S.InverseSquaredGravityModel(),SpaceAGORA.BoundSpacecraft((1,))),)))
    @test isnothing(SE._validate_ks_config!(SE._resolve_effector_scopes(public_cfg)))
    for selected in (Int[],[2])
        invalid=S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc],(SpaceAGORA.scoped(S.InverseSquaredGravityModel(),SpaceAGORA.BoundSpacecraft(Tuple(selected))),)))
        @test_throws ArgumentError SE._resolve_effector_scopes(invalid)
    end
end
