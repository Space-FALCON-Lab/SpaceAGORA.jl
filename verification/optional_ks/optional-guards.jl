ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
include("baseline-fixtures.jl")
using Test,KSPropagator
struct UnknownTrialForce <: S.AbstractTypes.AbstractForceTorqueModel end
function refused(config,fragment)
    err=try SE._validate_ks_config!(config);nothing catch e;e end
    @test err isa ArgumentError
    @test occursin(fragment,sprint(showerror,err))
end
@testset "bounded KS configuration guards" begin
    cfg=baseline_config(:two_body,:ks;duration=2.)
    @test SE._validate_ks_config!(cfg)===nothing
    @test !SpaceAGORA.trial_rhs_safe(UnknownTrialForce())
    @test SpaceAGORA.trial_rhs_safe(S.InverseSquaredGravityModel())
    @test SpaceAGORA.trial_rhs_safe(S.AerodynamicCoefficientConstant())
    mission=cfg.mission_configuration
    refused(S.SimConfig._with_configuration(cfg;mission_configuration=S.MissionConfiguration(
        mission_type=S.MissionTime,mission_time=2.,orientation_sim=true,keplerian=true)),"attitude")
    sc=only(cfg.dynamics_model.spacecraft)
    refused(S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc,deepcopy(sc)],cfg.dynamics_model.dynamic_effectors)),"exactly one")
    for forces in ((),(S.InverseSquaredGravityModel(),S.InverseSquaredJ2GravityModel()))
        refused(S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc],forces)),"exactly one configured")
    end
    refused(S.SimConfig._with_configuration(cfg;dynamics_model=S.DynamicsModel([sc],
        (S.InverseSquaredGravityModel(),UnknownTrialForce()))),"repeatable trial RHS")
    # Actual trial evaluation must reject a configured frozen sample contract.
    aero=baseline_config(:drag_only,:ks;duration=2.)
    ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="true"
    err=try run_simulation(aero;return_solution=true);nothing catch e;e end
    @test err isa ArgumentError
    @test occursin("DENSITY_FREEZE_PER_STEP=false",sprint(showerror,err))
    ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
    # Normal DP8 and KS algorithms must not reuse the other's integrator cache.
    normal=run_simulation(baseline_config(:two_body,:dp8;duration=2.);return_solution=true)
    ks=run_simulation(cfg;return_solution=true)
    again=run_simulation(baseline_config(:two_body,:dp8;duration=2.);return_solution=true)
    @test occursin("DP8",string(typeof(normal.alg)))
    @test occursin("CoupledKS",string(typeof(ks.alg)))
    @test typeof(again.alg)==typeof(normal.alg)
    @test all(s->s.t[end]==2.,(normal,ks,again))
end
