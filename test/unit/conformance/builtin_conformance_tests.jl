module BuiltinConformanceTests

using Test, SpaceAGORA
using SpaceAGORA.Conformance: check_force_torque_model, package_types, test_report

include(joinpath(@__DIR__, "..", "..", "fixtures", "conformance_instances.jl"))
using .ConformanceInstances: force_torque_instances, force_torque_exemptions, sensor_instances, sensor_exemptions,
    gnc_instances, gnc_exemptions, algorithm_instances, algorithm_exemptions, family_instances
using SpaceAGORA.Conformance: check_sensor, check_guidance_model, check_control_model

const SM = SpaceAGORA.SimulationModel
const C = SpaceAGORA.Conformance

# The status of a report's checkpoint round trip.
round_trip_status(report) = only(check.status for check in report.checks if check.id === :checkpoint_round_trip)

@testset "Built-in force/torque models pass the conformance checks" begin
    types = package_types(SpaceAGORA.AbstractForceTorqueModel)
    instances = force_torque_instances()
    exemptions = force_torque_exemptions()
    @testset "every concrete built-in type has an instance or a listed exemption" begin
        @test length(types) >= 15
        for T in types
            covered = count(entry -> typeof(entry.model) <: T, instances)
            exempt = haskey(exemptions, T)
            @testset let type = T, instance_count = covered, exempt = exempt
                @test (instance_count == 1) != exempt
            end
        end
        # Every instance is of a walked type, and every exemption names one, with a reason.
        @test all(entry -> any(T -> typeof(entry.model) <: T, types), instances)
        for (T, reason) in exemptions
            @test T in types
            @test !isempty(strip(reason))
        end
    end
    for entry in instances
        entry.runs || @test !isempty(strip(entry.reason))
        report = entry.runs ? check_force_torque_model(entry.model; fixture=entry.fixture) :
            SpaceAGORA.Conformance._static_force_report(entry.model)
        test_report(report)
    end
end

@testset "Built-in sensors pass the conformance checks" begin
    types = package_types(SpaceAGORA.AbstractSensor)
    instances = sensor_instances()
    exemptions = sensor_exemptions()
    @testset "every concrete built-in sensor kind has an instance or a listed exemption" begin
        @test length(types) >= 3
        for T in types
            covered = count(sensor -> typeof(sensor) <: T, instances)
            @testset let type = T, instance_count = covered, exempt = haskey(exemptions, T)
                @test (instance_count == 1) != exempt
            end
        end
        @test all(sensor -> any(T -> typeof(sensor) <: T, types), instances)
        for (T, reason) in exemptions
            @test T in types
            @test !isempty(strip(reason))
        end
    end
    for sensor in instances
        report = test_report(check_sensor(sensor))
        # Advisory in 0.4 (so test_report records a failure as broken), required here.
        @test round_trip_status(report) === :pass
    end
end

@testset "Built-in guidance and control models pass the conformance checks" begin
    abstract_types = (guidance_model=SpaceAGORA.AbstractGuidanceModel, navigation_model=SpaceAGORA.AbstractNavigationModel,
        control_model=SpaceAGORA.AbstractControlEffectorModel)
    instances = gnc_instances()
    exemptions = gnc_exemptions()
    @testset "every concrete built-in GNC type has an instance or a listed exemption" begin
        types = Dict(kind => C.package_types(T) for (kind, T) in pairs(abstract_types))
        @test length(types[:guidance_model]) >= 5 && length(types[:control_model]) >= 6
        for (kind, kind_types) in types, T in kind_types
            covered = count(entry -> entry.kind === kind && typeof(entry.model) <: T, instances)
            @testset let kind = kind, type = T, instance_count = covered, exempt = haskey(exemptions, T)
                @test (instance_count == 1) != exempt
            end
        end
        # The thrusters run in the control model without subtyping the control effector type.
        @test count(entry -> entry.model isa SM.ThrusterModels.BaseThrusterModel, instances) == 1
        for (T, reason) in exemptions
            @test any(kind_types -> T in kind_types, values(types))
            @test !isempty(strip(reason))
        end
        for entry in instances
            entry.runs || @test !isempty(strip(entry.reason))
        end
    end
    @testset "the exempt algorithm bindings have their guidance tick" begin
        fixture = C.conformance_fixture()
        for binding in (SpaceAGORA.RPOPlannerLifecycle.PlannerGuidance, SpaceAGORA.AlgorithmGuidance)
            @test hasmethod(SpaceAGORA.update_guidance!,
                Tuple{binding, fixture.call.state, fixture.call.params, Float64, Int})
        end
    end
    @testset "the energy depletion panel actuator has outputs but no tick of its own" begin
        aoa = SpaceAGORA.SolarPanelAngleOfAttackControlModel
        fixture = C.conformance_fixture(control=(aoa(),))
        @test !hasmethod(SpaceAGORA.update_control!, Tuple{aoa, fixture.call.state, fixture.call.params, Float64, Int})
        @test hasmethod(SpaceAGORA.control_wrench, Tuple{aoa, SpaceAGORA.StateSample, fixture.call.params, Float64})
    end
    for entry in instances
        if entry.runs
            check = entry.kind === :guidance_model ? check_guidance_model : check_control_model
            report = test_report(entry.fixture === nothing ? check(entry.model) : check(entry.model; fixture=entry.fixture))
            # Advisory in 0.4 (so test_report records a failure as broken), required here.
            @test round_trip_status(report) === :pass
        else
            test_report(C._static_gnc_report(entry.model, entry.kind))
        end
    end
end

@testset "Built-in planning algorithms and families pass the conformance checks" begin
    types = C.package_types(SpaceAGORA.AbstractGNCAlgorithm)
    instances = algorithm_instances()
    exemptions = algorithm_exemptions()
    families = family_instances()
    @testset "every concrete built-in algorithm and family has an instance or a listed exemption" begin
        @test SpaceAGORA.DirectRPOPlanner in types && SpaceAGORA.HYPRRPOPlanner in types
        for T in types
            covered = count(entry -> entry.algorithm isa T, instances)
            @testset let type = T, instance_count = covered, exempt = haskey(exemptions, T)
                @test (instance_count == 1) != exempt
            end
        end
        for (T, reason) in exemptions
            @test T in types
            @test !isempty(strip(reason))
        end
        family_types = C.package_types(SpaceAGORA.AbstractAlgorithmFamily)
        @test !isempty(family_types)
        for F in family_types
            @testset let family = F
                @test any(entry -> entry.family isa F, families)
            end
        end
        # Every rejection policy is exercised: :fail and :install_fallback by the RPO
        # family, :keep_previous by the burn-plan family, :fail by the joint-reference family.
        @test Set(SpaceAGORA.on_rejection(entry.family) for entry in families) ==
            Set((:fail, :install_fallback, :keep_previous))
    end
    for entry in instances
        test_report(C.check_algorithm(entry.algorithm, entry.requests; trusted=entry.trusted,
            simulation=entry.simulation))
    end
    for entry in families
        test_report(C.check_algorithm_family(entry.family, entry.algorithm, entry.requests; link=entry.link,
            trusted=entry.trusted))
    end
end

end # module BuiltinConformanceTests
