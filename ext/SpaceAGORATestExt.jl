module SpaceAGORATestExt

using SpaceAGORA
using SpaceAGORA.Conformance: ConformanceCheck, ConformanceReport, check_force_torque_model, check_service_facade,
    check_sensor, check_guidance_model, check_navigation_model, check_control_model, check_algorithm,
    check_algorithm_family
import SpaceAGORA.Conformance: test_report, test_force_torque_model, test_service_facade, test_sensor,
    test_guidance_model, test_navigation_model, test_control_model, test_algorithm, test_algorithm_family
using Test

function test_report(report::ConformanceReport)
    @testset "$(report.target) ($(report.kind))" begin
        for check in report.checks
            _record(check)
        end
    end
    SpaceAGORA.Conformance._write_report_file(report)
    return report
end

# One testset per check. A failure prints the check's detail as test context.
function _record(check::ConformanceCheck)
    @testset "$(check.id)" begin
        if check.status === :skip
            @test_skip check.detail
        elseif check.status === :advisory_fail
            @info "Advisory conformance check $(check.id) (added in $(check.since)) did not pass; " *
                "it becomes required in the next minor release." detail = check.detail
            @test_broken false
        else
            @testset let status = check.status, detail = check.detail
                @test status === :pass
            end
        end
    end
    return nothing
end

test_force_torque_model(model; kwargs...) = test_report(check_force_torque_model(model; kwargs...))
test_service_facade(facade::Module; kwargs...) = test_report(check_service_facade(facade; kwargs...))
test_sensor(sensor; kwargs...) = test_report(check_sensor(sensor; kwargs...))
test_guidance_model(model; kwargs...) = test_report(check_guidance_model(model; kwargs...))
test_navigation_model(model; kwargs...) = test_report(check_navigation_model(model; kwargs...))
test_control_model(model; kwargs...) = test_report(check_control_model(model; kwargs...))
test_algorithm(algorithm, requests; kwargs...) = test_report(check_algorithm(algorithm, requests; kwargs...))
test_algorithm_family(family, algorithm, requests; kwargs...) =
    test_report(check_algorithm_family(family, algorithm, requests; kwargs...))

end # module SpaceAGORATestExt
