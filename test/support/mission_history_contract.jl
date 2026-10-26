# Generated from SpaceAGORA.jl scripts/companion_support/mission_history_contract.jl.
# Source SHA-256: 7f536de3e962cbff083ee9507db1688ef6c3a8138a394ba5387259f5be259044; regenerate with scripts/sync_companion_support.jl.
module MissionHistoryContract
using Test, JSON
include("mission_history_cases.jl")

function check(validate)
    @test MISSION_HISTORY_FIXTURE_VERSION == 1
    @test length(unique(case.name for case in MISSION_HISTORY_CASES)) == length(MISSION_HISTORY_CASES)
    for case in MISSION_HISTORY_CASES
        @testset "$(case.name)" begin
            history = JSON.parse(case.json)
            original = deepcopy(history)
            if case.accepted
                @test validate(history) === history
                @test JSON.json(history) == JSON.json(original)
            else
                @test_throws ArgumentError validate(history)
            end
            @test isequal(history, original)
        end
    end
    # Non-JSON values must also fail, without coercing NaN/Inf into valid times.
    for bad in (NaN, Inf, -Inf)
        history = JSON.parse(MISSION_HISTORY_CASES[2].json)
        history["segments"][1]["end_s"] = bad
        @test_throws ArgumentError validate(history)
    end
end
end
