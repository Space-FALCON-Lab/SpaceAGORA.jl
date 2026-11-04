module GramPerturbationStudyScriptsTests
# The perturbed-density study scripts' outcome rules, without a run: the
# dispersed campaign's member classification (dispersed_campaign/member_status.jl)
# and the perturbed-density runs' attempt record (perturbed_density_modes/
# attempt_record.jl). Neither file loads SpaceAGORA.
using Test, SHA, TOML

const STUDY = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "telemetry_validation"))
include(joinpath(STUDY, "dispersed_campaign", "member_status.jl"))
include(joinpath(STUDY, "perturbed_density_modes", "attempt_record.jl"))
using .DispersedMemberStatus, .PerturbedModeAttempt

# A member row as run_campaign.jl builds it: the dispatcher's fields merged with
# odyssey_sample.jl's summary.
member(; kw...) = merge((dispatch_success=true, error="", retcode="Terminated", n_apo=41, n_peri=41), (; kw...))

@testset "dispersed campaign member status" begin
    H = 40
    @test member_status(member(); horizon_passes=H) == "complete"
    @test member_status(member(n_apo=42); horizon_passes=H) == "complete"
    @test member_status(member(retcode="Success"); horizon_passes=H) == "complete"
    # A solve that threw: odyssey_sample.jl catches it and returns retcode ERROR
    # with the error text.
    @test member_status(member(retcode="ERROR", error="DomainError(-1.0)", n_apo=0, n_peri=0);
                        horizon_passes=H) == "solver_failure"
    # An unsuccessful return without an exception.
    @test member_status(member(retcode="Unstable"); horizon_passes=H) == "solver_failure"
    @test member_status(member(retcode="MaxIters", n_apo=3, n_peri=2); horizon_passes=H) == "solver_failure"
    # Terminated short of the horizon (an impact), and the shape the old writer
    # produced for --orbits=41: not complete, although the retcode is Terminated.
    @test member_status(member(n_apo=2, n_peri=2); horizon_passes=H) == "early_termination"
    @test member_status(member(n_apo=40, n_peri=40); horizon_passes=H) == "early_termination"
    @test member_status(member(n_peri=39); horizon_passes=H) == "early_termination"
    # The dispatcher holds no value, so the row has no member fields.
    @test member_status((dispatch_success=false, error="ProcessExitedException(3)"); horizon_passes=H) ==
          "dispatch_failure"
    @test member_status((dispatch_success=false,); horizon_passes=H) == "dispatch_failure"
    @test_throws ArgumentError member_status(member(); horizon_passes=0)

    rows = [member(), member(n_apo=2, n_peri=2), member(retcode="ERROR", error="boom"),
            member(retcode="Unstable"), (dispatch_success=false, error="lost")]
    counts = status_counts(rows; horizon_passes=H)
    @test counts == Dict("complete" => 1, "early_termination" => 1, "solver_failure" => 2, "dispatch_failure" => 1)
    @test status_counts(rows[1:1]; horizon_passes=H) ==
          Dict("complete" => 1, "early_termination" => 0, "solver_failure" => 0, "dispatch_failure" => 0)
    @test Set(MEMBER_STATUSES) == Set(keys(counts))
end

@testset "campaign execution totals survive every dispatch failing" begin
    failed = [(dispatch_success=false, dispatch_elapsed_s=1.0, error="lost worker"),
              (dispatch_success=false, dispatch_elapsed_s=2.0, error="dispatch exception")]
    # These rows have neither solve_s nor pid, the shape that used to prevent
    # campaign.toml from being written when no member returned a value.
    totals = member_execution_totals(failed)
    @test totals == (sum_member_solve_s=0.0, distinct_member_pids=0)
    @test member_execution_totals(NamedTuple[]) == totals
    counts = status_counts(failed; horizon_passes=40)
    @test counts["dispatch_failure"] == 2
    summary = merge(Dict("n_dispatch_failed" => counts["dispatch_failure"]),
                    Dict(string(k) => v for (k,v) in pairs(totals)))
    io = IOBuffer()
    TOML.print(io, summary)
    @test TOML.parse(String(take!(io))) == summary
    returned = [(dispatch_success=true, solve_s=2.0, pid=101, retcode="Success"),
                (dispatch_success=true, solve_s=3.0, pid=101, retcode="ERROR"),
                (dispatch_success=true, solve_s=4.0, pid=202, retcode="Terminated")]
    @test member_execution_totals(vcat(failed, returned)) ==
          (sum_member_solve_s=9.0, distinct_member_pids=2)
end

@testset "perturbed-density run attempts" begin
    mktempdir() do dir
        tag = joinpath(dir, "A_s11_rep")
        # A completed attempt: its outputs, then the summary with their identities.
        id1 = begin_attempt!(tag)
        write(joinpath(tag, "simulation_results.csv"), "time,x\n0,1\n")
        write(joinpath(tag, "extrema.csv"), "event,index\napo,1\n")
        write(joinpath(tag, "perturbation_log.csv"), "kind\n1\n")
        status, reason = attempt_status(; error="", retcode="Terminated", have_trajectory=true,
                                        n_apo=119, n_peri=119, orbits=120,
                                        completed_orbits=120, termination_cause="orbit_count")
        @test (status, reason) == ("complete", "")
        finish_attempt!(tag, Dict{String, Any}("attempt_id" => id1, "status" => status))
        s = TOML.parsefile(joinpath(tag, "run_summary.toml"))
        @test s["attempt_id"] == id1
        @test s["output_sha256"]["simulation_results.csv"] == bytes2hex(sha256("time,x\n0,1\n"))
        @test sort(collect(keys(s["output_sha256"]))) == ["extrema.csv", "perturbation_log.csv", "simulation_results.csv"]
        @test !isfile(joinpath(tag, ".run_summary.toml.partial"))

        # A rerun that fails before writing a trajectory: the earlier attempt's
        # outputs are removed before it runs, and its summary names no output.
        id2 = begin_attempt!(tag)
        @test id2 != id1
        @test isempty(readdir(tag))
        status, reason = attempt_status(; error="solver failure", retcode="ERROR", have_trajectory=false,
                                        n_apo=0, n_peri=0, orbits=120)
        @test status == "failed"
        finish_attempt!(tag, Dict{String, Any}("attempt_id" => id2, "status" => status, "status_reason" => reason))
        s = TOML.parsefile(joinpath(tag, "run_summary.toml"))
        @test s["status"] == "failed"
        @test isempty(s["output_sha256"])
        @test readdir(tag) == ["run_summary.toml"]
    end

    ok = (error="", retcode="Terminated", have_trajectory=true, n_apo=119, n_peri=119, orbits=120,
          completed_orbits=120, termination_cause="orbit_count")
    status_of(; kw...) = first(attempt_status(; merge(ok, (; kw...))...))
    @test status_of() == "complete"
    @test status_of(n_apo=120) == "complete"
    @test status_of(retcode="Success", termination_cause="end_of_time_span") == "complete"
    @test status_of(error="boom") == "failed"
    @test status_of(retcode="Unstable") == "failed"
    @test status_of(have_trajectory=false) == "failed"
    @test status_of(n_apo=118) == "incomplete"
    @test status_of(n_peri=118) == "incomplete"
    # Near-end early termination has enough saved extrema to pass the old gate,
    # but has not reached the requested runtime orbit count.
    @test attempt_status(; error="", retcode="Terminated", have_trajectory=true,
                         n_apo=4, n_peri=4, orbits=5, completed_orbits=4,
                         termination_cause="terminated_before_orbit_count") ==
          ("incomplete", "4 completed orbit events for 5 requested")
    @test status_of(retcode="Success", termination_cause="end_of_time_span", completed_orbits=119) == "incomplete"
    # Missing, malformed or inconsistent metadata cannot produce a complete record.
    @test first(attempt_status(; error="", retcode="Terminated", have_trajectory=true,
                              n_apo=119, n_peri=119, orbits=120)) == "failed"
    for count in (nothing, missing, "120", 120.0, true, -1)
        @test status_of(completed_orbits=count) == "failed"
    end
    for cause in (nothing, missing, 1, "terminated_unknown", "end_of_time_span", "terminated_before_orbit_count")
        @test status_of(termination_cause=cause) == "failed"
    end
    @test status_of(completed_orbits=119) == "failed"  # inconsistent orbit_count cause
    @test status_of(retcode="Success") == "failed"  # inconsistent orbit_count cause
    @test status_of(orbits=0) == "failed"
    @test status_of(orbits=true) == "failed"
end
end
