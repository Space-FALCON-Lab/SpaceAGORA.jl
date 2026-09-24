# The paper harness's precompile workload must cover every (case, mode) the
# paper phases (P1-P7) run, at the size each case names.
#
# The workload's point list (benchmarks/studies/paper_parallelization_benchmarks/
# workload/SpaceAGORAPaperWorkload/src/points.jl) is derived from
# PAPER_BENCHMARK_PHASES and the case catalog rather than written out, so a case
# added to a phase is picked up without an edit. This gate pins that property:
# every case a P phase names must either be a workload point under every mode the
# phase runs it with, or be one of the documented native-GRAM exclusions, and
# every workload point must name a real catalog case and a real mode. It reads
# only the harness definitions; it does not build or load the image.

using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PPBW_STUDIES = joinpath(REPO_ROOT, "benchmarks", "studies")

# A sandbox module, so the harness's globals stay out of whatever includes this gate.
const PPBW_SANDBOX = Module(:PPBWorkloadCoverageGate)
for file in (
    joinpath(PPBW_STUDIES, "parallelization_performance", "cli.jl"),
    joinpath(PPBW_STUDIES, "parallelization_performance", "modes.jl"),
    joinpath(PPBW_STUDIES, "parallelization_performance", "cases.jl"),
    joinpath(PPBW_STUDIES, "paper_parallelization_benchmarks", "cli.jl"),
    joinpath(PPBW_STUDIES, "paper_parallelization_benchmarks", "workload",
             "SpaceAGORAPaperWorkload", "src", "points.jl"),
)
    Base.include(PPBW_SANDBOX, file)
end

@testset "paper precompile workload covers the paper phases" begin
    S = PPBW_SANDBOX
    phases = Base.invokelatest(S.ppb_workload_phases)
    ids = [p.id for p in phases]
    # The paper figure phases exist and are what the workload reads.
    @test issubset(["P1", "P2", "P3", "P4", "P5", "P6", "P6p", "P7"], ids)
    @test all(id -> startswith(id, "P"), ids)

    points = Base.invokelatest(S.ppb_workload_points)
    catalog = Base.invokelatest(S.ppc_case_catalog)
    modes = Base.invokelatest(S.ppc_mode_specs)
    have = Set((p.case, p.mode, p.kind) for p in points)
    excluded(case) = Base.invokelatest(S.ppb_workload_excluded, case)

    missing_points = String[]
    for phase in phases
        for case in phase.cases, mode in phase.modes
            excluded(case) && continue
            (case, mode, :perf) in have || push!(missing_points, "$(phase.id): $(case) / $(mode)")
        end
        for case in phase.parity_cases, mode in phase.modes
            (mode == "serial" || excluded(case)) && continue
            (case, mode, :parity) in have || push!(missing_points, "$(phase.id) parity: $(case) / $(mode)")
        end
    end
    @test isempty(missing_points)
    isempty(missing_points) || println("workload points missing:\n  ", join(missing_points, "\n  "))

    # Every point is buildable: a catalog case, a known mode, a sample count that
    # reaches the campaign path exactly when the case is a campaign.
    for p in points
        @test haskey(catalog, p.case)
        @test haskey(modes, p.mode)
        @test p.samples == (catalog[p.case].montecarlo ? 2 : 1)
    end

    # The only exclusions are the native-GRAM cases, and they are all P6 traces.
    skipped = Base.invokelatest(S.ppb_workload_skipped)
    @test all(((ph, case),) -> ph in ("P6", "P6p") && occursin("gram", case), skipped)

    # Every constellation size the P phases name appears among the points.
    size_of(case) = (m = match(r"_([0-9]+)sat_", case); m === nothing ? 1 : parse(Int, m.captures[1]))
    phase_sizes = Set(size_of(c) for ph in phases for c in ph.cases if !excluded(c))
    @test phase_sizes == Set(size_of(p.case) for p in points)
end

# P7 is P1's comparison on one spacecraft over one orbit: the same modes and
# thread axis as P1, one satellite per case, one mission length shared by every
# row and equal to the two-body period of P1's spacecraft's orbit.
@testset "P7 short single-satellite phase" begin
    S = PPBW_SANDBOX
    by_id = Dict(p.id => p for p in S.PAPER_BENCHMARK_PHASES)
    p1, p7 = by_id["P1"], by_id["P7"]
    @test p7.modes == p1.modes
    @test p7.thread_mode == p1.thread_mode == :max_only
    @test p7.mc_samples == [1]
    @test p7.repeats == 11
    @test p7.warmup == p1.warmup
    catalog = Base.invokelatest(S.ppc_case_catalog)
    @test length(p7.cases) == 3
    @test all(c -> haskey(catalog, c) && !catalog[c].montecarlo, p7.cases)
    @test all(c -> occursin(r"_1sat_", c), p7.cases)
    @test all(c -> endswith(c, "_$(S.PPC_P7_MISSION_S)s"), p7.cases)
    ra, rp = Base.invokelatest(S.ppc_constellation_member_alts_m, 1)
    earth = S.Earth()
    a = earth.Rp_e + (ra + rp) / 2
    @test S.PPC_P7_MISSION_S == round(Int, 2π * sqrt(a^3 / earth.μ))
end

# The workload stays opt-in: with the image loaded the timed repeats move by a
# few percent, by an amount that depends on the build of the image (code
# placement; workload/README.md, "Why the timed repeats move"). Pin that the
# harness and the remote runner both default to off, as the README says, so the
# default cannot be flipped without the documentation and this gate changing too.
@testset "paper precompile workload is off by default" begin
    Base.include(PPBW_SANDBOX, joinpath(PPBW_STUDIES, "parallelization_performance", "execution.jl"))
    withenv("SPACEAGORA_PPB_WORKLOAD" => nothing) do
        @test Base.invokelatest(PPBW_SANDBOX.ppc_workload_env) === nothing
    end
    withenv("SPACEAGORA_PPB_WORKLOAD" => "0") do
        @test Base.invokelatest(PPBW_SANDBOX.ppc_workload_env) === nothing
    end
    # Asking for it with no image built is an error, never a silent fallback.
    withenv("SPACEAGORA_PPB_WORKLOAD" => "1",
            "SPACEAGORA_PPB_WORKLOAD_ENV" => joinpath(mktempdir(), "no_workload_env")) do
        @test_throws ErrorException Base.invokelatest(PPBW_SANDBOX.ppc_workload_env)
    end
    # The unset case above returns nothing under auto too when no current image
    # happens to exist, so pin the default value itself as well.
    execution = read(joinpath(PPBW_STUDIES, "parallelization_performance", "execution.jl"), String)
    @test occursin("get(ENV, \"SPACEAGORA_PPB_WORKLOAD\", \"0\")", execution)
    remote = read(joinpath(REPO_ROOT, "scripts", "remote", "spaceagora-remote"), String)
    @test occursin(r"\bworkload=\"off\"", remote)
end

println("ppb_workload_coverage_gate_ok")
