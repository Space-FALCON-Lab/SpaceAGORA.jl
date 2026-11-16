include("CampaignPaths.jl")
include(joinpath(CampaignPaths.BASE_CODE, "StudyRunner.jl"))
include(joinpath(CampaignPaths.BASE_CODE, "VerifyOutput.jl"))
include("ComparisonResults.jl")
using Arrow, CSV, DataFrames, Dates, LinearAlgebra, StaticArrays, Test, TOML
using Plots

const SC = OracleStudy.StudyConfig
const SM = OracleStudy.SM
const SP = parentmodule(SM.SchedulingPolicyModel)
const CR = ComparisonResults
const RESULTS = normpath(joinpath(@__DIR__, "..", "..", "Tests",
    "Output_3_rangeOnly_vs_positiveAlongTrack", "validation", "production_campaign"))
mkpath(RESULTS)
SC.atomic_toml(joinpath(RESULTS, "validation.toml"),
    Dict("passed" => false, "status" => "running", "started_utc" => string(now(UTC))))

@testset "Exact production grid and unchanged Output 2 contract" begin
    @test length(SC.grid_cases()) == 270
    laser, baseline = SC.grid_cases("range_along"), SC.grid_cases("no_laser")
    @test length(laser) == 162
    @test length(baseline) == 81
    @test allunique(SC.case_signature.(vcat(laser, baseline)))
    @test Set(c.helpers for c in laser) == Set((1, 100, 200, 300, 400, 500, 600, 700, 800))
    @test Set(c.altitude_km for c in laser) == Set((1000.0, 1050.0, 1100.0))
    @test Set(c.inclination_deg for c in laser) == Set((0.0, 0.5, 1.0))
    @test all(c -> c.power_W == 0 && c.policy === :no_laser, baseline)
    @test all(c -> c.power_W == 10_000 && c.orbits == 100, laser)
    @test_throws ArgumentError SC.StudyCase(1000, 0, 1, :no_laser; power_W=10_000)
    for c in vcat(laser, baseline)
        settings = SC.case_settings(c)
        @test settings["helper_altitude_km"] == 1000
        @test settings["helper_inclination_deg"] == 0
        @test settings["target_true_anomaly_deg"] == 1
        @test settings["target_index"] == c.helpers + 1
        @test settings["schema_version"] == 3
        @test settings["dt_max_s"] == 10
        @test settings["reltol"] == settings["abstol"] == 1e-12
        @test length(SC.trajectory_times(c)) == 1001
        @test SC.trajectory_times(c)[801] ≈ 0.8SC.mission_duration(c)
        @test SC.trajectory_times(c)[1001] == SC.mission_duration(c)
    end
    old_case = SC.StudyCase(1050, 0.5, 1, :gve_sma)
    old_info = TOML.parsefile(joinpath(SC.case_directory(SC.DEFAULT_OUTPUT, old_case), "run_info.toml"))
    @test SC.case_settings(old_case) == old_info["settings"]
    @test SC.case_signature(old_case) == old_info["run_id"]
end

function scheduler_fixture(policy)
    spacecraft = fill((n_terminal=1,), 4)
    model = SM.ConstellationInterLinkModel()
    for helper in (3, 1, 2)
        SM.register_allowed_link!(model, spacecraft, (helper, 1), (4, 1);
            parameters=SM.InterLinkParameters(range=100.0))
    end
    selector = SM.SchedulingPolicyModel(policy; target_idx=4)
    function select(offsets; active=trues(4))
        target = SVector(7_000_000.0, 0.0, 0.0)
        state = (sc=[(pos=target + SVector{3,Float64}(offset),
            vel=SVector(0.0, 7500.0, 0.0), mass=227.0)
            for offset in (offsets..., (0.0, 0.0, 0.0))],)
        result = SP.schedule_interlinks!(model, selector, spacecraft, state, SC.MU;
            candidate_force=(args...) -> error("Stateful scheduling must not use GVE scores"),
            is_active=active)
        @test count(model.graph.active) == length(result) <= 1
        isempty(result) ? 0 : only(result)[1][1]
    end
    (; select, model, selector)
end

@testset "Stateful native scheduling and deterministic helper order" begin
    range_only = scheduler_fixture(:naive_next_entering)
    @test range_only.select(((40, 0, 0), (50, 0, 0), (200, 0, 0))) == 1
    @test range_only.select(((80, 0, 0), (10, 0, 0), (200, 0, 0))) == 1
    @test range_only.select(((150, 0, 0), (10, 0, 0), (200, 0, 0))) == 0
    @test range_only.select(((150, 0, 0), (10, 0, 0), (200, 0, 0))) == 0
    @test range_only.select(((150, 0, 0), (10, 0, 0), (100, 0, 0))) == 3
    along = scheduler_fixture(:positive_along_track)
    @test along.select(((80, -60, 0), (0, -50, 0), (200, 0, 0))) == 2
    @test along.select(((0, -50, 0), (80, -10, 0), (200, 0, 0))) == 2
    @test along.select(((0, -50, 0), (0, 50, 0), (200, 0, 0))) == 1
    @test along.select(((0, 50, 0), (0, 50, 0), (80, 0, 0))) == 0
    @test along.select(((0, -100, 0), (0, -100, 0), (0, 0, 0))) == 1
    @test along.select(((0, -100, 0), (0, -100, 0), (0, 0, 0));
        active=Bool[false, true, true, true]) == 2
    disabled = scheduler_fixture(:no_laser)
    disabled.model.graph.active .= true
    @test disabled.select(((0, -50, 0), (0, -20, 0), (0, -10, 0))) == 0
    @test_throws ArgumentError SM.SchedulingPolicyModel(:positive_along_track; target_idx=nothing)
end

@testset "Common-frame velocity subtraction and RAAN wrapping" begin
    expected = SVector(3.0, 2.0, 4.0, sqrt(29.0))
    @test CR.velocity_difference_in_baseline_rtn(
        SVector(3.0, 7.0, 4.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 5.0, 0.0)) ≈ expected
    @test CR.velocity_difference_in_baseline_rtn(
        SVector(-7.0, 3.0, 4.0), SVector(0.0, 1.0, 0.0), SVector(-5.0, 0.0, 0.0)) ≈ expected
    elements = [(raan=deg2rad(value),) for value in (170.0, -175.0, -160.0)]
    @test rad2deg.(CR.unwrap_raan(elements)) ≈ [170.0, 185.0, 200.0]
    @test all(isnan, CR.unwrap_raan([(raan=NaN,), (raan=NaN,)]))
    @test_throws ErrorException CR.unwrap_raan([(raan=0.0,), (raan=NaN,)])
end

function synthetic_rows(campaign; windows=(80, 100))
    rows = NamedTuple[]
    for case in CR.SC.grid_cases(campaign), window in windows
        scale = case.policy === :no_laser ? 0.0 :
            (case.policy === :naive_next_entering ? 1.0 : 0.8)
        value = scale * case.helpers / 800
        rtn = NamedTuple{CR.RTN_FIELDS}(ntuple(k -> k % 2 == 0 ? -value : value, length(CR.RTN_FIELDS)))
        velocity = NamedTuple{CR.VELOCITY_FIELDS}((value, -value, 2value, sqrt(6)value))
        push!(rows, merge((scheduler=string(case.policy), target_altitude_km=case.altitude_km,
            delta_i_deg=case.inclination_deg, helpers=case.helpers, orbits=window,
            time_s=window * CR.SC.mission_duration(case) / 100, row=10window + 1,
            run_id=CR.SC.case_signature(case), source_trajectory="synthetic_fixture_only",
            source_trajectory_sha256="synthetic_fixture_only",
            da_m=100value, de=-1e-5value, di_deg=case.inclination_deg * value / 10,
            draan_deg=case.inclination_deg == 0 ? NaN : -value - 0.01case.altitude_km,
            raan_defined=case.inclination_deg != 0), rtn, velocity))
    end
    rows
end

@testset "Paper figure grids, units, data ranges and zero baselines" begin
    mktempdir() do output
        for (campaign, difference) in (("range_along", false), ("no_laser", false), ("range_along", true))
            rows = synthetic_rows(campaign)
            root = joinpath(output, "$campaign-$difference")
            files = CR.write_figures(rows, root, (80, 100); difference)
            @test length(files) == (difference ? 12 : 8)
            @test all(path -> filesize(path) > 10_000, files)
            @test all(path -> endswith(path, ".png") ?
                read(path)[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] :
                String(read(path)[1:4]) == "%PDF", files)
            for window in (80, 100)
                directory = joinpath(root, "images", window == 100 ? "for_paper" : "80orbits")
                manifest = TOML.parsefile(joinpath(directory, "figure_manifest.toml"))
                @test manifest["cases"] == (campaign == "no_laser" ? 81 : 162)
                spec = CR.plot_spec(:oe, window; difference)
                table = CSV.read(joinpath(directory, spec.csv), DataFrame)
                @test all(table.row .== 10window + 1)
                @test all(isnan, table.draan_deg[table.delta_i_deg .== 0])
                @test Set(table.scheduler) == Set(string.(CR.SC.campaign_settings(campaign).policies))
                @test Set(table.target_altitude_km) == Set((1000, 1050, 1100))
            end
        end
        rows = synthetic_rows("range_along"; windows=(100,))
        directory = joinpath(output, "range_along-false", "images", "for_paper")
        spec = CR.plot_spec(:oe, 100)
        figure = CR.PP.make_plot(:oe; directory, output=joinpath(output, "inspected"),
            policies=("naive_next_entering", "positive_along_track"),
            colors=CR.COLORS, styles=CR.STYLES, policy_labels=CR.POLICY_LABELS,
            specification=spec, source_limits=CR.plot_limits(rows, spec))
        @test length(figure.subplots) == 13
        for (r, field) in enumerate(spec.fields), (col, altitude) in enumerate(CR.PP.ALTITUDES)
            panel = figure.subplots[(r - 1) * 3 + col]
            curves = filter(series -> series[:seriestype] == :path, panel.series_list)
            expected = [sort(filter(row -> row.scheduler == string(policy) &&
                row.target_altitude_km == altitude && row.delta_i_deg == inclination, rows);
                by=row -> row.helpers)
                for policy in CR.SC.campaign_settings("range_along").policies
                for inclination in CR.PP.INCLINATIONS if !(field == :draan_deg && inclination == 0)]
            @test length(curves) == length(expected)
            for (curve, data) in zip(curves, expected)
                @test collect(curve[:x]) == collect(CR.PP.HELPERS)
                @test collect(curve[:y]) == [getproperty(row, field) * spec.scales[r] for row in data]
                @test all(v -> ylims(panel)[1] <= v <= ylims(panel)[2], curve[:y])
            end
        end
    end
end

@testset "Real production-format propagation, saved-data verification and pairing" begin
    mkpath(RESULTS)
    mktempdir(RESULTS) do root
        provenance = CampaignPaths.campaign_sources(SC.source_hashes)
        tracks = Dict{Symbol,Any}()
        for policy in (:naive_next_entering, :positive_along_track, :no_laser, :gve_sma)
            case = SC.StudyCase(1050, 0.5, 1, policy)
            directory = OracleStudy.run_case(case; output_root=root, save_schedule=true, provenance)
            report = StudyVerification.verify_output(directory)
            @test report["passed"] && report["rows"] == 1001
            plot_case = CR.SC.StudyCase(1050, 0.5, 1, policy)
            track = CR.read_target(directory, plot_case)
            tracks[policy] = track
            @test CR.endpoint_indices(track, (80, 100)) == [801, 1001]
            if policy === :gve_sma
                original = Arrow.Table(joinpath(SC.case_directory(SC.DEFAULT_OUTPUT, case), "trajectory.feather"))
                reproduced = Arrow.Table(joinpath(directory, "trajectory.feather"))
                @test all(name -> original[name] == reproduced[name], propertynames(original))
            end
        end
        baseline = tracks[:no_laser]
        for policy in (:naive_next_entering, :positive_along_track)
            laser = tracks[policy]
            series = CR.difference_series(laser, baseline)
            @test series.rtn == laser.rtn
            @test maximum(series.velocity[4, :]) > 0
            @test series.elements[1, end] == laser.elements[end].a - baseline.elements[end].a
            @test series.elements[4, end] == rad2deg(laser.raan[end] - baseline.raan[end])
            @test all(isapprox.(series.velocity[4, :],
                [norm(laser.velocities[:, i] - baseline.velocities[:, i]) for i in eachindex(laser.time_s)]))
            path = CR.save_difference(laser, baseline, series, joinpath(root, "differences"))
            saved = Arrow.Table(path)
            @test length(saved.time_s) == 1001
            @test saved.velocity_difference_R_mps == series.velocity[1, :]
            @test saved.da_m == series.elements[1, :]
            @test saved.dv_T_mps == laser.rtn[findfirst(==(:dv_T_mps), CR.RTN_FIELDS), :]
            endpoints = CR.difference_endpoints(laser, baseline, series, path, (80, 100))
            @test [row.row for row in endpoints] == [801, 1001]
            @test_throws ErrorException CR.validate_pair(laser,
                merge(baseline, (time_s=baseline.time_s .+ 1.0,)))
            wrong_info = deepcopy(baseline.info)
            wrong_info["settings"]["J2"] = 0.0
            @test_throws ErrorException CR.validate_pair(laser, merge(baseline, (info=wrong_info,)))
            wrong_info = deepcopy(baseline.info)
            wrong_info["initial_states"][1]["position_m"][1] += 1.0
            @test_throws ErrorException CR.validate_pair(laser, merge(baseline, (info=wrong_info,)))
            wrong_info = deepcopy(baseline.info)
            wrong_info["source_hashes"]["changed_source"] = "different"
            @test_throws ErrorException CR.validate_pair(laser, merge(baseline, (info=wrong_info,)))
        end
        for (policy, helpers) in ((:naive_next_entering, 800), (:gve_inc, 1))
            case = SC.StudyCase(1100, 1, helpers, policy; orbits=0.02)
            directory = OracleStudy.run_case(case; output_root=root, provenance)
            @test StudyVerification.verify_output(directory)["passed"]
        end
    end
end

SC.atomic_toml(joinpath(RESULTS, "validation.toml"), Dict(
    "passed" => true, "completed_utc" => string(now(UTC)),
    "checks" => ["243-case grid", "range-only and positive-along-track transitions",
        "zero-power baseline", "common-frame velocity subtraction", "RAAN unwrapping",
        "80/100-orbit PNG/PDF and exact plotted data", "four full 100-orbit N1 cases",
        "bit-identical reproduction of saved Output 2 GVE-SMA trajectory",
        "800-helper production-format pilot", "GVE-INC regression pilot",
        "settings/source/time/initial-state mismatch rejection"]))
println("All production campaign checks passed.")
