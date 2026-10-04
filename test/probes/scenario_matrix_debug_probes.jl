# Bounded utility probes. Run with julia --project=. test/probes/scenario_matrix_debug_probes.jl.
# All reports use synthetic summaries or populated result caches. No propagation,
# private telemetry, plotting, reference export, or native GRAM is exercised.
using Test
using DataFrames
using CSV
using SHA

const _DEBUG_PROBE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _DEBUG_PROBE_SCRIPTS = joinpath(_DEBUG_PROBE_ROOT, "scripts")
const _DEBUG_PROBE_MATRIX = joinpath(_DEBUG_PROBE_ROOT, "test", "gmat_scenario_matrix.jl")

module DebugProbeDefinitions end
module DebugProbeJ2 end
module DebugProbeTB end

function _debug_source_snapshot()
    paths = [_DEBUG_PROBE_MATRIX; [joinpath(_DEBUG_PROBE_SCRIPTS, name) for name in
        ("scenario_matrix_debug_support.jl", "tb_matrix_debug_defs_only.jl",
         "j2_parity_debug.jl", "tb_matrix_debug.jl")]]
    return Dict(path => bytes2hex(sha256(read(path))) for path in paths)
end

function _debug_capture_stdout(f)
    mktemp() do path, io
        redirect_stdout(io) do
            f()
        end
        flush(io)
        seekstart(io)
        return read(io, String)
    end
end

function _debug_summary()
    return DataFrame(scenario=String[], event=String[], rmse_km=Float64[])
end

function _debug_add_axes!(summary, scenario, values; axes=("x", "y", "z"))
    for (axis, value) in zip(axes, values)
        push!(summary, (scenario, "state_$(axis)_time", value))
    end
    return summary
end

function _debug_result(support, summary)
    return support.TV.VerificationResult(; summary, errors=DataFrame(),
        summary_path="fixture-summary.csv", errors_path="fixture-errors.csv",
        plots_dir="", profile=:quick, enforce=false, total_runtime_s=0.0)
end

function _debug_with_cached_results(f, support, basilisk, stk)
    caches = (support._BASILISK_MATRIX_RESULT_CACHE, support._BASILISK_MATRIX_SUMMARY_CACHE,
        support._BASILISK_MATRIX_CACHE_KEY, support._STK_MATRIX_RESULT_CACHE,
        support._STK_MATRIX_SUMMARY_CACHE, support._STK_MATRIX_CACHE_KEY)
    original = map(ref -> ref[], caches)
    try
        withenv("SPACEAGORA_GMAT_SCENARIOS" => nothing) do
            key = support._basilisk_matrix_cache_key()
            for (ref, value) in zip(caches, (basilisk, basilisk.summary, key, stk, stk.summary, key))
                ref[] = value
            end
            f()
        end
    finally
        for (ref, value) in zip(caches, original)
            ref[] = value
        end
    end
end

@testset "Matrix definitions load without executing diagnostics" begin
    before = _debug_source_snapshot()
    script_names = readdir(_DEBUG_PROBE_SCRIPTS)
    output = _debug_capture_stdout() do
        Base.include(DebugProbeDefinitions, joinpath(_DEBUG_PROBE_SCRIPTS, "tb_matrix_debug_defs_only.jl"))
    end
    support = DebugProbeDefinitions.ScenarioMatrixDebugSupport
    @test isempty(output)
    @test support.SCENARIO_MATRIX_DEFINITIONS_ONLY
    @test all(getfield(support, name)[] === nothing for name in
        (:_BASILISK_MATRIX_RESULT_CACHE, :_STK_MATRIX_RESULT_CACHE,
         :_CYGNSS_48HR_RESULT_CACHE, :_CYGNSS_96HR_RESULT_CACHE,
         :_CYGNSS_CYG04_96HR_RESULT_CACHE, :_CYGNSS_GMAT_RESULT_CACHE))
    @test all(isdefined(support, name) for name in
        (:_run_scenario_matrix_testsets, :_run_basilisk_scenario_matrix_result_once,
         :_run_stk_scenario_matrix_result_once, :_plot_cygnss_drag_force_timeseries,
         :_plot_cygnss_drag_tangential_timeseries, :_export_spaceagora_examples))

    # Repeated compatibility includes reuse the same definitions module.
    Base.include(DebugProbeDefinitions, joinpath(_DEBUG_PROBE_SCRIPTS, "tb_matrix_debug_defs_only.jl"))
    @test DebugProbeDefinitions.ScenarioMatrixDebugSupport === support
    for (mod, script) in ((DebugProbeJ2, "j2_parity_debug.jl"), (DebugProbeTB, "tb_matrix_debug.jl"))
        Core.eval(mod, :(const ScenarioMatrixDebugSupport = $support))
        @test isempty(_debug_capture_stdout(() -> Base.include(mod, joinpath(_DEBUG_PROBE_SCRIPTS, script))))
    end
    @test _debug_source_snapshot() == before
    @test readdir(_DEBUG_PROBE_SCRIPTS) == script_names
end

@testset "Direct matrix includes retain their runner boundary" begin
    # Replace only the campaign function with a counter before evaluating the
    # file. This exercises the real automatic-call condition without permitting
    # an accidental campaign even on machines with all private assets present.
    for mode in (nothing, false, true)
        mod = Module(gensym(:MatrixBoundaryProbe))
        calls = Ref(0)
        Core.eval(mod, :(const _probe_campaign_calls = $calls))
        if mode !== nothing
            Core.eval(mod, :(const SCENARIO_MATRIX_DEFINITIONS_ONLY = $mode))
        end
        replacements = Ref(0)
        mapexpr = function (expr)
            if expr isa Expr && expr.head == :function &&
                    expr.args[1] isa Expr && expr.args[1].head == :call &&
                    expr.args[1].args[1] == :_run_scenario_matrix_testsets
                replacements[] += 1
                return :(function _run_scenario_matrix_testsets()
                    _probe_campaign_calls[] += 1
                    nothing
                end)
            end
            return expr
        end
        Base.include(mapexpr, mod, _DEBUG_PROBE_MATRIX)
        @test replacements[] == 1
        @test calls[] == (mode === true ? 0 : 1)
    end
end

@testset "J2 default dispatch and ranking" begin
    support = DebugProbeDefinitions.ScenarioMatrixDebugSupport
    summary = _debug_summary()
    names = ["earth_j2_tbfalse", "earth_j2_tbtrue", "mars_j2_tbfalse",
             "mars_j2_tbtrue", "venus_j2_tbfalse", "moon_j2_tbfalse"]
    for (i, name) in enumerate(names)
        _debug_add_axes!(summary, name, (3.0i, 4.0i, 0.0))
    end
    _debug_add_axes!(summary, "venus_j0_tbfalse", (300.0, 400.0, 0.0))
    push!(summary, (names[1], "non_position_diagnostic", 999.0))
    basilisk = _debug_result(support, summary)
    stk = _debug_result(support, _debug_summary())
    checks = Ref(0)
    _debug_with_cached_results(support, basilisk, stk) do
        # The default runner itself is exercised; its matching cache is the
        # safety boundary and contains a different result from the STK cache.
        output = _debug_capture_stdout() do
            withenv("SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT" => "planet_j2",
                    "SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS" => names[1],
                    "SPACEAGORA_DEBUG_COMPARE_J2" => "1") do
                DebugProbeJ2.run_j2_parity_debug(; check_inputs=() -> (checks[] += 1))
            end
        end
        @test checks[] == 1
        @test occursin("J2 source default: planet_j2", output)
        @test occursin("Planet override list: " * names[1], output)
        @test occursin("Compare J2 analytic/generic: 1", output)
        @test !occursin("venus_j0_tbfalse", output)
        @test !occursin("non_position_diagnostic", output)
        worst = strip.(split(strip(last(split(output, "Worst J2 cases:"))), '\n'))
        @test worst == ["$(names[i]): xyz_norm_rmse=$(5.0i) km" for i in 6:-1:2]
        @test DebugProbeJ2._j2_combined_xyz_rmse(summary, names[1]) == 5.0
    end
    calls = Symbol[]
    @test_throws ArgumentError DebugProbeJ2.run_j2_parity_debug(;
        check_inputs=() -> (push!(calls, :check); throw(ArgumentError("missing fixture"))),
        runner=() -> (push!(calls, :run); basilisk))
    @test calls == [:check]
end

@testset "TB default dispatch, CSV compatibility, and missing cases" begin
    support = DebugProbeDefinitions.ScenarioMatrixDebugSupport
    basilisk_summary, stk_summary = _debug_summary(), _debug_summary()
    scenarios = ["$(body)_$(gravity)_$(tb)" for body in ("earth", "mars", "venus", "moon")
        for gravity in ("j0", "j2", "j50") for tb in ("tbfalse", "tbtrue")]
    expected = Dict{Tuple{String,String},Union{Nothing,Float64}}()
    for (i, name) in enumerate(scenarios)
        for (target, summary, scale) in (("GMAT", basilisk_summary, 1.0), ("STK", stk_summary, 2.0))
            values = (scale*i, 2scale*i, 2scale*i)
            expected[(name, target)] = 3000.0scale*i
            if target == "GMAT" && name in ("earth_j0_tbfalse", "earth_j0_tbtrue")
                values = (3.0, 4.0, 0.0)
                expected[(name, target)] = 5000.0
            elseif target == "GMAT" && name == "earth_j2_tbfalse"
                values = (0.001, 0.002, 0.003)
                expected[(name, target)] = sqrt(14.0)
            elseif target == "GMAT" && name == "moon_j50_tbtrue"
                _debug_add_axes!(summary, name, (1.0, 2.0); axes=("x", "y"))
                expected[(name, target)] = nothing
                continue
            elseif target == "STK" && name == "mars_j2_tbtrue"
                expected[(name, target)] = nothing
                continue
            end
            _debug_add_axes!(summary, name, values)
        end
    end
    push!(basilisk_summary, ("earth_j0_tbfalse", "non_position_diagnostic", 999.0))
    basilisk = _debug_result(support, basilisk_summary)
    stk = _debug_result(support, stk_summary)
    mktempdir() do tmp
        csv_path = joinpath(tmp, "tb_matrix_rmse_m.csv")
        checks = Ref(0)
        output = _debug_with_cached_results(support, basilisk, stk) do
            _debug_capture_stdout() do
                DebugProbeTB.run_tb_matrix_debug(; csv_path, check_inputs=() -> (checks[] += 1))
            end
        end
        @test checks[] == 1
        @test first(readlines(csv_path)) == "planet,gravity,target,third_body,rmse_m"
        result = CSV.read(csv_path, DataFrame)
        @test nrow(result) == 48
        @test names(result) == ["planet", "gravity", "target", "third_body", "rmse_m"]
        @test collect(zip(result.target[1:4], result.third_body[1:4])) ==
            [("GMAT", "tbfalse"), ("STK", "tbfalse"), ("GMAT", "tbtrue"), ("STK", "tbtrue")]
        @test count(ismissing, result.rmse_m) == 2
        for (i, row) in enumerate(eachrow(result))
            scenario = "$(row.planet)_$(row.gravity)_$(row.third_body)"
            @test scenario == scenarios[cld(i, 2)]
            @test row.target == (isodd(i) ? "GMAT" : "STK")
            value = expected[(scenario, row.target)]
            if value === nothing
                @test ismissing(row.rmse_m)
            else
                @test row.rmse_m ≈ value rtol=2eps(Float64) atol=0
            end
        end
        fractional = only(filter(line -> startswith(line, "earth,j2,GMAT,tbfalse,"), readlines(csv_path)))
        @test parse(Float64, last(split(fractional, ','))) != round(sqrt(14.0); digits=3)
        @test occursin("moon_j50_tbtrue (GMAT)", output)
        @test occursin("mars_j2_tbtrue (STK)", output)
        @test occursin("earth_j0 (GMAT): both = 5000.0 m", output)
        @test count(line -> occursin(": both = ", line), split(output, '\n')) == 1
        @test occursin("=== SCRIPT DONE ===", output)

        calls = Symbol[]
        blocked_path = joinpath(tmp, "blocked.csv")
        @test_throws ArgumentError DebugProbeTB.run_tb_matrix_debug(; csv_path=blocked_path,
            check_inputs=() -> (push!(calls, :check); throw(ArgumentError("missing STK fixture"))),
            basilisk_runner=() -> (push!(calls, :basilisk); basilisk),
            stk_runner=() -> (push!(calls, :stk); stk))
        @test calls == [:check]
        @test !ispath(blocked_path)
    end
end

@testset "Optional matrix input preflight" begin
    # Evaluate the actual preflight function with file-only fixture dependencies.
    # This checks every failure stage without mocking or invoking a propagator.
    source = Meta.parse(read(joinpath(_DEBUG_PROBE_SCRIPTS, "scenario_matrix_debug_support.jl"), String))
    functions = filter(expr -> expr isa Expr && expr.head == :function, source.args[3].args)
    preflight = only(filter(expr -> expr.args[1].args[1] == :require_matrix_inputs, functions))
    mktempdir() do tmp
        mod = Module(gensym(:MatrixInputProbe))
        basilisk_dir = joinpath(tmp, "Basilisk_Examples_Full")
        stk_dir = joinpath(tmp, "stk_results")
        scenarios = Ref(Set(["earth_j0_tbfalse"]))
        Core.eval(mod, quote
            const _GMAT_REPO_ROOT = $tmp
            const _BASILISK_REFERENCE_DIR = $basilisk_dir
            const _STK_RESULTS_DIR = $stk_dir
            const _FETCH_REFERENCES_HINT = "fetch fixture references"
            const _GMAT_HARMONICS_EARTH_FILE = "earth.csv"
            const _GMAT_HARMONICS_MARS_FILE = "mars.csv"
            const _GMAT_HARMONICS_VENUS_FILE = "venus.csv"
            const _GMAT_HARMONICS_MOON_FILE = "moon.csv"
            const _probe_scenarios = $scenarios
            _basilisk_reference_available() = isdir(_BASILISK_REFERENCE_DIR)
            _stk_reference_available() = isdir(_STK_RESULTS_DIR)
            _active_basilisk_expected_scenario_names() = copy(_probe_scenarios[])
            _scenario_basilisk_path(name) = joinpath(_BASILISK_REFERENCE_DIR, name * ".feather")
            _scenario_stk_path(name) = joinpath(_STK_RESULTS_DIR, name * ".csv")
            function _gmat_planetary_kernel_relpath()
                isfile(joinpath(_GMAT_REPO_ROOT, "de430.bsp")) || throw(ArgumentError("missing fixture kernel"))
                return "de430.bsp"
            end
        end)
        Core.eval(mod, preflight)
        check(; kwargs...) = Base.invokelatest(() -> getfield(mod, :require_matrix_inputs)(; kwargs...))
        @test_throws ArgumentError check()
        mkpath(basilisk_dir)
        empty!(scenarios[])
        @test_throws ArgumentError check()
        push!(scenarios[], "earth_j0_tbfalse")
        @test_throws ArgumentError check() # selected reference missing
        touch(joinpath(basilisk_dir, "earth_j0_tbfalse.feather"))
        @test_throws ArgumentError check(; stk=true) # STK directory missing
        mkpath(stk_dir)
        @test_throws ArgumentError check(; stk=true) # STK scenario missing
        touch(joinpath(stk_dir, "earth_j0_tbfalse.csv"))
        @test_throws ArgumentError check() # gravity files missing
        for name in ("earth.csv", "mars.csv", "venus.csv", "moon.csv")
            touch(joinpath(tmp, name))
        end
        @test_throws ArgumentError check() # kernel missing
        touch(joinpath(tmp, "de430.bsp"))
        @test check() === nothing
        @test check(; stk=true) === nothing
        rm(joinpath(stk_dir, "earth_j0_tbfalse.csv"))
        @test check() === nothing # Basilisk-only diagnostic has no STK dependency
        @test_throws ArgumentError check(; stk=true)
    end
end

@testset "Full-arc entrypoint uses the maintained definition boundary" begin
    # A separate process checks the command entrypoint. Empty arguments must
    # reach its usage check before any campaign or output write.
    script = joinpath(_DEBUG_PROBE_SCRIPTS, "xval_fullarc.jl")
    code = """
        using Test
        err = try
            include($(repr(script)))
            main()
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("usage: scripts/xval_fullarc.jl", sprint(showerror, err))
        @test isempty(_BASILISK_MATRIX_CACHE_KEY[])
        @test _BASILISK_MATRIX_RESULT_CACHE[] === nothing
        @test _STK_MATRIX_RESULT_CACHE[] === nothing
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --project=$(_DEBUG_PROBE_ROOT) -e $code`
    @test success(pipeline(cmd; stdout=devnull))
    support = DebugProbeDefinitions.ScenarioMatrixDebugSupport
    @test !isdefined(support, :_LUNA_STK_ADJUSTED_HARMONICS_FILE)
    @test support._matrix_scenario_overrides("moon_j2_tbfalse", :stk)["gravity_harmonics_file"] ==
          support._STK_HARMONICS_MOON_FILE
    @test support._matrix_scenario_overrides("moon_j2_tbfalse", :gmat)["gravity_harmonics_file"] ==
          support._GMAT_HARMONICS_MOON_FILE
    @test support._matrix_scenario_overrides("mars_j2_tbfalse", :gmat)["gravity_harmonics_order"] == 0
    @test support._matrix_scenario_overrides("mars_j2_tbfalse", :basilisk)["gravity_harmonics_order"] == 2
end
