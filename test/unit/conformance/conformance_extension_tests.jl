module ConformanceExtensionTests

using Test, SpaceAGORA, TOML
using SpaceAGORA.Conformance
using SpaceAGORA.Conformance: ConformanceCheck, ConformanceReport

const C = SpaceAGORA.Conformance
const ROOT = pkgdir(SpaceAGORA)

# Collects what test_report records without passing failures on to the suite.
mutable struct RecordingTestSet <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
end
RecordingTestSet(description; kwargs...) = RecordingTestSet(description, Any[])
Test.record(ts::RecordingTestSet, result) = (push!(ts.results, result); result)
function Test.finish(ts::RecordingTestSet)
    parent = Test.get_testset()
    parent isa RecordingTestSet && parent !== ts && Test.record(parent, ts)
    return ts
end
leaves(ts::RecordingTestSet) = reduce(vcat, (r isa RecordingTestSet ? leaves(r) : Any[r] for r in ts.results); init=Any[])

# Loads SpaceAGORA and Test in the order given and prints what it finds.
const LOAD_PROBE = """
    order = ARGS[1]
    test_id = Base.PkgId(Base.UUID("8dfed614-e22c-5e08-85e1-65c5234f0b40"), "Test")
    order == "test-first" && Core.eval(Main, :(using Test))
    using SpaceAGORA
    with_spaceagora = haskey(Base.loaded_modules, test_id)
    order == "spaceagora-first" && Core.eval(Main, :(using Test))
    ext = Base.get_extension(SpaceAGORA, :SpaceAGORATestExt)
    report = SpaceAGORA.Conformance.ConformanceReport
    owner = ext === nothing ? nothing : parentmodule(which(SpaceAGORA.Conformance.test_report, Tuple{report}))
    print(join((Base.JLOptions().use_compiled_modules, with_spaceagora, ext !== nothing, owner === ext), ","))
    """

function load_probe(order::String, flags::Vector{String}=String[])
    command = `$(Base.julia_cmd()) $flags --startup-file=no --project=$(ROOT) -e $LOAD_PROBE $order`
    fields = split(read(command, String), ',')
    return (compiled_modules=parse(Int, fields[1]), test_with_spaceagora=fields[2] == "true",
        extension=fields[3] == "true", methods_from_extension=fields[4] == "true")
end

failing_report() = ConformanceReport("T", :force_torque_model, [
    ConformanceCheck(:a, :pass, v"0.3.0", ""),
    ConformanceCheck(:b, :skip, v"0.3.0", "why"),
    ConformanceCheck(:c, :advisory_fail, v"0.4.0", "new check"),
    ConformanceCheck(:d, :fail, v"0.3.0", "what differed"),
])

@testset "Conformance reports and the Test extension" begin
    @testset "reports" begin
        @test_throws ArgumentError ConformanceCheck(:x, :unknown, v"0.3.0", "")
        good = [ConformanceCheck(:a, :pass, v"0.3.0", ""), ConformanceCheck(:b, :skip, v"0.3.0", "why"),
            ConformanceCheck(:c, :advisory_fail, v"0.4.0", "new check")]
        @test passed(ConformanceReport("T", :force_torque_model, good))
        failing = ConformanceReport("T", :force_torque_model, [good; ConformanceCheck(:d, :fail, v"0.3.0", "what differed")])
        @test !passed(failing)
        @test passed(good[2]) && !passed(failing.checks[end])
        # K-4: a check is advisory in the minor release that adds it, then required;
        # the first release's checks are required from the start.
        @test C._failure_status(v"0.3.0", v"0.2.0") === :fail
        @test C._failure_status(v"0.3.0", v"0.3.4") === :fail
        @test C._failure_status(v"0.4.0", v"0.3.9") === :advisory_fail
        @test C._failure_status(v"0.4.0", v"0.4.2") === :advisory_fail
        @test C._failure_status(v"0.4.0", v"0.5.0") === :fail
        text = sprint(show, MIME"text/plain"(), failing)
        @test occursin("Conformance of T (force_torque_model): FAILED", text)
        @test occursin("what differed", text) && occursin("advisory_fail", text)
        @test occursin("scientific validation stays with the model's author", text)
    end

    @testset "the extension defines the test functions and records each check" begin
        ext = Base.get_extension(SpaceAGORA, :SpaceAGORATestExt)
        @test ext isa Module
        for f in (test_report, test_force_torque_model, test_service_facade, test_sensor, test_guidance_model,
                test_navigation_model, test_control_model, test_algorithm, test_algorithm_family)
            @test !isempty(methods(f)) && all(m -> m.module === ext, methods(f))
        end
        local recorded
        @test_logs (:info, r"Advisory conformance check c") match_mode=:any begin
            recorded = @testset RecordingTestSet "recorded" begin
                @test test_report(failing_report()) isa ConformanceReport
            end
        end
        results = leaves(recorded)
        @test count(r -> r isa Test.Pass, results) == 2   # check a, and the isa above
        @test count(r -> r isa Test.Broken, results) == 2  # the skip and the advisory failure
        fails = filter(r -> r isa Test.Fail, results)
        @test length(fails) == 1 && occursin("what differed", something(only(fails).context, ""))
    end

    @testset "test_report writes the report file the contributor workflow reads" begin
        @test C.REPORT_DIRECTORY_VARIABLE == "SPACEAGORA_CONFORMANCE_REPORT_DIR"
        @test C._write_report_file(failing_report(), "") === nothing
        mktempdir() do dir
            target = joinpath(dir, "reports")
            @test_logs (:info, r"Advisory conformance check c") match_mode=:any begin
                withenv(C.REPORT_DIRECTORY_VARIABLE => target) do
                    @testset RecordingTestSet "written" begin
                        test_report(failing_report())
                        test_report(ConformanceReport("Second model", :service_facade, ConformanceCheck[]))
                    end
                end
            end
            files = readdir(target; join=true)
            @test length(files) == 2 && all(endswith(".toml"), files)
            written = TOML.parsefile(first(files))
            @test written["target"] == "T" && written["kind"] == "force_torque_model"
            @test [check["id"] for check in written["checks"]] == ["a", "b", "c", "d"]
            @test written["checks"][4] == Dict("id" => "d", "status" => "fail", "since" => "0.3.0", "detail" => "what differed")
            @test TOML.parsefile(last(files))["target"] == "Second model"
            @test occursin("Second_model", basename(last(files)))
        end
    end

    @testset "the error hint names Test when the extension is missing" begin
        exc = MethodError(test_report, (1,))
        @test occursin("using Test", sprint(C._test_extension_hint, exc, false))
        @test isempty(sprint(C._test_extension_hint, exc, true))
        @test isempty(sprint(C._test_extension_hint, MethodError(check_service_facade, (1,)), false))
    end

    @testset "Test is a weak dependency with its extension" begin
        project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
        @test !haskey(project["deps"], "Test")
        @test project["weakdeps"]["Test"] == "8dfed614-e22c-5e08-85e1-65c5234f0b40"
        @test project["extensions"]["SpaceAGORATestExt"] == "Test"
        @test isfile(joinpath(ROOT, "ext", "SpaceAGORATestExt.jl"))
    end

    @testset "the public API gate holds the Conformance section to the module's exports" begin
        gate = read(joinpath(ROOT, "test", "gates", "ci_public_api_surface_gate.jl"), String)
        symbols = read(joinpath(ROOT, "docs", "public_api_symbols.jl"), String)
        entry = "            (owner = :Conformance, symbol = :passed, rendered = \"SpaceAGORA.Conformance.passed\"),\n"
        @test occursin(entry, symbols)
        function run_gate(symbols_text)
            mktempdir() do dir
                mkpath(joinpath(dir, "docs"))
                mkpath(joinpath(dir, "test", "gates"))
                write(joinpath(dir, "docs", "public_api_symbols.jl"), symbols_text)
                path = joinpath(dir, "test", "gates", "ci_public_api_surface_gate.jl")
                write(path, gate)
                isolated = Module(:PublicAPIGateControl)
                Core.eval(isolated, :(include(path) = Base.include(@__MODULE__, path)))
                try
                    redirect_stdout(devnull) do
                        Base.include(isolated, path)
                    end
                    return nothing
                catch err
                    return err
                end
            end
        end
        @test run_gate(symbols) === nothing
        failure = run_gate(replace(symbols, entry => ""))
        @test failure isa LoadError && failure.error isa ErrorException
        @test occursin("exports of SpaceAGORA.Conformance differ: undocumented [:passed]", failure.error.msg)
    end

    @testset "the extension loads in either order and from existing caches only" begin
        for order in ("spaceagora-first", "test-first")
            loaded = load_probe(order)
            @testset let order = order, loaded = loaded
                @test loaded.compiled_modules == Base.JLOptions().use_compiled_modules
                @test loaded.extension && loaded.methods_from_extension
            end
        end
        existing = load_probe("spaceagora-first", ["--compiled-modules=existing"])
        @test existing.compiled_modules == 2
        @test existing.extension && existing.methods_from_extension
    end
end

end # module ConformanceExtensionTests
