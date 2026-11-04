#!/usr/bin/env julia
# Opt-in, real-package loading diagnostic. No loader installation, workers,
# native constructors, or simulation calls. Package imports can run __init__.
module GRAMPackageReadiness
using TOML
using SHA
using Dates

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CASES = ("discoverable_space_first", "discoverable_gram_first",
               "fallback_space_first", "fallback_gram_first", "absent_gram", "absent_dependency")
const SPACE_UUID = "afbfb69f-5c0b-4832-b760-43725dff8540"
const GRAM_UUID = "b50455af-6a46-4eae-bf92-8039261dd674"
const USAGE = """
Usage: julia --startup-file=no --compiled-modules=existing \\
  test/smoke/gram_package_readiness.jl --gram-project /prepared/GRAMSuite.jl \\
  --output /new/report/directory [--space-project /prepared/SpaceAGORA.jl] \\
  [--timeout-seconds 300]

Uses only existing packages/depot. Each case runs in a fresh Julia process with
an explicit LOAD_PATH, offline Pkg settings, and automatic precompilation off.
The SpaceAGORA project must resolve its dependencies but must not already expose
GRAMSuite: fallback and absent-package cases require initial nondiscovery.
Six cases must pass; missing prerequisites, timeouts and import failures are
nonzero exits, never passing skips. Success covers this Julia version and these
observed package sources only. Native and worker readiness remain unverified.
"""

function options(args)
    result = Dict("space-project" => ROOT, "timeout-seconds" => "300")
    allowed = Set(("space-project", "gram-project", "output", "timeout-seconds", "child"))
    length(args) % 2 == 0 || error(USAGE)
    for i in 1:2:length(args)
        startswith(args[i], "--") || error(USAGE)
        key = args[i][3:end]
        key in allowed || error("Unknown option: $(args[i])")
        result[key] = args[i + 1]
    end
    for key in ("space-project", "gram-project", "output")
        haskey(result, key) || error("Missing --$key\n$USAGE")
        result[key] = abspath(result[key])
    end
    timeout = tryparse(Int, result["timeout-seconds"])
    timeout !== nothing && timeout > 0 || error("Timeout must be a positive integer")
    haskey(result, "child") && !(result["child"] in CASES) && error("Unknown child case")
    return result
end

sha(path) = bytes2hex(sha256(read(path)))
function file_metadata(root)
    result = Dict{String, Any}("path" => realpath(root), "files" => Dict{String, String}())
    for name in ("Project.toml", "Manifest.toml", "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
        path = joinpath(root, name)
        isfile(path) && (result["files"][name] = sha(path))
    end
    sources = String[]
    for dir in ("src", "ext")
        isdir(joinpath(root, dir)) || continue
        for (base, _, files) in walkdir(joinpath(root, dir))
            for file in sort(files)
                endswith(file, ".jl") || continue
                path = joinpath(base, file)
                push!(sources, relpath(path, root) * "\t" * sha(path))
            end
        end
    end
    result["julia_source_sha256"] = bytes2hex(sha256(join(sort(sources), "\n")))
    result["julia_source_file_count"] = length(sources)
    # An enclosing repository is not evidence of the package's own revision.
    git = Sys.which("git")
    if git !== nothing
        try
            top = readchomp(pipeline(`$git -C $root rev-parse --show-toplevel`, stderr=devnull))
            if realpath(top) == realpath(root)
                result["git_revision"] = readchomp(`$git -C $root rev-parse HEAD`)
                result["git_status"] = readchomp(`$git -C $root status --porcelain --untracked-files=no`)
            end
        catch
            # Hashes still identify the actual files in a source-only snapshot.
        end
    end
    result["revision_available"] = haskey(result, "git_revision")
    return result
end

function project_identity(path, name, uuid)
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    get(project, "name", "") == name || error("Expected $name project at $path")
    get(project, "uuid", "") == uuid || error("Unexpected $name UUID at $path")
    return project
end

function environment()
    Dict("active_project" => something(Base.active_project(), ""),
         "load_path" => copy(LOAD_PATH), "depot_path" => copy(DEPOT_PATH),
         "directory" => pwd(), "args" => copy(ARGS))
end

function require_fact(report, key, condition)
    report["facts"][key] = condition
    condition || error("Unmet loading contract: $key")
    return nothing
end

function package_metadata(mod)
    path = pathof(mod)
    path === nothing && error("Package has no source path: $mod")
    result = file_metadata(dirname(dirname(path)))
    result["entry_path"] = realpath(path)
    result["entry_sha256"] = sha(path)
    result["module"] = string(mod)
    result["package_id"] = string(Base.PkgId(mod))
    result["version"] = string(Base.pkgversion(mod))
    return result
end

# Julia 1.12 also ages module bindings introduced by imports. Observe those
# bindings in the latest world without moving the counter probe's calling frame.
binding(mod, name) = Base.invokelatest(getfield, mod, name)

function imported_pair!(caller, gram_first)
    if gram_first
        Core.eval(caller, :(import GRAMSuite))
        Core.eval(caller, :(import SpaceAGORA))
    else
        Core.eval(caller, :(import SpaceAGORA))
        Core.eval(caller, :(import GRAMSuite))
    end
    return binding(caller, :SpaceAGORA), binding(caller, :GRAMSuite)
end

# This frame begins before either import. The only executed extension function
# reads two Julia atomic counters; it does not construct or query native models.
function positive_case!(report, opts, caller)
    case = opts["child"]
    gram_first = endswith(case, "gram_first")
    fallback = startswith(case, "fallback")
    before = Base.find_package("GRAMSuite")
    report["initial_gram_discovery"] = something(before, "")
    require_fact(report, "initial_discovery_matches_case", fallback ? before === nothing : before !== nothing)
    if fallback
        pushfirst!(LOAD_PATH, opts["gram-project"])
    end
    expected_path = copy(LOAD_PATH)
    report["import_load_path"] = expected_path
    require_fact(report, "selected_gram_source", realpath(Base.find_package("GRAMSuite")) ==
                 realpath(joinpath(opts["gram-project"], "src", "GRAMSuite.jl")))
    sa, gram = imported_pair!(caller, gram_first)
    report["SpaceAGORA"] = package_metadata(sa)
    report["GRAMSuite"] = package_metadata(gram)
    report["StaticArrays"] = package_metadata(binding(gram, :StaticArrays))
    require_fact(report, "imported_gram_source", realpath(pathof(gram)) ==
                 realpath(joinpath(opts["gram-project"], "src", "GRAMSuite.jl")))
    require_fact(report, "selected_space_source", realpath(pathof(sa)) ==
                 realpath(joinpath(opts["space-project"], "src", "SpaceAGORA.jl")))
    require_fact(report, "space_uuid", string(Base.PkgId(sa).uuid) == SPACE_UUID)
    require_fact(report, "gram_uuid", string(Base.PkgId(gram).uuid) == GRAM_UUID)
    ext = Base.get_extension(sa, :SpaceAGORAGRAMSuiteExt)
    report["extension_attached"] = ext !== nothing
    require_fact(report, "extension_attached", ext !== nothing)
    report["extension_source"] = string(pathof(ext))
    sm = binding(sa, :SimulationModel)
    em = binding(sm, :EnvironmentModels)
    model = binding(em, :GRAMAtmosphereModel)
    gram_model = binding(gram, :GRAMAtmosphereModel)
    require_fact(report, "canonical_type_owner", parentmodule(model) === em)
    require_fact(report, "extension_model_module_identity", binding(ext, :EM) === em)
    require_fact(report, "wrapper_and_package_types_distinct", model !== gram_model)
    require_fact(report, "initial_time_identity", binding(sa, :InitialTime) === binding(sm, :InitialTime))
    method = Base.invokelatest(which, model, Tuple{})
    report["model_type"] = string(model)
    report["model_parent_module"] = string(parentmodule(model))
    report["wrapper_model_type"] = string(gram_model)
    report["constructor_method"] = Dict("owner" => string(method.module),
        "signature" => string(method.sig), "file" => string(method.file), "line" => method.line)
    require_fact(report, "constructor_owned_by_extension", method.module === ext)
    require_fact(report, "shared_native_lock_identity", binding(ext, :GRAM_LOCK) ===
                 binding(binding(sa, :RuntimeServices), :GRAM_LOCK))
    require_fact(report, "wrapper_lock_hook_identity", binding(gram, :_GRAM_DEFAULT_LOCK_HOOK)[] ===
                 binding(ext, :GRAM_LOCK))
    require_fact(report, "ephemeris_hook_identity", binding(gram, :_GRAM_EPHEMERIS_STATE_FN)[] ===
                 binding(ext, :_gram_spice_ephemeris_state))
    require_fact(report, "caller_aliases_preserved", binding(caller, :InitialTime) === :initial_time_sentinel &&
                 binding(caller, :GRAMAtmosphereModel) === :atmosphere_sentinel)
    repeated_sa, repeated_gram = imported_pair!(caller, gram_first)
    require_fact(report, "repeated_import_identity", repeated_sa === sa && repeated_gram === gram &&
                 Base.get_extension(sa, :SpaceAGORAGRAMSuiteExt) === ext)
    require_fact(report, "load_path_preserved_during_import", LOAD_PATH == expected_path)
    # Direct old-frame visibility is observational: importing from a compiled
    # frame can require invokelatest. Readiness requires the safe route to work.
    probe = binding(ext, :_native_atmospheres_live)
    direct = try
        Dict("outcome" => "returned", "value" => probe())
    catch err
        Dict("outcome" => "error", "error_type" => string(typeof(err)),
             "message" => sprint(showerror, err))
    end
    report["older_frame_direct_probe"] = direct
    result = Base.invokelatest(probe)
    report["invokelatest_counter_probe"] = result
    require_fact(report, "nonnative_probe_zero_live_models", result === 0)
    require_fact(report, "zero_native_models_created", binding(ext, :_NATIVE_ATMOSPHERES_CREATED)[] === 0)
    require_fact(report, "zero_native_models_finalized", binding(ext, :_NATIVE_ATMOSPHERES_FINALIZED)[] === 0)
    require_fact(report, "native_wrapper_not_loaded", binding(gram, :_GRAM_WRAPPER)[] === nothing)
    require_fact(report, "native_wrapper_file_unset", binding(gram, :_GRAM_WRAPPER_FILE)[] == "")
    report["package_loading_ready"] = true
end

function child(opts)
    report = Dict{String, Any}("case" => opts["child"], "status" => "unmet_acceptance",
        "julia_version" => string(VERSION), "julia_executable" => joinpath(Sys.BINDIR, Base.julia_exename()),
        "extension_attached" => false, "package_loading_ready" => false,
        "native_readiness" => "unverified", "worker_readiness" => "unverified",
        "facts" => Dict{String, Bool}(), "started_at_utc" => string(now(UTC)),
        "julia_settings" => Dict(key => get(ENV, key, "") for key in
            ("JULIA_LOAD_PATH", "JULIA_DEPOT_PATH", "JULIA_PKG_OFFLINE", "JULIA_PKG_PRECOMPILE_AUTO")))
    original = environment()
    original_env = Dict(ENV)
    report["environment_before"] = original
    original_path = copy(LOAD_PATH)
    try
        caller = Module(:GRAMReadinessCaller)
        Core.eval(caller, :(const InitialTime = :initial_time_sentinel))
        Core.eval(caller, :(const GRAMAtmosphereModel = :atmosphere_sentinel))
        if opts["child"] == "absent_gram"
            require_fact(report, "gram_initially_absent", Base.find_package("GRAMSuite") === nothing)
            Core.eval(caller, :(import SpaceAGORA))
            sa = binding(caller, :SpaceAGORA)
            report["SpaceAGORA"] = package_metadata(sa)
            require_fact(report, "space_available_without_gram", realpath(pathof(sa)) ==
                realpath(joinpath(opts["space-project"], "src", "SpaceAGORA.jl")))
            require_fact(report, "extension_absent_without_gram",
                Base.get_extension(sa, :SpaceAGORAGRAMSuiteExt) === nothing)
            report["boundary"] = "expected_missing_package_boundary"
            import_error = try
                Core.eval(caller, :(import GRAMSuite))
                nothing
            catch err
                err
            end
            require_fact(report, "missing_package_import_failed", import_error !== nothing)
            report["expected_import_error"] = sprint(showerror, import_error)
            require_fact(report, "missing_package_error_identified", occursin("Package GRAMSuite not found", report["expected_import_error"]))
            require_fact(report, "extension_absent_after_failed_gram_import",
                Base.get_extension(sa, :SpaceAGORAGRAMSuiteExt) === nothing)
            require_fact(report, "no_gram_binding", !isdefined(caller, :GRAMSuite))
            report["negative_control_passed"] = true
            require_fact(report, "absent_case_load_path_unchanged", LOAD_PATH == original_path)
        elseif opts["child"] == "absent_dependency"
            report["boundary"] = "expected_missing_dependency_boundary"
            require_fact(report, "gram_discoverable_without_dependency", realpath(Base.find_package("GRAMSuite")) ==
                realpath(joinpath(opts["gram-project"], "src", "GRAMSuite.jl")))
            require_fact(report, "staticarrays_initially_absent", Base.find_package("StaticArrays") === nothing)
            import_error = try
                Core.eval(caller, :(import GRAMSuite))
                nothing
            catch err
                err
            end
            require_fact(report, "missing_dependency_import_failed", import_error !== nothing)
            report["expected_import_error"] = sprint(showerror, import_error)
            message = report["expected_import_error"]
            require_fact(report, "missing_dependency_error_identified", occursin("StaticArrays", message) &&
                (occursin("not found", message) || occursin("does not seem to be installed", message)))
            require_fact(report, "no_gram_binding", !isdefined(caller, :GRAMSuite))
            require_fact(report, "absent_case_load_path_unchanged", LOAD_PATH == original_path)
            report["negative_control_passed"] = true
        else
            positive_case!(report, opts, caller)
        end
        report["status"] = "passed"
    catch err
        report["error_type"] = string(typeof(err))
        report["error"] = sprint(showerror, err, catch_backtrace())
    finally
        # Restoring our explicit test-only fallback does not exercise or change
        # any production loader's restoration or activation policy.
        report["environment_before_restoration"] = environment()
        empty!(LOAD_PATH)
        append!(LOAD_PATH, original_path)
        report["environment_after"] = environment()
        report["environment_variables_preserved"] = Dict(ENV) == original_env
        report["environment_preserved"] = report["environment_after"] == original && report["environment_variables_preserved"]
        if !report["environment_preserved"]
            report["status"] = "unmet_acceptance"
        end
        report["finished_at_utc"] = string(now(UTC))
        open(joinpath(opts["output"], opts["child"] * ".toml"), "w") do io
            TOML.print(io, report; sorted=true)
        end
    end
    return report["status"] == "passed" ? 0 : 1
end

function parent(opts)
    project = project_identity(opts["space-project"], "SpaceAGORA", SPACE_UUID)
    gram_project = project_identity(opts["gram-project"], "GRAMSuite", GRAM_UUID)
    # Use Julia's existing Pkg version-spec parser, not a hand-written compat
    # approximation; this only reads metadata and does not activate or resolve.
    pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    compat = project["compat"]["julia"]
    Base.invokelatest(in, VERSION, Base.invokelatest(pkg.Types.semver_spec, compat)) || error("Julia $VERSION does not satisfy $compat")
    gram_julia_compat = gram_project["compat"]["julia"]
    Base.invokelatest(in, VERSION, Base.invokelatest(pkg.Types.semver_spec, gram_julia_compat)) ||
        error("Julia $VERSION does not satisfy GRAMSuite's $gram_julia_compat compatibility")
    gram_version = VersionNumber(gram_project["version"])
    gram_compat = project["compat"]["GRAMSuite"]
    Base.invokelatest(in, gram_version, Base.invokelatest(pkg.Types.semver_spec, gram_compat)) ||
        error("GRAMSuite $gram_version does not satisfy SpaceAGORA's $gram_compat compatibility")
    out = opts["output"]
    ispath(out) && error("Output must be a new directory: $out")
    mkpath(out)
    summary = Dict{String, Any}("julia_version" => string(VERSION), "julia_compat" => compat,
        "coverage" => "this_runtime_only", "gram_version" => string(gram_version),
        "gram_compat" => gram_compat, "gram_julia_compat" => gram_julia_compat, "SpaceAGORA_source" => file_metadata(opts["space-project"]),
        "GRAMSuite_source" => file_metadata(opts["gram-project"]),
        "native_readiness" => "unverified", "worker_readiness" => "unverified",
        "cases" => Dict{String, Any}())
    timeout = parse(Int, opts["timeout-seconds"])
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    for case in CASES
        paths = [opts["space-project"], "@stdlib"]
        startswith(case, "discoverable") && insert!(paths, 2, opts["gram-project"])
        depots = copy(DEPOT_PATH)
        if case == "absent_dependency"
            paths = [opts["gram-project"], "@stdlib"]
            empty_depot = joinpath(out, "empty-depot")
            mkdir(empty_depot)
            depots = [empty_depot]
        end
        sep = Sys.iswindows() ? ';' : ':'
        cmd = `$julia --startup-file=no --history-file=no --compiled-modules=existing --depwarn=error --project=$(opts["space-project"]) $(@__FILE__) --child $case --space-project $(opts["space-project"]) --gram-project $(opts["gram-project"]) --output $out`
        cmd = addenv(cmd, "JULIA_LOAD_PATH" => join(paths, sep), "JULIA_DEPOT_PATH" => join(depots, sep),
                     "JULIA_PKG_OFFLINE" => "true", "JULIA_PKG_PRECOMPILE_AUTO" => "0")
        process = open(joinpath(out, case * ".log"), "w") do io
            run(pipeline(ignorestatus(cmd), stdout=io, stderr=io); wait=false)
        end
        timed_out = timedwait(() -> process_exited(process), timeout; pollint=0.1) === :timed_out
        if timed_out
            kill(process)
            timedwait(() -> process_exited(process), 5; pollint=0.1) === :timed_out && kill(process, Base.SIGKILL)
        end
        wait(process)
        path = joinpath(out, case * ".toml")
        result = isfile(path) ? TOML.parsefile(path) : Dict("status" => "missing_report")
        passed = !timed_out && success(process) && get(result, "status", "") == "passed"
        summary["cases"][case] = Dict("passed" => passed, "exit_code" => process.exitcode,
            "timed_out" => timed_out, "report_status" => result["status"],
            "package_loading_ready" => get(result, "package_loading_ready", false),
            "negative_control_passed" => get(result, "negative_control_passed", false))
        println(case, ": ", passed ? "passed" : "unmet acceptance")
    end
    summary["positive_cases_ready"] = all(summary["cases"][case]["passed"] &&
        summary["cases"][case]["package_loading_ready"] for case in CASES if !startswith(case, "absent"))
    summary["negative_boundaries_passed"] = all(summary["cases"][case]["passed"] &&
        summary["cases"][case]["negative_control_passed"] for case in CASES if startswith(case, "absent"))
    summary["passed"] = summary["positive_cases_ready"] && summary["negative_boundaries_passed"]
    open(joinpath(out, "summary.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    return summary["passed"] ? 0 : 1
end

function main(args=ARGS)
    args == ["--help"] && (println(USAGE); return 0)
    opts = options(args)
    return haskey(opts, "child") ? child(opts) : parent(opts)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(GRAMPackageReadiness.main())
end
