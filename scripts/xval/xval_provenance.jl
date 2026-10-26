# Provenance for full-arc runs. Kept independent of the simulation so failure and
# tampering checks can run without a numerical campaign or external datasets.
module XvalProvenance
using SHA, JSON

file_identity(path) = Dict("path" => abspath(path), "sha256" => bytes2hex(open(sha256, path)))
const PRIMARY_CONTROLS = Set(["XVAL_SCENARIOS", "SPACEAGORA_WARN_NORMALIZE",
                              "SPACEAGORA_WARN_DEPRECATED_CONFIG"])
controls(env=ENV) = Dict(k => v for (k, v) in env if
    (startswith(k, "XVAL_") || startswith(k, "SPACEAGORA_")) && !isempty(v))
function check_primary_controls(target, variant; env=ENV)
    overrides = setdiff(Set(keys(controls(env))), PRIMARY_CONTROLS)
    if target in ("gmat", "stk") && variant == "committed" && !isempty(overrides)
        throw(ArgumentError("Committed primary runs refuse scientific overrides: " *
            join(sort!(collect(overrides)), ", ") * ". Use an explicit sensitivity variant."))
    end
    nothing
end
function source_identity(root)
    try
        Dict("commit" => readchomp(`git -C $root rev-parse HEAD`),
             "tree" => readchomp(`git -C $root rev-parse 'HEAD^{tree}'`),
             "dirty" => !isempty(readchomp(`git -C $root status --porcelain --untracked-files=all`)))
    catch
        Dict("commit" => "unknown", "tree" => "unknown", "dirty" => true)
    end
end
function verify_inputs(inputs)
    for input in inputs
        file_identity(input["path"])["sha256"] == input["sha256"] ||
            error("Full-arc input changed during the run: $(input["path"])")
    end
    nothing
end
function write_record(rundir, record)
    path = joinpath(rundir, "run_info.json")
    open(path * ".tmp", "w") do io
        JSON.print(io, record, 2)
        println(io)
    end
    mv(path * ".tmp", path; force=true)
end
function start_record(root, rundir, target, variant, scenarios)
    check_primary_controls(target, variant)
    ispath(rundir) && error("Full-arc output already exists; choose a new run directory: $rundir")
    record = Dict{String,Any}(
        "schema_version" => 1, "status" => "running", "target" => target,
        "variant" => variant, "scenarios" => scenarios, "source" => source_identity(root),
        "controls" => controls(), "cases" => Dict{String,Any}(),
        "reference_provenance" => "unverified_generation_settings")
    mkpath(rundir)
    write_record(rundir, record)
    record
end
function finish_record!(record, root, rundir)
    record["source_end"] = source_identity(root)
    record["source_end"] == record["source"] || error("Source identity changed during the full-arc run")
    controls() == record["controls"] || error("Full-arc controls changed during the run")
    Set(keys(record["cases"])) == Set(record["scenarios"]) || error("Incomplete full-arc case records")
    for case in values(record["cases"])
        verify_inputs(case["inputs"])
        verify_inputs(case["kernels"])
    end
    record["results_sha256"] = file_identity(joinpath(rundir, "results.csv"))["sha256"]
    record["status"] = "complete"
    write_record(rundir, record)
    record
end
end
