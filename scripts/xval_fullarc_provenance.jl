# Metadata for full-arc exports. This records execution identity, not independent
# scientific validation of the external reference trajectories.
module XvalFullarcProvenance
using SHA
using TOML

const SCHEMA_VERSION = 1

file_digest(path::AbstractString) = bytes2hex(open(sha256, path))

# Include future model/solver overrides by default, so adding a sensitivity knob
# cannot silently make a committed run eligible for primary export.
function input_overrides(env=ENV)
    return Dict{String, String}(String(k) => String(v) for (k, v) in env
        if (startswith(k, "XVAL_") && k != "XVAL_SCENARIOS" ||
            startswith(k, "SPACEAGORA_SPICE_") || startswith(k, "SPACEAGORA_TELEMETRY_") ||
            startswith(k, "SPACEAGORA_SOLVER_") || k == "SPACEAGORA_GMAT_PARITY_SOLVER") && !isempty(v))
end

function run_info(target::String, variant::String, scenarios; env=ENV)
    overrides = input_overrides(env)
    primary = target in ("gmat", "stk") && variant == "committed"
    if primary && !isempty(overrides)
        throw(ArgumentError("committed primary runs do not accept model/reference overrides: " *
            join(sort!(collect(keys(overrides))), ", ") *
            ". Use an explicit sensitivity variant name."))
    end
    return Dict{String, Any}(
        "provenance_schema" => SCHEMA_VERSION,
        "status" => "running", "target" => target, "variant" => variant,
        "primary_eligible" => primary, "input_overrides" => overrides,
        "scenarios" => String.(scenarios))
end

# Replace the completion marker before touching any previous run artifacts.
function write_info(path::AbstractString, info)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        TOML.print(io, info)
    end
    mv(tmp, path; force=true)
    return nothing
end

function finish_run!(rundir::AbstractString, info)
    artifacts = Dict{String, String}("results.csv" => file_digest(joinpath(rundir, "results.csv")))
    for scenario in info["scenarios"], file in ("manifest.toml", "series.arrow")
        relative = scenario * "/" * file
        artifacts[relative] = file_digest(joinpath(rundir, relative))
    end
    info["artifacts_sha256"] = artifacts
    info["status"] = "complete"
    write_info(joinpath(rundir, "run_info.toml"), info)
    return nothing
end
end
