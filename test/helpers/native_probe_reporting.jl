# Test-only protocol shared by the native coverage child and its parent.
module NativeProbeReporting
const PREFIX = "spaceagora_native_gram_probes="
const COMPLETED = PREFIX * "completed"
const SKIPPED = PREFIX * "skipped"

function native_probe_required(value::AbstractString=get(ENV, "SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES", "0"))
    value == "1" && return true
    value == "0" && return false
    throw(ArgumentError("SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES must be 0 or 1"))
end

function run_native_probes(run::Function; available::Bool,
        required::Bool=native_probe_required(), io::IO=stdout)
    if !available
        println(io, SKIPPED)
        required && error("Native GRAM probes are required but native setup is unavailable")
        return :skipped
    end
    run() # A failed include/test must return no completion marker.
    println(io, COMPLETED)
    return :completed
end

function native_probe_status(output::AbstractString)
    markers = filter(line -> startswith(line, PREFIX), strip.(split(output, '\n')))
    isempty(markers) && return :missing
    length(markers) == 1 || return :invalid
    only(markers) == COMPLETED && return :completed
    only(markers) == SKIPPED && return :skipped
    return :invalid
end

function report_native_probes(io::IO, output::AbstractString, child_succeeded::Bool;
        required::Bool=native_probe_required())
    status = native_probe_status(output)
    println(io, PREFIX, status, " required=", required, " child_success=", child_succeeded)
    return child_succeeded &&
        (status === :completed || (!required && status === :skipped))
end
end
