# Shared options for local, reproducible viewer demonstrations.
# This helper never removes an existing directory or reuses old results.
function viewer_demo_options(name, default_duration; argv=ARGS, env=ENV)
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    duration = Float64(default_duration)
    output = joinpath(get(env, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(root, "output", "viewer_demos")), name)
    seen = Set{String}()
    i = 1
    while i <= length(argv)
        key = argv[i]
        key in ("--duration-s", "--output-dir") || throw(ArgumentError("unknown viewer demo option: $key"))
        key in seen && throw(ArgumentError("duplicate viewer demo option: $key"))
        push!(seen, key)
        i < length(argv) || throw(ArgumentError("$key requires a value"))
        value = argv[i + 1]
        if key == "--duration-s"
            parsed = tryparse(Float64, value)
            parsed === nothing && throw(ArgumentError("--duration-s must be a finite positive number"))
            duration = parsed
        else
            isempty(strip(value)) && throw(ArgumentError("--output-dir must not be empty"))
            output = value
        end
        i += 2
    end
    isfinite(duration) && duration > 0 || throw(ArgumentError("--duration-s must be a finite positive number"))
    output = abspath(output)
    ispath(output) && (!isdir(output) || !isempty(readdir(output))) &&
        throw(ArgumentError("output path must be absent or an empty directory: $output"))
    return (duration_s=duration, output_dir=output)
end
