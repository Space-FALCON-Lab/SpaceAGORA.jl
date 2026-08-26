# Shared helpers for the telemetry validation study (Odyssey + Venus Express
# aerobraking reconstructions of record). See README.md in this directory.

const STUDY_DIR = normpath(@__DIR__)
const REPO_ROOT = normpath(joinpath(STUDY_DIR, "..", "..", ".."))
const MANIFEST_DIR = joinpath(STUDY_DIR, "manifests")
const VENDORED_GRAMSUITE = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
const RECORD_DOC = "docs/spaceagora_aerobraking_reconstruction_record.md"

# One entry per run of the record. `needs_gram` runs require the native GRAM
# Suite build (data/GRAMSuite.jl/GRAM Suite 2.0 + scripts/ensure_gram_native.jl).
const RUN_SPECS = (
    (
        name="odyssey_tolson",
        manifest=joinpath(MANIFEST_DIR, "odyssey_tolson.toml"),
        needs_gram=false,
        description="Odyssey primary reconstruction (Tolson accelerometer-derived density)"
    ),
    (
        name="odyssey_marsgram",
        manifest=joinpath(MANIFEST_DIR, "odyssey_marsgram.toml"),
        needs_gram=true,
        description="Odyssey MarsGRAM TES Mapping Year 2 comparator (context run, ~4.3-5.8 h)"
    ),
    (
        name="odyssey_tolson_sigma_minus",
        manifest=joinpath(MANIFEST_DIR, "odyssey_tolson_sigma_minus.toml"),
        needs_gram=false,
        description="Density-uncertainty envelope, -1 sigma (informational)"
    ),
    (
        name="odyssey_tolson_sigma_plus",
        manifest=joinpath(MANIFEST_DIR, "odyssey_tolson_sigma_plus.toml"),
        needs_gram=false,
        description="Density-uncertainty envelope, +1 sigma (informational; impacts before campaign end)"
    ),
    (
        name="vex_venusgram",
        manifest=joinpath(MANIFEST_DIR, "vex_venusgram.toml"),
        needs_gram=true,
        description="Venus Express 2014 reconstruction with VenusGRAM (~18 min)"
    ),
)

const RUN_ALIASES = Dict(
    "core" => ["odyssey_tolson", "vex_venusgram"],
    "envelope" => ["odyssey_tolson_sigma_minus", "odyssey_tolson_sigma_plus"],
    "all" => [spec.name for spec in RUN_SPECS],
)

run_spec(name::AbstractString) = begin
    idx = findfirst(spec -> spec.name == name, RUN_SPECS)
    idx === nothing && error(
        "Unknown run '$name'. Available: $(join([spec.name for spec in RUN_SPECS], ", ")) " *
        "or aliases $(join(sort(collect(keys(RUN_ALIASES))), ", "))."
    )
    RUN_SPECS[idx]
end

function resolve_runs(raw::AbstractString)::Vector{String}
    names = String[]
    for token in split(raw, ",")
        token = strip(token)
        isempty(token) && continue
        expanded = get(RUN_ALIASES, token, nothing)
        for name in something(expanded, [String(token)])
            run_spec(name)  # validates
            name in names || push!(names, String(name))
        end
    end
    isempty(names) && error("No runs selected (got '$raw').")
    return names
end

results_root(out_root::AbstractString="") = isempty(out_root) ?
    joinpath(STUDY_DIR, "results", gethostname()) : abspath(out_root)

# Mirrors benchmarks/studies/telemetry_orbit_accuracy_study.jl: the vendored
# GRAMSuite submodule keeps its own project and must be instantiated against
# it before `import GRAMSuite` can succeed from the root environment.
function instantiate_vendored_gramsuite!()
    # Pkg may be loaded fresh here (newer world than this running function),
    # so every Pkg call goes through invokelatest.
    Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    previous_project = something(Base.active_project(), "")
    try
        Base.invokelatest(Pkg.activate, VENDORED_GRAMSUITE; io=devnull)
        Base.invokelatest(Pkg.instantiate; io=devnull)
    finally
        if !isempty(previous_project)
            Base.invokelatest(Pkg.activate, dirname(previous_project); io=devnull)
        else
            Base.invokelatest(Pkg.activate, REPO_ROOT; io=devnull)
        end
    end
    return nothing
end

function load_gramsuite!()
    try
        @eval import GRAMSuite
    catch err
        if isdir(VENDORED_GRAMSUITE)
            if Base.find_package("GRAMSuite") === nothing
                pushfirst!(LOAD_PATH, VENDORED_GRAMSUITE)
            end
            instantiate_vendored_gramsuite!()
            @eval import GRAMSuite
        else
            rethrow(err)
        end
    end
    return nothing
end

function parse_kv_args(argv::Vector{String}, positional_profiles=("quick", "full"))
    opts = Dict{String, String}()
    for arg in argv
        if arg in positional_profiles
            opts["profile"] = arg
        elseif startswith(arg, "--") && occursin("=", arg)
            k, v = split(arg[3:end], "="; limit=2)
            opts[String(k)] = String(v)
        else
            error("Unsupported argument '$arg'. Use [quick|full] and --key=value flags.")
        end
    end
    return opts
end

parse_bool_flag(raw::AbstractString) = lowercase(strip(raw)) in ("1", "true", "yes", "on")
