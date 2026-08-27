# The benchmark harness's mode table must agree with the shipped profile it names.
#
# benchmarks/studies/parallelization_performance/modes.jl duplicates the routing
# knobs that src/parallel/routing/profile_definitions.jl already declares. That
# duplicate has gone stale twice: R5's inner scheduler was reverted to "static"
# in 622ae2a0 and its thermal_mode to "auto" later, and in both cases the
# harness kept measuring the old value -- so the benchmark backing the paper's
# routing claims was reporting a profile SpaceAGORA does not ship.
#
# scripts/paired_profile_probe.jl solved this by deriving its environment from
# ParallelProfiles.profile_env_pairs, which is the shipped source of truth. This
# harness should do the same; until it does, this gate fails the build whenever
# a value in the table disagrees with the profile it claims to represent.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

push!(LOAD_PATH, REPO_ROOT)
using SpaceAGORA

const MODES_PATH = joinpath(REPO_ROOT, "benchmarks", "studies",
                            "parallelization_performance", "modes.jl")
isfile(MODES_PATH) || error("modes.jl not found at $(MODES_PATH)")
const MODES_SRC = read(MODES_PATH, String)

# (harness mode name => shipped profile). Only the modes that claim to BE a
# shipped profile are checked; the *_no* attribution arms deliberately deviate
# on the one knob they isolate, so they are checked on everything else via the
# base mode they are derived from.
const CHECKED = [
    ("serial", "R0"), ("outer_threads", "R1_a"), ("outer_process", "R1_b"),
    ("inner_only", "R2"), ("outer_inner_static", "R3"),
    ("outer_inner_adaptive", "R4"), ("full_smart", "R5"),
]

# Pull one `key="value"` out of the PPCModeSpec block for `mode`.
function _spec_field(mode::String, key::String)::Union{Nothing, String}
    anchor = "\"$(mode)\" => PPCModeSpec("
    i = findfirst(anchor, MODES_SRC)
    i === nothing && return nothing
    # The block ends at the next `" => PPCModeSpec(` or end of file.
    rest = MODES_SRC[last(i):end]
    nxt = findfirst("\" => PPCModeSpec(", rest[2:end])
    block = nxt === nothing ? rest : rest[1:first(nxt)]
    m = match(Regex("\\b" * key * "\\s*=\\s*\"([^\"]*)\""), block)
    return m === nothing ? nothing : String(m.captures[1])
end

violations = String[]
for (mode, profile) in CHECKED
    cfg = SpaceAGORA.ParallelProfiles.profile_config(profile)
    for (field, shipped) in (
        ("thermal", cfg.thermal_mode),
        ("density", cfg.density_mode),
        ("control", cfg.control_mode),
        ("multibody", cfg.multibody_mode),
        ("effector", cfg.effector_mode),
    )
        harness = _spec_field(mode, field)
        harness === nothing && continue
        if harness != shipped
            push!(violations,
                  "modes.jl mode \"$(mode)\" sets $(field)=\"$(harness)\" but " *
                  "profile $(profile) ships $(field)=\"$(shipped)\"")
        end
    end
end

if !isempty(violations)
    error("Benchmark mode table disagrees with the shipped profiles:\n" *
          join(violations, "\n") *
          "\n\nFix modes.jl (or the profile), do not silence this: the harness " *
          "is what backs the paper's routing numbers.")
end

println("ppc_mode_profile_parity_gate_ok")
