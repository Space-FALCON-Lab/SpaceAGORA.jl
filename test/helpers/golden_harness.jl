using TOML

"""
Root directory holding one subdirectory per reproducible-propagation
("golden") scenario. Each scenario directory contains:

- `config.jl`      — defines a zero-argument function `build_<scenario>_config()`
                     returning a `SimulationConfiguration`.
- `fixture.csv`     — the checked-in golden trajectory: one row per comparison
                     time, columns `time, pos_atol_m, vel_atol_mps`, then
                     `sc<k>_pos_<1|2|3>` / `sc<k>_vel_<1|2|3>` for each
                     spacecraft `k` (1-indexed, matching the order the
                     spacecraft were passed to `build_config`/`build_config_multi`).
- `metadata.toml`   — provenance: `git_sha`, `julia_version`, `generated_at`,
                     `description`, `tier` (`"pr"` or `"nightly"`).

See `test/README.md` for how to add a new scenario and
`scripts/regenerate_golden.jl` for how to (re)generate `fixture.csv`.
"""
const GOLDEN_ROOT = joinpath(REPO_ROOT, "test", "golden")

"""
    run_and_compare_golden(scenario_name, build_config_fn; spacecraft_count=1)

Run the `SimulationConfiguration` returned by `build_config_fn()` and compare
the resulting trajectory against `test/golden/<scenario_name>/fixture.csv` at
each fixture time (linearly interpolated), for each of `spacecraft_count`
spacecraft. Skips (via `@test_skip`) if the fixture is absent, matching the
convention used before this scenario had its own directory.
"""
function run_and_compare_golden(
    scenario_name::AbstractString,
    build_config_fn::Function;
    spacecraft_count::Int=1
)
    scenario_dir = joinpath(GOLDEN_ROOT, scenario_name)
    fixture_path = joinpath(scenario_dir, "fixture.csv")
    if !isfile(fixture_path)
        @test_skip "Golden fixture is not present for scenario '$scenario_name' ($fixture_path)"
        return
    end
    fixture = CSV.read(fixture_path, DataFrame)

    args = build_config_fn()
    df = run_case(args)
    @test nrow(df) > 1
    times = Vector{Float64}(df.time)

    for row in eachrow(fixture)
        t = Float64(row.time)
        pos_atol = Float64(row.pos_atol_m)
        vel_atol = Float64(row.vel_atol_mps)

        for k in 1:spacecraft_count
            for c in 1:3
                pos_col = Symbol("sc$(k)_pos_$(c)")
                vel_col = Symbol("sc$(k)_vel_$(c)")
                sim_pos = interp_linear(times, getproperty(df, pos_col), t)
                sim_vel = interp_linear(times, getproperty(df, vel_col), t)
                @test isapprox(sim_pos, Float64(row[pos_col]); atol=pos_atol, rtol=0.0)
                @test isapprox(sim_vel, Float64(row[vel_col]); atol=vel_atol, rtol=0.0)
            end
        end
    end
end

"""
    golden_scenario_tier(scenario_name) -> Symbol

Reads `tier` from `test/golden/<scenario_name>/metadata.toml` (`:pr` or
`:nightly`, defaulting to `:pr` if unspecified/absent).
"""
function golden_scenario_tier(scenario_name::AbstractString)
    metadata_path = joinpath(GOLDEN_ROOT, scenario_name, "metadata.toml")
    isfile(metadata_path) || return :pr
    metadata = TOML.parsefile(metadata_path)
    return Symbol(get(metadata, "tier", "pr"))
end

"""
    golden_scenario_enabled(scenario_name) -> Bool

Whether `scenario_name` should run given the `SPACEAGORA_GOLDEN_TIER`
environment variable (`"pr"` by default; set to `"all"` to also run
`:nightly`-tier scenarios, e.g. from the nightly-stress CI workflow).
"""
function golden_scenario_enabled(scenario_name::AbstractString)
    requested_tier = get(ENV, "SPACEAGORA_GOLDEN_TIER", "pr")
    requested_tier == "all" && return true
    return golden_scenario_tier(scenario_name) == Symbol(requested_tier)
end
