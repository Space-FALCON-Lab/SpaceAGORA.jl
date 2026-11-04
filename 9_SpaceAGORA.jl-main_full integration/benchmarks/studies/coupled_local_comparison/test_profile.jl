# Evaluate only the profiler's solve wrapper with synthetic configurations: no
# SpaceAGORA import, native library, allocation profiler, or simulation runs.
module ProfileConfigTests
using Test

function find_solve(ex)
    ex isa Expr || return nothing
    if ex.head in (:function, :(=)) && ex.args[1] isa Expr &&
            ex.args[1].head == :call && ex.args[1].args[1] === :solve
        return ex
    end
    for arg in ex.args
        found = find_solve(arg)
        found === nothing || return found
    end
    return nothing
end

Base.@kwdef struct Settings
    results::Bool = true
    verbose::Bool = false
    results_directory::String = "synthetic-output"
    save_csv::Bool = true
end

const SpaceAGORA = (SimulationModel=(SimConfig=(
    _with_configuration=(args; overrides...) -> merge(args, (; overrides...)),
),),)

case = ""
cfg = (profile="smoke", seed=123)

# Model the documented case builder's N-dependent durations; the wrapper must
# preserve whichever duration and other settings its original case supplies.
function ppc_single_config(name, cfg; seed)
    n = parse(Int, match(r"^stack([0-9]+)_", name)[1])
    saved = endswith(name, "_saved")
    duration = cfg.profile == "test" ? 10.0 : cfg.profile == "smoke" ?
        (n >= 128 && !saved ? 300.0 : 120.0) : (n >= 128 && !saved ? 3600.0 : 600.0)
    return (mission_configuration=(mission_time=duration, data_rate=10.0),
        simulation_settings=Settings(results=saved), spacecraft_count=n, seed=seed)
end
ppc_solve_once(args, cfg) = args

definition = find_solve(Meta.parseall(read(joinpath(@__DIR__, "profile_allocations.jl"), String)))
definition === nothing && error("profile solve wrapper not found")
Core.eval(@__MODULE__, definition)

@testset "saving-off profile preserves the original workload" begin
    for n in (32, 128, 256), profile in ("test", "smoke", "full")
        global case = "stack$(n)_e6_actuated_saved"
        global cfg = (; profile, seed=123)
        saved = solve(case)
        off = solve(replace(case, r"_saved$" => ""))
        @test saved.simulation_settings.results
        @test !off.simulation_settings.results
        @test saved.mission_configuration == off.mission_configuration
        @test saved.spacecraft_count == off.spacecraft_count && saved.seed == off.seed
        for name in fieldnames(Settings)
            name === :results && continue
            @test getfield(saved.simulation_settings, name) == getfield(off.simulation_settings, name)
        end
    end
end
end
