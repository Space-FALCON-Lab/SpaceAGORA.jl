const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const VENDORED_GRAMSUITE = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")

function instantiate_vendored_gramsuite!()
    Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    previous_project = something(Base.active_project(), "")
    try
        Pkg.activate(VENDORED_GRAMSUITE; io=devnull)
        Pkg.instantiate(; io=devnull)
    finally
        if !isempty(previous_project)
            Pkg.activate(dirname(previous_project); io=devnull)
        else
            Pkg.activate(REPO_ROOT; io=devnull)
        end
    end
    return nothing
end

if Base.find_package("GRAMSuite") === nothing && isdir(VENDORED_GRAMSUITE)
    pushfirst!(LOAD_PATH, VENDORED_GRAMSUITE)
end

try
    import GRAMSuite
catch err
    if isdir(VENDORED_GRAMSUITE)
        instantiate_vendored_gramsuite!()
        @eval import GRAMSuite
    else
        rethrow(err)
    end
end

using SpaceAGORA
using SpaceAGORA.TelemetryVerification

function _study_cli_args(argv::Vector{String})::Vector{String}
    args = copy(argv)
    has_plot_flag = any(arg -> startswith(arg, "--plots="), args)
    if !has_plot_flag && !haskey(ENV, "SPACEAGORA_TELEMETRY_PLOTS")
        push!(args, "--plots=false")
    end
    return args
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_verification_cli(_study_cli_args(ARGS))
end
