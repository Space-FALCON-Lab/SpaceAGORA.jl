const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPO_PROJECT = joinpath(REPO_ROOT, "Project.toml")

if something(Base.active_project(), "") != REPO_PROJECT
    import Pkg
    Pkg.activate(REPO_ROOT; io=devnull)
end

using SpaceAGORA

function _instantiate_vendored_gramsuite!(vendored_gramsuite::String)
    Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    previous_project = something(Base.active_project(), "")
    try
        Pkg.activate(vendored_gramsuite; io=devnull)
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

function _missing_dependency_error(err)::Bool
    msg = sprint(showerror, err)
    return occursin("is required but does not seem to be installed", msg) ||
        occursin("Run `Pkg.instantiate()` to install all recorded dependencies", msg)
end

function ensure_gramsuite_loaded!()
    if !isdefined(SpaceAGORA, :GRAMSuite)
        vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
        try
            if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
                pushfirst!(LOAD_PATH, vendored_gramsuite)
            end
            @eval import GRAMSuite
        catch err
            if isdir(vendored_gramsuite) && _missing_dependency_error(err)
                _instantiate_vendored_gramsuite!(vendored_gramsuite)
                try
                    @eval import GRAMSuite
                    return nothing
                catch retry_err
                    throw(ErrorException(
                        "GRAM-backed examples could not load `GRAMSuite` after instantiating its vendored project. " *
                        "Original load failure: $(sprint(showerror, err)) | " *
                        "Retry failure: $(sprint(showerror, retry_err))"
                    ))
                end
            end
            throw(ErrorException(
                "GRAM-backed examples require loading `GRAMSuite` so the SpaceAGORA GRAM extension can load. " *
                "Original load failure: $(sprint(showerror, err))"
            ))
        end
    end
    return nothing
end

function setup_gram_example!(mod::Module=@__MODULE__)
    ensure_gramsuite_loaded!()
    if !isdefined(mod, :InitialTime)
        Core.eval(mod, :(const InitialTime = $SM.InitialTime))
    end
    if !isdefined(mod, :GRAMAtmosphereModel)
        Core.eval(mod, :(const GRAMAtmosphereModel = $SM.GRAMAtmosphereModel))
    end
    return nothing
end

if !isdefined(@__MODULE__, :SimulationEngine)
    const SimulationEngine = SpaceAGORA.SimulationEngine
end
if !isdefined(@__MODULE__, :SimulationModel)
    const SimulationModel = SpaceAGORA.SimulationModel
end
if !isdefined(@__MODULE__, :RuntimeServices)
    const RuntimeServices = SpaceAGORA.RuntimeServices
end
using .SimulationModel

if !isdefined(@__MODULE__, :SM)
    const SM = SimulationModel
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SpaceAGORA.run_simulation
end
if !isdefined(@__MODULE__, :quat_mult)
    const quat_mult = SimulationModel.quat_mult
end
if !isdefined(@__MODULE__, :SPICE_PATH)
    const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
end

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft, run_and_report
