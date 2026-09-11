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

function print_thread_diagnostics(args; label::String="")
    PP  = SimulationModel.ParallelPolicy
    SE  = SimulationEngine
    sep = "─" ^ 60

    header = isempty(label) ? "Thread Diagnostics" : "Thread Diagnostics — $label"
    println("[SpaceAGORA] $header")
    println(sep)

    total_threads = Threads.nthreads()
    inner_budget  = PP.effective_inner_thread_budget()
    budget_note   = inner_budget < total_threads ? "  ← capped by SPACEAGORA_INNER_THREAD_BUDGET" : ""
    println("  Julia threads available   : $total_threads")
    println("  Inner thread budget       : $inner_budget$budget_note")

    effector_cap  = SE._effector_max_threads()
    effector_mode = SE._effector_parallel_mode()
    heavy_only    = SE._effector_heavy_only()
    work_threshold_µs = round(Int, SE._effector_work_ns_per_worker_threshold() / 1e3)
    println("  Effector thread cap       : $effector_cap  (SPACEAGORA_EFFECTOR_MAX_THREADS)")
    println("  Effector parallel mode    : $effector_mode")
    println("  Heavy-work gate           : $heavy_only  (SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY)")

    calib_raw  = lowercase(strip(get(ENV, "SPACEAGORA_RHS_CALIBRATE", "auto")))
    effs       = args.dynamics_model.dynamic_effectors
    n_sats     = length(args.dynamics_model.spacecraft)
    calib_note = n_sats < 2 ? "  ← skipped (requires ≥2 active satellites)" : ""
    println("  RHS plan calibration      : $calib_raw$calib_note")

    println()
    println("  Effectors ($(length(effs))):")
    for (i, eff) in enumerate(effs)
        safe   = SE._dynamic_effector_threadsafe(eff) ? "thread-safe" : "NOT thread-safe"
        tname  = string(nameof(typeof(eff)))
        detail = if eff isa GravitationalHarmonicsModel
            " ($(eff.L)×$(eff.M) deg/ord)"
        elseif eff isa NBodyGravityModel
            " ($(join(eff.body_names, ", ")))"
        else
            ""
        end
        println("    [$i] $tname$detail  —  $safe")
    end

    println()
    plan = SE._rhs_execution_plan(args, nothing, effs, n_sats)
    ed   = plan.effector_decision
    dom  = hasproperty(plan, :dominant_axis) ? string(plan.dominant_axis) : string(plan.mode)
    println("  Routing ($n_sats satellite$(n_sats == 1 ? "" : "s"), cold cost model):")
    println("    mode            : $(plan.mode)")
    println("    dominant axis   : $dom")
    println("    effector threads: use=$(ed.use_threads)  allotment=$(ed.allotment)  mode=$(ed.mode)")
    if !ed.use_threads && heavy_only
        println("    → heavy-only gate: threading activates once observed cost/worker ≥ $work_threshold_µs µs")
        harm_idx = findfirst(e -> e isa GravitationalHarmonicsModel, effs)
        if harm_idx !== nothing
            h = effs[harm_idx]
            println("      ($(h.L)×$(h.M) harmonics will exceed this threshold within the first few RHS calls)")
        end
    elseif ed.use_threads
        println("    → $(ed.allotment) effector thread(s) active from first RHS call")
    end
    println(sep)
    return nothing
end
