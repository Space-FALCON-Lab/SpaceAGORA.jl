if !isdefined(@__MODULE__, :_bridge_include_once!)
    function _bridge_include_once!(flag::Symbol, path::String)
        if !isdefined(@__MODULE__, flag)
            include(path)
            @eval const $(flag) = true
        end
        return nothing
    end
end

@inline function _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_density_models, joinpath(@__DIR__, "..", "..", "environment", "atmosphere", "density_models.jl"))
    _bridge_include_once!(:_bridge_inc_montecarlo_perturbations, joinpath(@__DIR__, "..", "..", "mission", "campaigns", "montecarlo_perturbations.jl"))
    _bridge_include_once!(:_bridge_inc_reference_system, joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    _bridge_include_once!(:_bridge_inc_misc_utils, joinpath(@__DIR__, "..", "..", "core", "utils", "misc.jl"))
    return nothing
end

@inline function _bridge_include_control_runtime_deps!()
    _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_gravity_models, joinpath(@__DIR__, "..", "..", "environment", "gravity", "gravity_models.jl"))
    _bridge_include_once!(:_bridge_inc_thermal_models, joinpath(@__DIR__, "..", "..", "vehicle", "thermal", "thermal_models.jl"))
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    return nothing
end

@inline function _bridge_include_control_entrypoint_deps!()
    _bridge_include_control_runtime_deps!()
    _bridge_include_once!(:_bridge_inc_time_switch_calcs, joinpath(@__DIR__, "..", "control", "heatload_control", "time_switch_calcs.jl"))
    _bridge_include_once!(:_bridge_inc_second_tsw_calcs, joinpath(@__DIR__, "..", "control", "heatload_control", "second_tsw_calcs.jl"))
    _bridge_include_once!(:_bridge_inc_security_mode, joinpath(@__DIR__, "..", "control", "heatload_control", "security_mode.jl"))
    return nothing
end

@inline function _bridge_include_heatload_utils_deps!()
    _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    return nothing
end

@inline function _bridge_include_heatload_switch_deps!()
    _bridge_include_once!(:_bridge_inc_utils_timeswitch, joinpath(@__DIR__, "..", "control", "heatload_control", "utils_timeswitch.jl"))
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    _bridge_include_once!(:_bridge_inc_eoms, joinpath(@__DIR__, "..", "control", "eoms.jl"))
    return nothing
end

@inline function _bridge_include_heatload_second_switch_deps!()
    _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    _bridge_include_once!(:_bridge_inc_eoms, joinpath(@__DIR__, "..", "control", "eoms.jl"))
    return nothing
end

@inline function _bridge_include_heatload_security_deps!()
    _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    return nothing
end

@inline function _bridge_include_targeting_runtime_deps!()
    _bridge_include_closed_form_core_deps!()
    _bridge_include_once!(:_bridge_inc_thermal_models, joinpath(@__DIR__, "..", "..", "vehicle", "thermal", "thermal_models.jl"))
    _bridge_include_once!(:_bridge_inc_closed_form_solution, joinpath(@__DIR__, "..", "control", "closed_form_solution.jl"))
    return nothing
end

if !isdefined(@__MODULE__, :CONTROL_BRIDGE_STATE_LOCK)
    const CONTROL_BRIDGE_STATE_LOCK = ReentrantLock()
end

if !isdefined(@__MODULE__, :_bridge_get_cnf)
    @inline function _bridge_get_cnf(args=nothing; cnf=nothing)
        if cnf !== nothing
            return cnf
        end
        if args isa NamedTuple && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        if !(args isa NamedTuple) && hasproperty(args, :cnf)
            return getproperty(args, :cnf)
        end
        if (@isdefined config) && isdefined(config, :cnf)
            return getproperty(config, :cnf)
        end
        throw(ArgumentError("Control state `cnf` not found. Pass `cnf=` or typed args.cnf field."))
    end
end

if !isdefined(@__MODULE__, :_bridge_get_solution)
    @inline function _bridge_get_solution(args=nothing; solution=nothing, cnf=nothing)
        if solution !== nothing
            return solution
        end
        if args isa NamedTuple && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if !(args isa NamedTuple) && hasproperty(args, :solution)
            return getproperty(args, :solution)
        end
        if cnf !== nothing && hasproperty(cnf, :solution)
            return getproperty(cnf, :solution)
        end
        if (@isdefined config) && isdefined(config, :solution)
            return getproperty(config, :solution)
        end
        throw(ArgumentError("Solution state `solution` not found. Pass `solution=` or typed args.solution field."))
    end
end
