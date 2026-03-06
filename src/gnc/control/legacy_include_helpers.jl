if !isdefined(@__MODULE__, :_sa_include_once!)
    function _sa_include_once!(flag::Symbol, path::String)
        if !isdefined(@__MODULE__, flag)
            include(path)
            @eval const $(flag) = true
        end
        return nothing
    end
end

@inline function _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_density_models, joinpath(@__DIR__, "..", "..", "environment", "atmosphere", "density_models.jl"))
    _sa_include_once!(:_sa_inc_montecarlo_perturbations, joinpath(@__DIR__, "..", "..", "mission", "campaigns", "montecarlo_perturbations.jl"))
    _sa_include_once!(:_sa_inc_reference_system, joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    _sa_include_once!(:_sa_inc_misc_utils, joinpath(@__DIR__, "..", "..", "core", "utils", "misc.jl"))
    return nothing
end

@inline function _sa_include_control_runtime_deps!()
    _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_gravity_models, joinpath(@__DIR__, "..", "..", "environment", "gravity", "gravity_models.jl"))
    _sa_include_once!(:_sa_inc_thermal_models, joinpath(@__DIR__, "..", "..", "vehicle", "thermal", "thermal_models.jl"))
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    return nothing
end

@inline function _sa_include_control_entrypoint_deps!()
    _sa_include_control_runtime_deps!()
    _sa_include_once!(:_sa_inc_time_switch_calcs, joinpath(@__DIR__, "heatload_control", "time_switch_calcs.jl"))
    _sa_include_once!(:_sa_inc_second_tsw_calcs, joinpath(@__DIR__, "heatload_control", "second_tsw_calcs.jl"))
    _sa_include_once!(:_sa_inc_security_mode, joinpath(@__DIR__, "heatload_control", "security_mode.jl"))
    return nothing
end

@inline function _sa_include_heatload_utils_deps!()
    _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    return nothing
end

@inline function _sa_include_heatload_switch_deps!()
    _sa_include_once!(:_sa_inc_utils_timeswitch, joinpath(@__DIR__, "heatload_control", "utils_timeswitch.jl"))
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    _sa_include_once!(:_sa_inc_eoms, joinpath(@__DIR__, "eoms.jl"))
    return nothing
end

@inline function _sa_include_heatload_second_switch_deps!()
    _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    _sa_include_once!(:_sa_inc_eoms, joinpath(@__DIR__, "eoms.jl"))
    return nothing
end

@inline function _sa_include_heatload_security_deps!()
    _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    return nothing
end

@inline function _sa_include_targeting_runtime_deps!()
    _sa_include_closed_form_core_deps!()
    _sa_include_once!(:_sa_inc_thermal_models, joinpath(@__DIR__, "..", "..", "vehicle", "thermal", "thermal_models.jl"))
    _sa_include_once!(:_sa_inc_closed_form_solution, joinpath(@__DIR__, "closed_form_solution.jl"))
    return nothing
end
