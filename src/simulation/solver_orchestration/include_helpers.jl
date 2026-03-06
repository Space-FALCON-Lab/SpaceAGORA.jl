if !isdefined(@__MODULE__, :_sa_solver_include_once!)
    function _sa_solver_include_once!(flag::Symbol, path::String)
        if !isdefined(@__MODULE__, flag)
            include(path)
            @eval const $(flag) = true
        end
        return nothing
    end
end

@inline function _sa_include_solver_orchestration_deps!()
    _sa_solver_include_once!(:_sa_inc_reference_system, joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
    _sa_solver_include_once!(:_sa_inc_reference_system_config, joinpath(@__DIR__, "..", "..", "core", "state", "reference_system_config.jl"))
    return nothing
end

@inline function _sa_include_solver_jacobian_deps!()
    _sa_include_solver_orchestration_deps!()
    _sa_solver_include_once!(:_sa_inc_quaternion_utils, joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
    _sa_solver_include_once!(:_sa_inc_gravity_models, joinpath(@__DIR__, "..", "..", "environment", "gravity", "gravity_models.jl"))
    _sa_solver_include_once!(:_sa_inc_density_models, joinpath(@__DIR__, "..", "..", "environment", "atmosphere", "density_models.jl"))
    _sa_solver_include_once!(:_sa_inc_aerodynamic_models, joinpath(@__DIR__, "..", "..", "dynamics", "translational", "aerodynamic_models.jl"))
    _sa_solver_include_once!(:_sa_inc_thermal_models, joinpath(@__DIR__, "..", "..", "vehicle", "thermal", "thermal_models.jl"))
    _sa_solver_include_once!(:_sa_inc_perturbations, joinpath(@__DIR__, "..", "..", "dynamics", "coupled", "perturbations.jl"))
    _sa_solver_include_once!(:_sa_inc_control_runtime, joinpath(@__DIR__, "..", "..", "gnc", "control", "control.jl"))
    _sa_solver_include_once!(:_sa_inc_propulsive_maneuvers, joinpath(@__DIR__, "..", "..", "gnc", "control", "propulsive_maneuvers.jl"))
    return nothing
end
