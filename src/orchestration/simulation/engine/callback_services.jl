# Implement the callback boundary through the established engine operations.
# The original engine functions retain their identities and numerical bodies.
import ..CallbackServices

@inline CallbackServices._effector_force_ii(args...; kwargs...) = _effector_force_ii(args...; kwargs...)
@inline CallbackServices._initialize_save_cache_buffers!(args...; kwargs...) = _initialize_save_cache_buffers!(args...; kwargs...)
@inline CallbackServices._inside_atmosphere_at_state(args...; kwargs...) = _inside_atmosphere_at_state(args...; kwargs...)
@inline CallbackServices._is_gravity_backbone_state(args...; kwargs...) = _is_gravity_backbone_state(args...; kwargs...)
@inline CallbackServices._rhs_plan_step_cache_enabled(args...; kwargs...) = _rhs_plan_step_cache_enabled(args...; kwargs...)
@inline CallbackServices._run_solver_config(args...; kwargs...) = _run_solver_config(args...; kwargs...)
@inline CallbackServices._set_state_position_velocity_ii!(args...; kwargs...) = _set_state_position_velocity_ii!(args...; kwargs...)
@inline CallbackServices._state_heat_loads(args...; kwargs...) = _state_heat_loads(args...; kwargs...)
@inline CallbackServices._state_mass_kg(args...; kwargs...) = _state_mass_kg(args...; kwargs...)
@inline CallbackServices._state_position_ii(args...; kwargs...) = _state_position_ii(args...; kwargs...)
@inline CallbackServices._state_quaternion(args...; kwargs...) = _state_quaternion(args...; kwargs...)
@inline CallbackServices._state_velocity_ii(args...; kwargs...) = _state_velocity_ii(args...; kwargs...)
@inline CallbackServices.build_state_sample(args...; kwargs...) = build_state_sample(args...; kwargs...)
@inline CallbackServices.sample_atmosphere(args...; kwargs...) = sample_atmosphere(args...; kwargs...)
@inline CallbackServices.sample_buffered_atmosphere(args...; kwargs...) = sample_buffered_atmosphere(args...; kwargs...)
@inline CallbackServices.sample_environment(args...; kwargs...) = sample_environment(args...; kwargs...)
@inline CallbackServices.sample_planet_frame(args...; kwargs...) = sample_planet_frame(args...; kwargs...)
