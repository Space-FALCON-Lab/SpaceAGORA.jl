module CallbackServices

# Internal operations used by callbacks. The engine installs implementations
# after its state and sampling code has loaded. Callbacks depend on these
# explicit operations, rather than discovering an engine module at runtime.
# Solver policies, frame conventions and mutation rules remain engine-owned.

function _initialize_save_cache_buffers! end
function _inside_atmosphere_at_state end
function _is_gravity_backbone_state end
function _rhs_plan_step_cache_enabled end
function _set_state_position_velocity_ii! end
function _solver_policy_mode end
function _state_heat_loads end
function _state_mass_kg end
function _state_position_ii end
function _state_quaternion end
function _state_velocity_ii end
function build_state_sample end
function sample_atmosphere end
function sample_buffered_atmosphere end
function sample_environment end
function sample_planet_frame end

end # module CallbackServices
