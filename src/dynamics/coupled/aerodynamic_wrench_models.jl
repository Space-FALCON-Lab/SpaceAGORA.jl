# include(joinpath(@__DIR__, "..", "..", "mission", "campaigns", "montecarlo_perturbations.jl"))

using SpecialFunctions
# using .SimulationModel

const sqrt_π = sqrt(π)
const inv_sqrt_π = 1 / sqrt(π)

@inline _parse_bool_env(name::String, default::Bool)::Bool = ParallelPolicy.parse_bool_env(name, default)

@inline function _multibody_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_MULTIBODY_PARALLEL")
end

@inline function _multibody_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD", 2)
end

@inline function _multibody_max_threads()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_MAX_THREADS", 4)
end

@inline function _threadid_capacity()::Int
    # Compatibility shim for legacy tests/callers after migrating scratch indexing
    # from thread ids to stable worker ids.
    return max(Threads.maxthreadid(), _multibody_max_threads())
end

@inline function _multibody_outer_parallel_hint()::Bool
    # Explicit hint for disabling intra-satellite threading under outer parallel execution.
    return ParallelPolicy.outer_parallel_active()
end

@inline function _multibody_use_threads(num_items::Int; heavy_work::Bool=true)::Bool
    return _multibody_thread_decision(num_items; heavy_work=heavy_work).use_threads
end

@inline function _multibody_thread_decision(num_items::Int; heavy_work::Bool=true)
    mode = _multibody_parallel_mode()
    threshold = _multibody_thread_threshold()
    outer_active = _multibody_outer_parallel_hint()
    allow_with_outer = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER", false)
    heavy_only = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY", true)
    policy = ParallelPolicy.thread_policy_decision(
        num_items;
        mode=mode,
        threshold=threshold,
        heavy_work=heavy_work,
        heavy_only=heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:multibody
    )
    allotment = min(policy.allotment, _multibody_max_threads())
    use_threads = policy.use_threads && allotment > 1
    return (use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode)
end

# Aerodynamic models
@kwdef struct AerodynamicCoefficientConstant <: AbstractForceTorqueModel

end

"""
    AerodynamicCoefficientfM()

Legacy Hart rectangular-prism free-molecular effector retained for backward
compatibility. It derives geometry from each spacecraft link's `dims` and does
not use the shared surface representation. New work should use
[`AerodynamicSurfaceModel`](@ref) with `flow_regime=:free_molecular` or
`:automatic`.
"""
@kwdef struct AerodynamicCoefficientfM <: AbstractForceTorqueModel

end

"""
    AerodynamicSurfaceModel(;
        geometry,
        flow_regime=:continuum,
        pressure_model=:regular_newtonian,
        fixed_alpha_rad=0.0,
        fixed_beta_rad=0.0,
    )

Unified surface aerodynamic effector. Construct a closed shared `geometry`
with [`sphere_cone_aerodynamic_geometry`](@ref) or
[`combine_aerodynamic_surfaces`](@ref), then select `:continuum`,
`:free_molecular`, `:transitional`, or `:automatic` flow. The default remains
Grant-Braun regular Newtonian continuum behavior for compatibility.

Automatic and transitional operation compute Knudsen number from the sampled
density and temperature, the planet's (or overridden) dynamic viscosity, and
the geometry reference length unless a separate characteristic length or
Knudsen override is supplied. Any regime with a nonzero free-molecular weight
requires an explicit positive `wall_temperature_k`.
"""
@kwdef struct AerodynamicSurfaceModel <: AbstractForceTorqueModel
    geometry::Union{Nothing, AerodynamicGeometry} = nothing
    flow_regime::Symbol = :continuum
    pressure_model::Symbol = :regular_newtonian
    fixed_alpha_rad::Float64 = 0.0
    fixed_beta_rad::Float64 = 0.0
    wall_temperature_k::Union{Nothing, Float64} = nothing
    normal_accommodation::Float64 = 1.0
    tangential_accommodation::Float64 = 1.0
    dynamic_viscosity_pa_s::Union{Nothing, Float64} = nothing
    knudsen_characteristic_length_m::Union{Nothing, Float64} = nothing
    knudsen_number_override::Union{Nothing, Float64} = nothing
    continuum_knudsen_limit::Float64 = 1.0e-3
    free_molecular_knudsen_limit::Float64 = 10.0
end

# Preserve the historical constructor and dispatch surface while giving the
# regime-neutral implementation an accurate public name.
const AerodynamicCoefficientNoBallisticFlight = AerodynamicSurfaceModel

@inline function _make_aero_scratch_workspace(n_threads::Int)::AeroScratchWorkspace
    n_threads >= 1 || throw(ArgumentError("Aerodynamic scratch workspace thread count must be >= 1, got $n_threads"))
    thread_force = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_threads]
    thread_cl = zeros(Float64, n_threads)
    thread_cd = zeros(Float64, n_threads)
    thread_area = zeros(Float64, n_threads)
    return AeroScratchWorkspace(thread_force, thread_cl, thread_cd, thread_area)
end

@inline function _ensure_aero_workspace_capacity!(
    workspace::AeroScratchWorkspace,
    n_threads::Int
)::AeroScratchWorkspace
    n_threads >= 1 || throw(ArgumentError("Aerodynamic scratch workspace thread count must be >= 1, got $n_threads"))
    if length(workspace.thread_force) < n_threads
        old_len = length(workspace.thread_force)
        resize!(workspace.thread_force, n_threads)
        @inbounds for tid in (old_len + 1):n_threads
            workspace.thread_force[tid] = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
    end
    if length(workspace.thread_cl) < n_threads
        resize!(workspace.thread_cl, n_threads)
        resize!(workspace.thread_cd, n_threads)
        resize!(workspace.thread_area, n_threads)
    end
    return workspace
end

@inline function _aero_workspace_for_sat!(
    param::ODEParams,
    sat_idx::Int,
    n_threads::Int
)::AeroScratchWorkspace
    workspaces = param.shared_buffers.aero_workspaces
    if sat_idx > length(workspaces)
        return _ensure_aero_workspace_capacity!(_make_aero_scratch_workspace(n_threads), n_threads)
    end

    sat_entry = @inbounds workspaces[sat_idx]
    workspace = if sat_entry === nothing
        created = _make_aero_scratch_workspace(n_threads)
        @inbounds workspaces[sat_idx] = created
        created
    else
        sat_entry::AeroScratchWorkspace
    end
    _ensure_aero_workspace_capacity!(workspace, n_threads)
    return workspace
end

@inline function _simulation_model_module_for_aero()
    mod = @__MODULE__
    while true
        if isdefined(mod, :planet_frame_lpi)
            return mod
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("SimulationModel.planet_frame_lpi not found in module ancestry for aerodynamic_wrench_models.jl")
end

@inline function _store_vector_cache!(
    cache::Vector{SVector{3, Float64}},
    sat_idx::Int,
    value::SVector{3, Float64}
)::Nothing
    if length(cache) < sat_idx
        resize!(cache, sat_idx)
        @inbounds for idx in eachindex(cache)
            if !isassigned(cache, idx)
                cache[idx] = SVector{3, Float64}(0.0, 0.0, 0.0)
            end
        end
    end
    @inbounds cache[sat_idx] = value
    return nothing
end

@inline function _store_aero_caches!(
    param::ODEParams,
    sat_idx::Int,
    drag_ii::SVector{3, Float64},
    lift_ii::SVector{3, Float64},
    cross_ii::SVector{3, Float64}
)::Nothing
    _store_vector_cache!(param.save_cache.drag_cache, sat_idx, drag_ii)
    _store_vector_cache!(param.save_cache.lift_cache, sat_idx, lift_ii)
    _store_vector_cache!(param.save_cache.cross_cache, sat_idx, cross_ii)
    return nothing
end

function collect_and_reset_link_wrenches!(bodies)
    # Collect into fresh vectors to avoid aliasing when there is only one link.
    force_acc = MVector{3, Float64}(0.0, 0.0, 0.0)
    torque_acc = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for b in bodies
        force_acc .+= SVector{3, Float64}(b.net_force)
        torque_acc .+= SVector{3, Float64}(b.net_torque)
        b.net_force .= SVector{3, Float64}(0.0, 0.0, 0.0)
        b.net_torque .= SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return SVector{3, Float64}(force_acc), SVector{3, Float64}(torque_acc)
end

@inline environment_requirements(::AerodynamicCoefficientConstant) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientfM) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientNoBallisticFlight) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline solver_partition(::AerodynamicCoefficientConstant) = :implicit
@inline solver_partition(::AerodynamicCoefficientfM) = :implicit
@inline solver_partition(::AerodynamicCoefficientNoBallisticFlight) = :implicit

@inline function _constant_drag_coefficient(alpha_rad::Float64)::Float64
    return 2 * (2.2 - 0.8) / pi * alpha_rad + 0.8
end

@inline function _aero_link_angles(
    spacecraft,
    body,
    root_index::Int,
    orientation_sim::Bool,
    vel_pi::SVector{3, Float64},
    θ_body::Float64,
)::Tuple{Float64, Float64, Union{Nothing, SMatrix{3, 3, Float64, 9}}}
    if orientation_sim
        R = rotate_to_inertial(spacecraft, body, root_index)
        body_frame_velocity = R' * vel_pi
        α_body = atan(body_frame_velocity[1], body_frame_velocity[3])
        β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3]))
        return α_body, β_body, SMatrix{3, 3, Float64, 9}(R)
    end
    α_body = body.root ? (pi / 2) : atan((rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0))[1], (rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0))[3])
    return α_body, 0.0, nothing
end

function _aero_pure_wrench(
    coefficient_mode::Symbol,
    x::StateSample,
    env::EnvironmentSample,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    x.spacecraft === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires StateSample.spacecraft."))
    planet_frame = env.planet_frame
    atmosphere = env.atmosphere
    planet_frame === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires env.planet_frame."))
    atmosphere === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires env.atmosphere."))

    rho = atmosphere.rho_kg_m3
    T = atmosphere.temperature_k
    wind = atmosphere.wind_pp
    if !isfinite(rho) || rho <= eps(Float64) || !isfinite(T) || T <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    spacecraft = x.spacecraft
    bodies = spacecraft.links
    root_index = 1
    orientation_sim = x.q_ib !== nothing
    planet = env.planet

    sound_velocity = sqrt(planet.γ * planet.R * T)
    vel_pp = planet_frame.vel_pp
    vel_pp_mag = norm(vel_pp)
    mach = vel_pp_mag / sound_velocity
    S = sqrt(planet.γ * 0.5) * mach
    h_pp = cross(planet_frame.pos_pp, vel_pp)
    h_pp_mag = norm(h_pp)
    if !isfinite(h_pp_mag) || h_pp_mag <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    bank_angle = deg2rad(0.0)
    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = vel_pp + wind_pp
    vel_pp_rw_mag = norm(vel_pp_rw)
    if vel_pp_rw_mag <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag
    h_pp_hat = h_pp / h_pp_mag
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    drag_pp_hat = -vel_pp_rw_hat
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat)
    q = 0.5 * rho * vel_pp_mag^2
    l_pi_t = planet_frame.l_pi'
    vel_pi = orientation_sim ? l_pi_t * vel_pp_rw : SVector{3, Float64}(0.0, 0.0, 0.0)
    lift_scale = q * cos(bank_angle)
    θ_body = acos(clamp(vel_pp_rw[1] / vel_pp_rw_mag, -1.0, 1.0))

    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for body in bodies
        α_body, β_body, R_body_to_inertial = _aero_link_angles(spacecraft, body, root_index, orientation_sim, vel_pi, θ_body)
        CL_body, CD_body, CS_body = if coefficient_mode == :fm
            coeffs = aerodynamic_coefficient_fM(body, T, S, α_body, β_body, θ_body)
            coeffs[1], coeffs[2], coeffs[3]
        else
            0.0, _constant_drag_coefficient(α_body), 0.0
        end

        drag_pp_body = q * CD_body * body.ref_area * drag_pp_hat
        lift_pp_body = lift_scale * CL_body * body.ref_area * lift_pp_hat
        cross_pp_body = orientation_sim ? (q * CS_body * body.ref_area * cross_pp_hat) : SVector{3, Float64}(0.0, 0.0, 0.0)
        force_body = l_pi_t * (drag_pp_body + lift_pp_body + cross_pp_body)
        force_ii .+= force_body
        if orientation_sim && !(R_body_to_inertial === nothing)
            torque_body .+= cross(SVector{3, Float64}(body.r), R_body_to_inertial' * force_body)
        end
    end

    return SVector{3, Float64}(force_ii), SVector{3, Float64}(torque_body)
end

function _surface_aerodynamic_pure_wrench(
    model::AerodynamicCoefficientNoBallisticFlight,
    x::StateSample,
    env::EnvironmentSample,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet_frame = env.planet_frame
    atmosphere = env.atmosphere
    planet = env.planet
    planet_frame === nothing && throw(ArgumentError(
        "Surface aerodynamic wrench evaluation requires env.planet_frame.",
    ))
    atmosphere === nothing && throw(ArgumentError(
        "Surface aerodynamic wrench evaluation requires env.atmosphere.",
    ))
    planet === nothing && throw(ArgumentError(
        "Surface aerodynamic wrench evaluation requires env.planet.",
    ))

    rho = atmosphere.rho_kg_m3
    temperature = atmosphere.temperature_k
    if !isfinite(rho) || rho <= 0.0 ||
       !isfinite(temperature) || temperature <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = atmosphere.wind_pp
    wind_pp = wN * uN + wE * uE - wU * uD
    relative_velocity_pp = planet_frame.vel_pp + wind_pp
    relative_speed = norm(relative_velocity_pp)
    if !isfinite(relative_speed) || relative_speed <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    geometry = model.geometry
    geometry === nothing && throw(ArgumentError(
        "AerodynamicCoefficientNoBallisticFlight requires an aerodynamic geometry; " *
        "construct one with `sphere_cone_aerodynamic_geometry` or `combine_aerodynamic_surfaces`.",
    ))
    isfinite(model.fixed_alpha_rad) || throw(DomainError(
        model.fixed_alpha_rad,
        "fixed angle of attack must be finite.",
    ))
    isfinite(model.fixed_beta_rad) || throw(DomainError(
        model.fixed_beta_rad,
        "fixed sideslip must be finite.",
    ))

    gamma = Float64(planet.γ)
    gas_constant = Float64(planet.R)
    sound_speed = sqrt(gamma * gas_constant * temperature)
    mach = relative_speed / sound_speed
    speed_ratio = sqrt(0.5 * gamma) * mach
    knudsen_number = if model.knudsen_number_override !== nothing
        model.knudsen_number_override
    elseif model.flow_regime === :automatic || model.flow_regime === :transitional
        dynamic_viscosity = model.dynamic_viscosity_pa_s === nothing ?
            Float64(planet.μ_fluid) : model.dynamic_viscosity_pa_s
        characteristic_length = model.knudsen_characteristic_length_m === nothing ?
            geometry.reference_length_m : model.knudsen_characteristic_length_m
        gas_knudsen_number(
            rho,
            temperature,
            gas_constant,
            dynamic_viscosity,
            characteristic_length,
        )
    else
        nothing
    end
    dynamic_pressure = 0.5 * rho * relative_speed^2
    force_scale = dynamic_pressure * geometry.reference_area_m2
    moment_scale = force_scale * geometry.reference_length_m
    l_pi_transpose = planet_frame.l_pi'

    if x.q_ib !== nothing
        body_to_inertial = rot(x.q_ib)'
        relative_velocity_ii = l_pi_transpose * relative_velocity_pp
        freestream_body = body_to_inertial' * (-relative_velocity_ii / relative_speed)
        regime_result = aerodynamic_regime_coefficients(
            geometry,
            freestream_body;
            flow_regime=model.flow_regime,
            knudsen_number=knudsen_number,
            continuum_limit=model.continuum_knudsen_limit,
            free_molecular_limit=model.free_molecular_knudsen_limit,
            pressure_model=model.pressure_model,
            mach=mach,
            gamma=gamma,
            speed_ratio=speed_ratio,
            temperature_inf_k=temperature,
            wall_temperature_k=model.wall_temperature_k,
            normal_accommodation=model.normal_accommodation,
            tangential_accommodation=model.tangential_accommodation,
        )
        coefficients = regime_result.coefficients
        force_ii = body_to_inertial * (force_scale * coefficients.force_body)
        torque_body = moment_scale * coefficients.moment_body
        return SVector{3, Float64}(force_ii), SVector{3, Float64}(torque_body)
    end

    regime_result = aerodynamic_regime_coefficients(
        geometry,
        model.fixed_alpha_rad,
        model.fixed_beta_rad;
        flow_regime=model.flow_regime,
        knudsen_number=knudsen_number,
        continuum_limit=model.continuum_knudsen_limit,
        free_molecular_limit=model.free_molecular_knudsen_limit,
        pressure_model=model.pressure_model,
        mach=mach,
        gamma=gamma,
        speed_ratio=speed_ratio,
        temperature_inf_k=temperature,
        wall_temperature_k=model.wall_temperature_k,
        normal_accommodation=model.normal_accommodation,
        tangential_accommodation=model.tangential_accommodation,
    )
    coefficients = regime_result.coefficients
    velocity_unit_pp = relative_velocity_pp / relative_speed
    drag_unit_pp = -velocity_unit_pp
    angular_momentum_pp = cross(planet_frame.pos_pp, relative_velocity_pp)
    angular_momentum_magnitude = norm(angular_momentum_pp)
    angular_momentum_magnitude > eps(Float64) || throw(DomainError(
        angular_momentum_pp,
        "a nondegenerate planet-relative trajectory plane is required for 3-DOF aerodynamic lift.",
    ))
    lift_unit_pp = normalize(cross(angular_momentum_pp / angular_momentum_magnitude, velocity_unit_pp))
    side_unit_pp = cross(drag_unit_pp, lift_unit_pp)
    force_pp = force_scale * (
        coefficients.CD * drag_unit_pp +
        coefficients.CL * lift_unit_pp +
        coefficients.CY_wind * side_unit_pp
    )
    force_ii = l_pi_transpose * force_pp
    torque_body = moment_scale * coefficients.moment_body
    return SVector{3, Float64}(force_ii), SVector{3, Float64}(torque_body)
end

@inline function wrench(
    model::AerodynamicCoefficientConstant,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _aero_pure_wrench(:constant, x, env)
end

@inline function wrench(
    model::AerodynamicCoefficientfM,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _aero_pure_wrench(:fm, x, env)
end

@inline function wrench(
    model::AerodynamicCoefficientNoBallisticFlight,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _surface_aerodynamic_pure_wrench(model, x, env)
end

# Calculate force/torque functions
function calcForceTorque(model::AerodynamicCoefficientConstant, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    m = param.m
    cnf = param.cnf
    orientation_sim = param.orientation_sim

    bodies, root_index = traverse_bodies(m.body, m.body.roots[1]) # Get all bodies in the simulation

    pos_ii = SVector{3, Float64}(x[1:3])
    vel_ii = SVector{3, Float64}(x[4:6])

    h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

    h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]

    # Inertial to planet relative transformation
    pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, cnf.et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
    pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

    vel_pp_mag = norm(vel_pp)

    h_pp = cross(pos_pp, vel_pp)
    
    h_pp_mag = norm(h_pp)
    h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
    
    bank_angle = deg2rad(0.0)

    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
    drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame

    if orientation_sim
        Rot = [MMatrix{3,3,Float64}(zeros(3, 3)) for i in eachindex(bodies)] # Rotation matrix from the root body to the spacecraft link
        @inbounds for (i, b) in enumerate(bodies)
            Rot[i] .= rotate_to_inertial(m.body, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
        end
    end
    
    CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
    total_area = 0.0 # Initialize total area

    α = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize angle of attack vector
    β = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize sideslip angle vector
    R = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from the root body to the spacecraft link
    # Determine angle of attack (α) and sideslip angle (β)
    # Vehicle Aerodynamic Forces
    # CL and CD
    @inbounds for (i, b) in enumerate(bodies)
        if orientation_sim
            R .= Rot[i] # Rotation matrix from the spacecraft link to the inertial frame
            body_frame_velocity = R' * m.planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
            
            α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
            β_body = atan(body_frame_velocity[2], norm([body_frame_velocity[1], body_frame_velocity[3]])) # Sideslip angle in radians
            α[i] = α_body # Angle of attack for the spacecraft link
            β[i] = β_body # Sideslip angle for the spacecraft link
            b.α = α_body
            b.β = β_body
            b.θ = acos(clamp(vel_pp_rw[1]/norm(vel_pp_rw), -1.0, 1.0)) # Elevation angle for the spacecraft link
        else
            # TODO: Change this so that it just uses above code even with orientation_sim = false
            if b.root
                # if the body is the root body, then the angle of attack is 90 degrees
                α[i] = pi/2
                b.α = pi/2 # Angle of attack for the root body
            else
                body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                α[i] = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                # α[i] = pi/2 # Angle of attack for the spacecraft link, temporary hard code for testing
                b.α = α[i] # Angle of attack for the spacecraft link
            end
        end

        CL_body = 0.0
        CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8

        # if montecarlo == true
        #     CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
        # end

        drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = q * CL_body * b.ref_area * lift_pp_hat * cos(bank_angle)     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
            cross_body = m.planet.L_PI' * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Planet relative cross force vector
            cross_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Inertial cross force vector
        end

        drag_body = m.planet.L_PI' * drag_pp_body   # Inertial drag force vector
        lift_body = m.planet.L_PI' * lift_pp_body   # Inertial lift force vector

        # Update the force on the spacecraft link
        b.net_force .+= drag_body + lift_body + cross_body # Update the force on the spacecraft link, inertial frame
        # aero_torque = q * Cl_body * b.ref_area * b.dims[1] * SVector{3, Float64}(1.0, 0.0, 0.0) + # Aerodynamic roll torque, body frame
        #               q * Cm_body * b.ref_area * b.dims[2] * SVector{3, Float64}(0.0, 1.0, 0.0) + # Aerodynamic pitch torque, body frame
        #               q * Cn_body * b.ref_area * b.dims[3] * SVector{3, Float64}(0.0, 0.0, 1.0)   # Aerodynamic yaw torque, body frame
        # b.net_torque .+= aero_torque # Update the torque on the spacecraft link, body frame
        b.net_torque .+= cross(b.r, rot_body_to_inertial' * (drag_body + lift_body + cross_body)) # Update the torque on the spacecraft link, body frame
        # Update the total CL/CD
        CL += CL_body * b.ref_area
        CD += CD_body * b.ref_area
        total_area += b.ref_area # Update the total area
        drag_ii += drag_body # Update the total drag force
        lift_ii += lift_body # Update the total lift force
        drag_pp += drag_pp_body # Update the total drag force in planet relative frame
        lift_pp += lift_pp_body # Update the total lift force in planet relative frame
    end
    
    # Normalize the aerodynamic coefficients
    CL = CL / total_area
    CD = CD / total_area

    force_ii, torque_ii = collect_and_reset_link_wrenches!(bodies)

    return force_ii, torque_ii
end

function calcForceTorque(model::AerodynamicCoefficientfM, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # m = param.m
    # cnf = param.cnf
    # orientation_sim = param.orientation_sim
    planet = param.args.environment_model.planet
    ephemerides_model = param.args.environment_model.ephemerides_model
    orientation_sim = param.args.mission_configuration.orientation_sim
    # bodies, root_index = traverse_bodies(m.body, m.body.roots[1]) # Get all bodies in the simulation
    spacecraft = param.args.dynamics_model.spacecraft[i]
    bodies = spacecraft.links # Include the root body of the spacecraft
    root_index = 1
    env_state = SimulationModel.SimulationCallbacks._stage_environment_state(
        x,
        param,
        i,
        param.shared_buffers.current_time[];
        write_buffers=true,
    )
    rho = env_state.rho
    T = env_state.T
    wind = env_state.wind

    # Skip expensive aerodynamic geometry when the flow is effectively vacuum.
    if !isfinite(rho) || rho <= eps(Float64)
        _store_aero_caches!(param, i, SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    sound_velocity = sqrt(planet.γ * planet.R * T) # Speed of sound at the current altitude based on the temperature from the shared buffer updated by the callback function
    pos_ii = env_state.pos_ii
    vel_ii = env_state.vel_ii

    h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

    h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]

    # Inertial to planet relative transformation
    # println("pos_ii: ", pos_ii, " vel_ii: ", vel_ii, " alt: ", norm(pos_ii) - planet.Rp_e, " m, vel_mag: ", norm(vel_ii), " m/s")
    # println("planet.L_PI: ", planet.L_PI)
    pos_pp = env_state.pos_pp
    vel_pp = env_state.vel_pp
    pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

    vel_pp_mag = norm(vel_pp)
    mach = vel_pp_mag / sound_velocity # Mach number of the flow at the current altitude
    S = sqrt(planet.γ * 0.5) * mach # Molecular speed ratio
    h_pp = cross(pos_pp, vel_pp)
    
    h_pp_mag = norm(h_pp)
    h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
    
    bank_angle = deg2rad(0.0)
    alt,lat,lon = env_state.alt, env_state.lat, env_state.lon

    uD, uN, uE = latlongtoNED((alt, lat, lon))
    wE, wN, wU = wind # positive to the east , m / s
    wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s 
    vel_pp_rw = vel_pp + wind_pp                  # relative wind vector, m / s
    vel_pp_rw_mag = norm(vel_pp_rw)
    if vel_pp_rw_mag <= eps(Float64)
        _store_aero_caches!(param, i, SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    mach = vel_pp_rw_mag / sound_velocity # Mach number of the flow at the current altitude
    S = sqrt(planet.γ * 0.5) * mach # Molecular speed ratio
    vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag     # relative wind unit vector
    vel_pp_rw_inv_mag = inv(vel_pp_rw_mag)
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
    drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame
    q = 0.5 * rho * vel_pp_rw_mag^2 # Dynamic pressure in planet relative frame, using the stage-consistent density value
    L_PI_t = env_state.l_pi'
    vel_pi = orientation_sim ? L_PI_t * vel_pp_rw : SVector{3, Float64}(0.0, 0.0, 0.0)
    lift_scale = q * cos(bank_angle)

    CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
    total_area = 0.0 # Initialize total area
    θ_body = acos(clamp(vel_pp_rw[1] * vel_pp_rw_inv_mag, -1.0, 1.0))
    decision = _multibody_thread_decision(length(bodies); heavy_work=true)
    use_threads = decision.use_threads
    n_workers = use_threads ? ParallelPolicy.thread_worker_count(length(bodies), decision.allotment) : 1

    zero_vec = SVector{3, Float64}(0.0, 0.0, 0.0)

    function compute_link_wrench!(idx::Int)
        b = bodies[idx]
        if orientation_sim
            R = rotate_to_inertial(spacecraft, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
            body_frame_velocity = R' * vel_pi # Velocity of the spacecraft link in inertial frame

            α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
            β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3])) # Sideslip angle in radians
            b.α = α_body
            b.β = β_body
            b.θ = θ_body
        else
            # TODO: Change this so that it just uses above code even with orientation_sim = false
            if b.root
                # if the body is the root body, then the angle of attack is 90 degrees
                b.α = pi / 2 # Angle of attack for the root body
            else
                body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                b.α = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                # b.α = pi/2 # Angle of attack for the spacecraft link, temporary hard code for testing
            end
        end

        # CL_body, CD_body = aerodynamic_coefficient_fM(α[i], m.body, T_p, S, m.aerodynamics, MonteCarlo)
        CL_body, CD_body, CS_body, _, _, _ = aerodynamic_coefficient_fM(b, T, S)

        drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = lift_scale * CL_body * b.ref_area * lift_pp_hat     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
            cross_body = L_PI_t * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = zero_vec # Planet relative cross force vector
            cross_body = zero_vec # Inertial cross force vector
        end

        drag_body = L_PI_t * drag_pp_body   # Inertial drag force vector
        lift_body = L_PI_t * lift_pp_body   # Inertial lift force vector
        force_body = drag_body + lift_body + cross_body

        return force_body, drag_body, lift_body, cross_body, CL_body * b.ref_area, CD_body * b.ref_area, b.ref_area
    end

    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    drag_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    lift_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    cross_ii = MVector{3, Float64}(0.0, 0.0, 0.0)

    # Determine angle of attack (α) and sideslip angle (β)
    # Vehicle Aerodynamic Forces
    # CL and CD
    started_ns = time_ns()
    if use_threads
        workspace = _aero_workspace_for_sat!(param, i, n_workers)
        thread_force = workspace.thread_force
        thread_drag = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_lift = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_cross = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_cl = workspace.thread_cl
        thread_cd = workspace.thread_cd
        thread_area = workspace.thread_area
        @inbounds for worker_id in 1:n_workers
            thread_force[worker_id] .= 0.0
            thread_drag[worker_id] .= 0.0
            thread_lift[worker_id] .= 0.0
            thread_cross[worker_id] .= 0.0
            thread_cl[worker_id] = 0.0
            thread_cd[worker_id] = 0.0
            thread_area[worker_id] = 0.0
        end

        ParallelPolicy.threaded_foreach_worker(length(bodies), decision.allotment) do worker_id, idx
            force_body, drag_body, lift_body, cross_body, cl_area, cd_area, area = compute_link_wrench!(idx)
            thread_force[worker_id] .+= force_body
            thread_drag[worker_id] .+= drag_body
            thread_lift[worker_id] .+= lift_body
            thread_cross[worker_id] .+= cross_body
            thread_cl[worker_id] += cl_area
            thread_cd[worker_id] += cd_area
            thread_area[worker_id] += area
        end

        @inbounds for worker_id in 1:n_workers
            force_ii .+= thread_force[worker_id]
            drag_ii .+= thread_drag[worker_id]
            lift_ii .+= thread_lift[worker_id]
            cross_ii .+= thread_cross[worker_id]
            CL += thread_cl[worker_id]
            CD += thread_cd[worker_id]
            total_area += thread_area[worker_id]
        end
    else
        @inbounds for idx in eachindex(bodies)
            force_body, drag_body, lift_body, cross_body, cl_area, cd_area, area = compute_link_wrench!(idx)
            force_ii .+= force_body
            drag_ii .+= drag_body
            lift_ii .+= lift_body
            cross_ii .+= cross_body
            CL += cl_area
            CD += cd_area
            total_area += area
        end
    end
    ParallelPolicy.record_policy_observation!(
        :multibody;
        mode=decision.mode,
        num_items=length(bodies),
        use_threads=use_threads,
        elapsed_ns=(time_ns() - started_ns)
    )

    # cnf.drag_pp = drag_pp
    # cnf.lift_pp = lift_pp

    # cnf.drag_ii = drag_ii
    # cnf.lift_ii = lift_ii
    
    # Normalize the aerodynamic coefficients
    if total_area > 0.0
        CL = CL / total_area
        CD = CD / total_area
    end

    # cnf.CL_current = CL
    # cnf.CD_current = CD

    # cnf.β_body = β
    # cnf.α_body = α

    _store_aero_caches!(param, i, SVector{3, Float64}(drag_ii), SVector{3, Float64}(lift_ii), SVector{3, Float64}(cross_ii))
    # println("Force_ii: ", force_ii, " Drag_ii: ", drag_ii, " vel_ii: ", vel_ii, " CL: ", CL, " CD: ", CD, " Alt: ", alt, " Lat: ", lat, " Lon: ", lon)
    return SVector{3, Float64}(force_ii), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function aerodynamic_coefficient_constant(α, body, T, S, args, montecarlo=false)
    """

    """

    CL_body = 0.0
    CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

# function aerodynamic_coefficient_fM(α, body, T, S, args, montecarlo=false)
#     """

#     """

#     σ = args.reflection_coefficient
#     Tw = T

#     function pressure(S, α, ρ_inf, vel, σ)
#         """

#         """

#         p = (ρ_inf*vel^2) / (2*S^2) * ((((2 - σ) / sqrt(pi))*S*sin(α) + sqrt(T/Tw)*σ/2) * exp(-(S*sin(α))^2) + ((2-σ)*((S*sin(α))^2 + 0.5) + (σ/2)*sqrt(pi)*(S*sin(α))) * (1 + erf(S*sin(α))))

#         return p
#     end

#     function τ(S, α, ρ_inf, vel, σ)
#         """

#         """

#         t = ((σ*cos(α)*ρ_inf*vel^2) / (sqrt(pi)*2*S)) * (exp(-(S*sin(α))^2) + sqrt(pi)*(S*sin(α)) * (1 + erf(S*sin(α))))

#         return t
#     end

#     function normalcoefficient(S, aoa, sigma)
#         CN = 1 / (S^2) * ((((2 - sigma) / sqrt(pi)) * S * sin(aoa) + sigma / 2) * exp(-(S * sin(aoa))^2) + ((2 - sigma) * ((S * sin(aoa))^2 + 0.5) + sigma / 2 * sqrt(pi) * (S * sin(aoa))) * (1 + erf(S * sin(aoa))))
#         return CN
#     end

#     function axialcoefficient(S, aoa, sigma)
#         CA = ((sigma * cos(aoa)) / (sqrt(pi) * S)) * (exp(-(S * sin(aoa))^2) + sqrt(pi) * (S * sin(aoa)) * (1 + erf(S * sin(aoa))))
#         return CA
#     end

#     # println("α: ", α)

#     # Solar Panels
#     CN_sa = normalcoefficient(S, α, σ)
#     CA_sa = axialcoefficient(S, α, σ)
#     CL_sa = CN_sa*cos(α) - CA_sa*sin(α)
#     CD_sa = CA_sa*cos(α) + CN_sa*sin(α)
#     # println("CL_sa: ", CL_sa)
#     # println("CD_sa: ", CD_sa)
#     # Spacecraft
#     # CN_sc = normalcoefficient(S, pi*0.5, σ)
#     # CA_sc = axialcoefficient(S, pi*0.5, σ)
#     # CL_sc = CN_sc*cos(pi*0.5) - CA_sc*sin(pi*0.5)
#     # CD_sc = CA_sc*cos(pi*0.5) + CN_sc*sin(pi*0.5)

#     # area_SA = config.get_SA_area(body, body.roots[1])
#     # area_SC = config.get_SC_area(body, body.roots[1])
    
#     # CD_body = (CD_sa*area_SA + CD_sc*area_SC) / (area_SA + area_SC)
#     # CL_body = (CL_sa*area_SA + CL_sc*area_SC) / (area_SA + area_SC)
#     if montecarlo == true
#         CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
#     end

#     return CL_sa, CD_sa
# end

function aerodynamic_coefficient_fM(body, T::Float64, S::Float64)
    return aerodynamic_coefficient_fM(body, T, S, body.α, body.β, body.θ)
end

function aerodynamic_coefficient_fM(
    body,
    T::Float64,
    S::Float64,
    α_input::Float64,
    β_input::Float64,
    θ_input::Float64,
)
    """
    Calculate the aerodynamic coefficients for a blunt body in ballistic flight using the F.M. model.
        Angle of attack and side slip are calculated as angles between CA (normal to flat plate), and wind-relative velocity vector.
    # Arguments
    - `body`: Body object containing dimensions and properties
    - `T`: Temperature, K
    - `S`: Molecular speed ratio
    - `args`: Dictionary containing additional parameters like reflection coefficient
    - `montecarlo`: Boolean flag to apply Monte Carlo perturbations
    # Returns
    - `CL`: Lift coefficient
    - `CD`: Drag coefficient
    - `CS`: Sideslip coefficient
    # Notes
    - Equations are from Hart et al. (2017, https://doi.org/10.2514/1.A33606), for a rectangular prism. 
    """

    α = α_input
    β = β_input
    θ = θ_input
    # Adjust angles to match model assumptions
    α -= π/2
    σ = body.reflection_coefficient
    Tw = T
    lx = body.dims[1]
    ly = body.dims[2]
    lz = body.dims[3]
    lx_over_ly = lx / ly
    lx_over_lz = lx / lz
    σN = σ
    σT = σ
    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)
    # println("cosα: ", cosα)
    # println("cosβ: ", cosβ)
    # println("sinα: ", sinα)
    # println("sinβ: ", sinβ)
    # println("S: ", S)
    # println("T: ", T)
    # println("sqrt_π: ", sqrt_π)
    # println("inv_sqrt_π: ", inv_sqrt_π)
    # println("σN: ", σN)
    # println("lx: ", lx, " ly: ", ly, " lz: ", lz)
    # println("α (deg): ", rad2deg(α))
    # println("β (deg): ", rad2deg(β))
    # println("θ (deg): ", rad2deg(θ))
    # Calculate the aerodynamic coefficients in the body frame (flat plate)
    # Axial
    CA = ((2-σN)/(S*sqrt_π)*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
            (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
            σT*cosα*cosβ*(lx_over_ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
            lx_over_lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Crossflow
    CS = lx_over_ly*(((2-σN)/(S*sqrt_π)*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
                (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
                (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
                σT*sinβ*(1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
                lx_over_lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Normal
    CN = lx_over_lz*((((2-σN)/(S*sqrt_π)*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
                (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
                (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
                σT*sinα*cosβ*(lx_over_ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
                (1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))

    # println("CA: ", CA)
    # println("CS: ", CS)
    # println("CN: ", CN)
    # Calculate the aerodynamic coefficients in the body frame (box)
    # Axial
    # CA = calculate_CA(lx, ly, lz, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T, θ)
    # # Crossflow
    # CS = calculate_CS(lx, ly, lz, θ, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T)
    # # Normal
    # CN = calculate_CN(lx, ly, lz, θ, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T)

    # # Calculate moment coefficients
    # Cl = calculate_Cl(σT, S, sinα, cosα, sinβ, cosβ, θ)
    # Cm = calculate_Cm(σT, S, sinα, cosα, sinβ, cosβ, θ)
    # Cn = calculate_Cn(σT, S, sinα, cosα, sinβ, cosβ, θ)

    # Calculate the aerodynamic coefficients in the wind frame
    CL = -sinα*CA + cosα*CN
    CD = cosα*cosβ*CA + sinβ*CS + sinα*cosβ*CN
    CS = -cosα*sinβ*CA + cosβ*CS - sinα*sinβ*CN

    # If doing Monte Carlo simulations, apply perturbations to the coefficients
    # if montecarlo == true
    #     CL, CD = monte_carlo_aerodynamics(CL, CD, args)
    # end

    return CL, CD, CS, 0.0, 0.0, 0.0 #, Cl, Cm, Cn
    # return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

"""
    aerodynamic_coefficient_no_ballistic_flight(
        alpha, body, args, T=0, S=0, aero=0, montecarlo=false;
        pressure_model=:regular_newtonian, mach=nothing, gamma=nothing,
        beta=0, geometry=nothing,
    ) -> (CL, CD)

Compatibility wrapper for the zero-sideslip Grant-Braun sphere-cone model.
`body.nose_radius`, `body.base_radius`, and `body.δ` are required, with `body.δ`
in radians. The full surface evaluator supports sideslip and shadowing; pass a
prebuilt `geometry` to avoid reconstructing quadrature in repeated calls.
Existing positional calls use regular Newtonian theory. Select
`:modified_newtonian` explicitly and provide `mach` and `gamma` to apply the
stagnation-pressure correction.
"""
function aerodynamic_coefficient_no_ballistic_flight(
    α,
    body,
    args,
    T=0,
    S=0,
    a=0,
    montecarlo=false;
    pressure_model::Symbol=:regular_newtonian,
    mach::Union{Nothing, Real}=nothing,
    gamma::Union{Nothing, Real}=nothing,
    beta::Real=0.0,
    geometry::Union{Nothing, NewtonianAerodynamicGeometry}=nothing,
)
    newtonian_geometry = if geometry !== nothing
        geometry
    elseif hasproperty(body, :newtonian_geometry)
        getproperty(body, :newtonian_geometry)
    else
        sphere_cone_newtonian_geometry(
            body.nose_radius,
            body.base_radius,
            body.δ,
        )
    end
    coefficients = newtonian_aerodynamic_coefficients(
        newtonian_geometry,
        α,
        beta;
        pressure_model=pressure_model,
        mach=mach,
        gamma=gamma,
    )
    CL_body = coefficients.CL
    CD_body = coefficients.CD

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

"""
    calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)

Calculates the expression C_A from the provided image.
This version incorporates the clarification that θ is a scalar variable (an angle).

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- cos_α, sin_α, cos_β, sin_β: trig functions of angle of attack and sideslip angles (angles)
- θ: theta (flow angle)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2 = t1 * cos_α * cos_β * σ_T * θ
    
    term1 = (common_factor_1_2 * (-sin_β) * paren_term_minus_1) / den1_2
    term2_frac = (common_factor_1_2 * sin_β * paren_term_plus_1) / den1_2

    # --- Group 2: Terms 3 and 4 (no denominators) ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    s_sq_ca_sq_cb_sq = s_sq * cos_α^2 * cos_β^2
    
    # Term 3
    exp_bracket3 = (s_ca_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket3 = (1 - erf_s_cacb) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) - 
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) + 
                    exp_s_cacb_sq * exp_bracket3
    term3 = 1/s_sq * θ * (-cos_α) * cos_β * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_ca_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_cacb + 1) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) +
                    exp_s_cacb_sq * exp_bracket4
    term4 = 1/s_sq * θ * cos_α * cos_β * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6 = t1 * cos_α * cos_β * σ_T * θ
    
    term5 = (common_factor_5_6 * (-cos_β * sin_α) * paren_term5) / den5_6
    term6 = (common_factor_5_6 * (cos_β * sin_α) * paren_term6) / den5_6

    # --- Final Combination ---
    # Based on the signs shown in the image:
    # T1 + T2_frac - 1/s² - T3 + T4 + T5 + T6
    Ca_result = term1 + term2_frac - term3 + term4 + term5 + term6
    
    return Ca_result
end

"""
    calculate_CS(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_S from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CS(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)
    s_sq_sin_β_sq = s_sq * sin_β^2

    # Term 1
    exp_bracket1 = (s_sin_β * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket1 = (1 - erf_s_sin_β) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) -
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket1
    num1 = t1 * θ * (-sin_β) * main_bracket1
    den1_2 = s_sq * t2
    term1 = num1 / den1_2

    # Term 2
    exp_bracket2 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sin_β * (σ_N - 2)) / sqrt_π
    main_bracket2 = (erf_s_sin_β + 1) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) +
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket2
    num2 = t1 * θ * sin_β * main_bracket2
    term2 = num2 / den1_2

    # --- Group 2: Terms 3 and 4 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    
    # Term 3
    paren_term3 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    num3 = sin_β * σ_T * θ * (-cos_α * cos_β) * paren_term3
    den3_4 = sqrt_π * s
    term3 = num3 / den3_4

    # Term 4
    paren_term4 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    num4 = sin_β * σ_T * θ * (cos_α * cos_β) * paren_term4
    term4 = num4 / den3_4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6_num = t1 * sin_β * σ_T * θ
    
    num5 = common_factor_5_6_num * (-cos_β * sin_α) * paren_term5
    term5 = num5 / den5_6

    num6 = common_factor_5_6_num * (cos_β * sin_α) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown at the start and end of each line in the image
    Cs_result = -term1 + term2 + term3 + term4 + term5 + term6
    
    return Cs_result
end

"""
    calculate_CN(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_N from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CN(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2_num = t1 * sin_α * cos_β * σ_T * θ
    
    num1 = common_factor_1_2_num * (-sin_β) * paren_term_minus_1
    term1 = num1 / den1_2

    num2_frac = common_factor_1_2_num * sin_β * paren_term_plus_1
    term2_frac = num2_frac / den1_2
    term2_sub = 1 / (s_sq * t3)

    # --- Group 2: Terms 3 and 4 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)
    s_sq_sa_sq_cb_sq = s_sq * sin_α^2 * cos_β^2

    # Term 3
    exp_bracket3 = (s_sa_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf + 1 / (s_sq * t3)
    main_bracket3 = (1 - erf_s_sacb) * ((2 - σ_N) * s_sq_sa_sq_cb_sq - 0.5) -
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket3
    term3 = t1 * θ * (-cos_β * sin_α) * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sa_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_sacb + 1) * ((2 - σ_N) * s_sq_sa_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket4
    term4 = t1 * θ * (cos_β * sin_α) * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)

    paren_term5 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    paren_term6 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    
    den5_6 = sqrt_π * s
    common_factor_5_6_num = sin_α * cos_β * σ_T * θ

    num5 = common_factor_5_6_num * (-cos_α * cos_β) * paren_term5
    term5 = num5 / den5_6
    
    num6 = common_factor_5_6_num * (cos_α * cos_β) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown in the image: 
    # T1 + T2_frac - 1/(s²t₃) - T3 + T4 + T5 + T6
    Cn_result = term1 + term2_frac - term2_sub - term3 + term4 + term5 + term6
    
    return Cn_result
end

"""
    calculate_Cl(σ_T, s, α, β, θ)

Calculates the expression C_L from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cl(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Calculation for the "Top Part" (first three lines) ---
    # This part is structured as: ( sin(β) * [ ... ] ) - 1/(sqrt(π)*|cos(β)|)
    
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    # The two θ terms can be combined and simplified algebraically:
    # θ(-x)A + θ(x)B = -θx A + θx B = θx(B - A)
    paren_A = s_sa_cb * erf_s_sacb + exp_s_sacb_sq / sqrt_π - s_sa_cb
    paren_B = -s_sa_cb * erf_s_sacb - exp_s_sacb_sq / sqrt_π - s_sa_cb
    
    term_in_paren_top = paren_B - paren_A # This simplifies to -2*s_sa_cb*erf - 2*exp/sqrt_π
    
    # Note: θ(-cos(β)sin(α)) becomes -θ*cos(β)sin(α) etc.
    sub_part_1 = sin_β * (θ * cos_β * sin_α * term_in_paren_top)
    
    sub_part_2 = 1 / (sqrt_π * abs(cos_β))
    
    top_part = sub_part_1 #- sub_part_2

    # --- Calculation for the "Bottom Part" (last three lines) ---
    # This is structured as: [prefactors] * ( θ(-sinβ)[...] - θ(sinβ)[...] )
    s_sin_β = s * sin_β
    exp_s_sin_β_sq_neg = exp(-s_sin_β^2)
    
    # Note: sqrt(sec(β)²) simplifies to 1/abs(cos(β))
    prefactor_bottom = sin_α * cos_β^3 * (1 / abs(cos_β)) * exp_s_sin_β_sq_neg

    # The bracketed terms can also be simplified:
    exp_s_sin_β_sq_pos = exp(s_sin_β^2)
    erf_s_sin_β = erf(s_sin_β)
    common_in_brackets = sqrt_π * s_sin_β * exp_s_sin_β_sq_pos * erf_s_sin_β
    
    bracket_A = common_in_brackets + sqrt_π * s_sin_β * (-exp_s_sin_β_sq_neg + 1)
    bracket_B = common_in_brackets + sqrt_π * s_sin_β * (exp_s_sin_β_sq_neg + 1)
    
    # θ(-sinβ)A - θ(sinβ)B = -θ*sinβ*A - θ*sinβ*B = -θ*sinβ*(A + B)
    sum_of_brackets = bracket_A + bracket_B # Simplifies to 2*common + 2*sqrt(π)*s_sin_β

    main_paren_bottom = -θ * sin_β * sum_of_brackets

    bottom_part = prefactor_bottom * main_paren_bottom
    
    # --- Final Combination ---
    total_inside_braces = top_part - sub_part_2 * bottom_part
    
    Cl_result = (1 / (2 * s)) * σ_T * total_inside_braces
    
    return Cl_result
end

"""
    calculate_Cm(σ_T, s, α, β, θ)

Calculates the expression C_m from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cm(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Define arguments for the two main patterns ---
    arg1 = s * sin_α * cos_β
    arg2 = s * cos_α * cos_β

    # --- Term 1 ---
    erf_arg1 = erf(arg1)
    exp_arg1_sq = exp(-arg1^2)
    
    # NOTE: The structure of this parenthesis is unusual compared to the others.
    # It is translated literally as it appears in the image.
    paren1 = sin_α * cos_β - sin_α * cos_β * erf_arg1 - exp_arg1_sq / (sqrt_π * s)
    # term1 = cos_α * (θ * (-cos_β * sin_α)) * paren1

    # --- Term 2 ---
    paren2 = arg1 * erf_arg1 + exp_arg1_sq / sqrt_π + arg1
    # term2 = (1/s) * cos_α * (θ * (cos_β * sin_α)) * paren2

    # --- Term 3 ---
    erf_arg2 = erf(arg2)
    exp_arg2_sq = exp(-arg2^2)
    
    paren3 = arg2 * erf_arg2 + exp_arg2_sq / sqrt_π - arg2
    # term3 = sin_α * (θ * (-cos_α * cos_β)) * paren3

    # --- Term 4 ---
    # NOTE: This term is not prefixed by sin(α), unlike its counterpart, Term 3.
    # This asymmetry is preserved from the image.
    paren4 = -arg2 * erf_arg2 - exp_arg2_sq / sqrt_π - arg2
    term4 = (θ * (cos_α * cos_β)) * paren4

    # --- Final Combination ---
    total_sum = cos_α*θ*(-cos_β*sin_α)*paren1 + 1/s*(cos_α*θ*cos_β*cos_α*paren2 + sin_α*(θ * (-cos_α * cos_β) * paren3 + term4))
    
    Cm_result = 0.5 * cos_β * σ_T * total_sum
    
    return Cm_result
end

"""
    calculate_Cn(σ_T, s, α, β, θ)

Calculates the expression C_n from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cn(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    abs_cos_β = abs(cos_β)

    # --- Define common exponential terms ---
    exp_arg_cos_term = s_sq * cos_α^2 * cos_β^2
    exp_arg_sin_term = s_sq * sin_β^2
    
    exp_cos_pos = exp(exp_arg_cos_term)
    exp_cos_neg = exp(-exp_arg_cos_term)
    exp_sin_pos = exp(exp_arg_sin_term)
    
    exp_main_pos = exp(exp_arg_cos_term + exp_arg_sin_term)
    exp_main_neg = exp(-(exp_arg_cos_term + exp_arg_sin_term))

    # --- Outer Factor ---
    outer_factor = (1 / (2 * sqrt_π * s * abs_cos_β)) * σ_T * exp_main_neg

    # --- Part 1: Terms factored by (-sin(β)|cos(β)|) ---
    s_ca_cb = s * cos_α * cos_β
    
    # The two θ terms inside the first main parenthesis can be combined and simplified.
    # The structure is prop. to θ*cosαcosβ * (Paren_B - Paren_A)
    # (Paren_B - Paren_A) simplifies to: 2 * sqrt(π) * s_ca_cb * exp(-s²cos²αcos²β)
    diff_paren_AB = 2 * sqrt_π * s_ca_cb * exp_cos_neg
    combo_AB = θ * cos_α * cos_β * diff_paren_AB
    
    part1 = (-sin_β * abs_cos_β) * combo_AB
    
    # --- Part 2: Terms factored by (cos³(β)sqrt(sec²(β))) ---
    s_sin_β = s * sin_β
    erf_s_sinb = erf(s_sin_β)

    # The two θ terms inside the second main parenthesis can also be simplified.
    # The structure is prop. to -θ*sinβ * (Paren_C + Paren_D)
    # (Paren_C + Paren_D) simplifies to: 
    # 2 * (sqrt(π)*s*sinβ*erf(s*sinβ)*exp_main_pos + exp(s²cos²αcos²β))
    sum_paren_CD = 2 * (sqrt_π * s_sin_β * erf_s_sinb * exp_main_pos + exp_cos_pos)
    combo_CD = θ * (-sin_β) * sum_paren_CD
    
    # Prefactor for Part 2 simplifies: cos³(β)sqrt(sec²(β)) = cos³(β)/|cos(β)|
    prefactor2 = cos_β^3 / abs_cos_β
    part2 = prefactor2 * combo_CD

    # --- Final Combination ---
    total_sum_in_brackets = part1 + part2
    
    Cn_result = outer_factor * total_sum_in_brackets
    
    return Cn_result
end
