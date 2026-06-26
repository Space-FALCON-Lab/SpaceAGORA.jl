module Stage0p5GeneratorPropagation

using ..ConstellationUtils
using LinearAlgebra
using DataFrames
using Statistics

"""
    _resolve_disturbance_config(opt_params::AbstractDict) -> Dict{String,Any}

Resolve disturbance configuration from preset or individual toggles.
"""
function _resolve_disturbance_config(opt_params::AbstractDict)
    preset = Int(get(opt_params, "high_fidelity_preset", 0))
    
    # If preset is specified, use it to set individual toggles
    if preset >= 0
        config = _get_preset_disturbance_config(preset)
        
        # Allow individual toggles to override preset
        if haskey(opt_params, "disturbance_enable_j2")
            config["enable_j2"] = Bool(opt_params["disturbance_enable_j2"])
        end
        if haskey(opt_params, "disturbance_enable_harmonics")
            config["enable_harmonics"] = Bool(opt_params["disturbance_enable_harmonics"])
        end
        if haskey(opt_params, "disturbance_enable_atmosphere")
            config["enable_atmosphere"] = Bool(opt_params["disturbance_enable_atmosphere"])
        end
        if haskey(opt_params, "disturbance_enable_srp")
            config["enable_srp"] = Bool(opt_params["disturbance_enable_srp"])
        end
        if haskey(opt_params, "disturbance_enable_third_body")
            config["enable_third_body"] = Bool(opt_params["disturbance_enable_third_body"])
        end
        
        return config
    else
        # Use individual toggles only
        return Dict{String,Any}(
            "enable_j2" => Bool(get(opt_params, "disturbance_enable_j2", false)),
            "enable_harmonics" => Bool(get(opt_params, "disturbance_enable_harmonics", false)),
            "harmonics_degree" => Int(get(opt_params, "disturbance_harmonics_degree", 4)),
            "harmonics_order" => Int(get(opt_params, "disturbance_harmonics_order", 4)),
            "enable_atmosphere" => Bool(get(opt_params, "disturbance_enable_atmosphere", false)),
            "atmosphere_model" => String(get(opt_params, "disturbance_atmosphere_model", "nrlmsise00")),
            "enable_srp" => Bool(get(opt_params, "disturbance_enable_srp", false)),
            "enable_albedo" => Bool(get(opt_params, "disturbance_enable_albedo", false)),
            "enable_ir" => Bool(get(opt_params, "disturbance_enable_ir", false)),
            "enable_third_body" => Bool(get(opt_params, "disturbance_enable_third_body", false)),
            "third_body_bodies" => get(opt_params, "disturbance_third_body_bodies", ["sun", "moon"]),
        )
    end
end

"""
    _get_preset_disturbance_config(preset::Int) -> Dict{String,Any}

Get disturbance configuration for a preset fidelity level.
"""
function _get_preset_disturbance_config(preset::Int)
    presets = Dict{Int, Dict{String,Any}}(
        0 => Dict{String,Any}(
            "enable_j2" => false,
            "enable_harmonics" => false,
            "enable_atmosphere" => false,
            "enable_srp" => false,
            "enable_third_body" => false,
        ),
        1 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => false,
            "enable_atmosphere" => false,
            "enable_srp" => false,
            "enable_third_body" => false,
        ),
        2 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => false,
            "enable_atmosphere" => true,
            "atmosphere_model" => "exponential",
            "enable_srp" => false,
            "enable_third_body" => false,
        ),
        3 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => false,
            "enable_atmosphere" => true,
            "atmosphere_model" => "nrlmsise00",
            "enable_srp" => false,
            "enable_third_body" => false,
        ),
        4 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => false,
            "enable_atmosphere" => true,
            "atmosphere_model" => "nrlmsise00",
            "enable_srp" => true,
            "enable_albedo" => false,
            "enable_ir" => false,
            "enable_third_body" => false,
        ),
        5 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => false,
            "enable_atmosphere" => true,
            "atmosphere_model" => "nrlmsise00",
            "enable_srp" => true,
            "enable_albedo" => false,
            "enable_ir" => false,
            "enable_third_body" => true,
            "third_body_bodies" => ["sun", "moon"],
        ),
        6 => Dict{String,Any}(
            "enable_j2" => true,
            "enable_harmonics" => true,
            "harmonics_degree" => 8,
            "harmonics_order" => 8,
            "enable_atmosphere" => true,
            "atmosphere_model" => "gram",
            "enable_srp" => true,
            "enable_albedo" => true,
            "enable_ir" => true,
            "enable_third_body" => true,
            "third_body_bodies" => ["sun", "moon", "jupiter"],
        ),
    )
    
    return get(presets, preset, presets[0])
end

"""
    run_high_fidelity_generator_propagation(config_dict::AbstractDict) -> Dict{String,Any}

Run high-fidelity multithreaded propagation to compute time-varying generator fields.
"""
function run_high_fidelity_generator_propagation(config_dict::AbstractDict)
    # Extract Stage 0 results
    candidate_bank = config_dict["stage0_candidate_bank"]
    seed_indices = config_dict["stage0_seed_indices"]
    
    # Get propagation parameters
    opt_params = config_dict["optimizer_params"]
    num_threads = Int(get(opt_params, "high_fidelity_num_threads", Threads.nthreads()))
    
    # Resolve disturbance configuration
    disturbance_config = _resolve_disturbance_config(opt_params)
    
    # Get simulation parameters
    dt = Float64(config_dict["sim_params"]["dt"])
    N = Int(config_dict["sim_params"]["N"])
    H = Int(config_dict["mission"]["n_horizons"])
    P = Int(config_dict["client_bounds"]["n_clients"])
    Kd = Int(opt_params["controllable_set_dirs"])
    
    # Get planet for propagation
    planet = get_planet_from_config(config_dict)
    
    # Log disturbance configuration
    constellation_log("stage0p5", "Disturbance configuration"; 
        preset=get(opt_params, "high_fidelity_preset", -1),
        j2=disturbance_config["enable_j2"],
        harmonics=disturbance_config["enable_harmonics"],
        atmosphere=disturbance_config["enable_atmosphere"],
        srp=disturbance_config["enable_srp"],
        third_body=disturbance_config["enable_third_body"])
    
    # Propagate candidate satellites in parallel
    M = length(seed_indices)
    sat_trajectories = _propagate_candidates_parallel(
        candidate_bank, seed_indices, dt, N, H, planet, disturbance_config;
        num_threads=num_threads
    )
    
    # Compute time-varying generator fields
    h_fwd_time_varying = _compute_time_varying_h_fwd(
        sat_trajectories, config_dict, Kd, P, H, N
    )
    h_Wcorr_time_varying = _compute_time_varying_h_Wcorr(
        sat_trajectories, config_dict, Kd, P, H, N
    )
    
    return Dict{String,Any}(
        "h_fwd_time_varying" => h_fwd_time_varying,  # [Kd, M, P, H]
        "h_Wcorr_time_varying" => h_Wcorr_time_varying,  # [Kd, M, P, H]
        "sat_trajectories" => sat_trajectories,  # [6, N+1, M]
        "propagation_mode" => "high_fidelity",
        "disturbance_config" => disturbance_config,
    )
end

"""
    _propagate_candidates_parallel(candidate_bank::DataFrame, seed_indices::Vector{Int},
                                     dt::Float64, N::Int, H::Int, planet, disturbance_config::AbstractDict;
                                     num_threads::Int=Threads.nthreads()) -> Array{Float64,3}

Propagate candidate satellites in parallel using Threads.@threads with configurable disturbances.
"""
function _propagate_candidates_parallel(candidate_bank::DataFrame, seed_indices::Vector{Int},
                                         dt::Float64, N::Int, H::Int, planet, disturbance_config::AbstractDict;
                                         num_threads::Int=Threads.nthreads())
    M = length(seed_indices)
    trajectories = zeros(Float64, 6, N+1, M)
    
    # Set thread count if specified
    if num_threads > 0 && num_threads != Threads.nthreads()
        BLAS.set_num_threads(num_threads)
    end
    
    # Build disturbance model based on configuration
    disturbance_model = _build_disturbance_model(planet, disturbance_config)
    
    Threads.@threads for i in 1:M
        sat_idx = seed_indices[i]
        a = candidate_bank.a_km[sat_idx] * 1000.0  # Convert km to m
        e = candidate_bank.e[sat_idx]
        inc = deg2rad(candidate_bank.inc_deg[sat_idx])
        raan = deg2rad(candidate_bank.raan_deg[sat_idx])
        arg_p = deg2rad(candidate_bank.arg_p_deg[sat_idx])
        ta = deg2rad(candidate_bank.ta_deg[sat_idx])
        
        # Initialize state
        state = oe2cart(a, e, inc, raan, arg_p, ta, planet)
        
        for k in 0:N
            if k > 0
                # Propagate with disturbances
                state = _propagate_with_disturbances(
                    state, dt, planet, disturbance_model, disturbance_config
                )
            end
            trajectories[:, k+1, i] .= state
        end
    end
    
    return trajectories
end

"""
    _build_disturbance_model(planet, disturbance_config::AbstractDict)

Build disturbance model based on configuration.
"""
function _build_disturbance_model(planet, disturbance_config::AbstractDict)
    # This would integrate with SpaceAGORA's full disturbance models
    # For now, return a placeholder that the propagation function can use
    return Dict{String,Any}(
        "j2" => disturbance_config["enable_j2"] ? planet.J2 : 0.0,
        "atmosphere" => disturbance_config["enable_atmosphere"],
        "atmosphere_model" => get(disturbance_config, "atmosphere_model", "nrlmsise00"),
        "srp" => disturbance_config["enable_srp"],
        "third_body" => disturbance_config["enable_third_body"],
        "third_body_bodies" => get(disturbance_config, "third_body_bodies", []),
    )
end

"""
    _propagate_with_disturbances(state::Vector{Float64}, dt::Float64, planet, 
                                  disturbance_model::AbstractDict, disturbance_config::AbstractDict)

Propagate state with disturbances using SpaceAGORA's full dynamics.
"""
function _propagate_with_disturbances(state::Vector{Float64}, dt::Float64, planet, 
                                      disturbance_model::AbstractDict, disturbance_config::AbstractDict)
    # For now, use Keplerian propagator with J2 if enabled
    # Full integration with SpaceAGORA's disturbance models would go here
    
    r = state[1:3]
    v = state[4:6]
    
    # Compute acceleration
    a_grav = _compute_gravity_acceleration(r, planet, disturbance_config)
    a_total = a_grav
    
    # Add atmospheric drag if enabled
    if disturbance_config["enable_atmosphere"]
        a_drag = _compute_drag_acceleration(r, v, planet, disturbance_model)
        a_total += a_drag
    end
    
    # Add SRP if enabled
    if disturbance_config["enable_srp"]
        a_srp = _compute_srp_acceleration(r, planet, disturbance_model)
        a_total += a_srp
    end
    
    # Add third-body if enabled
    if disturbance_config["enable_third_body"]
        a_third = _compute_third_body_acceleration(r, planet, disturbance_model)
        a_total += a_third
    end
    
    # Simple Euler integration (replace with RK4 or SpaceAGORA's integrator)
    v_new = v + a_total * dt
    r_new = r + v * dt
    
    return vcat(r_new, v_new)
end

"""
    _compute_gravity_acceleration(r::Vector{Float64}, planet, disturbance_config::AbstractDict)

Compute gravitational acceleration with J2 if enabled.
"""
function _compute_gravity_acceleration(r::Vector{Float64}, planet, disturbance_config::AbstractDict)
    r_norm = norm(r)
    μ = planet.μ
    
    if disturbance_config["enable_j2"]
        J2 = planet.J2
        Rp = planet.Rp_e
        x, y, z = r
        r2 = r_norm^2
        
        # J2 acceleration
        a_j2 = 1.5 * J2 * μ * Rp^2 / r_norm^5
        ax = -μ * x / r_norm^3 + a_j2 * x / r_norm * (5 * z^2 / r2 - 1)
        ay = -μ * y / r_norm^3 + a_j2 * y / r_norm * (5 * z^2 / r2 - 1)
        az = -μ * z / r_norm^3 + a_j2 * z / r_norm * (5 * z^2 / r2 - 3)
        
        return [ax, ay, az]
    else
        # Inverse-square only
        return -μ * r / r_norm^3
    end
end

"""
    _compute_drag_acceleration(r::Vector{Float64}, v::Vector{Float64}, planet, disturbance_model::AbstractDict)

Compute atmospheric drag acceleration.
"""
function _compute_drag_acceleration(r::Vector{Float64}, v::Vector{Float64}, planet, disturbance_model::AbstractDict)
    # Placeholder for atmospheric drag
    # Would integrate with SpaceAGORA's atmosphere models
    return zeros(3)
end

"""
    _compute_srp_acceleration(r::Vector{Float64}, planet, disturbance_model::AbstractDict)

Compute solar radiation pressure acceleration.
"""
function _compute_srp_acceleration(r::Vector{Float64}, planet, disturbance_model::AbstractDict)
    # Placeholder for SRP
    # Would integrate with SpaceAGORA's SRP models
    return zeros(3)
end

"""
    _compute_third_body_acceleration(r::Vector{Float64}, planet, disturbance_model::AbstractDict)

Compute third-body gravitational acceleration.
"""
function _compute_third_body_acceleration(r::Vector{Float64}, planet, disturbance_model::AbstractDict)
    # Placeholder for third-body gravity
    # Would integrate with SpaceAGORA's N-body models
    return zeros(3)
end

"""
    oe2cart(a::Real, e::Real, i::Real, Ω::Real, ω::Real, ν::Real, planet)

Convert orbital elements to Cartesian state vector.
"""
function oe2cart(a::Real, e::Real, i::Real, Ω::Real, ω::Real, ν::Real, planet)
    μ = planet.μ
    
    # Semi-latus rectum
    p = a * (1 - e^2)
    
    # Distance
    r = p / (1 + e * cos(ν))
    
    # Position in perifocal frame
    x_p = r * cos(ν)
    y_p = r * sin(ν)
    
    # Velocity in perifocal frame
    v_fac = sqrt(μ / p)
    vx_p = -v_fac * sin(ν)
    vy_p = v_fac * (e + cos(ν))
    
    # Rotation matrices
    cO = cos(Ω)
    sO = sin(Ω)
    co = cos(ω)
    so = sin(ω)
    ci = cos(i)
    si = sin(i)
    
    # Rotate to ECI
    x = (cO * co - sO * so * ci) * x_p + (-cO * so - sO * co * ci) * y_p
    y = (sO * co + cO * so * ci) * x_p + (-sO * so + cO * co * ci) * y_p
    z = (so * si) * x_p + (co * si) * y_p
    
    vx = (cO * co - sO * so * ci) * vx_p + (-cO * so - sO * co * ci) * vy_p
    vy = (sO * co + cO * so * ci) * vx_p + (-sO * so + cO * co * ci) * vy_p
    vz = (so * si) * vx_p + (co * si) * vy_p
    
    return [x, y, z, vx, vy, vz]
end

"""
    get_planet_from_config(config_dict::AbstractDict)

Get planet parameters from configuration.
"""
function get_planet_from_config(config_dict::AbstractDict)
    phys = config_dict["physical_constants"]
    return Dict{String,Any}(
        "μ" => phys["mu"],
        "J2" => phys["J2"],
        "Rp_e" => phys["radius"],
    )
end

"""
    _compute_time_varying_h_fwd(sat_trajectories::Array{Float64,3}, config_dict::AbstractDict,
                                Kd::Int, P::Int, H::Int, N::Int) -> Array{Float64,4}

Compute time-varying forward reachable-set support coefficients.
"""
function _compute_time_varying_h_fwd(sat_trajectories::Array{Float64,3}, config_dict::AbstractDict,
                                     Kd::Int, P::Int, H::Int, N::Int)
    M = size(sat_trajectories, 3)
    h_fwd_tv = zeros(Float64, Kd, M, P, H)
    
    # Get client trajectories
    client_trajectories = config_dict["client_trajectories"]  # [6, N+1, P]
    
    # Build direction bank
    dirs_mat = _build_polyhedral_dirs(Kd)
    
    # For each horizon, compute support coefficients
    horizon_steps = div(N, H)
    
    for n in 1:H
        t_idx = n * horizon_steps
        for p in 1:P
            client_pos = client_trajectories[1:3, t_idx, p]
            
            for i in 1:M
                sat_pos = sat_trajectories[1:3, t_idx, i]
                
                # Compute support in each direction
                for q in 1:Kd
                    dir = dirs_mat[:, q]
                    # Support function: h = d^T * (sat_pos - client_pos)
                    h_fwd_tv[q, i, p, n] = dot(dir, sat_pos - client_pos)
                end
            end
        end
    end
    
    return h_fwd_tv
end

"""
    _compute_time_varying_h_Wcorr(sat_trajectories::Array{Float64,3}, config_dict::AbstractDict,
                                   Kd::Int, P::Int, H::Int, N::Int) -> Array{Float64,4}

Compute time-varying correction authority support coefficients.
"""
function _compute_time_varying_h_Wcorr(sat_trajectories::Array{Float64,3}, config_dict::AbstractDict,
                                      Kd::Int, P::Int, H::Int, N::Int)
    M = size(sat_trajectories, 3)
    h_Wcorr_tv = zeros(Float64, Kd, M, P, H)
    
    # Get effector parameters
    u_max = config_dict["effector_params"]["max_thrust"] / config_dict["effector_params"]["sc_mass"]
    
    # Build direction bank
    dirs_mat = _build_polyhedral_dirs(Kd)
    
    # For each horizon, compute correction authority
    horizon_steps = div(N, H)
    
    for n in 1:H
        t_idx = n * horizon_steps
        
        for i in 1:M, p in 1:P, q in 1:Kd
            # Correction authority: h_Wcorr = u_max * dt * ||d||
            dir = dirs_mat[:, q]
            h_Wcorr_tv[q, i, p, n] = u_max * config_dict["sim_params"]["dt"] * norm(dir)
        end
    end
    
    return h_Wcorr_tv
end

"""
    _build_polyhedral_dirs(Kd::Int) -> Array{Float64,2}

Build a bank of polyhedral direction vectors.
"""
function _build_polyhedral_dirs(Kd::Int)
    # Simple uniform distribution on sphere
    dirs = zeros(Float64, 3, Kd)
    
    for i in 1:Kd
        # Golden spiral for uniform distribution
        theta = 2π * i / (1 + sqrt(5))
        phi = acos(1 - 2 * (i - 0.5) / Kd)
        
        dirs[1, i] = sin(phi) * cos(theta)
        dirs[2, i] = sin(phi) * sin(theta)
        dirs[3, i] = cos(phi)
    end
    
    return dirs
end

export run_high_fidelity_generator_propagation

end # module Stage0p5GeneratorPropagation
