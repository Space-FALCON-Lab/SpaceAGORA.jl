module AccessMatrix

using LinearAlgebra
using Base.Threads
using ..PhysicsAdapter
using ..ConstellationConfig

"""
    logistic_access(d::Real; Rmax::Real, k::Real, x0::Real) -> Float64

Logistic function to model access probability based on distance.

# Arguments
- `d::Real`: Distance between satellite and client [m]
- `Rmax::Real`: Maximum effective range [m]
- `k::Real`: Steepness of the logistic curve
- `x0::Real`: Midpoint of the logistic curve (in units of d/Rmax)

# Returns
- `Float64`: Access probability in [0, 1]
"""
function logistic_access(d::Real; Rmax::Real, k::Real, x0::Real)
    x = d / Rmax
    access = 1 / (1 + exp(k * (x - x0)))
    return access
end

"""
    get_D_values(config_dict::AbstractDict, client_pos::AbstractVector{<:Real}, sat_pos::AbstractVector{<:Real}) -> Vector{Float64}

Compute the template thrust vector D from satellite to client.

D = u_max * (r_sat - r_client) / ||r_sat - r_client||

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `client_pos::AbstractVector{<:Real}`: Client position [x, y, z]
- `sat_pos::AbstractVector{<:Real}`: Satellite position [x, y, z]

# Returns
- `Vector{Float64}`: Thrust direction vector scaled by max acceleration
"""
function get_D_values(config_dict::AbstractDict, client_pos::AbstractVector{<:Real}, sat_pos::AbstractVector{<:Real})
    Δr1 = sat_pos[1] - client_pos[1]
    Δr2 = sat_pos[2] - client_pos[2]
    Δr3 = sat_pos[3] - client_pos[3]
    
    Δr = [Δr1, Δr2, Δr3]
    
    m_sc = config_dict["effector_params"]["sc_mass"]
    thrust = config_dict["effector_params"]["max_thrust"]
    u_max = thrust / m_sc
    
    D = u_max * (Δr / norm(Δr))
    
    return D
end

"""
    build_discrete_access_direction_tables(
        config_dict::AbstractDict,
        xbar::AbstractArray{<:Real,3},
        h::Real;
        time_offset::Real=0.0,
    ) -> Tuple{Array{Float64,3}, Array{Float64,4}}

Build discrete access direction tables for constellation optimization.

Computes:
- F[i_idx, p, k]: Binary access indicator (1 if within range, 0 otherwise)
- D[3, i_idx, p, k]: Thrust direction vectors

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `xbar::AbstractArray{<:Real,3}`: Client trajectories [3, N_knot, P]
- `h::Real`: Time step between knots
- `time_offset::Real=0.0`: Time offset for propagation

# Returns
- `Tuple{Array{Float64,3}, Array{Float64,4}}`: (F, D) access matrices
"""
function build_discrete_access_direction_tables(
    config_dict::AbstractDict,
    xbar::AbstractArray{<:Real,3},
    h::Real;
    time_offset::Real = 0.0,
)
    M = Int(config_dict["constellation_params"]["max_sats"])
    P = Int(config_dict["client_bounds"]["n_clients"])
    N_knot = size(xbar, 2)
    dt_access = Float64(config_dict["sim_params"]["dt"])
    active_ics = config_dict["active_ics"]
    const_ic = config_dict["constellation_params"]["initial_conditions"]
    μ = Float64(config_dict["physical_constants"]["mu"])
    Rmax = Float64(config_dict["effector_params"]["range"])
    m_sc = Float64(config_dict["effector_params"]["sc_mass"])
    thrust = Float64(config_dict["effector_params"]["max_thrust"])
    u_max = thrust / m_sc
    
    F = Array{Float64}(undef, M, P, N_knot)
    D = Array{Float64}(undef, 3, M, P, N_knot)
    
    # Get planet for Keplerian propagator
    planet = get_planet_from_config(config_dict)
    
    Threads.@threads for i_idx in eachindex(active_ics)
        sat_idx = active_ics[i_idx]
        a_i = const_ic[:a][sat_idx]
        e_i = haskey(const_ic, :e) ? const_ic[:e][sat_idx] :
              haskey(const_ic, :ecc) ? const_ic[:ecc][sat_idx] : 0.0
        inc_i = const_ic[:inc][sat_idx]
        raan_i = const_ic[:raan][sat_idx]
        argp_i = const_ic[:arg_p][sat_idx]
        ta_i = const_ic[:ta][sat_idx]
        
        for k in 1:N_knot
            t = Float64(time_offset) + (k - 1) * Float64(h)
            
            # Use SpaceAGORA's Keplerian propagator
            sat_state = kepler_satellite_state_adapter(
                t / dt_access, dt_access, a_i, e_i, inc_i, raan_i, argp_i, ta_i, planet
            )
            sat_pos = sat_state[1:3]
            
            for p in 1:P
                xc1 = Float64(xbar[1, k, p])
                xc2 = Float64(xbar[2, k, p])
                xc3 = Float64(xbar[3, k, p])
                
                Δr1 = sat_pos[1] - xc1
                Δr2 = sat_pos[2] - xc2
                Δr3 = sat_pos[3] - xc3
                
                dist = sqrt(Δr1^2 + Δr2^2 + Δr3^2)
                F[i_idx, p, k] = Float64(dist <= Rmax)
                
                if dist > 0.0
                    scale = u_max / dist
                    D[1, i_idx, p, k] = scale * Δr1
                    D[2, i_idx, p, k] = scale * Δr2
                    D[3, i_idx, p, k] = scale * Δr3
                else
                    D[1, i_idx, p, k] = 0.0
                    D[2, i_idx, p, k] = 0.0
                    D[3, i_idx, p, k] = 0.0
                end
            end
        end
    end
    
    return F, D
end

"""
    compute_access_matrix(config_dict::AbstractDict, client_trajectories::AbstractArray{<:Real,3}) -> Dict{String,Any}

Compute access matrix for constellation optimization.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary
- `client_trajectories::AbstractArray{<:Real,3}`: Client trajectories [3, N, P]

# Returns
- `Dict{String,Any}`: Dictionary containing F (access matrix) and D (direction vectors)
"""
function compute_access_matrix(config_dict::AbstractDict, client_trajectories::AbstractArray{<:Real,3})
    h = config_dict["sim_params"]["dt"]
    F, D = build_discrete_access_direction_tables(config_dict, client_trajectories, h)
    
    return Dict{String,Any}(
        "F" => F,
        "D" => D,
    )
end

export logistic_access, get_D_values, build_discrete_access_direction_tables, compute_access_matrix

end # module AccessMatrix
