module ROEBoundsCalculator

using LinearAlgebra
using Statistics

"""
    compute_orbital_bounds_from_cluster(client_orbitals::Matrix{Float64}; sigma_mult::Float64=3.0, margin::Float64=0.05) -> Matrix{Float64}

Compute orbital element bounds from debris cluster using statistical analysis.
Returns a (6, 2) matrix with [min, max] for each orbital element.

# Arguments
- `client_orbitals`: Matrix of client orbital elements (n_clients × 6) [a, e, inc, raan, arg_p, ta]
- `sigma_mult`: Number of standard deviations for statistical bounds (default: 3.0)
- `margin`: Additional safety margin as fraction of range (default: 0.05)

# Returns
Matrix (6, 2) with orbital element bounds:
- Column 1: minimum values
- Column 2: maximum values
"""
function compute_orbital_bounds_from_cluster(client_orbitals::Matrix{Float64}; sigma_mult::Float64=3.0, margin::Float64=0.05)
    n_clients, n_elements = size(client_orbitals)
    
    if n_clients == 0
        return zeros(6, 2)
    end
    
    bounds = zeros(6, 2)
    
    for i in 1:6
        elements = client_orbitals[:, i]
        
        # Compute statistics
        mean_val = mean(elements)
        std_val = std(elements)
        
        # Statistical bounds (mean ± kσ)
        if std_val > 1e-8
            min_val = mean_val - sigma_mult * std_val
            max_val = mean_val + sigma_mult * std_val
        else
            # If no variance, use min/max with small margin
            min_val = minimum(elements)
            max_val = maximum(elements)
        end
        
        # Add safety margin
        range_val = max_val - min_val
        min_val -= margin * range_val
        max_val += margin * range_val
        
        # Element-specific constraints
        if i == 1  # Semi-major axis (must be positive)
            min_val = max(min_val, 1000.0)  # Minimum 1 km
        elseif i == 2  # Eccentricity (must be [0, 1))
            min_val = max(min_val, 0.0)
            max_val = min(max_val, 0.99)
        elseif i in [3, 4, 5, 6]  # Angles (inc, raan, arg_p, ta) - normalize to [0, 2π]
            min_val = mod(min_val, 2π)
            max_val = mod(max_val, 2π)
            if max_val < min_val
                max_val += 2π
            end
        end
        
        bounds[i, 1] = min_val
        bounds[i, 2] = max_val
    end
    
    return bounds
end

"""
    roe_to_orbital_elements(roe::Vector{Float64}, reference_orbit::Vector{Float64}) -> Vector{Float64}

Convert relative orbital elements (ROE) to absolute orbital elements.

# Arguments
- `roe`: Relative orbital elements [δa, δex, δey, δix, δiy, δλ]
- `reference_orbit`: Reference orbit [a_ref, e_ref, inc_ref, raan_ref, arg_p_ref, ta_ref]

# Returns
Absolute orbital elements [a, e, inc, raan, arg_p, ta]

# Notes
This implements the standard ROE to orbital element conversion used in CAPOConstellation.
The conversion follows the Schaub and Alfriend formulation.
"""
function roe_to_orbital_elements(roe::Vector{Float64}, reference_orbit::Vector{Float64})
    # Extract ROE components
    δa = roe[1]
    δex = roe[2]
    δey = roe[3]
    δix = roe[4]
    δiy = roe[5]
    δλ = roe[6]
    
    # Extract reference orbit
    a_ref = reference_orbit[1]
    e_ref = reference_orbit[2]
    inc_ref = reference_orbit[3]
    raan_ref = reference_orbit[4]
    arg_p_ref = reference_orbit[5]
    ta_ref = reference_orbit[6]
    
    # Compute mean motion
    μ = 3.986004418e14  # Earth gravitational parameter (m^3/s^2)
    n = sqrt(μ / a_ref^3)
    
    # Convert ROE to orbital element differences
    # Semi-major axis
    a = a_ref + δa
    
    # Eccentricity vector components
    ex_ref = e_ref * cos(arg_p_ref)
    ey_ref = e_ref * sin(arg_p_ref)
    
    ex = ex_ref + δex
    ey = ey_ref + δey
    
    e = sqrt(ex^2 + ey^2)
    arg_p = atan(ey, ex)
    
    # Inclination vector components
    ix_ref = tan(inc_ref / 2) * cos(raan_ref)
    iy_ref = tan(inc_ref / 2) * sin(raan_ref)
    
    ix = ix_ref + δix
    iy = iy_ref + δiy
    
    inc = 2 * atan(sqrt(ix^2 + iy^2))
    raan = atan(iy, ix)
    
    # True anomaly
    ta = ta_ref + δλ
    
    # Normalize angles to [0, 2π]
    arg_p = mod(arg_p, 2π)
    inc = mod(inc, 2π)
    raan = mod(raan, 2π)
    ta = mod(ta, 2π)
    
    return [a, e, inc, raan, arg_p, ta]
end

"""
    orbital_elements_to_roe(orbit::Vector{Float64}, reference_orbit::Vector{Float64}) -> Vector{Float64}

Convert absolute orbital elements to relative orbital elements (ROE).

# Arguments
- `orbit`: Absolute orbital elements [a, e, inc, raan, arg_p, ta]
- `reference_orbit`: Reference orbit [a_ref, e_ref, inc_ref, raan_ref, arg_p_ref, ta_ref]

# Returns
Relative orbital elements [δa, δex, δey, δix, δiy, δλ]
"""
function orbital_elements_to_roe(orbit::Vector{Float64}, reference_orbit::Vector{Float64})
    # Extract orbit
    a = orbit[1]
    e = orbit[2]
    inc = orbit[3]
    raan = orbit[4]
    arg_p = orbit[5]
    ta = orbit[6]
    
    # Extract reference orbit
    a_ref = reference_orbit[1]
    e_ref = reference_orbit[2]
    inc_ref = reference_orbit[3]
    raan_ref = reference_orbit[4]
    arg_p_ref = reference_orbit[5]
    ta_ref = reference_orbit[6]
    
    # Compute mean motion
    μ = 3.986004418e14  # Earth gravitational parameter (m^3/s^2)
    n = sqrt(μ / a_ref^3)
    
    # ROE: semi-major axis difference
    δa = a - a_ref
    
    # ROE: eccentricity vector differences
    ex = e * cos(arg_p)
    ey = e * sin(arg_p)
    ex_ref = e_ref * cos(arg_p_ref)
    ey_ref = e_ref * sin(arg_p_ref)
    
    δex = ex - ex_ref
    δey = ey - ey_ref
    
    # ROE: inclination vector differences
    ix = tan(inc / 2) * cos(raan)
    iy = tan(inc / 2) * sin(raan)
    ix_ref = tan(inc_ref / 2) * cos(raan_ref)
    iy_ref = tan(inc_ref / 2) * sin(raan_ref)
    
    δix = ix - ix_ref
    δiy = iy - iy_ref
    
    # ROE: true anomaly difference
    δλ = ta - ta_ref
    
    return [δa, δex, δey, δix, δiy, δλ]
end

"""
    compute_reference_orbit(client_orbitals::Matrix{Float64}) -> Vector{Float64}

Compute reference orbit from cluster mean orbital elements.

# Arguments
- `client_orbitals`: Matrix of client orbital elements (n_clients × 6)

# Returns
Reference orbit [a, e, inc, raan, arg_p, ta]
"""
function compute_reference_orbit(client_orbitals::Matrix{Float64})
    n_clients = size(client_orbitals, 1)
    
    if n_clients == 0
        return zeros(6)
    end
    
    # Compute mean orbital elements
    mean_orbit = mean(client_orbitals, dims=1)[:]
    
    # Normalize angles to [0, 2π]
    for i in [3, 4, 5, 6]
        mean_orbit[i] = mod(mean_orbit[i], 2π)
    end
    
    return mean_orbit
end

"""
    compute_shell_bounds_from_roe(roe_bounds::Matrix{Float64}, reference_orbit::Vector{Float64}) -> Matrix{Float64}

Convert ROE bounds to orbital element bounds.

# Arguments
- `roe_bounds`: ROE bounds (6, 2) with [min, max] for each ROE component
- `reference_orbit`: Reference orbit [a, e, inc, raan, arg_p, ta]

# Returns
Orbital element bounds (6, 2) with [min, max] for each orbital element
"""
function compute_shell_bounds_from_roe(roe_bounds::Matrix{Float64}, reference_orbit::Vector{Float64})
    orbital_bounds = zeros(6, 2)
    
    for i in 1:6
        # Convert ROE min/max to orbital element min/max
        roe_min = roe_bounds[i, 1]
        roe_max = roe_bounds[i, 2]
        
        roe_min_vec = zeros(6)
        roe_max_vec = zeros(6)
        roe_min_vec[i] = roe_min
        roe_max_vec[i] = roe_max
        
        orbital_min = roe_to_orbital_elements(roe_min_vec, reference_orbit)
        orbital_max = roe_to_orbital_elements(roe_max_vec, reference_orbit)
        
        orbital_bounds[i, 1] = orbital_min[i]
        orbital_bounds[i, 2] = orbital_max[i]
    end
    
    return orbital_bounds
end

export compute_orbital_bounds_from_cluster, roe_to_orbital_elements, orbital_elements_to_roe, compute_reference_orbit, compute_shell_bounds_from_roe

end # module ROEBoundsCalculator
