"""
This module provides diagnostic functions to evaluate the state of the satellite system.
"""

"""
    This function computes the total linear momentum vector, its magnitude,
        and the magnitude of linear momentum for each satellite.

    Inputs:
        u: state vector (6N)
        masses: vector of each body’s mass

    Returns:
        P: total momentum vector
        Pmag: magnitude of total momentum
"""
function total_momentum(u, masses)

    N = length(masses);
    #println("==================== Number of satellites for total momentum: $N =========================")
    P = zeros(3); pmag = zeros(N) # P is total momentum vector, pmag is the magnitude of momentum for each satellite
    @inbounds for i in 1:N
        v = @SVector [u[idx(i,4)], u[idx(i,5)], u[idx(i,6)]]
        p = masses[i]*v
        P += p; pmag[i] = norm(p) # norm(p) is the magnitude of momentum for satellite i
    end
    return P, norm(P), pmag
end

"""
    This function computes the specific angular momentum magnitude for each satellite 
        and the total angular momentum magnitude of the system.

    Inputs:
        u: state vector (6N)   
        masses: vector of each body’s mass

    Returns:
        h_each: vector of specific angular momentum magnitude for each satellite
        H: total angular momentum magnitude of the system
"""
function angular_momentum(u, masses)

    N = length(masses) # number of satellites
    h_each = zeros(N); # specific angular momentum magnitude for each satellite
    H = zeros(3) # total angular momentum vector
    @inbounds for i in 1:N # @inbounds means to skip bounds checking for performance
        r = @SVector [u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]] # position vector of satellite i
        p = @SVector [masses[i]*u[idx(i,4)], masses[i]*u[idx(i,5)], masses[i]*u[idx(i,6)]] # momentum vector of satellite i
        h = cross(r,p) # angular momentum vector of satellite i
        h_each[i] = norm(h); H += h
    end
    return h_each, norm(H)
end

"""
    This function computes the total mechanical energy (kinetic + potential) for each satellite.

    Inputs:
        u: state vector (6N)   
        masses: vector of each body’s mass
        mu: gravitational parameter

    Returns:
        E: vector of total mechanical energy for each satellite
"""
function orbital_energy(u, masses, mu)

    N = length(masses); E = zeros(N)
    @inbounds for i in 1:N
        x,y,z = u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]
        r = sqrt(x^2+y^2+z^2) + 1e-12
        vx,vy,vz = u[idx(i,4)], u[idx(i,5)], u[idx(i,6)]
        K = 0.5*masses[i]*(vx^2+vy^2+vz^2)
        U = -mu*masses[i]/r
        E[i] = K+U
    end
    return E
end

"""
    This function computes the total mechanical impulse and work done by laser/cavity forces over the entire simulation.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            plus all keys required by laser_forces()

    Returns:
        ΔP: dictionary of total mechanical impulse keyed by (i,j) tuples, where ΔP[(i,j)] is the impulse on satellite j due to satellite i
        ΔE: dictionary of total work done keyed by (i,j) tuples, where ΔE[(i,j)] is the work done on satellite j due to satellite i
"""
function evaluate_laser_exchanges(sol, p) # sol is the ODE solution object, p is the parameters dictionary
    N = p[:N]; # number of satellites
    t = sol.t; u = sol.u # sol – ODE solution object # sol.t is the time vector, sol.u is the state vector at each time
    ΔP = Dict{Tuple{Int,Int},MVector{3,Float64}}() # ΔP is the change in momentum # MVector is a mutable static vector, mutable means you can change its elements
    ΔE = Dict{Tuple{Int,Int},Float64}() # ΔE is sum of the change in energy for all satellites

    @inline vel(u_, j) = @SVector [u_[idx(j,4)], u_[idx(j,5)], u_[idx(j,6)]] # trailing underscore "_" to avoid naming conflicts with other u

    @inbounds for k in 1:length(t)-1
        dt  = t[k+1]-t[k]
        Fk, current_helpers  = laser_forces(u[k],   p)
        Fk1, current_helpers1 = laser_forces(u[k+1], p)
        keys_union = union(keys(Fk), keys(Fk1)) # union of keys from Fk and Fk1, an example of the keys is (i,j) meaning force on j due to i

        vL = [vel(u[k],   j) for j in 1:N] # vL is the velocity at the left time step
        vR = [vel(u[k+1], j) for j in 1:N]

        for key in keys_union
            i, j = key
            FL = haskey(Fk,  key) ? Fk[key]  : SVector(0.0,0.0,0.0) # if key (i, j) exists in Fk, then FL = Fk[key], else FL = SVector(0.0,0.0,0.0)
            FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0,0.0,0.0)

            Favg = 0.5*(FL + FR) # average force over the time step
            vj   = 0.5*(vL[j] + vR[j]) # average velocity of satellite j over the time step

            # Here, 0.5*(Left+Right) is the trapezoidal rule
            # This applies trapezoidal integration of mechanical impulse and work from *all* forces.
            # Keys (i,j) mean "force on j due to i".

            if !haskey(ΔP, key); ΔP[key] = MVector(0.0,0.0,0.0); end # if key (i, j) does not exist in ΔP, then initialize it to MVector(0.0,0.0,0.0)
            if !haskey(ΔE, key); ΔE[key] = 0.0; end
            ΔP[key] .+= Favg*dt
            ΔE[key]  += dot(Favg, vj)*dt
        end # this loop goes through each key in the union of keys from Fk and Fk1 and calculates the average force and velocity, then updates the change in momentum and energy for that key
    end
    return ΔP, ΔE
end

"""
    saturation_number(h_helper_km, i_helper_deg, h_target_km, i_target_deg, L_max_m) -> N_sat

Compute the helper-satellite saturation number: the minimum number of helpers
needed to maintain continuous laser coverage of the target orbit.

Formula:
    α* = acos( (rh² + rt² − L²) / (2·rh·rt·cos(Δi)) )
    N_sat = ⌈π / α*⌉

# Arguments
- `h_helper_m`   : helper orbit altitude [m]
- `i_helper_deg` : helper orbit inclination [deg]
- `h_target_m`   : target orbit altitude [m]
- `i_target_deg` : target orbit inclination [deg]
- `L_max_m`      : maximum laser range [m]

# Returns
- `N_sat`     : saturation number (Int)
- `alpha_star`: half-angle α* [deg] (for reference)
"""
function saturation_number(h_helper_m, i_helper_deg, h_target_m, i_target_deg, L_max_m)
    r_h  = R_EARTH + h_helper_m   # helper orbit radius [m]
    r_t  = R_EARTH + h_target_m   # target orbit radius [m]
    Δi   = abs(i_helper_deg - i_target_deg) * π / 180.0  # inclination difference [rad]

    arg = (r_h^2 + r_t^2 - L_max_m^2) / (2 * r_h * r_t * cos(Δi))

    if arg < -1.0 || arg > 1.0
        error("L_max_m=$L_max_m m is too small to bridge the two orbits: acos argument = $arg (must be in [-1, 1])")
    end

    α_star = acos(arg)                   # [rad]
    N_sat  = ceil(Int, π / α_star)

    return N_sat, rad2deg(α_star)
end
