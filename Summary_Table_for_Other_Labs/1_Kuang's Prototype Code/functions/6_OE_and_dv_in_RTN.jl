"""
This module provides functions to compute orbital elements and delta-v in the RTN frame.
"""

const IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))

"""
    Compute orbital elements time series for each satellite from solution object.

    Inputs:
        sol: solution object from ODE solver
        μ: gravitational parameter

    Returns:
        elems: vector of vectors of named tuples containing orbital elements for each satellite over time
"""
# Elements for every sat and time sample
function elements_time_series(sol, μ, coe_fn = rv2coe) #time seriies is a sequence of data points, typically consisting of successive measurements made over a time interval.
    # coe_fn: converter function to use; defaults to rv2coe (angles in (-π,π]).
    #         Pass rv2coe_2pi for [0,2π) angles (avoids ±180° jump in difference plots).
    N = Int(length(sol.u[1]) ÷ 6) # number of satellites
    elems = [Vector{NamedTuple}(undef, length(sol.t)) for _ in 1:N] # elems is a vector of named tuples, each named tuple contains the orbital elements for a satellite at a time sample
    # _ in for _ in 1:N is a throwaway variable.It means "I don't care about the value, just repeat this N times."
    for k in eachindex(sol.t) # the index (e.g. 1, 2, 3) for time value in sol.t
        uk = sol.u[k] # at this time, index is k, and state is uk      
        for i in 1:N
            r = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]]
            v = @SVector [uk[idx(i,4)], uk[idx(i,5)], uk[idx(i,6)]] #r and v at thiss time for satellite i
            elems[i][k] = coe_fn(r, v, μ) #orbital elements at this time for satellite i
        end
    end
    return elems
end

"""
    Print initial and final orbital elements for each satellite in a formatted manner.
    Inputs:
        sol: solution object from ODE solver
        μ: gravitational parameter
        degrees: whether to print angles in degrees (default true)

    Returns:
        None (prints to console)
"""
# Pretty print initial/final elements for each sat
function print_initial_final_elements(sol, μ; degrees=true)
    elems = elements_time_series(sol, μ)
    todeg(x) = degrees ? rad2deg(x) : x
    for i in 1:length(elems)
        e0 = elems[i][1]; ef = elems[i][end]
        println("\nSatellite $i:")
        @printf("  a0 = %12.3f m   af = %12.3f m   Δa = %+10.3f m\n", e0.a, ef.a, ef.a-e0.a)
        @printf("  e0 = %12.6g     ef = %12.6g     Δe = %+12.6g\n", e0.e, ef.e, ef.e-e0.e)
        @printf("  i0 = %9.6f°     if = %9.6f°     Δi = %+9.6f°\n", todeg(e0.i), todeg(ef.i), todeg(ef.i-e0.i))
        @printf("  Ω0 = %9.6f°     Ωf = %9.6f°     ΔΩ = %+9.6f°\n", todeg(e0.Ω), todeg(ef.Ω), todeg(ef.Ω-e0.Ω))
        if !isnan(e0.ν)
            @printf("  ω0 = %9.6f°     ωf = %9.6f°\n", todeg(e0.ω), todeg(ef.ω))
            @printf("  ν0 = %9.6f°     νf = %9.6f°\n", todeg(e0.ν), todeg(ef.ν))
        else
            @printf("  u0 = %9.6f°     uf = %9.6f°     Δu = %+9.6f°\n",
                    todeg(e0.u), todeg(ef.u), todeg(ef.u-e0.u))
        end
    end
end

"""
    Compute RTN (Radial, Tangential, Normal) basis vectors from position and velocity.
    
    Inputs:
        r: position vector
        v: velocity vector

    Returns:
        (er, et, en): tuple of RTN basis unit vectors
"""
# RTN basis and Δv time series # RTN means Radial, Tangential, Normal
function rtn_basis(r::SVector{3,Float64}, v::SVector{3,Float64}) #input position and velocity vectors, output RTN basis vectors
    er = r / (norm(r) + 1e-12) # radial unit vector
    h  = cross(r, v); en = h / (norm(h) + 1e-12) # normal unit vector
    et = cross(en, er) # tangential unit vector
    return (er, et, en)
end

"""
    Compute cumulative delta-v in RTN frame for each satellite over time due to laser/cavity forces.
    # this dv is already cumulative with ΔP[j] .+= Fav[j]*dt!!!
    
    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
    Returns:
        t: time vector
        Δv_RTN_hist: vector of RTN delta-v time series for each satellite  
"""
# cumulative Δv(t) per satellite due to lasers/cavity; returns RTN time series
function delta_v_RTN_time_series(sol, p)
    N = p[:N]
    masses = p[:masses]
    t = sol.t; u = sol.u

    ΔP_RTN       = [zeros(3) for _ in 1:N]
    ΔP_RTN_hist  = [zeros(3, length(t)) for _ in 1:N]
    Δv_RTN_hist  = [zeros(3, length(t)) for _ in 1:N]

    for k in 1:length(t)-1
        dt  = t[k+1] - t[k]
        uk  = u[k]
        uk1 = u[k+1]

        # Endpoint forces (Cartesian)
        Fk, current_helpers  = laser_forces(uk,  p)
        Fk1, current_helpers1 = laser_forces(uk1, p)
        keys_union = union(keys(Fk), keys(Fk1))

        # Net Cartesian endpoint forces
        Fnet_k  = [SVector(0.0,0.0,0.0) for _ in 1:N]
        Fnet_k1 = [SVector(0.0,0.0,0.0) for _ in 1:N]
        for key in keys_union
            (_, j) = key
            FL = haskey(Fk,  key) ? Fk[key]  : SVector(0.0,0.0,0.0)
            FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0,0.0,0.0)
            Fnet_k[j]  += FL
            Fnet_k1[j] += FR
        end

        # Project to RTN at endpoints and integrate ΔP in RTN via trapezoid
        for j in 1:N
            r  = @SVector [uk[idx(j,1)],  uk[idx(j,2)],  uk[idx(j,3)]]
            v  = @SVector [uk[idx(j,4)],  uk[idx(j,5)],  uk[idx(j,6)]]
            r1 = @SVector [uk1[idx(j,1)], uk1[idx(j,2)], uk1[idx(j,3)]]
            v1 = @SVector [uk1[idx(j,4)], uk1[idx(j,5)], uk1[idx(j,6)]]
            er, et, en   = rtn_basis(r,  v)
            er1,et1,en1  = rtn_basis(r1, v1)

            F_R_k  = SVector(dot(Fnet_k[j],  er),  dot(Fnet_k[j],  et),  dot(Fnet_k[j],  en))
            F_R_k1 = SVector(dot(Fnet_k1[j], er1), dot(Fnet_k1[j], et1), dot(Fnet_k1[j], en1))

            Fav_RTN = 0.5*(F_R_k + F_R_k1)
            ΔP_RTN[j] += Fav_RTN * dt
            ΔP_RTN_hist[j][:, k+1] = ΔP_RTN[j]
            Δv_RTN_hist[j][:, k+1] = ΔP_RTN_hist[j][:, k+1] ./ masses[j]
        end
    end

    return t, Δv_RTN_hist
end

# # cumulative Δv(t) per satellite due to lasers/cavity; returns RTN time series
# function delta_v_RTN_time_series(sol, p)

#     N = p[:N]; # number of satellites
#     masses = p[:masses]; # vector of masses for all satellites
#     t = sol.t; u = sol.u
#     ΔP = [zeros(3) for _ in 1:N] # initialize change in momentum for each satellite
#     ΔP_hist = [zeros(3, length(t)) for _ in 1:N] # history of change in momentum for each satellite over time

#     for k in 1:length(t)-1
#         dt  = t[k+1]-t[k]
#         Fk  = laser_forces(u[k],   p)
#         Fk1 = laser_forces(u[k+1], p)
#         keys_union = union(keys(Fk), keys(Fk1))

#         Fav = [SVector(0.0,0.0,0.0) for _ in 1:N]
#         for key in keys_union
#             (_, j) = key
#             FL = haskey(Fk,  key) ? Fk[key]  : SVector(0.0,0.0,0.0) # the frame of the force is the same as the frame of the satellite
#             FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0,0.0,0.0)
#             Fav[j] += 0.5*(FL + FR)
#         end

#         for j in 1:N
#             ΔP[j] .+= Fav[j]*dt #. in .+= means element-wise addition # ΔP is cumulative change in momentum
#             ΔP_hist[j][:,k+1] = ΔP[j]
#         end
#     end

#     # Project to RTN at each time
#     Δv_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
#     for k in eachindex(t) # for every time instant
#         uk = u[k]
#         for j in 1:N # for sat j
#             r = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
#             v = @SVector [uk[idx(j,4)], uk[idx(j,5)], uk[idx(j,6)]]
#             er, et, en = rtn_basis(r, v) # RTN basis vectors at this time for sat j
#             dv = ΔP_hist[j][:,k] ./ masses[j] # [:, k] means do this for every row in column k
#             Δv_RTN_hist[j][:,k] = [dot(dv, er), dot(dv, et), dot(dv, en)]
#         end
#     end
#     return t, Δv_RTN_hist
# end

"""
    laser_force_RTN_time_series(sol, p)

Compute net laser/cavity force time series for each satellite, projected into the RTN frame.

Returns:
    t, F_RTN_hist  where F_RTN_hist[j] is 3×T (rows: R,T,N).

Side effects:
    Stores in p:
      :F_hist_cartesian      => averaged Cartesian F history (3×T per sat, stored at k+1)
      :ΔP_hist_cartesian     => cumulative ΔP history (3×T per sat, stored at k+1)
      :ΔP_RTN_hist           => cumulative ΔP projected to RTN (3×T per sat)
"""
function laser_force_RTN_time_series(sol, p)
    N = p[:N]
    t = sol.t; u = sol.u

    # Project to RTN and integrate ΔP in RTN
    #F_hist = [zeros(3, length(t)) for _ in 1:N]
    F_RTN_hist   = [zeros(3, length(t)) for _ in 1:N]
    ΔP_RTN       = [zeros(3) for _ in 1:N]
    ΔP_RTN_hist  = [zeros(3, length(t)) for _ in 1:N]

    # Trapezoidal integrate ΔP in RTN over segments
    for k in 1:length(t)-1
        dt  = t[k+1] - t[k]
        uk  = u[k]
        uk1 = u[k+1]

        # Endpoint forces (Cartesian)
        Fk, current_helpers  = laser_forces(uk,  p)
        Fk1, current_helpers1 = laser_forces(uk1, p)
        keys_union = union(keys(Fk), keys(Fk1))

        # Net Cartesian forces at endpoints
        Fnet_k  = [SVector(0.0,0.0,0.0) for _ in 1:N]
        Fnet_k1 = [SVector(0.0,0.0,0.0) for _ in 1:N]
        for key in keys_union
            (_, j) = key
            FL = haskey(Fk,  key) ? Fk[key]  : SVector(0.0,0.0,0.0)
            FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0,0.0,0.0)
            Fnet_k[j]  += FL
            Fnet_k1[j] += FR
        end

        # Bases at endpoints
        for j in 1:N
            r  = @SVector [uk[idx(j,1)],  uk[idx(j,2)],  uk[idx(j,3)]]
            v  = @SVector [uk[idx(j,4)],  uk[idx(j,5)],  uk[idx(j,6)]]
            r1 = @SVector [uk1[idx(j,1)], uk1[idx(j,2)], uk1[idx(j,3)]]
            v1 = @SVector [uk1[idx(j,4)], uk1[idx(j,5)], uk1[idx(j,6)]]
            er, et, en   = rtn_basis(r,  v)
            er1,et1,en1  = rtn_basis(r1, v1)

            # Endpoint RTN forces
            F_R_k  = SVector(dot(Fnet_k[j],  er),  dot(Fnet_k[j],  et),  dot(Fnet_k[j],  en))
            F_R_k1 = SVector(dot(Fnet_k1[j], er1), dot(Fnet_k1[j], et1), dot(Fnet_k1[j], en1))

            # Trapezoidal step in RTN
            Fav_RTN = 0.5*(F_R_k + F_R_k1)
            F_RTN_hist[j][:, k+1] = Fav_RTN
            ΔP_RTN[j] += Fav_RTN * dt
            ΔP_RTN_hist[j][:, k+1] = ΔP_RTN[j]
        end
    end

    return t, F_RTN_hist, ΔP_RTN_hist
end

# function laser_force_RTN_time_series(sol, p)
#     N = p[:N]
#     t = sol.t; u = sol.u

#     # Build Cartesian averaged force history, aligned like ΔP_hist (store at k+1)
#     F_hist = [zeros(3, length(t)) for _ in 1:N]

#     # Also accumulate ΔP in Cartesian, same method as delta_v_RTN_time_series
#     ΔP     = [zeros(3) for _ in 1:N]
#     ΔP_hist = [zeros(3, length(t)) for _ in 1:N]

#     for k in 1:length(t)-1
#         dt  = t[k+1] - t[k]
#         Fk  = laser_forces(u[k],   p)
#         Fk1 = laser_forces(u[k+1], p)
#         keys_union = union(keys(Fk), keys(Fk1))

#         Fav = [SVector(0.0,0.0,0.0) for _ in 1:N]
#         for key in keys_union
#             (_, j) = key
#             FL = haskey(Fk,  key) ? Fk[key]  : SVector(0.0,0.0,0.0)
#             FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0,0.0,0.0)
#             Fav[j] += 0.5*(FL + FR)
#         end

#         for j in 1:N
#             F_hist[j][:, k+1] = Fav[j]              # store averaged Cartesian force
#             ΔP[j] .+= Fav[j] * dt                   # accumulate momentum change
#             # HERE, ΔP and Fav are in Cartesian frame, although force is always positive in the direction of the beam
#             # but since the direction of the beam is rotating, so, Fav is also rotating, so ΔP will not continue goring, bu instead will sometimes goes up and sometime goes down 
#             ΔP_hist[j][:, k+1] = ΔP[j]              # record cumulative ΔP
#         end
#     end

#     # Project stored Cartesian forces to RTN at each time
#     F_RTN_hist   = [zeros(3, length(t)) for _ in 1:N]
#     ΔP_RTN_hist  = [zeros(3, length(t)) for _ in 1:N] # cumulative ΔP projected to RTN, which does have much meaning
#     for k in eachindex(t)
#         uk = u[k]
#         for j in 1:N
#             r = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
#             v = @SVector [uk[idx(j,4)], uk[idx(j,5)], uk[idx(j,6)]]
#             er, et, en = rtn_basis(r, v)

#             Fc = F_hist[j][:, k]
#             Pc = ΔP_hist[j][:, k]

#             F_RTN_hist[j][:, k]  = [dot(Fc, er), dot(Fc, et), dot(Fc, en)]
#             ΔP_RTN_hist[j][:, k] = [dot(Pc, er), dot(Pc, et), dot(Pc, en)]

#         end
#     end

#     return t, F_RTN_hist, ΔP_RTN_hist
# end

"""
    Compute RTN position, velocity, and acceleration time series for each satellite from solution object.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites

    Returns:
        t: time vector
        r_RTN_hist: vector of RTN position time series for each satellite
        v_RTN_hist: vector of RTN velocity time series for each satellite
        a_RTN_hist: vector of RTN acceleration time series for each satellite
"""
function r_v_a_RTN_time_series(sol, p)
    N = p[:N]  # number of satellites
    t = sol.t  # time vector
    u = sol.u  # state vector (position and velocity)

    # Initialize RTN position, velocity, and acceleration histories
    r_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
    v_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
    a_RTN_hist = [zeros(3, length(t)) for _ in 1:N]

    # Loop through each time step
    for k in eachindex(t)
        uk = u[k]  # state at time k
        duk = similar(uk)  # Allocate space for derivatives
        nbody_photon!(duk, uk, p, t[k])  # Compute derivatives manually
        for j in 1:N  # for each satellite
            # Extract position, velocity, and acceleration for satellite j
            r = @SVector [uk[idx(j, 1)], uk[idx(j, 2)], uk[idx(j, 3)]]
            v = @SVector [duk[idx(j, 1)], duk[idx(j, 2)], duk[idx(j, 3)]]
            a = @SVector [duk[idx(j, 4)], duk[idx(j, 5)], duk[idx(j, 6)]]

            # Compute RTN basis vectors
            er, et, en = rtn_basis(r, v)

            # Project position, velocity, and acceleration into RTN frame
            r_RTN_hist[j][:, k] = [dot(r, er), dot(r, et), dot(r, en)]
            v_RTN_hist[j][:, k] = [dot(v, er), dot(v, et), dot(v, en)]
            a_RTN_hist[j][:, k] = [dot(a, er), dot(a, et), dot(a, en)]
        end
    end

    return t, r_RTN_hist, v_RTN_hist, a_RTN_hist
end

"""
    Compute RTN acceleration time series for each satellite due to gravity and laser forces separately.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites

    Returns:
        t: time vector
        a_gravity_RTN_hist: RTN acceleration time series due to gravity for each satellite
        a_laser_RTN_hist: RTN acceleration time series due to laser forces for each satellite
"""
function a_gravity_laser_RTN_time_series(sol, p)
    N = p[:N]  # number of satellites
    t = sol.t  # time vector
    u = sol.u  # state vector (position and velocity)

    # Initialize RTN acceleration histories for gravity and laser forces
    a_gravity_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
    a_laser_RTN_hist = [zeros(3, length(t)) for _ in 1:N]

    # Loop through each time step
    for k in 1:length(t)-1
        dt = t[k+1] - t[k]  # Time step duration
        uk = u[k]           # State at time k
        uk1 = u[k+1]        # State at time k+1

        # Initialize position and velocity arrays
        r = Array{Float64}(undef, 3, N)  # Position array
        v = similar(r)                   # Velocity array

        # Unpack positions and velocities
        @inbounds for i in 1:N
            r[1, i] = uk[idx(i, 1)]
            r[2, i] = uk[idx(i, 2)]
            r[3, i] = uk[idx(i, 3)]
            v[1, i] = uk[idx(i, 4)]
            v[2, i] = uk[idx(i, 5)]
            v[3, i] = uk[idx(i, 6)]
        end

        # Compute accelerations due to gravity and laser forces
        a_gravity = similar(r)  # Acceleration due to gravity
        a_laser = similar(r)    # Acceleration due to laser forces

        # Gravity
        @inbounds for i in 1:N
            rx, ry, rz = r[1, i], r[2, i], r[3, i]
            ρ = sqrt(rx^2 + ry^2 + rz^2) + 1e-12  # Distance from central body
            fac = -p[:mu] / ρ^3  # Gravitational factor
            a_gravity[1, i] = fac * rx
            a_gravity[2, i] = fac * ry
            a_gravity[3, i] = fac * rz
        end

        # Laser forces (averaging over time step)
        Fk, current_helpers = laser_forces(uk, p)       # Forces at time k
        Fk1, current_helpers1 = laser_forces(uk1, p)    # Forces at time k+1
        keys_union = union(keys(Fk), keys(Fk1))  # Combine keys from both time steps

        @inbounds for j in 1:N
            a_laser[:, j] .= 0.0  # Reset laser acceleration for satellite j
        end

        @inbounds for key in keys_union
            (_, j) = key
            FL = haskey(Fk, key)  ? Fk[key]  : SVector(0.0, 0.0, 0.0)  # Force at time k
            FR = haskey(Fk1, key) ? Fk1[key] : SVector(0.0, 0.0, 0.0)  # Force at time k+1
            Fav = 0.5 * (FL + FR)  # Average force over the time step

            # Compute laser acceleration
            a_laser[1, j] += Fav[1] / p[:masses][j]
            a_laser[2, j] += Fav[2] / p[:masses][j]
            a_laser[3, j] += Fav[3] / p[:masses][j]
        end

        # Compute RTN basis and project accelerations into RTN frame
        for j in 1:N
            r_sat = @SVector [r[1, j], r[2, j], r[3, j]]
            v_sat = @SVector [v[1, j], v[2, j], v[3, j]]
            a_grav_sat = @SVector [a_gravity[1, j], a_gravity[2, j], a_gravity[3, j]]
            a_laser_sat = @SVector [a_laser[1, j], a_laser[2, j], a_laser[3, j]]

            # Compute RTN basis vectors
            er, et, en = rtn_basis(r_sat, v_sat)

            # Project accelerations into RTN frame
            a_gravity_RTN_hist[j][:, k] = [dot(a_grav_sat, er), dot(a_grav_sat, et), dot(a_grav_sat, en)]
            a_laser_RTN_hist[j][:, k] = [dot(a_laser_sat, er), dot(a_laser_sat, et), dot(a_laser_sat, en)]
        end
    end

    # Handle the last time step (k = length(t))
    k = length(t)
    for j in 1:N
        a_gravity_RTN_hist[j][:, k] = a_gravity_RTN_hist[j][:, k-1]
        a_laser_RTN_hist[j][:, k] = a_laser_RTN_hist[j][:, k-1]
    end

    return t, a_gravity_RTN_hist, a_laser_RTN_hist
end


"""
    Compute link status (range, clearance, in-range) time series for a satellite pair (i,j).
    
    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            plus all keys required by los_metrics()
        i: index of first satellite (default 1)
        j: index of second satellite (default 2)

    Returns:
        t: time vector
        rng: slant range time series [m]
        clr: atmosphere clearance time series [m]
        inrng: boolean time series indicating if link is enabled
"""
# Compute range/clearance/status over time for a pair (i,j)
function link_status_time_series(sol, p; i::Int=1, j::Int=2)

    t = sol.t
    use_los = get(p, :use_los, false) # whether to use line-of-sight checks
    R_atm   = get(p, :R_atm, R_ATMDEF)
    atm_cl  = get(p, :atm_clearance, 0.0) # minimum clearance [m] above atmosphere
    minR    = get(p, :min_range, 0.0) # minimum range [m] for link
    maxR    = get(p, :max_range, Inf)

    rng   = zeros(length(t))     # slant range [m], which is the distance between two points # (length(t)) meaning an array of zeros with the same length as t
    clr   = zeros(length(t))     # atmosphere clearance [m] (>=0 passes)
    inrng = falses(length(t))    # link enabled (both gates satisfied)

    for k in eachindex(t) # at each time instant
        uk = sol.u[k]
        ri = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]]
        rj = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
        met = los_metrics(ri, rj; R_atm=R_atm) # get line-of-sight metrics between the two satellites
        los_ok   = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl))
        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
        rng[k]   = met.slant_range  # rng is the slant range [m]
        clr[k]   = met.clearance # clr is the clearance above atmosphere
        inrng[k] = (los_ok && range_ok) # inrng is true if both line-of-sight and range conditions are satisfied
    end
    return t, rng, clr, inrng # return time, range, clearance, and in-range status at each time step
end

"""
    Estimate link duty cycle, which is the fraction of time the link is enabled, and upper bound on delta-v for a satellite pair (i,j).

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by los_metrics() and laser_forces()
        i: index of first satellite (default 1)
        j: index of second satellite (default 2)

    Returns:
        Named tuple with fields:
            duty: duty cycle (fraction of time link is enabled)
            on_time: total time link is enabled [s]
            a_nom: nominal acceleration if always on [m/s²]
            dv_upper_bound: upper bound on delta-v assuming always-on thrust [m/s]
"""

function link_duty_and_estimate(sol, p; i::Int=1, j::Int=2)
    # Duty (original "either-endpoint" rule)
    t, _, _, inrng = link_status_time_series(sol, p; i=i, j=j)
    dt = diff(t)
    onseg = (inrng[1:end-1] .| inrng[2:end])
    on_time = sum(dt[onseg])
    T = t[end] - t[1]
    duty = on_time/T # avoid div by zero

    # Nominal acceleration from cavity or single-pass
    c       = get(p, :c, C)
    masses  = p[:masses]
    cavmap  = get(p, :cavity, Dict{Tuple{Int,Int},Dict{Symbol,Any}}())
    Pm      = get(p, :Pmatrix, nothing)

    a_nom = 0.0
    if haskey(cavmap, (i,j)) # haskey(cavmap, (i,j)) means check if the key (i,j) exists in the dictionary cavmap
                             # cavmap is the :cavity dictionary in p, which contains cavity properties for satellite pairs
        cav = cavmap[(i,j)]  # cavity info (:B=>100, :Pin=>10, ...) between sats i and j
        B, Pin = cav[:B], cav[:Pin]
        a_nom = (2*B*Pin/c) / masses[j]      # accel on receiver j
    elseif haskey(cavmap, (j,i))
        cav = cavmap[(j,i)]
        B, Pin = cav[:B], cav[:Pin]
        a_nom = (2*B*Pin/c) / masses[i]      # flipped: accel on i
    elseif Pm !== nothing && Pm[i,j] > 0
        a_nom = (Pm[i,j]/c) / masses[j]      # single-pass i→j
    elseif Pm !== nothing && Pm[j,i] > 0
        a_nom = (Pm[j,i]/c) / masses[i]      # single-pass j→i
    end

    dv_upper = a_nom * on_time
    return (duty=duty, on_time=on_time, a_nom=a_nom, dv_upper_bound=dv_upper)
end

# function link_duty_and_estimate(sol, p; i::Int=1, j::Int=2)
    
#     t, _, _, inrng = link_status_time_series(sol, p; i=i, j=j)
#     dt = diff(t)
#     # Treat a segment as "on" if either endpoint is on
#     onseg = (inrng[1:end-1] .| inrng[2:end]) # onseg is a boolean array that is true if either endpoint of the segment is in range
#     on_time = sum(dt[onseg])
#     T = t[end] - t[1]

#     # Pull Pin,B and mass from p/cavity and masses
#     cav = p[:cavity][(1,2)] # cavity between sats i and j
#     Pin = cav[:Pin]; B = cav[:B]; m = p[:masses][1]  # assume both same mass
#     a_nom = (2*B*Pin/C) / m # nominal acceleration if always on
#     dv_upper = a_nom * on_time  # ignores thrust-direction rotation

#     return (duty = on_time/T, on_time = on_time, a_nom = a_nom, dv_upper_bound = dv_upper)
# end
