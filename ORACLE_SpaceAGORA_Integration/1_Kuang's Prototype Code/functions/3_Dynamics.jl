"""
This modules tells the ODE solver how the state u of satellites evolves.
"""

# --- Basic ODE: gravity + laser/cavity ---
"""
    This function computes du from U,
        where u = [ r1, v1, r2, v2, ..., rN, vN ] (6N vector)
        and du = [ v1, a1, v2, a2, ..., vN, aN ] (6N vector)

    This function tells the ODE solver how the state u evolves.

    Inputs:
        du: derivative vector to be filled (6N) (initialized to be zeros in the beginning of the function)
        u: state vector (6N)
        p: parameter dictionary with keys:
            :mu       - gravitational parameter
            :N        - number of satellites
            :masses   - vector of each body’s mass
            (optional J2 terms)
            :use_J2   - Bool to enable J2 acceleration (default false)
            :J2       - dimensionless J2 coefficient (default 1.08262668e-3)
            :Re       - equatorial radius of central body [m] (default 6378137.0)
            (optional atmospheric drag)
            :use_drag - Bool to enable atmospheric drag (default false)
            :Cd       - drag coefficient (scalar or length-N vector, default 2.2)
            :A        - cross-sectional area [m^2] (scalar or length-N vector, default 0.0)
            :rho0     - atmospheric density at reference altitude [kg/m^3] (default 3.614e-11)
            :H        - scale height [m] (default 88667.0)
            :h_ref    - reference altitude above equator [m] (default 700000.0)
            :omega_E  - planetary rotation rate [rad/s] (default 7.2921159e-5)
            plus all keys required by laser_forces()
        t: time (not used here, but required by ODE solvers)

    Returns:
        nothing (du is modified in place)
"""
function nbody_photon!(du, u, p, t) #u is the state vector (contains all the states), p is the parameters dictionary, t is time
    mu = p[:mu]; # gravitational parameter
    N = p[:N]; # number of satellites
    masses = p[:masses] # Vector of each body’s mass
    useJ2 = get(p, :use_J2, false)
    J2 = get(p, :J2, 1.08262668e-3)
    Re = get(p, :Re, 6378137.0)
    useDrag = get(p, :use_drag, false)
    Cd_def = get(p, :Cd, 2.2)
    A_def  = get(p, :A, 0.0)
    rho0   = get(p, :rho0, 3.614e-11)     # ~ density near 700 km (US Standard, rough)
    H      = get(p, :H, 88667.0)          # scale height ~ 88.7 km
    h_ref  = get(p, :h_ref, 700000.0)     # reference altitude 700 km
    omegaE = get(p, :omega_E, 7.2921159e-5)
    fill!(du, 0.0) #initializing the derivative vector to zeros

    # unpack
    r = Array{Float64}(undef, 3, N); v = similar(r) #similar() creates an array of the same type and size as r
    @inbounds for i in 1:N
        r[1,i]=u[idx(i,1)]; r[2,i]=u[idx(i,2)]; r[3,i]=u[idx(i,3)] # @inline idx(i, off) = 6*(i-1) + off  # state indexing helper
        v[1,i]=u[idx(i,4)]; v[2,i]=u[idx(i,5)]; v[3,i]=u[idx(i,6)] #idx is a function that calculates the index in the state vector for satellite i and offset off (1 to 6)
    end

    # gravity
    @inbounds for i in 1:N
        rx,ry,rz = r[1,i], r[2,i], r[3,i]
        ρ = sqrt(rx^2 + ry^2 + rz^2) + 1e-12 # ρ = distance from central body + small term to avoid div by 0
        fac = -mu/ρ^3
        du[idx(i,4)] += fac*rx # acceleration in x
        du[idx(i,5)] += fac*ry # acceleration in y
        du[idx(i,6)] += fac*rz # acceleration in z
    end

    # optional J2 perturbation (Earth oblateness)
    if useJ2
        @inbounds for i in 1:N
            rx, ry, rz = r[1,i], r[2,i], r[3,i]
            r2 = rx^2 + ry^2 + rz^2
            ρ = sqrt(r2) + 1e-12 #Small number avoids division by zero.
            z2_over_r2 = (rz^2) / (r2 + 1e-24)
            facJ2 = 1.5 * J2 * mu * Re^2 / (ρ^5)
            cxy = (5.0*z2_over_r2 - 1.0)
            cz  = (5.0*z2_over_r2 - 3.0)
            du[idx(i,4)] += facJ2 * rx * cxy
            du[idx(i,5)] += facJ2 * ry * cxy
            du[idx(i,6)] += facJ2 * rz * cz
        end
    end #https://docs.poliastro.space/en/stable/autoapi/poliastro/core/perturbations/index.html

    # atmospheric drag (optional)
    if useDrag
        # atmosphere co-rotates with planet about +Z with rate omegaE
        ωx, ωy, ωz = 0.0, 0.0, omegaE
        @inbounds for i in 1:N
            rx, ry, rz = r[1,i], r[2,i], r[3,i]
            vx, vy, vz = v[1,i], v[2,i], v[3,i]
            # altitude above reference sphere
            h = sqrt(rx^2 + ry^2 + rz^2) - Re
            ρ = rho0 * exp(-(h - h_ref)/H)
            # relative velocity w.r.t. rotating atmosphere: v_rel = v - ω × r
            vrelx = vx - (ωy*rz - ωz*ry)
            vrely = vy - (ωz*rx - ωx*rz)
            vrelz = vz - (ωx*ry - ωy*rx)
            vrel_mag = sqrt(vrelx^2 + vrely^2 + vrelz^2) + 1e-12
            # per-satellite Cd and area (scalar or vector)
            Cd_i = isa(Cd_def, AbstractVector) ? Cd_def[i] : Cd_def
            A_i  = isa(A_def, AbstractVector)  ? A_def[i]  : A_def
            if A_i != 0.0 && ρ > 0.0
                coeff = -0.5 * Cd_i * A_i * ρ / masses[i]
                du[idx(i,4)] += coeff * vrel_mag * vrelx
                du[idx(i,5)] += coeff * vrel_mag * vrely
                du[idx(i,6)] += coeff * vrel_mag * vrelz
            end
        end
    end

    # laser/cavity forces
    F, current_helpers = laser_forces(u, p)  # Dict((i,j)=>F_on_j)
    #all_helpers[t, :, :] .= current_helpers

    for ((i,j), Fj) in F
        du[idx(j,4)] += Fj[1]/masses[j]
        du[idx(j,5)] += Fj[2]/masses[j]
        du[idx(j,6)] += Fj[3]/masses[j]
    end

    # kinematics
    @inbounds for i in 1:N
        du[idx(i,1)] = v[1,i]
        du[idx(i,2)] = v[2,i]
        du[idx(i,3)] = v[3,i]
    end
    return nothing
end

# --- Augmented ODE: gravity + laser/cavity + Δv state for tracked sat ---
#augumented track of delta-v for cetain satellites j=sat_trk
"""
    This function computes du from U,
        where u = [ 6N orbital states ; (optional) 3 Δv states ]
            = [ r1, v1, r2, v2, ..., rN, vN ; Δvx, Δvy, Δvz ] (6N or 6N+3 vector)
        and du = [ 6N orbital states ; (optional) 3 Δv states ]
            = [ v1, a1, v2, a2, ..., vN, aN ; a_trx, a_try, a_trz ] (6N or 6N+3 vector)

    This function tells the ODE solver how the state u evolves.

    Inputs:
        du: derivative vector to be filled (6N or 6N+3) (initialized to be zeros in the beginning of the function)
        u: state vector (6N or 6N+3)
        p: parameter dictionary with keys:
            :mu       - gravitational parameter
            :N        - number of satellites
            :masses   - vector of each body’s mass
            :track_dv_sat - index of satellite to track Δv for (1 to N), or 0 to disable
            (optional J2 terms)
            :use_J2   - Bool to enable J2 acceleration (default false)
            :J2       - dimensionless J2 coefficient (default 1.08262668e-3)
            :Re       - equatorial radius of central body [m] (default 6378137.0)
            plus all keys required by laser_forces()
        t: time (not used here, but required by ODE solvers)

    Returns:
        nothing (du is modified in place)
"""
function nbody_photon_aug!(du, u, p, t)
    mu      = p[:mu]
    N       = p[:N]
    masses  = p[:masses]
    sat_trk = get(p, :track_dv_sat, 0)  # 0 => no Δv tracking # track_dv_sat – Index of satellite to track Δv for (1 to N), or 0 to disable.
    mode = get(p, :dv_target_mode, :magnitude) 
    useJ2 = get(p, :use_J2, false)
    J2 = get(p, :J2, 1.08262668e-3)
    Re = get(p, :Re, 6378137.0)
    useDrag = get(p, :use_drag, false)
    Cd_def = get(p, :Cd, 2.2)
    A_def  = get(p, :A, 0.0)
    rho0   = get(p, :rho0, 3.614e-11)
    H      = get(p, :H, 88667.0)
    h_ref  = get(p, :h_ref, 700000.0)
    omegaE = get(p, :omega_E, 7.2921159e-5)
    #println("DV tracking mode: ", mode)
    # layout: [ 6N orbital states ; (optional) 3 Δv states ]
    n6 = 6N # n6 is the number of orbital states #6N is the number of satellites times 6 (3 position + 3 velocity)
    # In Julia, a numeric literal (like 2 or 3.14) placed immediately before a variable (an identifier) or a parenthesized expression implies multiplication, so you do not need the *
    # However, the * is required for multiplication between two variables or identifiers. 
    fill!(du, 0.0)

    # unpack r,v
    r = Array{Float64}(undef, 3, N)
    v = similar(r)
    @inbounds for i in 1:N
        r[1,i]=u[idx(i,1)]; r[2,i]=u[idx(i,2)]; r[3,i]=u[idx(i,3)]
        v[1,i]=u[idx(i,4)]; v[2,i]=u[idx(i,5)]; v[3,i]=u[idx(i,6)]
    end

    # gravity
    @inbounds for i in 1:N
        rx,ry,rz = r[1,i], r[2,i], r[3,i]
        ρ = sqrt(rx^2 + ry^2 + rz^2) + 1e-12
        fac = -mu/ρ^3
        du[idx(i,4)] += fac*rx
        du[idx(i,5)] += fac*ry
        du[idx(i,6)] += fac*rz
    end

    # optional J2 perturbation (Earth oblateness)
    if useJ2
        @inbounds for i in 1:N
            rx, ry, rz = r[1,i], r[2,i], r[3,i]
            r2 = rx^2 + ry^2 + rz^2
            ρ = sqrt(r2) + 1e-12
            z2_over_r2 = (rz^2) / (r2 + 1e-24) 
            facJ2 = 1.5 * J2 * mu * Re^2 / (ρ^5)
            cxy = (5.0*z2_over_r2 - 1.0)
            cz  = (5.0*z2_over_r2 - 3.0)
            du[idx(i,4)] += facJ2 * rx * cxy
            du[idx(i,5)] += facJ2 * ry * cxy
            du[idx(i,6)] += facJ2 * rz * cz
        end
    end

    # initialize non-grav acceleration accumulator for tracked sat
    a_tr = SVector(0.0, 0.0, 0.0)

    # atmospheric drag (optional)
    if useDrag
        ωx, ωy, ωz = 0.0, 0.0, omegaE
        @inbounds for i in 1:N
            rx, ry, rz = r[1,i], r[2,i], r[3,i]
            vx, vy, vz = v[1,i], v[2,i], v[3,i]
            h = sqrt(rx^2 + ry^2 + rz^2) - Re # altitude above reference sphere
            ρ = rho0 * exp(-(h - h_ref)/H) # atmospheric density
            # relative velocity w.r.t. rotating atmosphere: v_rel = v - ω × r
            vrelx = vx - (ωy*rz - ωz*ry)
            vrely = vy - (ωz*rx - ωx*rz)
            vrelz = vz - (ωx*ry - ωy*rx)
            vrel_mag = sqrt(vrelx^2 + vrely^2 + vrelz^2) + 1e-12 # avoid div by 0
            # per-satellite Cd and area (scalar or vector)
            Cd_i = isa(Cd_def, AbstractVector) ? Cd_def[i] : Cd_def # get drag coefficient for satellite i #isa() checks if Cd_def is a vector
            A_i  = isa(A_def, AbstractVector)  ? A_def[i]  : A_def
            if A_i != 0.0 && ρ > 0.0
                coeff = -0.5 * Cd_i * A_i * ρ / masses[i]
                ax = coeff * vrel_mag * vrelx
                ay = coeff * vrel_mag * vrely
                az = coeff * vrel_mag * vrelz
                du[idx(i,4)] += ax
                du[idx(i,5)] += ay
                du[idx(i,6)] += az
                if sat_trk == i
                    # include drag in tracked non-grav acceleration
                    a_tr += SVector(ax, ay, az)
                end
            end
        end
    end

    # lasers/cavities
    F, current_helpers = laser_forces(u, p)  # Dict((i,j)=>F_on_j)
    for ((i,j), Fj) in F
        du[idx(j,4)] += Fj[1]/masses[j] # @inline idx(i, off) = 6*(i-1) + off
        du[idx(j,5)] += Fj[2]/masses[j]
        du[idx(j,6)] += Fj[3]/masses[j]
        if sat_trk == j
            a_tr += Fj ./ masses[j]  # non-grav accel of tracked sat #./ is element-wise division
            # Since Fj is a vector and masses[j] is a scalar, this element-wise division is the same as:
            # a_tr[1] += Fj[1] / masses[j]
            # a_tr[2] += Fj[2] / masses[j]
            # a_tr[3] += Fj[3] / masses[j]
        end
    end

    # kinematics
    @inbounds for i in 1:N
        du[idx(i,1)] = v[1,i]
        du[idx(i,2)] = v[2,i]
        du[idx(i,3)] = v[3,i]
    end

    # Δv dynamics (only if tracking enabled)
    if sat_trk != 0
        if mode == "component" # mode == "component" mean tracking RTN components of dv,
                              # so need to first convert everyting into RTN frome at here for acceleration
            # RTN triad from current r,v of the tracked satellite
            r_tr = @SVector [r[1,sat_trk], r[2,sat_trk], r[3,sat_trk]]
            v_tr = @SVector [v[1,sat_trk], v[2,sat_trk], v[3,sat_trk]]
            er = r_tr / (norm(r_tr) + 1e-12)
            en = cross(r_tr, v_tr); en /= (norm(en) + 1e-12)
            et = cross(en, er)

            # project non-grav accel into RTN
            a_RTN = SVector(dot(a_tr, er), dot(a_tr, et), dot(a_tr, en))

            # integrate Δv in RTN coordinates
            du[n6+1] = a_RTN[1]  # a_R
            du[n6+2] = a_RTN[2]  # a_T
            du[n6+3] = a_RTN[3]  # a_N
            #println("a_RTN_R: ", a_RTN[1])
  
        else
            du[n6+1] = a_tr[1]
            du[n6+2] = a_tr[2]
            du[n6+3] = a_tr[3]
            #println("a_tr_x: ", a_tr[1])
        end
    end
    return nothing
end
