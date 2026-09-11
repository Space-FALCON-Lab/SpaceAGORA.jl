"""
This module provides functions to compute laser forces on satellites, including cavity effects.
Forces are in Earth Centered Inertial (ECI) coordinates.
"""

# helper reservation utility
function reserve_link(i::Int, j::Int, current_helpers)
    if current_helpers[i, j] === 0 && current_helpers[j, i] === 0
        return true
    else
        return false
    end
end


"""
    Calculate laser radiation pressure forces on each satellite due to all others, including cavity forces where applicable.

    Inputs:
        u: state vector (6N)
        p: parameter dictionary with keys:
            :N               - number of satellites
            :Pmatrix        - power matrix (Pmatrix[i,j] = laser power from i to j)
            :c              - speed of light
            :use_los        - If true, apply line-of-sight and clearance checks before applying forces.
            :cavity         - dictionary of cavity properties, keyed by (i,j) tuples
                            (e.g., cavity[(1,2)] = Dict())
            :R_atm          - atmosphere radius for line-of-sight blocking (default R_ATMDEF)
            :atm_clearance  - minimum required clearance above atmosphere for line-of-sight
            :min_range      - minimum range for laser force application
            :max_range      - maximum range for laser force application

    Returns:
        forces: dictionary of forces keyed by (i,j) tuples, where forces[(i,j)] is the force on satellite j due to satellite i
"""
function laser_forces(u, p)
    N       = p[:N]
    Pm      = p[:Pmatrix]
    c       = p[:c]
    mu      = p[:mu]
    use_los = get(p, :use_los, false) # :use_los is the key you are looking up in p #use_los – If true, apply line-of-sight and clearance checks before applying forces.
    cavmap  = get(p, :cavity, Dict{Tuple{Int,Int},Dict{Symbol,Any}}()) #T{Iuplent,Int} is the key, Dict{Symbol,Any} is the value/content type
    R_atm   = get(p, :R_atm, R_ATMDEF)
    atm_cl  = get(p, :atm_clearance, 0.0)
    minR    = get(p, :min_range, 0.0)
    maxR    = get(p, :max_range, Inf)

    current_helpers = fill(0, (N, N))

    # positions at this state
    r = Array{Float64}(undef, 3, N)
    @inbounds for i in 1:N
        r[1,i]=u[idx(i,1)]; r[2,i]=u[idx(i,2)]; r[3,i]=u[idx(i,3)] # @inline idx(i, off) = 6*(i-1) + off  # state indexing helper
    end

    # --- GVE-optimal helper selection (optional) ---
    # p[:gve_schedule] is stored as a String (not Symbol) so that SciMLBase's
    # symbolic-map interpreter does not mistake it for a state-variable name.
    # Supported values: "gve_sma", "gve_ecc", "gve_inc", "gve_raan", "gve_argp"
    # Set p[:gve_target_idx] to specify the target satellite index
    # (defaults to the first entry of p[:target_ids], or N if not set).
    gve_sched = Symbol(get(p, :gve_schedule, "none"))  # convert String → Symbol internally
    if gve_sched in (:gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp)
        target_ids_p = get(p, :target_ids, [N])
        target_idx   = get(p, :gve_target_idx, isempty(target_ids_p) ? N : target_ids_p[1])
        cavmap = _gve_select_best_cavity(
            cavmap, target_idx, u, r,
            gve_sched, use_los, Float64(R_atm), Float64(atm_cl),
            Float64(minR), Float64(maxR), Float64(mu),
        )
    end

    forces = Dict{Tuple{Int,Int},SVector{3,Float64}}() # SVector type is for convenient 3D vector notation #https://juliaarrays.github.io/StaticArrays.jl/stable/api/#SVector
    @inline function add!(D, key::Tuple{Int,Int}, f::SVector{3,Float64}) # define a custom function to add forces to the dictionary
        if haskey(D, key); D[key] = D[key] + f else D[key] = f end
    end

    # mark pairs handled as cavities
    in_cavity = fill(false, N, N)
    for (pair, _) in cavmap
        i, j = pair
        in_cavity[i,j] = true; in_cavity[j,i] = true
    end

    # --- Single-pass beams for old sats ---
    # For each pair of satellites (i,j), if i != j and Pmatrix[i,j] > 0 and not in a cavity, calculate the force on j due to i
    @inbounds for i in 1:N, j in 1:N # F_on_Satj_due_to_Sati
        P = Pm[i,j] # power from i to j #Pm – Power matrix; Pmatrix[i,j] = laser power from i to j.
        if i!=j && P > 0 && !in_cavity[i,j] # if i and j are different satellites, power is positive, and the pair is not a cavity
            ri = @SVector [r[1,i], r[2,i], r[3,i]]
            rj = @SVector [r[1,j], r[2,j], r[3,j]]
            met = los_metrics(ri, rj; R_atm=R_atm) # get line-of-sight metrics between the two satellites
            los_ok   = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl)) # check if line-of-sight is okay based on the use_los flag and clearance, OR (||) if use_los is false, then don't need to check, so it's automatically okay
            range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
            link_available_ok = reserve_link(i, j, current_helpers)
            if los_ok && range_ok && link_available_ok
            #if los_ok && range_ok
                F = (P/c) * met.direction # Force calculations use ECI positions directly！！！！！！！！
                add!(forces, (i,j),  F)   # on j
                add!(forces, (j,i), -F)   # recoil on i

                current_helpers[i, :] .= 1 # 1 means single_pass # i is busy and not available to be helper
                current_helpers[:, i] .= 1 # 1 means single_pass # i is busy and not available to be target
                current_helpers[j, :] .= 1 # 1 means single_pass # j is busy and not available to be helper
                current_helpers[:, j] .= 1 # 1 means single_pass # j is busy and not available to be target
            end
        end
    end

    # --- Open cavity forces fror new sats ---
    # For each cavity pair (i,j), calculate the forces due to the cavity properties
    # cavity means a laser system between two satellites that enhances the laser power through constructive interference, allowing for stronger forces to be applied.
    for (pair, cav) in cavmap
        i, j = pair
        ri = @SVector [r[1,i], r[2,i], r[3,i]]
        rj = @SVector [r[1,j], r[2,j], r[3,j]]
        met = los_metrics(ri, rj; R_atm=R_atm)
        los_ok   = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl)) # check if line-of-sight is okay based on the use_los flag and clearance, OR (||) if use_los is false, then don't need to check, so it's automatically okay
        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
        link_available_ok = reserve_link(i, j, current_helpers)
        if los_ok && range_ok && link_available_ok
        #if los_ok && range_ok
            B   = cav[:B] # B – Cavity power-buildup factor (circulating vs input).
            Pin = cav[:Pin] # Pin – Input power for the cavity.
            û   = met.direction #type SVector{3,Float64} # unit vector from i to j # Force calculations use ECI positions directly！！！！！！！！

            # internal equal-opposite forces
            Fint = (B*Pin / c) * û 
            add!(forces, (i,j),  Fint)   # on j
            add!(forces, (j,i), -Fint)   # on i

            # leakage thrusts (fractions of circulating power B*Pin)
            fi  = get(cav, :leak_i_frac, 0.0) # fraction of B*Pin that leaks from sat i # cav is the cavity property pointed by key i,j in cavmap # 0.0 means no leakage by default
            fj  = get(cav, :leak_j_frac, 0.0)
            dir_i_along = get(cav, :leak_i_along,  true)   # along +û or -û # this determines the direction of leakage for sat i # true means along +û, false means along -û
            dir_j_along = get(cav, :leak_j_along, false)  # this determines the direction of leakage for sat j # true means along +û, false means along -û

            if fi > 0
                Fleak_i = (fi*B*Pin / c) * (dir_i_along ?  û : -û) # leakage thrust on sat i # dir_i_along ?  û : -û is a ternary operator that returns û if dir_i_along is true, and -û if dir_i_along is false
                add!(forces, (i,i), Fleak_i)   # thrust on sat i # (i,i) key means force on i due to i (leakage) #add!: if haskey(D, key); D[key] = D[key] + f else D[key] = f end
            end
            if fj > 0
                Fleak_j = (fj*B*Pin / c) * (dir_j_along ?  û : -û) # type SVector{3,Float64} # leakage thrust on sat j
                add!(forces, (j,j), Fleak_j)   # thrust on sat j
            end
            current_helpers[i, :] .= 2 # 2 means open_cavity
            current_helpers[:, i] .= 2 # 2 means open_cavity
            current_helpers[j, :] .= 2 # 2 means open_cavity
            current_helpers[:, j] .= 2 # 2 means open_cavity
        end
    end

    return forces, current_helpers
end

# =============================================================================
# GVE (Gauss Variational Equations) helper-selection scheduler
# =============================================================================

"""
    _gve_score_laser(elem, target_idx, helper_idx, u, mu) -> Float64

Compute the instantaneous GVE rate-of-change score for orbital element `elem`
that would result from firing the open-cavity laser link between `helper_idx`
(the helper) and `target_idx` (the target) given state vector `u` and
gravitational parameter `mu`.

The force magnitude is identical for every helper, so it cancels when comparing
scores, and only the direction-dependent scalar factor is returned.

Positive score means the chosen element would increase; negative means decrease.
The scheduler activates the helper with the highest positive score.

Supported `elem` values:
  :gve_sma   — maximise semi-major axis rate  ȧ
  :gve_ecc   — maximise eccentricity rate     ė
  :gve_inc   — maximise inclination rate      i̇
  :gve_raan  — maximise RAAN rate             Ω̇
  :gve_argp  — maximise arg-of-periapsis rate ω̇

Physics reference: Gauss Variational Equations (Battin, Schaub & Junkins).
"""
function _gve_score_laser(
    elem::Symbol,
    target_idx::Int,
    helper_idx::Int,
    u::AbstractVector{Float64},
    mu::Float64,
)::Float64
    # --- target & helper positions ---
    tgt_pos = SVector{3,Float64}(u[idx(target_idx,1)], u[idx(target_idx,2)], u[idx(target_idx,3)])
    tgt_vel = SVector{3,Float64}(u[idx(target_idx,4)], u[idx(target_idx,5)], u[idx(target_idx,6)])
    hlp_pos = SVector{3,Float64}(u[idx(helper_idx,1)], u[idx(helper_idx,2)], u[idx(helper_idx,3)])

    # Unit laser-force direction on the target: helper → target
    rel = tgt_pos - hlp_pos
    rho = norm(rel)
    rho > 0.0 || return 0.0
    f̂ = rel / rho

    # --- RTN basis of the target ---
    r = norm(tgt_pos)
    r > 0.0 || return 0.0
    er = tgt_pos / r
    h_vec = cross(tgt_pos, tgt_vel)
    hn    = norm(h_vec)
    hn > 0.0 || return 0.0
    en = h_vec / hn
    et = cross(en, er)

    aR = dot(f̂, er)   # Radial component of unit force
    aT = dot(f̂, et)   # Along-track component
    aN = dot(f̂, en)   # Cross-track component

    # --- Orbital elements of the target (vis-viva + angular momentum) ---
    v2    = dot(tgt_vel, tgt_vel)
    a_orb = -mu / (v2 - 2.0*mu/r)    # semi-major axis via vis-viva
    a_orb > 0.0 || return 0.0         # skip hyperbolic trajectories

    h_sq      = dot(h_vec, h_vec)
    p_slr     = h_sq / mu             # semi-latus rectum

    e_vec     = cross(tgt_vel, h_vec) / mu - tgt_pos / r
    e         = norm(e_vec)
    e_sq      = clamp(e * e, 0.0, 1.0 - 1e-12)
    sqrt_1me2 = sqrt(1.0 - e_sq)
    n_mean    = sqrt(mu / (a_orb^3))  # mean motion

    # True anomaly ν
    ν = acos(clamp(dot(e_vec / max(e, 1e-12), tgt_pos / r), -1.0, 1.0))
    dot(tgt_pos, tgt_vel) < 0.0 && (ν = 2π - ν)

    # --- GVE score for the requested element ---
    if elem === :gve_sma
        # ȧ = 2/(n√(1-e²)) * (e sinν · aR  +  p/r · aT)
        return (2.0 / (n_mean * sqrt_1me2)) * (e*sin(ν)*aR + (p_slr/r)*aT)

    elseif elem === :gve_ecc
        # ė = √(1-e²)/(na) * [sinν · aR + (cosν + (e+cosν)/(1+e cosν)) · aT]
        coeff_T = cos(ν) + (e + cos(ν)) / (1.0 + e*cos(ν))
        return (sqrt_1me2 / (n_mean * a_orb)) * (sin(ν)*aR + coeff_T*aT)

    else
        # :gve_inc, :gve_raan, :gve_argp — all need inclination i and arg-of-latitude u
        i_rad  = acos(clamp(h_vec[3] / hn, -1.0, 1.0))
        sin_i  = sin(i_rad)
        cos_i  = cos(i_rad)

        # Ascending-node vector: n_asc = ẑ × h
        n_asc = cross(SVector(0.0, 0.0, 1.0), h_vec)
        n_mag = norm(n_asc)

        # Argument of latitude u = ν + ω (robust to circular / equatorial)
        u_lat = if n_mag > 1e-12 && e > 1e-12
            # General case
            ω = acos(clamp(dot(n_asc/n_mag, e_vec/e), -1.0, 1.0))
            e_vec[3] < 0.0 && (ω = 2π - ω)
            ν + ω
        elseif n_mag > 1e-12
            # Circular (e ≈ 0): use angle from ascending node to position
            u_tmp = acos(clamp(dot(n_asc/n_mag, tgt_pos/r), -1.0, 1.0))
            tgt_pos[3] < 0.0 ? 2π - u_tmp : u_tmp
        else
            # Equatorial: use true longitude from x-axis
            atan(tgt_pos[2], tgt_pos[1])
        end

        denom = n_mean * a_orb^2 * sqrt_1me2

        if elem === :gve_inc
            # i̇ = r cos(u) / (na²√(1-e²)) · aN
            return (r * cos(u_lat) / denom) * aN

        elseif elem === :gve_raan
            # Ω̇ = r sin(u) / (na²√(1-e²) sin i) · aN  [singular at i = 0]
            abs(sin_i) < 1e-6 && return 0.0
            return (r * sin(u_lat) / (denom * sin_i)) * aN

        else  # :gve_argp
            # ω̇ = √(1-e²)/(nae) [-cosν · aR + (1+r/p) sinν · aT]
            #       - r sin(u) cos(i) / (na²√(1-e²) sin i) · aN
            abs(e) < 1e-6     && return 0.0  # circular: ω undefined
            abs(sin_i) < 1e-6 && return 0.0  # equatorial: ω undefined
            term_RT = (sqrt_1me2 / (n_mean * a_orb * e)) *
                      (-cos(ν)*aR + (1.0 + r/p_slr)*sin(ν)*aT)
            term_N  = -(r * sin(u_lat) * cos_i / (denom * sin_i)) * aN
            return term_RT + term_N
        end
    end
end

"""
    _gve_select_best_cavity(cavmap, target_idx, u, p, r, elem,
                             use_los, R_atm, atm_cl, minR, maxR)
                             -> Dict{Tuple{Int,Int}, Dict{Symbol,Any}}

Given the full `cavmap` (all helper↔target cavity pairs), return a new Dict
containing only the single cavity pair whose helper produces the highest
positive GVE score for `elem` from the target's current state.
Returns an empty Dict if no in-range helper has a positive score.
"""
function _gve_select_best_cavity(
    cavmap, target_idx::Int,
    u,
    r::Array{Float64,2},
    elem::Symbol,
    use_los::Bool, R_atm::Float64, atm_cl::Float64,
    minR::Float64, maxR::Float64,
    mu::Float64,
)
    best_score = 0.0
    best_pair  = nothing

    for (pair, _) in cavmap
        i, j = pair
        # Determine which end is the helper (the other end must be the target)
        if j == target_idx
            helper_idx = i
        elseif i == target_idx
            helper_idx = j
        else
            continue   # neither end is the designated target — skip
        end

        # LOS and range check
        ri  = @SVector [r[1,i], r[2,i], r[3,i]]
        rj  = @SVector [r[1,j], r[2,j], r[3,j]]
        met = los_metrics(ri, rj; R_atm=R_atm)
        los_ok   = (!use_los) || ((!met.blocked) && (met.clearance >= atm_cl))
        range_ok = (met.slant_range >= minR) && (met.slant_range <= maxR)
        (los_ok && range_ok) || continue

        score = _gve_score_laser(elem, target_idx, helper_idx, u, mu)
        if score > best_score
            best_score = score
            best_pair  = pair
        end
    end

    best_pair === nothing && return Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
    return Dict(best_pair => cavmap[best_pair])
end
