"""
This module provides functions to compute laser forces on satellites, including cavity effects.
Forces are in Earth Centered Inertial (ECI) coordinates.
"""

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
    use_los = get(p, :use_los, false) # :use_los is the key you are looking up in p #use_los – If true, apply line-of-sight and clearance checks before applying forces.
    cavmap  = get(p, :cavity, Dict{Tuple{Int,Int},Dict{Symbol,Any}}()) #T{Iuplent,Int} is the key, Dict{Symbol,Any} is the value/content type
    R_atm   = get(p, :R_atm, R_ATMDEF)
    atm_cl  = get(p, :atm_clearance, 0.0)
    minR    = get(p, :min_range, 0.0)
    maxR    = get(p, :max_range, Inf)

    # positions at this state
    r = Array{Float64}(undef, 3, N)
    @inbounds for i in 1:N
        r[1,i]=u[idx(i,1)]; r[2,i]=u[idx(i,2)]; r[3,i]=u[idx(i,3)] # @inline idx(i, off) = 6*(i-1) + off  # state indexing helper
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
            if los_ok && range_ok
                F = (P/c) * met.direction # Force calculations use ECI positions directly！！！！！！！！
                add!(forces, (i,j),  F)   # on j
                add!(forces, (j,i), -F)   # recoil on i
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
        if los_ok && range_ok
            B   = cav[:B] # B – Cavity power-buildup factor (circulating vs input).
            Pin = cav[:Pin] # Pin – Input power for the cavity.
            û   = met.direction #type SVector{3,Float64} # unit vector from i to j # Force calculations use ECI positions directly！！！！！！！！

            # internal equal-opposite forces
            Fint = (2*B*Pin / c) * û 
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
        end
    end

    return forces
end
