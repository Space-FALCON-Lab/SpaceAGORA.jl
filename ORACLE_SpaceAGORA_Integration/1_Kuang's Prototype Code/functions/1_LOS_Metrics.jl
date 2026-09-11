"""
This module provides functions to compute line-of-sight metrics between satellites.
"""

"""        
    Calculate the geometry between positions ri (emitter) and rj (receiver).

    Inputs:
        ri: position of emitter (3-vector)
        rj: position of receiver (3-vector)
        R_atm: atmosphere radius for line-of-sight blocking (default R_ATMDEF)

    Returns:
        slant_range, rclosest, rclosest_norm, clearance (= rclosest_norm - R_atm),
        blocked (rclosest_norm ≤ R_atm), direction (unit i→j).
"""
function los_metrics(ri::SVector{3,Float64}, rj::SVector{3,Float64}; R_atm=R_ATMDEF)

    d = rj - ri #vector
    d2 = dot(d,d); #distance squared
    dnorm = sqrt(d2) + 1e-12 #distance
    # parameter t on segment [0,1] where chord is closest to origin
    t = clamp(-dot(ri, d) / d2, 0.0, 1.0)
    rcl = ri + t*d # closest point on line segment to origin
    rcln = norm(rcl)
    return (slant_range = dnorm,
            rclosest = rcl,
            rclosest_norm = rcln,
            clearance = rcln - R_atm,
            blocked = (rcln <= R_atm), # returns a boolean, true if blocked
            direction = d / dnorm)
end
