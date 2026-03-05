Base.@kwdef struct RuntimePolicyConfig
    warn_normalize::Bool = true
    allow_typed_normalize::Bool = false
    gram_per_sat_instances::Bool = false
    srp_ephemeris_cache::Bool = true
    nbody_ephemeris_cache::Bool = true
    planet_frame_cache::Bool = true
    spice_rhs_memo::Bool = true
end
