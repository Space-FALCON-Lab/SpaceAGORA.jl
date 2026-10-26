"""
    get_magnetic_field_dipole(r_ecef::AbstractVector)

Calculates the Earth's magnetic field using a fast tilted dipole approximation.

This model is significantly faster than the full WMM and is suitable for use
inside performance-critical code like numerical integrators.

# Args

- `r_ecef`: The position vector of the spacecraft in an Earth-Centered,
            Earth-Fixed (ECEF) frame [meters].

# Returns

- A 3-element `SVector` representing the magnetic field `[Bx, By, Bz]` in the
  ECEF frame [Tesla].
"""
function get_magnetic_field_dipole(r_ecef::AbstractVector, L_PI::MMatrix{3, 3, Float64})
    r_norm = norm(r_ecef)
    r_hat = r_ecef / r_norm

    # Cosine of the magnetic colatitude
    cos_colat = dot(r_hat, M_HAT_ECEF)

    # Dipole field equation with Earth's dipole MOMENT pointing toward the
    # geomagnetic SOUTH pole (m_hat = -M_HAT_ECEF):
    #   B = B0 (R/r)^3 [3 (m_hat.r_hat) r_hat - m_hat]
    #     = B0 (R/r)^3 [M_HAT_ECEF - 3 cos_colat r_hat]
    # Checks: equator (cos_colat = 0) gives B = B0 * M_HAT_ECEF (north, B0);
    # north pole gives -2 B0 M_HAT_ECEF (down). The pre-2026-07 sign returned
    # the ANTIPARALLEL field (170 deg from IGRF at LEO).
    B_ecef = B0_2020 * (R_EARTH_MODEL / r_norm)^3 * (M_HAT_ECEF - 3 * cos_colat * r_hat)

    return L_PI' * B_ecef
end


"""
    calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)

Calculates the magnetic torque exerted on a magnetic dipole by an external
magnetic field.

**τ** = **m** × **B**

The magnetic dipole moment `m` and the magnetic field `B` vectors **must** be
expressed in the same reference frame. The resulting torque vector `τ` will be
in that same frame.

# Args

- `m`: The magnetic dipole moment vector of the spacecraft [A·m²].
- `B`: The external magnetic field vector [Tesla].

# Returns

- A 3-element `SVector` representing the magnetic torque `[τ_x, τ_y, τ_z]` [N·m].
"""
function calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)
    # Ensure inputs are StaticArrays for performance
    m_svector = SVector{3, Float64}(m)
    B_svector = SVector{3, Float64}(B)

    # Calculate the torque using the cross product
    # τ = m × B
    τ = cross(m_svector, B_svector)

    return τ
end

"""
    EddyCurrentDampingModel <: AbstractForceTorqueModel

Passive magnetic damping from eddy currents induced in a conductive
spacecraft structure rotating in the geomagnetic field:

    τ_body = k_e · B × (B × ω)

which is dissipative for `k_e > 0` (τ·ω = k_e[(B·ω)² − |B|²|ω|²] ≤ 0 by
Cauchy–Schwarz) and vanishes when ω is parallel to B — so the angular-velocity
component about the field line is untouched while the perpendicular component
decays. This is the spin-down and field-alignment mechanism for passively
magnet-stabilized spacecraft that carry no hysteresis material.

`k_e` has units N·m·s/T² and lumps the structure's conductivity and geometry.
Field source options are identical to [`MagneticTorqueRodModel`](@ref)
(`:dipole` fast tilted-dipole default, `:igrf` with a fixed `igrf_year`).
Produces zero wrench when `orientation_sim=false` — existing translational
configurations are unaffected.
"""
struct EddyCurrentDampingModel <: AbstractForceTorqueModel
    k_e::Float64
    field_model::Symbol
    igrf_year::Float64

    function EddyCurrentDampingModel(; k_e::Real, field_model::Symbol=:dipole, igrf_year::Real=NaN)
        isfinite(k_e) && k_e > 0.0 ||
            throw(ArgumentError("EddyCurrentDampingModel k_e must be finite and > 0, got $k_e"))
        field_model in (:dipole, :igrf) ||
            throw(ArgumentError("field_model must be :dipole or :igrf, got $(repr(field_model))"))
        if field_model === :igrf && !(isfinite(igrf_year) && 1900.0 <= igrf_year < 2035.0)
            throw(ArgumentError(
                "field_model=:igrf requires igrf_year in [1900, 2035) (decimal year, e.g. 2025.4)"))
        end
        if field_model === :igrf && igrf_year > 2030.0
            @warn "IGRF field accuracy is reduced for epochs past 2030 (igrf_year = $(Float64(igrf_year))); warning once at construction, per-evaluation library warnings are suppressed."
        end
        return new(Float64(k_e), field_model, Float64(igrf_year))
    end
end

"""
    eddy_damping_torque(k_e, B_body, ω_body)

Body-frame eddy-current damping torque `k_e · B × (B × ω)` [N·m].
"""
@inline function eddy_damping_torque(
    k_e::Float64, B_body::SVector{3, Float64}, ω_body::SVector{3, Float64}
)::SVector{3, Float64}
    return k_e .* cross(B_body, cross(B_body, ω_body))
end

@inline environment_requirements(::EddyCurrentDampingModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::EddyCurrentDampingModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if x.q_ib === nothing || x.ω_body === nothing
        return zero_wrench
    end
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("EddyCurrentDampingModel wrench requires env.planet_frame."))
    B_ii = _magnetic_field_inertial(
        model,
        planet_frame.l_pi,
        planet_frame.pos_pp,
        planet_frame.lat_rad,
        planet_frame.lon_rad,
        planet_frame.alt_m,
    )
    B_body = rot(x.q_ib) * B_ii
    return SVector{3, Float64}(0.0, 0.0, 0.0), eddy_damping_torque(model.k_e, B_body, x.ω_body)
end

"""
    MagneticTorqueRodModel <: AbstractForceTorqueModel

Dynamic effector for spacecraft magnetic torque rods / magnetorquers.

Sums the body-frame dipole moment `m` over every [`Magnet`](@ref) attached to
any link of the spacecraft, samples the Earth magnetic field at the
spacecraft's current position, and returns the resulting body-frame torque
`τ = m × B` ([`calculate_magnetic_torque`](@ref)). Produces zero force and
zero torque when the spacecraft carries no magnets or when attitude state is
unavailable (`orientation_sim=false`).

Field source (`field_model`):

- `:dipole` (default): fast tilted-dipole approximation
  ([`get_magnetic_field_dipole`](@ref)), epoch-2020 constants, typically good
  to 10-20% at LEO. Bit-identical to the pre-option behavior.
- `:igrf`: the International Geomagnetic Reference Field evaluated at the
  spacecraft's geodetic position. Requires `igrf_year` (decimal year, e.g.
  `2025.4`) — the secular field drift over a single simulation is negligible,
  so one epoch per model is deliberate and keeps the integrator hot loop free
  of calendar conversions.

Both sources are Earth models; neither is meaningful at another central body.
"""
struct MagneticTorqueRodModel <: AbstractForceTorqueModel
    field_model::Symbol
    igrf_year::Float64

    function MagneticTorqueRodModel(; field_model::Symbol=:dipole, igrf_year::Real=NaN)
        field_model in (:dipole, :igrf) ||
            throw(ArgumentError("field_model must be :dipole or :igrf, got $(repr(field_model))"))
        # SatelliteToolboxGeomagneticField's IGRF hard-rejects epochs outside
        # [1900, 2035), so reject unsupported epochs here at configuration time
        # instead of at the first wrench evaluation. The library also warns
        # about reduced accuracy past 2030 on every evaluation (no maxlog);
        # _magnetic_field_inertial suppresses that and we warn once here.
        if field_model === :igrf && !(isfinite(igrf_year) && 1900.0 <= igrf_year < 2035.0)
            throw(ArgumentError(
                "field_model=:igrf requires igrf_year in [1900, 2035) (decimal year, e.g. 2025.4)"))
        end
        if field_model === :igrf && igrf_year > 2030.0
            @warn "IGRF field accuracy is reduced for epochs past 2030 (igrf_year = $(Float64(igrf_year))); warning once at construction, per-evaluation library warnings are suppressed."
        end
        return new(field_model, Float64(igrf_year))
    end
end

"""
    _magnetic_field_inertial(model, l_pi, pos_pp, lat_rad, lon_rad, alt_m)

Inertial-frame magnetic field [Tesla] for the model's configured field source.
The IGRF branch converts from the library's nT to Tesla here so both branches
share one unit contract.
"""
@inline function _magnetic_field_inertial(
    model::Union{MagneticTorqueRodModel, EddyCurrentDampingModel},
    l_pi::SMatrix{3, 3, Float64, 9},
    pos_pp::SVector{3, Float64},
    lat_rad::Float64,
    lon_rad::Float64,
    alt_m::Float64,
)::SVector{3, Float64}
    if model.field_model === :igrf
        # show_warnings=Val(false): the library's reduced-accuracy warning for
        # epochs past 2030 has no maxlog and this runs once per RHS call;
        # the model constructors emit it once instead.
        B_ned_nT = igrf(model.igrf_year, alt_m, lat_rad, lon_rad, Val(:geodetic); show_warnings=Val(false))
        B_pp_nT = ned_to_ecef(B_ned_nT, lat_rad, lon_rad, alt_m)
        return SVector{3, Float64}(l_pi' * B_pp_nT) .* 1e-9
    end
    return get_magnetic_field_dipole(pos_pp, MMatrix{3, 3, Float64}(l_pi))
end

@inline function _total_dipole_moment_body(spacecraft)::SVector{3, Float64}
    m_total = SVector{3, Float64}(0.0, 0.0, 0.0)
    for link in spacecraft.links
        for magnet in link.magnets
            m_total = m_total + SVector{3, Float64}(magnet.m)
        end
    end
    return m_total
end

@inline environment_requirements(::MagneticTorqueRodModel) = EffectorEnvironmentRequirements(planet_frame=true)

"""
    rod_field_body(model::MagneticTorqueRodModel, pos_ii, q_ib, param) -> SVector{3}

The magnetic field (T) in the body frame at inertial position `pos_ii` and
attitude `q_ib`, from `model`'s field source, with the same arithmetic as the
model's torque. Mounted torque rods act in this field.
"""
function rod_field_body(model::MagneticTorqueRodModel, pos_ii::SVector{3, Float64}, q_ib::SVector{4, Float64}, param::ODEParams)::SVector{3, Float64}
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    planet = param.args.environment_model.planet
    l_pi = ephemerides_frame_lpi(param.args.environment_model.ephemerides_model, planet, et)
    pos_pp = l_pi * pos_ii
    if model.field_model === :igrf
        alt_m, lat_rad, lon_rad = rtolatlong(pos_pp, planet)
        B_ii = _magnetic_field_inertial(model, SMatrix{3, 3, Float64, 9}(l_pi), pos_pp, lat_rad, lon_rad, alt_m)
    else
        B_ii = get_magnetic_field_dipole(pos_pp, MMatrix{3, 3, Float64}(l_pi))
    end
    return rot(q_ib) * B_ii
end

@inline function wrench(
    model::MagneticTorqueRodModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if x.q_ib === nothing || x.spacecraft === nothing
        return zero_wrench
    end

    m_body = _total_dipole_moment_body(x.spacecraft)
    iszero(m_body) && return zero_wrench

    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("MagneticTorqueRodModel wrench requires env.planet_frame."))

    B_ii = _magnetic_field_inertial(
        model,
        planet_frame.l_pi,
        planet_frame.pos_pp,
        planet_frame.lat_rad,
        planet_frame.lon_rad,
        planet_frame.alt_m,
    )
    B_body = rot(x.q_ib) * B_ii
    torque_body = calculate_magnetic_torque(m_body, B_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque_body
end
