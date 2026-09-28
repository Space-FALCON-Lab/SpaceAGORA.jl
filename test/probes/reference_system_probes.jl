using Test
using Dates
using LinearAlgebra
using StaticArrays
using SPICE

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REFSYS_PATH = joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl")

function _ancestry_has_runtime_services(mod::Module)::Bool
    while true
        isdefined(mod, :RuntimeServices) && return true
        parent = parentmodule(mod)
        parent === mod && return false
        mod = parent
    end
end

# Probe the "RuntimeServices not found" error branch of _spice_lock. This is only
# reachable while no module in the ancestry defines RuntimeServices, i.e. before
# simulation_model.jl is included (it installs RuntimeServices in its parent).
module _SpiceLockErrorProbe end
Base.include(_SpiceLockErrorProbe, REFSYS_PATH)
const PROBE_ANCESTRY_HAS_RS = _ancestry_has_runtime_services(_SpiceLockErrorProbe)
const SPICE_LOCK_PROBE_RESULT = try
    _SpiceLockErrorProbe._spice_lock()
catch err
    err
end

using SpaceAGORA
const SimulationModel = SpaceAGORA.SimulationModel
using .SimulationModel
# Exercise the production owner; the isolated error-path sandbox above remains raw.
const RuntimeServices = SpaceAGORA.RuntimeServices
using .SimulationModel.FrameTransforms

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

_rz(θ::Float64) = SMatrix{3, 3, Float64}(
    cos(θ), -sin(θ), 0.0,
    sin(θ), cos(θ), 0.0,
    0.0, 0.0, 1.0
)' # rows [c s 0; -s c 0; 0 0 1]: J2000 -> body-fixed spin rotation

function _refsys_state_from_oe(a, e, i, Ω, ω, ν)
    R, V = orbitalelemtorv(SVector{7, Float64}(a, e, i, Ω, ω, ν, 0.0), EARTH)
    return SVector{3, Float64}(R), SVector{3, Float64}(V)
end

@testset "Reference System Coverage Probes" begin
    # ------------------------------------------------------------------
    # _spice_lock: success path and (fresh-process only) error path.
    # ------------------------------------------------------------------
    @test _spice_lock() === RuntimeServices.SPICE_LOCK
    if PROBE_ANCESTRY_HAS_RS
        # Another suite already installed RuntimeServices before this file ran;
        # the probe then resolves the lock instead of erroring.
        @test SPICE_LOCK_PROBE_RESULT === RuntimeServices.SPICE_LOCK
    else
        @test SPICE_LOCK_PROBE_RESULT isa ErrorException
        @test occursin("SPICE_LOCK not found", SPICE_LOCK_PROBE_RESULT.msg)
    end

    # ------------------------------------------------------------------
    # Small pure helpers.
    # ------------------------------------------------------------------
    @test _wrap_2pi(7.0) ≈ 7.0 - 2pi
    @test _wrap_2pi(-0.5) ≈ 2pi - 0.5
    @test _wrap_2pi(1.25) == 1.25
    @test _safe_acos(1.0 + 1e-12) == 0.0
    @test _safe_acos(-1.0 - 1e-12) ≈ pi
    @test _safe_acos(0.0) ≈ pi / 2
    @test _planet_flattening(EARTH) ≈ (EARTH.Rp_e - EARTH.Rp_p) / EARTH.Rp_e
    @test 0.0 < _planet_flattening(EARTH) < 0.01

    @test _spice_body_fixed_frame((name = "Moon",)) == "MOON_PA_DE421"
    @test _spice_body_fixed_frame((name = "Earth",)) == "ITRF93"
    @test _spice_body_fixed_frame((name = "Mars",)) == "IAU_MARS"

    # ------------------------------------------------------------------
    # Legacy (no ephemeris time) frame conversions using planet.L_PI.
    # ------------------------------------------------------------------
    EARTH.L_PI .= _rz(0.4)
    r_i = SVector{3, Float64}(EARTH.Rp_e + 500e3, 1.0e5, 2.0e5)
    v_i = SVector{3, Float64}(10.0, 7.6e3, 12.0)

    r_p_legacy, v_p_legacy = r_intor_p!(r_i, v_i, EARTH)
    @test norm(r_p_legacy) ≈ norm(r_i) rtol = 1e-12
    @test r_p_legacy ≈ SVector{3, Float64}(EARTH.L_PI * r_i) rtol = 1e-12
    r_i_back, v_i_back = r_pintor_i(r_p_legacy, v_p_legacy, EARTH)
    @test r_i_back ≈ r_i rtol = 1e-12
    @test v_i_back ≈ v_i rtol = 1e-12

    # ------------------------------------------------------------------
    # SPICE-backed transforms: Earth high-precision path (ITRF93 available).
    # ------------------------------------------------------------------
    et = 0.0 # J2000 epoch, inside earth_latest_high_prec.bpc coverage
    r_p_hp, v_p_hp = r_intor_p!(r_i, v_i, EARTH, et)
    itrf = SMatrix{3, 3, Float64}(pxform("J2000", "ITRF93", et))
    @test r_p_hp ≈ itrf * r_i rtol = 1e-12
    @test norm(r_p_hp) ≈ norm(r_i) rtol = 1e-12
    r_i_hp, v_i_hp = r_pintor_i(r_p_hp, v_p_hp, EARTH, et)
    @test r_i_hp ≈ r_i rtol = 1e-9
    @test v_i_hp ≈ v_i rtol = 1e-9

    # Non-Earth branch (IAU_<NAME> frame from the text PCK).
    mars_nt = (name = "Mars",)
    state_mars = _j2000_to_body_fixed_state(r_i, v_i, mars_nt, et)
    iau_mars = SMatrix{3, 3, Float64}(pxform("J2000", "IAU_MARS", et))
    @test SVector{3, Float64}(state_mars[1], state_mars[2], state_mars[3]) ≈ iau_mars * r_i rtol = 1e-12
    state_mars_back = _body_fixed_to_j2000_state(
        SVector{3, Float64}(state_mars[1], state_mars[2], state_mars[3]),
        SVector{3, Float64}(state_mars[4], state_mars[5], state_mars[6]),
        mars_nt,
        et
    )
    @test SVector{3, Float64}(state_mars_back[1], state_mars_back[2], state_mars_back[3]) ≈ r_i rtol = 1e-9
    @test SVector{3, Float64}(state_mars_back[4], state_mars_back[5], state_mars_back[6]) ≈ v_i rtol = 1e-9

    # Moon special-case frame: MOON_PA_DE421 has no loaded frame kernel/orientation
    # data in the vendored kernel set, so the (non-Earth) SPICE call must throw.
    moon_nt = (name = "Moon",)
    @test_throws SPICE.SpiceError _j2000_to_body_fixed_state(r_i, v_i, moon_nt, et)
    @test_throws SPICE.SpiceError _body_fixed_to_j2000_state(r_i, v_i, moon_nt, et)

    # ------------------------------------------------------------------
    # Earth fallback paths: strip kernels so ITRF93 (and then everything) fails.
    # ------------------------------------------------------------------
    kclear()
    # kclear() wipes CSPICE's kernel pool but not the furnished-kernel /
    # planet-instance caches in src/environment/ephemerides/planets.jl; without
    # this, the Earth("", SPICE_PATH) restore below would hit the cache from
    # the `const EARTH = ...` call above and skip re-furnishing entirely.
    SimulationModel.Planets._reset_furnished_kernels!()
    # No kernels at all: high-precision path throws, fallback throws too, and the
    # error from the fallback propagates out of the catch block.
    @test_throws SPICE.SpiceError r_intor_p!(r_i, v_i, EARTH, et)
    @test_throws SPICE.SpiceError r_pintor_i(r_p_hp, v_p_hp, EARTH, et)

    # Text PCK + LSK only: ITRF93 still unavailable -> IAU_EARTH fallback engages.
    furnsh(joinpath(SPICE_PATH, "pck", "pck00011.tpc"))
    furnsh(joinpath(SPICE_PATH, "lsk", "naif0012.tls"))
    @test_throws SPICE.SpiceError pxform("J2000", "ITRF93", et)
    r_p_fb, v_p_fb = r_intor_p!(r_i, v_i, EARTH, et)
    iau_earth = SMatrix{3, 3, Float64}(pxform("J2000", "IAU_EARTH", et))
    @test r_p_fb ≈ iau_earth * r_i rtol = 1e-12
    @test r_p_fb != r_p_hp # fallback frame differs from the high-precision one
    r_i_fb, v_i_fb = r_pintor_i(r_p_fb, v_p_fb, EARTH, et)
    @test r_i_fb ≈ r_i rtol = 1e-9
    @test v_i_fb ≈ v_i rtol = 1e-9

    # Restore the full kernel set for the remainder of the process.
    kclear()
    SimulationModel.Planets._reset_furnished_kernels!()
    Earth("", SPICE_PATH)
    @test SMatrix{3, 3, Float64}(pxform("J2000", "ITRF93", et)) ≈ itrf rtol = 1e-14

    # ------------------------------------------------------------------
    # Ephemerides-model dispatch of r_intor_p!.
    # ------------------------------------------------------------------
    spice_model = SimulationModel.SpiceEphemeridesModel()
    r_p_spice, v_p_spice = r_intor_p!(r_i, v_i, EARTH, et, spice_model)
    @test r_p_spice ≈ r_p_hp rtol = 1e-12
    @test v_p_spice ≈ v_p_hp rtol = 1e-12

    simple_model = SimulationModel.SimpleEphemeridesModel(
        reference_epoch_seconds = 0.0,
        prime_meridian_at_reference_rad = 0.3
    )
    et_simple = 250.0
    r_p_simple, v_p_simple = r_intor_p!(r_i, v_i, EARTH, et_simple, simple_model)
    l_pi = SimulationModel.planet_frame_lpi(EARTH, et_simple, simple_model)
    @test r_p_simple ≈ SVector{3, Float64}(l_pi * r_i) rtol = 1e-12
    @test v_p_simple ≈ SVector{3, Float64}(l_pi * (v_i - cross(EARTH.ω, r_i))) rtol = 1e-12

    # ------------------------------------------------------------------
    # Orbital element <-> Cartesian conversions and singular-case branches.
    # ------------------------------------------------------------------
    a0, e0, i0, Ω0, ω0, ν0 = 7.0e6, 0.1, 0.5, 1.0, 2.0, 1.0
    R0, V0 = orbitalelemtorv(SVector{7, Float64}(a0, e0, i0, Ω0, ω0, ν0, 0.0), EARTH)
    p0 = a0 * (1 - e0^2)
    @test norm(R0) ≈ p0 / (1 + e0 * cos(ν0)) rtol = 1e-12
    @test dot(V0, V0) / 2 - EARTH.μ / norm(R0) ≈ -EARTH.μ / (2 * a0) rtol = 1e-12
    @test norm(cross(SVector{3}(R0), SVector{3}(V0))) ≈ sqrt(EARTH.μ * p0) rtol = 1e-12

    # Generic struct/NamedTuple overload dispatches to the SVector method.
    R0_nt, V0_nt = orbitalelemtorv((a = a0, e = e0, i = i0, Ω = Ω0, ω = ω0, ν = ν0), EARTH)
    @test R0_nt ≈ R0 atol = 0.0
    @test V0_nt ≈ V0 atol = 0.0

    # Roundtrip: generic inclined elliptical, ascending (dot(r,v) > 0).
    oe_rt = rvtoorbitalelement(SVector{3, Float64}(R0), SVector{3, Float64}(V0), 123.5, EARTH)
    @test oe_rt[1] ≈ a0 rtol = 1e-9
    @test oe_rt[2] ≈ e0 rtol = 1e-9
    @test oe_rt[3] ≈ i0 rtol = 1e-9
    @test oe_rt[4] ≈ Ω0 rtol = 1e-9
    @test oe_rt[5] ≈ ω0 rtol = 1e-9
    @test oe_rt[6] ≈ ν0 rtol = 1e-9
    @test oe_rt[7] == 123.5

    # 6-element method + descending branch (ν > π so dot(r,v) < 0)
    # + e_vec[3] < 0 branch (sin(ω) < 0 for ω = 4.0).
    r_desc, v_desc = _refsys_state_from_oe(7.2e6, 0.15, 0.9, 2.5, 4.0, 4.2)
    oe_desc = rvtoorbitalelement(r_desc, v_desc, EARTH)
    @test length(oe_desc) == 6
    @test oe_desc[1] ≈ 7.2e6 rtol = 1e-9
    @test oe_desc[2] ≈ 0.15 rtol = 1e-9
    @test oe_desc[5] ≈ 4.0 rtol = 1e-9
    @test oe_desc[6] ≈ 4.2 rtol = 1e-9

    # Equatorial elliptical: Ω forced to 0, ω recovered as longitude of periapsis.
    r_eq, v_eq = _refsys_state_from_oe(7.0e6, 0.05, 0.0, 0.3, 0.4, 1.0)
    oe_eq = rvtoorbitalelement(r_eq, v_eq, 0.0, EARTH)
    @test oe_eq[3] ≈ 0.0 atol = 1e-9
    @test oe_eq[4] == 0.0
    @test oe_eq[5] ≈ 0.7 rtol = 1e-9 # Ω + ω collapses into ω
    @test oe_eq[6] ≈ 1.0 rtol = 1e-9

    # Circular inclined: ω forced to 0, ν is the argument of latitude.
    r_ci, v_ci = _refsys_state_from_oe(6.9e6, 0.0, 0.6, 1.0, 0.0, 0.8)
    oe_ci = rvtoorbitalelement(r_ci, v_ci, 0.0, EARTH)
    @test oe_ci[2] ≈ 0.0 atol = 1e-11
    @test oe_ci[5] == 0.0
    @test oe_ci[4] ≈ 1.0 rtol = 1e-9
    @test oe_ci[6] ≈ 0.8 rtol = 1e-8

    # Circular inclined, below the equator (r_z < 0 -> ν = 2π - acos(...)).
    r_cib, v_cib = _refsys_state_from_oe(6.9e6, 0.0, 0.6, 1.0, 0.0, 4.0)
    oe_cib = rvtoorbitalelement(r_cib, v_cib, 0.0, EARTH)
    @test oe_cib[6] ≈ 4.0 rtol = 1e-8

    # Circular equatorial: ν is the true longitude.
    r_ce, v_ce = _refsys_state_from_oe(7.1e6, 0.0, 0.0, 0.0, 0.0, 1.2)
    oe_ce = rvtoorbitalelement(r_ce, v_ce, 0.0, EARTH)
    @test oe_ce[2] ≈ 0.0 atol = 1e-11
    @test oe_ce[3] ≈ 0.0 atol = 1e-9
    @test oe_ce[4] == 0.0
    @test oe_ce[5] == 0.0
    @test oe_ce[6] ≈ 1.2 rtol = 1e-8

    # ------------------------------------------------------------------
    # Spherical / geodetic coordinate helpers.
    # ------------------------------------------------------------------
    rad1 = rtoalfadeltar([1.0, 1.0, 1.0]) # m > 0 branch
    @test rad1[1] ≈ sqrt(3.0) rtol = 1e-12
    @test rad1[2] ≈ pi / 4 rtol = 1e-12
    @test rad1[3] ≈ asin(1 / sqrt(3.0)) rtol = 1e-12
    rad2 = rtoalfadeltar([1.0, -1.0, 0.5]) # m <= 0 branch: RA = 2π - acos(...)
    @test rad2[2] > pi
    @test alfadeltartor(rad2) ≈ [1.0, -1.0, 0.5] rtol = 1e-12
    @test alfadeltartor(rad1) ≈ [1.0, 1.0, 1.0] rtol = 1e-12

    # latlongtor: at ϕ = 0, h = 0, α_g0 = 0, t = t0 the point sits on the
    # equatorial x-axis at one equatorial radius.
    r_eq0 = latlongtor([0.0, 0.0, 0.0], EARTH, 0.0, 0.0, 0.0)
    @test r_eq0 ≈ [EARTH.Rp_e, 0.0, 0.0] atol = 1e-6
    # Planet rotation advances the inertial longitude by ω3 * Δt.
    Δt = 1000.0
    r_eq_rot = latlongtor([0.0, 0.0, 0.0], EARTH, 0.0, Δt, 0.0)
    α_rot = EARTH.ω[3] * Δt
    @test r_eq_rot ≈ [EARTH.Rp_e * cos(α_rot), EARTH.Rp_e * sin(α_rot), 0.0] atol = 1e-6
    # Nonzero latitude exercises the trig terms; norm stays physical.
    r_lat = latlongtor([0.7, 1.1, 5000.0], EARTH, 0.2, 10.0, 0.0)
    @test all(isfinite, r_lat)
    @test EARTH.Rp_p < norm(r_lat) < 1.02 * EARTH.Rp_e
    # Geodetic construction must match the explicit prime-vertical formulas
    # (the same ones latlongtoOE uses): N = a/sqrt(1 - e2 sin^2 ϕ) and a z
    # component of ((1 - e2) N + h) sin ϕ.
    ϕ_g, λ_g, h_g = 0.7, 1.1, 5000.0
    e2_g = 1 - (EARTH.Rp_p / EARTH.Rp_e)^2
    N_g = EARTH.Rp_e / sqrt(1 - e2_g * sin(ϕ_g)^2)
    r_geo = latlongtor([ϕ_g, λ_g, h_g], EARTH, 0.0, 0.0, 0.0)
    @test r_geo ≈ [
        (N_g + h_g) * cos(ϕ_g) * cos(λ_g),
        (N_g + h_g) * cos(ϕ_g) * sin(λ_g),
        ((1 - e2_g) * N_g + h_g) * sin(ϕ_g)
    ] rtol = 1e-12
    # Bowring roundtrip recovers the geodetic lat/lon (PCI == PCPF at t = t0,
    # α_g0 = 0).
    lla_rt = rtolatlong(SVector{3, Float64}(r_geo), EARTH)
    @test lla_rt[2] ≈ ϕ_g atol = 1e-9
    @test lla_rt[3] ≈ λ_g atol = 1e-12
    # The ALTITUDE must round-trip too: the closed form previously used
    # e²·N·sin³φ instead of e²·N·sin²φ, costing ~6 km at mid-latitudes.
    @test lla_rt[1] ≈ h_g atol = 1e-3
    for (ϕ_chk, h_chk) in ((deg2rad(45.0), 125e3), (deg2rad(-60.0), 250e3), (deg2rad(20.0), 90e3))
        r_chk = latlongtor([ϕ_chk, 0.9, h_chk], EARTH, 0.0, 0.0, 0.0)
        lla_chk = rtolatlong(SVector{3, Float64}(r_chk), EARTH)
        @test lla_chk[1] ≈ h_chk atol = 1e-3
        @test lla_chk[2] ≈ ϕ_chk atol = 1e-9
    end

    # latlongtoOE (prints diagnostics; silence them). L_PI was set to a rotation
    # above, so the legacy PCPF->J2000 conversion is well-defined.
    ϕ_oe, λ_oe, h_oe, v_mag = 0.3, 0.5, 400e3, 7.7e3
    OE_geo = redirect_stdout(devnull) do
        latlongtoOE([ϕ_oe, λ_oe, h_oe], EARTH, 0.0, pi / 2, v_mag)
    end
    @test length(OE_geo) == 6
    @test all(isfinite, OE_geo)
    # Geodetic radius of the launch point (same construction as latlongtor).
    N_oe = EARTH.Rp_e / sqrt(1 - e2_g * sin(ϕ_oe)^2)
    r_geo_oe = norm([
        (N_oe + h_oe) * cos(ϕ_oe) * cos(λ_oe),
        (N_oe + h_oe) * cos(ϕ_oe) * sin(λ_oe),
        ((1 - e2_g) * N_oe + h_oe) * sin(ϕ_oe)
    ])
    # Vis-viva: semi-major axis from the point radius and speed.
    @test OE_geo[1] ≈ -EARTH.μ / (2 * (v_mag^2 / 2 - EARTH.μ / r_geo_oe)) rtol = 1e-9
    # γ = 0 with α = π/2 gives dot(r, v) = 0, so the state sits at an apsis
    # and e = |v^2 r / μ - 1|; v exceeds circular speed, so it is periapsis.
    @test OE_geo[2] ≈ abs(v_mag^2 * r_geo_oe / EARTH.μ - 1.0) rtol = 1e-9
    @test min(OE_geo[6], 2pi - OE_geo[6]) ≈ 0.0 atol = 1e-6
    # Roundtrip: the elements reproduce a Cartesian state with the original
    # radius/speed, and rvtoorbitalelement/orbitalelemtorv self-agree on it.
    R_geo_oe, V_geo_oe = orbitalelemtorv(SVector{7, Float64}(vcat(OE_geo, 0.0)), EARTH)
    @test norm(R_geo_oe) ≈ r_geo_oe rtol = 1e-9
    @test norm(V_geo_oe) ≈ v_mag rtol = 1e-9
    @test abs(dot(R_geo_oe, V_geo_oe)) < 1e-3 * r_geo_oe
    oe_back = rvtoorbitalelement(SVector{3, Float64}(R_geo_oe), SVector{3, Float64}(V_geo_oe), EARTH)
    R_back, V_back = orbitalelemtorv(SVector{7, Float64}(vcat(oe_back, 0.0)), EARTH)
    @test R_back ≈ R_geo_oe rtol = 1e-9
    @test V_back ≈ V_geo_oe rtol = 1e-9

    # rtolatlong: equator and pole, plus default-argument method.
    lla_eq = rtolatlong(SVector{3, Float64}(EARTH.Rp_e, 0.0, 0.0), EARTH)
    @test lla_eq[1] ≈ 0.0 atol = 1e-6
    @test lla_eq[2] ≈ 0.0 atol = 1e-12
    @test lla_eq[3] ≈ 0.0 atol = 1e-12
    lla_pole = rtolatlong(SVector{3, Float64}(0.0, 0.0, EARTH.Rp_p + 1000.0), EARTH, false)
    @test lla_pole[1] ≈ 1000.0 atol = 1e-6
    @test lla_pole[2] ≈ pi / 2 atol = 1e-9
    lla_west = rtolatlong(SVector{3, Float64}(EARTH.Rp_e / sqrt(2), -EARTH.Rp_e / sqrt(2), 0.0), EARTH)
    @test lla_west[3] ≈ -pi / 4 atol = 1e-12

    # Ephemerides-model overload just forwards to the two-argument method.
    lla_fwd = rtolatlong(SVector{3, Float64}(EARTH.Rp_e, 0.0, 0.0), EARTH, simple_model)
    @test lla_fwd == lla_eq

    # Topography branch: legacy duck-typed contract (mirrors the suite-03
    # sandbox test) — a planet exposing a 6-arg topography_function plus
    # Clm_topo/Slm_topo/A_topo gets altitude = |r| − elevation. Real Planet
    # structs lack those fields, so the branch throws a field error for them
    # instead of the pre-fix UndefVarError on the includer `args` global.
    planet_topo_probe = (
        Rp_e=10.0,
        Rp_p=9.0,
        topography_function=(a, Clm, Slm, lat, lon, A) -> 9.5,
        Clm_topo=zeros(1, 1),
        Slm_topo=zeros(1, 1),
        A_topo=0.0
    )
    rp_topo_probe = SVector{3, Float64}(0.1, 0.5, -8.944723618090451)
    lla_topo_probe = rtolatlong(rp_topo_probe, planet_topo_probe, true)
    @test isapprox(lla_topo_probe[1], norm(rp_topo_probe) - 9.5; atol=1e-12, rtol=0.0)
    @test_throws Exception rtolatlong(SVector{3, Float64}(EARTH.Rp_e, 0.0, 0.0), EARTH, true)

    # rtolatlongrad: planetocentric radius/lat/lon.
    rll = rtolatlongrad(SVector{3, Float64}(EARTH.Rp_e, 0.0, 0.0), EARTH)
    @test rll[1] ≈ EARTH.Rp_e rtol = 1e-12
    @test rll[2] ≈ 0.0 atol = 1e-12
    @test rll[3] ≈ 0.0 atol = 1e-12
    rll2 = rtolatlongrad(SVector{3, Float64}(1.0e6, -1.0e6, sqrt(2.0) * 1.0e6), EARTH)
    @test rll2[1] ≈ 2.0e6 rtol = 1e-12
    @test rll2[2] ≈ pi / 4 rtol = 1e-12
    @test rll2[3] ≈ -pi / 4 rtol = 1e-12

    # latlongtoNED: orthonormal right-handed triad; canonical values at (0, 0).
    uD, uN, uE = latlongtoNED([0.0, 0.0, 0.0])
    @test uD ≈ SVector{3, Float64}(-1.0, 0.0, 0.0) atol = 1e-14
    @test uN ≈ SVector{3, Float64}(0.0, 0.0, 1.0) atol = 1e-14
    @test uE ≈ SVector{3, Float64}(0.0, 1.0, 0.0) atol = 1e-14
    uD2, uN2, uE2 = latlongtoNED([0.0, 0.7, -1.2])
    for u in (uD2, uN2, uE2)
        @test norm(u) ≈ 1.0 rtol = 1e-12
    end
    @test abs(dot(uD2, uN2)) < 1e-12
    @test abs(dot(uD2, uE2)) < 1e-12
    @test abs(dot(uN2, uE2)) < 1e-12
    @test cross(uN2, uE2) ≈ uD2 rtol = 1e-12

    # ------------------------------------------------------------------
    # LVLH quaternion (all four DCM-to-quaternion extraction branches) and
    # quaternion vector rotation.
    # ------------------------------------------------------------------
    for raan in 0.0:0.9:5.5, inc in (0.05, 1.2, 2.2, 3.05), u in 0.0:1.1:6.0
        aop = 0.6 * u
        nu = u - aop
        q = orbital_elements_to_lvlh_quaternion(raan, inc, aop, nu)
        @test norm(q) ≈ 1.0 rtol = 1e-10

        su, cu = sincos(u)
        si, ci = sincos(inc)
        sO, cO = sincos(raan)
        r_hat = [cO * cu - sO * su * ci, sO * cu + cO * su * ci, su * si]
        h_hat = [sO * si, -cO * si, ci]

        # q rotates inertial vectors into LVLH: nadir is -Z, momentum is -Y.
        @test rotate_vector_by_quaternion(r_hat, Vector{Float64}(q)) ≈ [0.0, 0.0, -1.0] atol = 1e-9
        @test rotate_vector_by_quaternion(h_hat, Vector{Float64}(q)) ≈ [0.0, -1.0, 0.0] atol = 1e-9
    end

    # Identity rotation leaves vectors unchanged.
    @test rotate_vector_by_quaternion([1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 1.0]) ≈ [1.0, 2.0, 3.0] atol = 1e-14

    # ------------------------------------------------------------------
    # RTN / HCW utilities.
    # ------------------------------------------------------------------
    r_t = SVector{3, Float64}(7.0e6, 0.0, 0.0)
    v_t = SVector{3, Float64}(0.0, 7.5e3, 0.0)
    C = rtn_dcm_from_inertial(r_t, v_t)
    @test C[:, 1] ≈ SVector{3, Float64}(1.0, 0.0, 0.0) atol = 1e-14
    @test C[:, 2] ≈ SVector{3, Float64}(0.0, 1.0, 0.0) atol = 1e-14
    @test C[:, 3] ≈ SVector{3, Float64}(0.0, 0.0, 1.0) atol = 1e-14
    @test C' * C ≈ SMatrix{3, 3, Float64}(I) atol = 1e-12
    @test det(C) ≈ 1.0 rtol = 1e-12

    @test_throws ArgumentError rtn_dcm_from_inertial(SVector{3, Float64}(0.0, 0.0, 0.0), v_t)
    # Parallel r and v: zero angular momentum.
    @test_throws ArgumentError rtn_dcm_from_inertial(r_t, 2.0 * r_t)

    @test _rtn_rate_rad_s(r_t, v_t) ≈ 7.5e3 / 7.0e6 rtol = 1e-12
    @test_throws ArgumentError _rtn_rate_rad_s(SVector{3, Float64}(0.0, 0.0, 0.0), v_t)

    # Relative-state roundtrip: inertial -> RTN -> inertial.
    r_c = r_t + SVector{3, Float64}(120.0, -80.0, 45.0)
    v_c = v_t + SVector{3, Float64}(-0.4, 0.9, 0.2)
    rel = inertial_to_rtn_relative_state(r_c, v_c, r_t, v_t)
    @test rel[1] ≈ 120.0 atol = 1e-9 # radial offset in this axis-aligned geometry
    @test rel[2] ≈ -80.0 atol = 1e-9
    @test rel[3] ≈ 45.0 atol = 1e-9
    r_c_back, v_c_back = rtn_to_inertial_relative_state(
        SVector{3, Float64}(rel[1], rel[2], rel[3]),
        SVector{3, Float64}(rel[4], rel[5], rel[6]),
        r_t,
        v_t
    )
    @test r_c_back ≈ r_c rtol = 1e-12
    @test v_c_back ≈ v_c rtol = 1e-12

    # Co-located chaser has an exactly zero relative state.
    rel0 = inertial_to_rtn_relative_state(r_t, v_t, r_t, v_t)
    @test rel0 == SVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Acceleration rotation: pure radial command maps onto r_hat.
    a_inertial = rtn_accel_to_inertial(SVector{3, Float64}(2.0, 0.0, 0.0), r_t, v_t)
    @test a_inertial ≈ 2.0 * r_t / norm(r_t) atol = 1e-12
end

# ----------------------------------------------------------------------
# SimpleEphemeridesModel Earth rotation: pinned against an independent GMST
# computation and cross-checked against the SPICE-backed frame. Guards the
# geographic registration of every lat/lon-keyed model (IGRF, tilted dipole,
# atmospheres, tesseral harmonics) evaluated under the kernel-free backend —
# a zero prime meridian at J2000 would park Greenwich on the vernal equinox,
# ~280 deg away from where it belongs.
# ----------------------------------------------------------------------
@testset "SimpleEphemeridesModel Earth rotation pins" begin
    # Independent GMST formulation (IAU-82 in degree form — not a transcription
    # of the seconds-form series in simple_ephemerides.jl).
    function gmst82_deg(jd_ut1::Float64)
        d = jd_ut1 - 2451545.0
        T = d / 36525.0
        return mod(280.46061837 + 360.98564736629 * d + 3.87933e-4 * T^2 - T^3 / 3.871e7, 360.0)
    end
    zangle_deg(R) = rad2deg(atan(R[1, 2], R[1, 1]))  # rows [c s 0; -s c 0; 0 0 1]
    rotangle_deg(R) = rad2deg(acos(clamp((tr(R) - 1) / 2, -1.0, 1.0)))
    wrapdiff_deg(a, b) = mod(a - b + 180.0, 360.0) - 180.0

    simple_default = SimulationModel.SimpleEphemeridesModel()
    spice_backend = SimulationModel.SpiceEphemeridesModel()
    j2000_utc = DateTime(2000, 1, 1, 12, 0, 0)

    for utc in (DateTime(2000, 1, 1, 12), DateTime(2007, 6, 15, 3, 30),
                DateTime(2015, 2, 24), DateTime(2024, 12, 31, 23, 59))
        # The simple backend's timeline: leap-second-free UTC seconds past J2000.
        et_simple = Dates.value(utc - j2000_utc) / 1000.0
        l_pi = SimulationModel.planet_frame_lpi(EARTH, et_simple, simple_default)

        # Pure spin about +z: exact invariant z block, orthonormal to roundoff.
        @test l_pi[3, 3] == 1.0 && l_pi[1, 3] == 0.0 && l_pi[2, 3] == 0.0
        @test l_pi' * l_pi ≈ SMatrix{3, 3, Float64}(I) atol = 1e-14

        # Pin against the independent GMST computation (same UT1 ≈ UTC reading).
        θ_expected = gmst82_deg(2451545.0 + et_simple / 86400.0)
        @test wrapdiff_deg(zangle_deg(l_pi), θ_expected) ≈ 0.0 atol = 1e-5

        # Cross-check against the SPICE frame chain: GMST vs the full
        # precession/nutation/polar-motion transform stays sub-half-degree
        # across the vendored kernels' 2000-2025 span.
        et_spice = utc2et(Dates.format(utc, "yyyy-mm-ddTHH:MM:SS"))
        l_pi_spice = SimulationModel.planet_frame_lpi(EARTH, et_spice, spice_backend)
        @test rotangle_deg(l_pi * l_pi_spice') < 0.5
    end

    # Explicit prime-meridian override keeps the legacy linear model exactly.
    override = SimulationModel.SimpleEphemeridesModel(
        reference_epoch_seconds = 100.0,
        prime_meridian_at_reference_rad = 0.3
    )
    l_pi_o = SimulationModel.planet_frame_lpi(EARTH, 350.0, override)
    @test zangle_deg(l_pi_o) ≈ rad2deg(0.3 + EARTH.ω[3] * 250.0) rtol = 1e-12

    # Non-Earth planets keep the zero-at-reference default convention.
    mars_like = (name = "Mars", ω = SVector{3, Float64}(0.0, 0.0, 7.08823596e-5))
    l_pi_m = SimulationModel.planet_frame_lpi(mars_like, 1.0e4, simple_default)
    @test zangle_deg(l_pi_m) ≈ rad2deg(7.08823596e-5 * 1.0e4) rtol = 1e-12

    # The NaN planet-true default must key the ephemerides cache without
    # throwing, and distinctly from an explicit zero override.
    key_default = SimulationModel.ephemerides_cache_key(simple_default)
    key_zero = SimulationModel.ephemerides_cache_key(
        SimulationModel.SimpleEphemeridesModel(prime_meridian_at_reference_rad = 0.0))
    @test key_default != key_zero
end

# Keep this last: the earlier error-path checks deliberately have no lunar
# orientation kernels loaded. Suite 09 already runs this file in its own process.
@testset "Planet constructors preserve a compatible shared PCK pool" begin
    lock(RuntimeServices.SPICE_LOCK) do
        planets_owner = SimulationModel.Planets
        reset_pool!() = (kclear(); planets_owner._reset_furnished_kernels!())
        pck10 = joinpath(SPICE_PATH, "pck", "pck00010.tpc")
        pck11 = joinpath(SPICE_PATH, "pck", "pck00011.tpc")
        lsk = joinpath(SPICE_PATH, "lsk", "naif0012.tls")
        try
            # Independent old-Titan baseline, before loading any newer PCK.
            reset_pool!()
            furnsh(pck10)
            furnsh(lsk)
            titan_epochs = (utc2et("2004-10-26T15:30:04.608"),
                            utc2et("2005-04-16T19:11:45.772"))
            titan_frames = [pxform("J2000", "IAU_TITAN", et) for et in titan_epochs]
            titan_radii = bodvrd("TITAN", "RADII") .* 1e3
            moon_radii = bodvrd("MOON", "RADII") .* 1e3
            furnsh(joinpath(SPICE_PATH, "spk", "satellites", "SPICELunaCurrentKernel.bpc"))
            furnsh(joinpath(SPICE_PATH, "tf", "SPICELunaFrameKernel.tf"))
            moon_epochs = (0.0, 86400.0)
            moon_frames = Dict(frame => [pxform("J2000", frame, et) for et in moon_epochs]
                               for frame in ("IAU_MOON", "MOON_PA_DE421"))

            # Mars needs an internally consistent newer PCK. Loading 00010 after
            # 00011 used to overwrite only part of its nutation/precession model.
            reset_pool!()
            furnsh(pck11)
            furnsh(lsk)
            mars_et = utc2et("2001-11-06T18:44:35.434")
            mars_epochs = (mars_et - 180.0, mars_et, mars_et + 180.0)
            mars_frames = [pxform("J2000", "IAU_MARS", et) for et in mars_epochs]

            constructors = (mars=Mars, titan=Titan, moon=Moon)
            for order in ((:mars, :titan, :moon), (:titan, :mars, :moon),
                          (:moon, :titan, :mars))
                @testset "Constructor order $order" begin
                    reset_pool!()
                    made = Dict{Symbol, Any}()
                    for name in order
                        made[name] = getproperty(constructors, name)("", SPICE_PATH)
                        # Check immediately after every constructor. A later
                        # constructor must not conceal a transient broken pool.
                        if haskey(made, :mars)
                            @test Mars("", SPICE_PATH) === made[:mars]
                            for (et, expected) in zip(mars_epochs, mars_frames)
                                @test pxform("J2000", "IAU_MARS", et) == expected
                                @test SimulationModel.planet_frame_lpi(
                                    made[:mars], et, SimulationModel.SpiceEphemeridesModel()
                                ) == expected
                            end
                        end
                        if haskey(made, :titan)
                            @test Titan("", SPICE_PATH) === made[:titan]
                            for (et, expected) in zip(titan_epochs, titan_frames)
                                @test pxform("J2000", "IAU_TITAN", et) == expected
                            end
                        end
                        if haskey(made, :moon)
                            @test Moon("", SPICE_PATH) === made[:moon]
                            for frame in ("IAU_MOON", "MOON_PA_DE421")
                                for (et, expected) in zip(moon_epochs, moon_frames[frame])
                                    @test pxform("J2000", frame, et) == expected
                                end
                            end
                            for (et, expected) in zip(moon_epochs, moon_frames["MOON_PA_DE421"])
                                @test SimulationModel.planet_frame_lpi(
                                    made[:moon], et, SimulationModel.SpiceEphemeridesModel()
                                ) == expected
                            end
                        end
                    end
                    @test made[:mars].μ == 0.4282837285418775e5 * 1e9
                    @test made[:titan].Rp_e == titan_radii[1]
                    @test made[:titan].Rp_p == titan_radii[3]
                    @test made[:titan].Rp_m == sum(titan_radii) / 3.0
                    @test made[:moon].Rp_e == moon_radii[1]
                    @test made[:moon].Rp_p == moon_radii[3]
                    @test abspath(pck11) in planets_owner._FURNISHED_KERNELS
                    @test !(abspath(pck10) in planets_owner._FURNISHED_KERNELS)
                    kernel_count = ktotal("ALL")
                    for _ in 1:3, name in (:mars, :titan, :moon)
                        @test getproperty(constructors, name)("", SPICE_PATH) === made[name]
                    end
                    @test ktotal("ALL") == kernel_count
                end
            end

            # A 00010-only bundle must fail at the required PCK, before changing
            # the pool or publishing a cached Titan instance.
            mktempdir() do old_bundle
                mkpath(joinpath(old_bundle, "pck"))
                cp(pck10, joinpath(old_bundle, "pck", "pck00010.tpc"))
                kernel_count = ktotal("ALL")
                furnished = copy(planets_owner._FURNISHED_KERNELS)
                err = try
                    Titan("", old_bundle)
                    nothing
                catch failure
                    failure
                end
                @test err isa ArgumentError
                @test occursin("Required SPICE kernel not found", sprint(showerror, err))
                @test occursin(joinpath(old_bundle, "pck", "pck00011.tpc"), sprint(showerror, err))
                @test !haskey(planets_owner._TITAN_CACHE, ("", old_bundle))
                @test ktotal("ALL") == kernel_count
                @test planets_owner._FURNISHED_KERNELS == furnished
            end
        finally
            reset_pool!()
            Earth("", SPICE_PATH)
        end
    end
end
