using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SM = SpaceAGORA.SimulationModel
const QM = SM.QuaternionMath
const FT = SM.FrameTransforms
const RS = SpaceAGORA.RuntimeServices

# Raw inclusion remains supported for the legacy reference-system sandbox.
# Evaluate the unchanged arithmetic in a separate namespace to detect changes
# caused by ownership, imports, or module-scoped configuration.
module RawFrameSandbox
using SpaceAGORA.SimulationModel.EphemeridesModels: ephemerides_requires_spice, planet_frame_lpi
const args = :legacy_includer
include(joinpath(@__DIR__, "..", "..", "src", "core", "interfaces", "reference_system.jl"))
end

module QuaternionSandbox
include(joinpath(@__DIR__, "..", "..", "src", "core", "numerics", "quaternion_utils.jl"))
end

_bits(x::Float64) = reinterpret(UInt64, x)
_bits(x::AbstractArray) = map(_bits, x)
_bits(x::Tuple) = map(_bits, x)
_same_result(a, b) = typeof(a) === typeof(b) && _bits(a) == _bits(b)

@testset "Shared math ownership" begin
    quaternion_names = (
        :IDENTITY_QUATERNION, :quat_mult, :project_unit_quaternion, :hat,
        :rot, :error_quaternion, :qToEulerAngles, :dcm_to_quaternion,
    )
    quaternion_consumers = (SM, SM.Kinematics, SM.ControlHooks, SM.DynamicEffectors.PerturbationEffectors)
    for consumer in quaternion_consumers, name in quaternion_names
        @test getfield(consumer, name) === getfield(QM, name)
    end
    @test parentmodule(QM.rot) === QM

    frame_names = (
        :_EARTH_HIGH_PREC_BODY_FIXED_FRAME, :_EARTH_FALLBACK_BODY_FIXED_FRAME,
        :_spice_lock, :_spice_frame_lock, :r_intor_p!, :r_pintor_i,
        :_spice_body_fixed_frame, :_body_fixed_state_xform,
        :_j2000_to_body_fixed_state, :_body_fixed_to_j2000_state,
        :_planet_flattening, :orbitalelemtorv, :_wrap_2pi, :_safe_acos,
        :_rvtoorbitalelement_core, :rvtoorbitalelement, :rtoalfadeltar,
        :alfadeltartor, :latlongtor, :latlongtoOE, :rtolatlong,
        :rtolatlongrad, :latlongtoNED, :orbital_elements_to_lvlh_quaternion,
        :rotate_vector_by_quaternion, :rtn_dcm_from_inertial, :_rtn_rate_rad_s,
        :inertial_to_rtn_relative_state, :rtn_to_inertial_relative_state,
        :rtn_accel_to_inertial,
    )
    frame_consumers = (
        SpaceAGORA.SimulationEngine, SM.SimulationCallbacks, SM.GuidanceHooks,
        SM.ControlHooks, SM.DynamicEffectors.AerodynamicEffectors,
        SpaceAGORA.TelemetryVerification,
    )
    for consumer in frame_consumers, name in frame_names
        @test getfield(consumer, name) === getfield(FT, name)
    end
    @test SM.DynamicEffectors.PerturbationEffectors.rtolatlong === FT.rtolatlong
    @test parentmodule(FT.rtolatlong) === FT
    for consumer in frame_consumers
        @test consumer._spice_lock() === RS.SPICE_LOCK
        @test consumer._spice_frame_lock() === RS.tracked_lock(:spice_frame)
    end
end

@testset "Quaternion results and direction conventions" begin
    baseline = QuaternionSandbox.QuaternionMath
    qs = (
        SVector(0.0, 0.0, 0.0, 1.0),
        SVector(1.0, 0.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0, 0.0),
        SVector(0.0, 0.0, 1.0, 0.0),
        normalize(SVector(0.2, -0.4, 0.1, 0.8)),
    )
    for q in qs
        @test _same_result(QM.rot(q), baseline.rot(q))
        @test _same_result(QM.project_unit_quaternion(q), baseline.project_unit_quaternion(q))
        @test _same_result(QM.qToEulerAngles(q), baseline.qToEulerAngles(q))
        @test _same_result(QM.dcm_to_quaternion(QM.rot(q)), baseline.dcm_to_quaternion(QM.rot(q)))
        for p in qs
            @test _same_result(QM.quat_mult(q, p), baseline.quat_mult(q, p))
            @test _same_result(QM.error_quaternion(q, p), baseline.error_quaternion(q, p))
        end
    end
    for q in (zeros(4), [NaN, 0.0, 0.0, 1.0], [Inf, 0.0, 0.0, 1.0])
        @test QM.project_unit_quaternion(q) === QM.IDENTITY_QUATERNION
    end
    v = SVector(2.0, -3.0, 5.0)
    w = SVector(-1.0, 4.0, 2.0)
    @test _same_result(QM.hat(v), baseline.hat(v))
    @test QM.hat(v) * w ≈ cross(v, w)
    qz = SVector(0.0, 0.0, sin(pi / 4), cos(pi / 4))
    x = SVector(1.0, 0.0, 0.0)
    # rot is passive (inertial to body); the vector helper is active.
    @test QM.rot(qz) * x ≈ SVector(0.0, -1.0, 0.0) atol=1e-15
    @test FT.rotate_vector_by_quaternion(collect(x), collect(qz)) ≈ [0.0, 1.0, 0.0] atol=1e-15
    @test QM.rot(QM.dcm_to_quaternion(QM.rot(qz))) ≈ QM.rot(qz) atol=1e-15
end

@testset "Frame results across ownership contexts" begin
    # No native kernels or external assets are required for these probes.
    planet = (
        name="Earth", Rp_e=6378137.0, Rp_p=6356752.314245,
        μ=3.986004418e14, ω=SVector(0.0, 0.0, 7.292115e-5),
        L_PI=QM.rot(normalize(SVector(0.1, -0.2, 0.3, 0.9))),
    )
    ephemerides = SM.SimpleEphemeridesModel()
    oes = (
        SVector(7.0e6, 0.0, 0.0, 0.0, 0.0, 0.2, 100.0),
        SVector(8.0e6, 0.2, 0.7, 0.4, 1.2, 2.1, 200.0),
        SVector(9.0e6, 0.1, 2.4, 3.1, 0.8, 4.7, 300.0),
    )
    for oe in oes
        @test _same_result(FT.orbitalelemtorv(oe, planet), RawFrameSandbox.orbitalelemtorv(oe, planet))
        r_vec, v_vec = FT.orbitalelemtorv(oe, planet)
        r, v = SVector{3, Float64}(r_vec), SVector{3, Float64}(v_vec)
        for (name, arguments) in (
            (:rvtoorbitalelement, (r, v, planet)),
            (:rvtoorbitalelement, (r, v, oe[7], planet)),
            (:r_intor_p!, (r, v, planet)),
            (:r_intor_p!, (r, v, planet, 4321.0, ephemerides)),
            (:r_pintor_i, (r, v, planet)),
            (:rtoalfadeltar, (r,)),
            (:rtolatlong, (r, planet)),
            (:rtolatlongrad, (r, planet)),
            (:rtn_dcm_from_inertial, (r, v)),
            (:inertial_to_rtn_relative_state, (r .+ SVector(3.0, -2.0, 1.0), v .+ SVector(0.1, 0.2, -0.1), r, v)),
            (:rtn_to_inertial_relative_state, (SVector(3.0, -2.0, 1.0), SVector(0.1, 0.2, -0.1), r, v)),
            (:rtn_accel_to_inertial, (SVector(0.1, 0.2, -0.1), r, v)),
            (:orbital_elements_to_lvlh_quaternion, (oe[4], oe[3], oe[5], oe[6])),
        )
            @test _same_result(getfield(FT, name)(arguments...), getfield(RawFrameSandbox, name)(arguments...))
        end
    end
    for latitude in (-pi / 2, -pi / 4, 0.0, pi / 4, pi / 2)
        for (name, arguments) in (
            (:latlongtor, (SVector(latitude, 0.4, 125000.0), planet, 0.0, 0.0, 0.0)),
            (:latlongtoNED, (SVector(125000.0, latitude, 0.4),)),
        )
            @test _same_result(getfield(FT, name)(arguments...), getfield(RawFrameSandbox, name)(arguments...))
        end
        r = SVector{3, Float64}(FT.latlongtor(SVector(latitude, 0.4, 125000.0), planet, 0.0, 0.0, 0.0))
        @test _same_result(FT.rtolatlong(r, planet), RawFrameSandbox.rtolatlong(r, planet))
    end
    @test_throws ArgumentError FT.rtn_dcm_from_inertial(zeros(3), zeros(3))

    seen = Ref{Any}(:unset)
    topo = merge(planet, (
        topography_function=(received, clm, slm, lat, lon, a) -> (seen[] = received; a),
        Clm_topo=nothing, Slm_topo=nothing, A_topo=planet.Rp_e,
    ))
    r = SVector(planet.Rp_e + 125000.0, 0.0, 0.0)
    @test FT.rtolatlong(r, topo, true)[1] == 125000.0
    @test seen[] === nothing
    @test RawFrameSandbox.rtolatlong(r, topo, true)[1] == 125000.0
    @test seen[] === RawFrameSandbox.args
end
