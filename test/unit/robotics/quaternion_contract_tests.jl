module QuaternionContractTests

using Test
using LinearAlgebra
using StaticArrays
using ComponentArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const QM = SM.QuaternionMath
const CM = SM.ClothMultibody
const CA = SM.ClothRobotArmDynamics
const RB = SM.Robotics
const RP = SM.RobotArmPlanning
const DR = SM.DynamicsRotational
const QI = SVector(0.0, 0.0, 0.0, 1.0)
const Z3 = SVector(0.0, 0.0, 0.0)
const Z4 = SVector(0.0, 0.0, 0.0, 0.0)

# A stationary one-link reference: actual plan and state types, no planner or solver.
function one_link_plan()
    model = RB.default_cloth_arm_model(
        link_lengths_m=(0.2,), link_radii_m=(0.01,),
        link_masses_kg=(1.0,), joint_axes=((0.0, 0.0, 1.0),),
    )
    base = RB.ClothArmBasePose([0.0, 0.0, 0.0])
    tip = SVector(0.2, 0.0, 0.0)
    return RP.RobotArmPlan(
        model, base, [0.0, 1.0], zeros(1, 2), zeros(1, 2), zeros(1, 2),
        repeat(reshape(collect(tip), 3, 1), 1, 2), [0.0], [0.0], tip, 0.0,
        :cloth_quintic,
    )
end

function cloth_qdot(model, q, omega)
    x = CM.compliant_state_vector([Z3], [QI])
    # Write the raw state after construction to exercise normalization in the RHS reader.
    x[4:7] .= q
    x[11:13] .= omega
    before = copy(x)
    dx = CM.compliant_multibody_dynamics(model, x)
    @test x == before
    return SVector{4, Float64}(dx[4:7])
end

function arm_qdot(plan, q, omega)
    sc = ComponentVector(merge(
        (pos=zeros(3), vel=zeros(3), q=collect(QI), ω=zeros(3)),
        CA.coupled_cloth_robot_arm_state_shape(plan),
    ))
    CA.initialize_coupled_cloth_robot_arm_state!(sc, plan)
    sc.arm_q[:, 1] .= q
    sc.arm_ω[:, 1] .= omega
    before = copy(sc)
    du = zero(sc)
    CA.assign_coupled_cloth_robot_arm_rhs!(
        du, sc, plan, 0.0, zeros(3), zeros(3);
        k_translation_n_m=0.0, c_translation_n_s_m=0.0,
        k_rotation_n_m_rad=0.0, c_rotation_n_m_s_rad=0.0,
    )
    @test sc == before
    return SVector{4, Float64}(du.arm_q[:, 1])
end

@testset "Quaternion contracts across current owners" begin
    s = sqrt(0.5)
    qx = SVector(s, 0.0, 0.0, s)
    qy = SVector(0.0, s, 0.0, s)
    qz = SVector(0.0, 0.0, s, s)

    @testset "Raw and normalized products retain different scale contracts" begin
        xy = SVector(0.5, 0.5, 0.5, 0.5)
        yx = SVector(0.5, 0.5, -0.5, 0.5)
        for product in (QM.quat_mult, RB._quat_mul, CM._quat_raw_mul, CA._quat_raw_mul)
            @test product(qx, qy) ≈ xy atol=1e-15
            @test product(qy, qx) ≈ yx atol=1e-15
            @test product(2qx, 3qy) ≈ 6xy atol=2e-15
            @test product(-2qx, 3qy) ≈ -6xy atol=2e-15
            @test product(Z4, qy) == Z4
            @test product(qx, Z4) == Z4
            @test product(2QI, 3QI) == 6QI
            @test product(qx, qy) isa SVector{4, Float64}
        end
        for product in (CM._quat_mul, CA._quat_mul)
            @test product(qx, qy) ≈ xy atol=1e-15
            @test product(qy, qx) ≈ yx atol=1e-15
            @test product(2qx, 3qy) ≈ xy atol=1e-15
            @test product(-2qx, 3qy) ≈ -xy atol=1e-15
            @test product(Z4, qy) ≈ qy atol=1e-15
            @test product(qx, Z4) ≈ qx atol=1e-15
            @test product(2QI, 3QI) == QI
            @test product(qx, qy) isa SVector{4, Float64}
        end
        for owner in (CM, CA)
            @test owner._quat_conj(3qz) ≈ SVector(0.0, 0.0, -s, s) atol=1e-15
            @test owner._quat_mul(3qz, owner._quat_conj(3qz)) ≈ QI atol=1e-15
        end
    end

    @testset "Normalization preserves cutoffs and invalid-input fallback" begin
        ex = SVector(1.0, 0.0, 0.0, 0.0)
        local_cutoff = eps(Float64)
        canonical_cutoff = sqrt(eps(Float64))
        for project in (CM._unit_quat, CA._unit_quat, RB._unit_quat)
            for magnitude in (0.0, prevfloat(local_cutoff), local_cutoff)
                @test project(magnitude * ex) == QI
            end
            for magnitude in (nextfloat(local_cutoff), 1.0e-10, canonical_cutoff)
                @test project(magnitude * ex) == ex
                @test project(-magnitude * ex) == -ex
            end
        end
        for magnitude in (0.0, local_cutoff, 1.0e-10, prevfloat(canonical_cutoff), canonical_cutoff)
            @test QM.project_unit_quaternion(magnitude * ex) == QI
        end
        @test QM.project_unit_quaternion(nextfloat(canonical_cutoff) * ex) == ex
        @test QM.project_unit_quaternion(-nextfloat(canonical_cutoff) * ex) == -ex
        for project in (QM.project_unit_quaternion, CM._unit_quat, CA._unit_quat, RB._unit_quat)
            @test project(3qz) ≈ qz atol=1e-15
            @test project(-3qz) ≈ -qz atol=1e-15
            @test project(collect(3qz)) ≈ qz atol=1e-15
            for bad in (NaN, Inf, -Inf), component in 1:4
                input = collect(QI)
                input[component] = bad
                @test project(input) == QI
            end
        end
    end

    @testset "Rotation direction and input scale are observable" begin
        ex, ey, ez = SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)
        # Independent quarter-turn vector expectations, not a copied rotation formula.
        for (q, v, active) in ((qz, ex, ey), (qx, ey, ez), (qy, ez, ex))
            @test QM.rot(q) * v ≈ -active atol=1e-15
            @test QM.rot(2q) ≈ 4QM.rot(q) atol=1e-15
            for owner in (CM, CA, RB)
                @test owner._rot(q) * v ≈ active atol=1e-15
                @test owner._rot(2q) ≈ owner._rot(q) atol=1e-15
                @test owner._rot(-q) ≈ owner._rot(q) atol=1e-15
                @test owner._rot(q) ≈ transpose(QM.rot(q)) atol=1e-15
                @test owner._rot(q) isa SMatrix{3, 3, Float64}
            end
        end
        @test QM.rot(Z4) == zeros(3, 3)
        for owner in (CM, CA, RB)
            @test owner._rot(Z4) == Matrix{Float64}(I, 3, 3)
        end
        plan = one_link_plan()
        # Forward kinematics consumes the active local rotation and normalizes its product.
        pose = RB.cloth_fk(plan.model, plan.base_pose, [pi / 2])
        @test pose.end_effector_position ≈ SVector(0.0, 0.2, 0.0) atol=1e-15
        @test pose.link_quaternions[1] ≈ qz atol=1e-15
        # The exact SVector constructor stores a non-unit base; FK still projects link attitudes.
        scaled_base = RB.ClothArmBasePose(Z3, 2qz)
        scaled_pose = RB.cloth_fk(plan.model, scaled_base, [0.0])
        @test scaled_pose.end_effector_position ≈ pose.end_effector_position atol=1e-15
        @test norm(scaled_pose.link_quaternions[1]) ≈ 1.0 atol=1e-15
    end

    @testset "Direct RHS callers retain raw body-rate magnitude" begin
        model = CM.CompliantMultibodyModel(
            [CM.CompliantBody(:free_body, 1.0, SMatrix{3, 3, Float64}(I))],
            CM.CompliantJoint[], Z3, QI,
        )
        plan = one_link_plan()
        omega = SVector(2.0, -3.0, 4.0)
        # qz ⊗ [omega, 0] / 2, with a nonparallel body rate to distinguish product order.
        expected = SVector(2.5s, -0.5s, 2.0s, -2.0s)
        for derivative in (q -> cloth_qdot(model, q, omega), q -> arm_qdot(plan, q, omega))
            @test derivative(qz) ≈ expected atol=2e-15
            @test derivative(3qz) ≈ expected atol=2e-15
            @test derivative(-qz) ≈ -expected atol=2e-15
            @test derivative(Z4) ≈ SVector(1.0, -1.5, 2.0, 0.0) atol=1e-15
            # Local readers retain the norm cutoff, not the canonical squared-norm cutoff.
            @test derivative(SVector(1e-10, 0.0, 0.0, 0.0)) ≈ SVector(0.0, -2.0, -1.5, -1.0) atol=1e-15
        end
        for derivative in ((q, w) -> cloth_qdot(model, q, w), (q, w) -> arm_qdot(plan, q, w))
            @test derivative(qz, Z3) == Z4
            @test derivative(qz, 2omega) ≈ 2expected atol=4e-15
            @test derivative(qz, -omega) ≈ -expected atol=2e-15
            @test dot(qz, derivative(qz, omega)) ≈ 0.0 atol=1e-15
            @test norm(derivative(qz, omega)) ≈ norm(omega) / 2 atol=1e-15
        end
        # Canonical rotational kinematics is raw in q as well as omega, unlike the two readers.
        @test DR.quaternion_derivative(omega, qz) ≈ expected atol=2e-15
        @test DR.quaternion_derivative(omega, 3qz) ≈ 3expected atol=4e-15
        @test DR.quaternion_derivative(Z3, qz) == Z4
    end
end

end # module QuaternionContractTests
