# Link-frame kinematics helpers: both branches of the frame rotations and both
# rotate_link methods, which no simulation path exercises (links are rotated
# through the multibody joints), so the file sat below the coverage floor.
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const KM = SpaceAGORA.SimulationModel.Kinematics

@testset "Link-frame kinematics" begin
    root = Link(root=true, m=140.0, ref_area=1.2)
    child = Link(root=false, m=10.0, ref_area=0.5)
    model = SpacecraftModel(
        joints=Joint[], links=Link[root, child], root=root, instant_actuation=true,
        prop_mass=15.0, inertia_tensor=root.inertia, n_reaction_wheels=0, n_thrusters=0,
        initial_condition=InitialCondition(ra=6.9e6, rp=6.88e6, i=28.0, ω=15.0, Ω=20.0, ν=0.0), id=1,
    )

    # Rotate the child about z by 90 degrees through the DCM method.
    theta = pi / 2
    dcm = SMatrix{3, 3, Float64}(cos(theta), sin(theta), 0.0, -sin(theta), cos(theta), 0.0, 0.0, 0.0, 1.0)
    rotate_link(child, dcm)
    @test norm(child.q) ≈ 1.0
    @test isapprox(KM.rot(SVector{4, Float64}(child.q)), dcm; atol=1e-12) ||
          isapprox(KM.rot(SVector{4, Float64}(child.q)), dcm'; atol=1e-12)

    # The quaternion method normalises what it is given.
    q_raw = SVector{4, Float64}(2.0, 0.0, 0.0, 0.0)
    rotate_link(child, q_raw)
    @test child.q ≈ [1.0, 0.0, 0.0, 0.0]

    # Root links are rotated by the attitude state, never directly.
    @test_throws AssertionError rotate_link(root, q_raw)
    @test_throws AssertionError rotate_link(root, dcm)

    # Frame rotations: the root's link frame is the body frame; a child's
    # inertial rotation composes the root attitude with its own.
    @test rotate_to_body(root) == I(3)
    rotate_link(child, dcm)
    @test isapprox(rotate_to_body(child), KM.rot(SVector{4, Float64}(child.q))'; atol=1e-12)
    R_root = rotate_to_inertial(model, root, 1)
    @test isapprox(R_root, KM.rot(SVector{4, Float64}(root.q))'; atol=1e-12)
    R_child = rotate_to_inertial(model, child, 1)
    @test isapprox(R_child, KM.rot(SVector{4, Float64}(root.q))' * KM.rot(SVector{4, Float64}(child.q))'; atol=1e-12)
    @test isapprox(R_child' * R_child, I(3); atol=1e-12)
end
println("kinematics_probes_ok")
