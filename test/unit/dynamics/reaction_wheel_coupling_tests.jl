"""A control effector representing an idealized reaction wheel: the commanded body
torque is entirely wheel-attributable (reacts against stored wheel momentum), so
`calcControlForceTorque` and `calcReactionWheelTorque` return the same torque."""
struct TestReactionWheelEffector
    torque_body::SVector{3, Float64}
end

function SimulationModel.calcControlForceTorque(m::TestReactionWheelEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)
    return SVector{3, Float64}(0.0, 0.0, 0.0), m.torque_body
end

function SimulationModel.calcReactionWheelTorque(m::TestReactionWheelEffector, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)
    return m.torque_body
end

function _make_reaction_wheel_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    q0::SVector{4, Float64},
    ω0::SVector{3, Float64},
)
    root = Link{3}(
        root=true,
        m=500.0,
        ref_area=12.0,
        J_rw=SMatrix{3, 3, Float64}(1.0I),
    )
    a = (EARTH.Rp_e + ra_alt_m + EARTH.Rp_e + rp_alt_m) / 2.0
    e = (ra_alt_m - rp_alt_m) / (ra_alt_m + rp_alt_m + 2 * EARTH.Rp_e)
    ic = InitialCondition(a, e, 35.0, 40.0, 10.0, 175.0, q0, ω0)
    return SpacecraftModel(Joint[], [root], root, true, root.m, 0.0, root.inertia, 0, 0, ic, 1)
end

@testset "Reaction Wheel Momentum State Couples Into Rigid-Body RHS" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    ω0 = SVector{3, Float64}(0.05, -0.03, 0.02)
    h_wheels0 = SVector{3, Float64}(0.02, -0.01, 0.015)
    torque_body = SVector{3, Float64}(0.01, -0.005, 0.008)

    sc = _make_reaction_wheel_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, q0=q0, ω0=ω0)
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    sc.inertia_tensor = inertia_tensor

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(gravity_gradient=false),),
        control_effectors=(TestReactionWheelEffector(torque_body),),
        control_rates=[1.0],
        keplerian=true
    )

    u0 = build_initial_conditions(args)
    u0.sc[1].h_wheels .= h_wheels0
    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    rw_assembly = sc.root.rw_assembly
    J_rw = rw_assembly.J_rw
    J_rw_pinv = rw_assembly.J_rw_pinv

    # du.h_wheels must match ḣ_wheel = J_rw_pinv * (-rw_torque_body).
    dh_wheels_expected = J_rw_pinv * (-torque_body)
    @test isapprox(SVector{3, Float64}(du0.sc[1].h_wheels), dh_wheels_expected; atol=1e-12, rtol=1e-10)

    # du.ω must include the ω×h_wheel gyroscopic term in addition to the reaction torque.
    h_wheel_body = J_rw * h_wheels0
    ωdot_expected = SimulationModel.DynamicsRotational.angular_acceleration(
        ω0, inertia_tensor, torque_body; include_gyroscopic=true, h_wheel_body=h_wheel_body,
    )
    @test isapprox(SVector{3, Float64}(du0.sc[1].ω), ωdot_expected; atol=1e-12, rtol=1e-10)

    # Momentum-conservation regression: with zero external torque, the total system
    # angular momentum H_body = I*ω + J_rw*h_wheels obeys the torque-free transport
    # theorem Ḣ_body = -ω×H_body exactly — the idealized reaction-wheel torque command
    # must cancel out of the total-system budget regardless of its value.
    H_body = inertia_tensor * ω0 + J_rw * h_wheels0
    Ḣ_body = inertia_tensor * ωdot_expected + J_rw * dh_wheels_expected
    @test isapprox(Ḣ_body, -cross(ω0, H_body); atol=1e-10, rtol=1e-9)

    # Short-horizon integration check: norm(H_body) is rotation-invariant, so for
    # torque-free motion it must stay constant even though H_body itself precesses.
    function rhs!(du, u, p, t)
        du .= 0.0
        spacecraft_dynamics!(du, u, p, t)
        return nothing
    end
    function rk4_step(u, p, t, dt)
        k1 = zero(u); rhs!(k1, u, p, t)
        k2 = zero(u); rhs!(k2, u .+ (dt / 2) .* k1, p, t + dt / 2)
        k3 = zero(u); rhs!(k3, u .+ (dt / 2) .* k2, p, t + dt / 2)
        k4 = zero(u); rhs!(k4, u .+ dt .* k3, p, t + dt)
        return u .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end

    u = deepcopy(u0)
    t = 0.0
    dt = 0.05
    H_norm_0 = norm(inertia_tensor * SVector{3, Float64}(u.sc[1].ω) + J_rw * SVector{3, Float64}(u.sc[1].h_wheels))
    for _ in 1:40
        u = rk4_step(u, p, t, dt)
        t += dt
    end
    H_norm_T = norm(inertia_tensor * SVector{3, Float64}(u.sc[1].ω) + J_rw * SVector{3, Float64}(u.sc[1].h_wheels))
    @test isapprox(H_norm_T, H_norm_0; atol=1e-9, rtol=1e-6)
end

@testset "Zero-Wheel Spacecraft Unaffected By Reaction Wheel Wiring" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    ω0 = SVector{3, Float64}(0.02, -0.03, 0.015)
    sc = make_spacecraft(
        ra_alt_m=500e3, rp_alt_m=500e3, orientation_state=(q0, ω0),
    )
    inertia_tensor = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    sc.inertia_tensor = inertia_tensor
    applied_torque = SVector{3, Float64}(0.12, -0.08, 0.05)

    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=10.0,
        EI_km=120.0,
        dynamic_effectors=(ConstantTorqueModel(applied_torque),),
        keplerian=true
    )

    @test sc.root.rw_assembly.n_wheels == 0
    u0 = build_initial_conditions(args)
    @test !hasproperty(u0.sc[1], :h_wheels)

    du0 = copy(u0)
    du0 .= 0.0
    p = ODEParams(n_sats=1, args=args)
    spacecraft_dynamics!(du0, u0, p, 0.0)

    ω = SVector{3, Float64}(u0.sc[1].ω)
    ωdot_expected = SimulationModel.DynamicsRotational.angular_acceleration(
        ω, inertia_tensor, applied_torque; include_gyroscopic=true,
    )
    @test isapprox(SVector{3, Float64}(du0.sc[1].ω), ωdot_expected; atol=1e-12, rtol=1e-10)
end
