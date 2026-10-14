module RobotArmGravityTests

using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using ComponentArrays

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const CRA = SM.ClothRobotArmDynamics

const MASSES = (0.5, 0.4, 0.3)
const K_TRANS = 5.0e3

function static_plan()
    arm = SM.default_cloth_arm_model(link_lengths_m=(0.9, 0.8, 0.6), link_radii_m=(0.06, 0.05, 0.04),
        link_masses_kg=MASSES, mount_offset_body=(0.5, 0.0, 0.6))
    base = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    q0 = [0.1, 0.2, -0.1]
    target = SM.cloth_fk(arm, base, q0).end_effector_position
    plan = SM.plan_robot_arm_motion(arm, base, q0, target;
        config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=10.0))
    # The arm must be commanded to hold still, otherwise the test measures the plan.
    @test maximum(abs.(plan.q_ref .- plan.q_ref[:, 1])) < 1.0e-6
    return plan
end

# Parent-minus-child attachment point separation of joint i (the spring deflection).
function joint_deflections(sc_view, plan)
    links = plan.model.links
    return map(eachindex(links)) do i
        parent_idx, parent_point = i == 1 ? (0, plan.model.mount_offset_body) :
            (i - 1, links[i - 1].vector_parent - links[i - 1].com_offset_parent)
        parent = CRA._coupled_parent_kinematics(sc_view, parent_idx, parent_point)
        child = CRA._coupled_parent_kinematics(sc_view, i, -links[i].com_offset_parent)
        norm(parent.point - child.point)
    end
end

@testset "Arm links feel gravity at their own positions" begin
    plan = static_plan()
    planet = make_no_gram_planet(:earth)
    sc = make_three_body_spacecraft(
        bus_dims=(1.0, 1.0, 1.2), panel_dims=(0.01, 1.0, 0.6), bus_mass=120.0, panel_mass_each=4.0, panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=planet.Rp_e + 410e3, rp=planet.Rp_e + 400e3, i=51.6, ω=0.0, Ω=30.0, ν=0.0),
        prop_mass=0.0, id=1)
    base_cfg = make_example_config(
        planet=planet, spacecraft=sc, mission_time=120.0,
        initial_time=SM.InitialTime(year=2024, month=3, day=1, hour=12, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true, keplerian=true, EI_km=120.0,
        verbose=false, results=false)
    args = SM.SimConfig._with_configuration(base_cfg;
        mission_configuration=SM.MissionConfiguration(mission_type=base_cfg.mission_configuration.mission_type,
            keplerian=true, number_of_orbits=1, mission_time=120.0, orientation_sim=true, num_steps_to_save=1000),
        control_model=SM.ControlModel(
            control_effectors=(SM.RobotArmControlEffector(plan=plan, spacecraft_idx=1, k_translation_n_m=K_TRANS),),
            control_rates=[0.1]))

    sol = run_simulation(args; return_solution=true)
    max_defl = 0.0
    for u in sol.u
        max_defl = max(max_defl, maximum(joint_deflections(u.sc[1], plan)))
    end

    r_orbit = planet.Rp_e + 410e3
    g = planet.μ / r_orbit^2
    prefix_static = maximum(MASSES) * g / K_TRANS
    # Physical residual: the tidal difference between a link and the bus, 3 mu/r^3 times the
    # link offset (<= 3 m) times m/k, about 1e-9 m here. The measured value also carries the
    # integrator's floor: positions are ~7e6 m, so solver tolerance and rounding put ~1e-8 m of
    # noise into the difference of two absolute positions. The tolerance is 1e-4 of the pre-fix
    # static deflection m*g/k (~9e-4 m), which clears that floor and is still far below the defect.
    tidal = 3 * planet.μ / r_orbit^3 * 3.0 * maximum(MASSES) / K_TRANS
    tol = 1.0e-4 * prefix_static
    @info "arm gravity engine run" max_defl prefix_static tidal tol
    @test max_defl < tol
end

@testset "Pre-fix behavior (no link gravity) loads the springs with m*g" begin
    # One link, bus and link in free fall in a point-mass field. Without link gravity the
    # link does not fall, so its spring must hold the bus-to-link separation at m*g/k.
    arm = SM.default_cloth_arm_model(link_lengths_m=(0.5,), link_radii_m=(0.05,), link_masses_kg=(0.5,), joint_axes=((0.0, 0.0, 1.0),),
        mount_offset_body=(0.0, 0.0, 0.0))
    base = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    plan = SM.plan_robot_arm_motion(arm, base, [0.0], SM.cloth_fk(arm, base, [0.0]).end_effector_position;
        config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=1.0))
    m = plan.model.links[1].mass_kg
    mu = 3.986004418e14
    r0 = 6.778e6
    gfun(r) = -mu .* r ./ norm(r)^3
    k = 100.0
    function settle(link_gravity_ii)
        shape = merge((pos=zeros(3), vel=zeros(3), mass=120.0, heat_loads=zeros(1),
            q=Float64[0, 0, 0, 1], ω=zeros(3)), CRA.coupled_cloth_robot_arm_state_shape(plan))
        u = ComponentVector(shape)
        u.pos .= [r0, 0.0, 0.0]
        u.vel .= [0.0, 7.67e3, 0.0]
        CRA.initialize_coupled_cloth_robot_arm_state!(u, plan)
        function f(x)
            du = zero(x)
            forces = MVector{3, Float64}(0.0, 0.0, 0.0)
            torques = MVector{3, Float64}(0.0, 0.0, 0.0)
            CRA.assign_coupled_cloth_robot_arm_rhs!(du, x, plan, 0.0, forces, torques;
                k_translation_n_m=k, c_translation_n_s_m=30.0, k_rotation_n_m_rad=0.0,
                c_rotation_n_m_s_rad=0.0, link_gravity_ii=link_gravity_ii)
            du.pos .= x.vel
            du.vel .= forces ./ x.mass .+ gfun(SVector{3, Float64}(x.pos))
            return du
        end
        dt = 2.0e-3
        for _ in 1:Int(round(20.0 / dt))   # RK4, damped joint: settles well inside 20 s
            k1 = f(u); k2 = f(u .+ (dt / 2) .* k1); k3 = f(u .+ (dt / 2) .* k2); k4 = f(u .+ dt .* k3)
            u .+= (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        end
        return joint_deflections(u, plan)[1]
    end
    expected = m * mu / r0^2 / k
    defl_old = settle(nothing)
    defl_new = settle(gfun)
    @info "pre-fix vs fixed static deflection" expected defl_old defl_new
    @test defl_old ≈ expected rtol=0.10
    # Tidal residual only: far below the static value (link sits at the bus's own offset).
    @test defl_new < 1.0e-3 * expected
end

end # module
