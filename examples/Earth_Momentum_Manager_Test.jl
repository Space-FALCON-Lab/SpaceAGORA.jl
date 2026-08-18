## Closed-loop magnetic momentum management demo + regression.
#
# An LVLH-cascade attitude controller holds a 10-degree off-nadir roll, so the
# gravity-gradient torque at the setpoint is persistent and the (implicit)
# reaction wheels absorb it continuously. The MagneticMomentumManagerModel
# runs as a 1 Hz discrete control effector: it advances its wheel-momentum
# accumulator from the re-evaluated cascade command and commands the
# cross-product rod law tau = m x B against a tilted-dipole field.
#
# Regression asserts (CI examples smoke). The window spans >1 orbit so the
# field direction rotates and the field-parallel part of the settling
# transient becomes reachable by the rods:
#   1. boundedness: with the manager on, accumulated wheel momentum stays far
#      below the mu=0 (accounting-only) drift;
#   2. solver invariance: the discrete accumulator is step-size independent
#      (tight vs loose tolerances agree) — the property a continuous-effector
#      accumulator would violate under adaptive stepping.
include(joinpath(@__DIR__, "common.jl"))
using SPICE
using StaticArrays
using LinearAlgebra

const PE = SimulationModel.DynamicEffectors.PerturbationEffectors

planet = Earth("", SPICE_PATH)

# 700 km circular-ish, generic attitude ICs
q0 = normalize(SVector{4, Float64}(0.05, -0.03, 0.02, 0.998))
w0 = SVector{3, Float64}(0.0, -1.06e-3, 0.0)
a = planet.Rp_e + 700e3
ic = InitialCondition(a, 1.0e-4, 51.6, 0.0, 30.0, 170.0, q0, w0)

spacecraft = make_three_body_spacecraft(
    bus_dims=(1.8, 2.0, 1.5),
    panel_dims=(0.01, 1.2, 0.6),
    bus_mass=320.0,
    panel_mass_each=4.0,
    panel_offset_y=1.0,
    ic=ic,
    prop_mass=0.0,
    id=106
)

# attitude controller: hold 10 deg roll off nadir (generic round-number gains)
roll = deg2rad(10.0)
q_cmd = SVector(sin(roll / 2), 0.0, 0.0, cos(roll / 2))
cascade = PE.LVLHCascadeAttitudeControlModel(
    q_cmd_lb=q_cmd, k_out=[0.01, 0.01, 0.01], w_max=2.0e-3,
    k_rate=[0.6, 0.6, 0.6], tau_max=5.0e-3)

# tilted static dipole field (inertial), Earth-like magnitude
const M_DIP = 8.0e15 .* normalize(SVector(sin(deg2rad(11.0)), 0.0, cos(deg2rad(11.0))))
b_dipole(t, r_ii) = begin
    rn = norm(r_ii)
    rhat = r_ii / rn
    (3.0 * dot(M_DIP, rhat) .* rhat .- M_DIP) ./ rn^3
end
tau_cascade(t, r, v, q, w) = PE._lvlh_cascade_torque(cascade, r, v, q, w)

function run_managed(mu; reltol=1e-9)
    manager = MagneticMomentumManagerModel(
        mu_gain=mu, m_max_am2=10.0,
        commanded_torque=tau_cascade, b_field_ii=b_dipole)
    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=8_000.0,
        initial_time=InitialTime(year=2023, month=4, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredGravityModel(gravity_gradient=true), cascade),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        keplerian=true,
        EI_km=140.0
    )
    args = SimulationConfiguration(
        file_paths=base.file_paths,
        simulation_settings=base.simulation_settings,
        mission_configuration=base.mission_configuration,
        environment_model=base.environment_model,
        dynamics_model=base.dynamics_model,
        guidance_model=base.guidance_model,
        navigation_model=base.navigation_model,
        control_model=ControlModel(control_effectors=(manager,), control_rates=[1.0]),
        initial_time=base.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol, abstol_orbit=reltol,
            reltol_quaternion=reltol, abstol_quaternion=reltol,
            dt_max_orbit=3.0
        )
    )
    # isolate_state=false: run_simulation deep-copies args by default, so the
    # in-sim manager would be a copy and this handle would stay untouched.
    run_simulation(args; isolate_state=false)
    return manager
end

println("=== momentum manager on (mu=0.01) ===")
mgr_on = run_managed(3.0e-3)
println("  |H_w| end: ", norm(mgr_on.h_wheels), " N m s; |m| held: ", norm(mgr_on.held_dipole_am2), " A m^2; ticks: ", mgr_on.ticks)

println("=== accounting only (mu=0) ===")
mgr_off = run_managed(0.0)
println("  |H_w| end: ", norm(mgr_off.h_wheels), " N m s (unmanaged drift)")

@assert norm(mgr_off.h_wheels) > 0.0 "expected nonzero wheel-momentum drift at mu=0"
@assert norm(mgr_on.h_wheels) < 0.5 * norm(mgr_off.h_wheels) "manager failed to bound wheel momentum"

println("=== solver invariance (tight tolerances) ===")
mgr_tight = run_managed(3.0e-3; reltol=1e-11)
dev = norm(mgr_tight.h_wheels - mgr_on.h_wheels)
println("  |dH| loose vs tight: ", dev, " N m s")
@assert dev < 1.0e-3 * max(norm(mgr_off.h_wheels), 1e-12) "discrete accumulator is not step-size invariant"

println("Earth_Momentum_Manager_Test ok")
