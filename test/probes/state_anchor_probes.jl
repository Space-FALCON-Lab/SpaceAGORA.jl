# Scheduled state anchors: the callback lands on each anchor time exactly,
# overwrites the satellite's inertial position and velocity, reports the drift
# it removed, and the propagation continues from the anchored state. Checked
# on a single-satellite Earth orbit by anchoring onto a different orbit
# mid-run and evaluating the solution on both sides of the anchor.
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const REPO_ROOT_ANCHOR = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH_ANCHOR = joinpath(REPO_ROOT_ANCHOR, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const SE_ANCHOR = SpaceAGORA.SimulationEngine

@testset "Scheduled state anchors re-initialise the propagated state" begin
    planet = Earth("", SPICE_PATH_ANCHOR)
    root = Link(root=true, m=140.0, ref_area=1.2)
    ic = InitialCondition(ra=planet.Rp_e + 900e3, rp=planet.Rp_e + 800e3, i=28.0, ω=15.0, Ω=20.0, ν=0.0)
    sc = SpacecraftModel(
        joints=Joint[], links=Link[root], root=root, instant_actuation=true, prop_mass=15.0,
        inertia_tensor=root.inertia, n_reaction_wheels=0, n_thrusters=0, initial_condition=ic, id=1,
    )
    cfg = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1, mission_time=3000.0, orientation_sim=false, num_steps_to_save=50),
        environment_model=EnvironmentModel(planet=planet, EI=120.0, density_model=ExponentialAtmosphereModel(planet), thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0),
    )

    # Anchor onto a circular orbit 200 km above the initial ellipse at t=1500 s.
    r_anchor = planet.Rp_e + 1100e3
    v_anchor = sqrt(planet.μ / r_anchor)
    pos = SVector{3, Float64}(r_anchor, 0.0, 0.0)
    vel = SVector{3, Float64}(0.0, v_anchor * cosd(28.0), v_anchor * sind(28.0))
    anchor = StateAnchor(1500.0, 1, vcat(pos, vel))
    @test anchor.elapsed_s == 1500.0
    @test_throws ArgumentError StateAnchor(1500.0, 1, [1.0, 2.0, 3.0])
    @test_throws ArgumentError StateAnchor(NaN, 1, vcat(pos, vel))
    @test_throws ArgumentError StateAnchor(1500.0, 0, vcat(pos, vel))
    @test_throws ArgumentError StateAnchor(1500.0, 1, vcat(pos, vel); orbit_count=0)
    @test StateAnchor(1500.0, 1, vcat(pos, vel)).orbit_count === nothing
    @test StateAnchor(1500.0, 1, vcat(pos, vel); orbit_count=8).orbit_count == 8
    @test_throws ArgumentError get_state_anchor_callback(StateAnchor[])

    # Capture the anchor report lines (redirect_stdout needs a file, not an IOBuffer).
    function run_captured(anchors)
        withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
            mktempdir() do tmp
                cd(tmp) do
                    log_path = joinpath(tmp, "anchor_log.txt")
                    sol = open(log_path, "w") do f
                        redirect_stdout(f) do
                            run_simulation(cfg; return_solution=true, extra_callbacks=(get_state_anchor_callback(anchors),))
                        end
                    end
                    return sol, read(log_path, String)
                end
            end
        end
    end
    sol, text = run_captured([anchor])
    @test occursin("state_anchor sat=1 index=1 t_s=1500.0", text)
    @test sol !== nothing
    @test 1500.0 in sol.t                       # the anchor is a tstop, hit exactly

    radius_at(t) = norm(SE_ANCHOR._state_position_ii(sol(t), 1))
    speed_at(t) = norm(SE_ANCHOR._state_velocity_ii(sol(t), 1))
    # Before the anchor: the initial ellipse (800 to 900 km altitude).
    @test all(t -> planet.Rp_e + 795e3 <= radius_at(t) <= planet.Rp_e + 905e3, 0.0:100.0:1400.0)
    # After it: the anchored circle, to within the integrator over half an
    # orbit. The solution at the anchor time itself is the pre-anchor state
    # (the save precedes the affect) and the interpolant of the step that
    # straddles the anchor mixes both, so the checks start one save interval
    # (10 s) after it.
    @test all(t -> abs(radius_at(t) - r_anchor) < 2e3, 1510.0:100.0:3000.0)
    @test all(t -> abs(speed_at(t) - v_anchor) < 2.0, 1510.0:100.0:3000.0)
    # The reported drift is the distance between the propagated position and
    # the anchor: two points on different orbits, of the order of the radius.
    m = match(r"drift_pos_m=([0-9.eE+-]+)", text)
    @test m !== nothing
    @test 1e6 < parse(Float64, m.captures[1]) < 3e7

    # The gravity-backbone split keeps a second-order state whose symplectic
    # stepper ends the solve at the next tstop after an in-place edit, so the
    # callback refuses that solver mode at initialisation instead of
    # truncating the run silently.
    @test_throws ArgumentError withenv("SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split") do
        run_captured([anchor])
    end

    # Several anchors at one instant are all applied in that instant's single
    # callback invocation, and the ones after it still fire.
    pos_alt = SVector{3, Float64}(0.0, r_anchor, 0.0)
    vel_alt = SVector{3, Float64}(-v_anchor * cosd(28.0), 0.0, v_anchor * sind(28.0))
    r_late = planet.Rp_e + 1300e3
    v_late = sqrt(planet.μ / r_late)
    late = StateAnchor(2500.0, 1, vcat(SVector{3, Float64}(r_late, 0.0, 0.0), SVector{3, Float64}(0.0, v_late * cosd(28.0), v_late * sind(28.0))))
    sol_multi, text_multi = run_captured([anchor, StateAnchor(1500.0, 1, vcat(pos_alt, vel_alt)), late])
    @test occursin("index=1 t_s=1500.0", text_multi)
    @test occursin("index=2 t_s=1500.0", text_multi)
    @test occursin("index=3 t_s=2500.0", text_multi)
    radius_multi(t) = norm(SE_ANCHOR._state_position_ii(sol_multi(t), 1))
    # The last of the simultaneous anchors wins (same radius here, different
    # phase), and the late anchor moves the orbit up by 200 km.
    @test all(t -> abs(radius_multi(t) - r_anchor) < 2e3, 1510.0:100.0:2400.0)
    @test all(t -> abs(radius_multi(t) - r_late) < 2e3, 2510.0:100.0:3000.0)
    pos_at_1600 = SE_ANCHOR._state_position_ii(sol_multi(1600.0), 1)
    @test pos_at_1600[2] > 0.5 * r_anchor    # on the second anchor's phase, not the first's

    # An anchor can reset the orbit counter to the anchored trajectory's count:
    # a three-orbit mission anchored at 1500 s with orbit_count=3 completes at
    # the first apoapsis after the anchor instead of two orbits later.
    cfg_orbits = SimulationConfiguration(
        simulation_settings=cfg.simulation_settings,
        mission_configuration=MissionConfiguration(mission_type=MissionOrbits, keplerian=true, number_of_orbits=3, mission_time=30000.0, orientation_sim=false, num_steps_to_save=50),
        environment_model=cfg.environment_model, dynamics_model=cfg.dynamics_model,
        guidance_model=cfg.guidance_model, navigation_model=cfg.navigation_model, control_model=cfg.control_model,
        initial_time=cfg.initial_time, integration_tolerances=cfg.integration_tolerances,
    )
    sol_count, text_count = withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        mktempdir() do tmp
            cd(tmp) do
                log_path = joinpath(tmp, "anchor_log.txt")
                sol = open(log_path, "w") do f
                    redirect_stdout(f) do
                        run_simulation(cfg_orbits; return_solution=true, extra_callbacks=(get_state_anchor_callback([StateAnchor(1500.0, 1, vcat(pos, vel); orbit_count=3)]),))
                    end
                end
                return sol, read(log_path, String)
            end
        end
    end
    @test occursin("orbit_count=3 was=1", text_count)
    @test occursin("termination_cause=orbit_count sat=1 orbits=3", text_count)
    period_anchor_s = 2π * sqrt(r_anchor^3 / planet.μ)
    @test 1500.0 < sol_count.t[end] < 1500.0 + period_anchor_s

    # Anchors at or before the initial time are skipped, later ones still apply.
    sol2, text2 = run_captured([StateAnchor(0.0, 1, vcat(pos, vel)), anchor])
    @test occursin("index=2 t_s=1500.0", text2)
    @test !occursin("index=1 t_s=0.0", text2)
    @test abs(norm(SE_ANCHOR._state_position_ii(sol2(100.0), 1)) - radius_at(100.0)) < 1.0
end
println("state_anchor_probes_ok")
