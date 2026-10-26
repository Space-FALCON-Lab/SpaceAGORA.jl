module PowerBudgetTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA, TOML
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel

# A controller that holds levels for the mounted actuators of spacecraft 1.
mutable struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = i == 1 ? c.levels : nothing
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

const PLANET = SM.Earth()
const EPHEM = SM.SimpleEphemeridesModel()
const EPOCH = SM.InitialTime(year=2024, month=3, day=20, hour=0, minute=0, second=0.0)
const ORBIT_R = PLANET.Rp_e + 500e3

# A 10 kg bus in a 500 km equatorial orbit, attitude held at the identity, with
# an array facing the Sun at the start (24.5 W lit), a 5 Wh battery, a 5 W load
# and a laser terminal drawing 10 W at its 2 W command.
function power_spacecraft(; soc=0.5, battery=true, laser=true, wheel=false, facing=1.0, area=0.06, recovery_soc=0.05)
    sun0 = SM.ephemerides_sun_ii(EPHEM, PLANET, SM.ephemerides_epoch_seconds(EPHEM, EPOCH))
    ic = SM.CartesianInitialCondition(SVector(ORBIT_R, 0.0, 0.0), SVector(0.0, sqrt(PLANET.μ / ORBIT_R), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    bus = SM.Link(root=true, m=10.0, dims=MVector(0.1, 0.2, 0.3))
    sc = SM.SpacecraftModel(links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=ic, id=1)
    add_power!(sc, bus, SolarArray(normal=facing * sun0, area_m2=area, efficiency=0.3))
    wheel && add_actuator!(sc, bus, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
    battery && add_power!(sc, bus, Battery(capacity_wh=5.0, initial_soc=soc, recovery_soc=recovery_soc))
    add_power!(sc, bus, PowerLoad(power_w=5.0, name="avionics"))
    laser && add_actuator!(sc, bus, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=2.0, efficiency=0.2))
    return sc
end

function power_run(sc; mission_time=2 * 5677.0, solver=nothing, levels=[2.0], effectors=(HoldLevels(collect(Float64, levels)),),
                   settings=(;), kw...)
    base = make_example_config(planet=PLANET, spacecraft=sc, mission_time=mission_time, initial_time=EPOCH,
        dynamic_effectors=(SM.InverseSquaredGravityModel(),), ephemerides_model=EPHEM,
        orientation_sim=true, keplerian=false, EI_km=100.0, verbose=false, results=false)
    extra = solver === nothing ? (;) : (; solver_config=SM.SolverConfig(solver_mode=solver))
    if !isempty(settings)
        st = base.simulation_settings
        fields = NamedTuple{fieldnames(typeof(st))}(map(n -> getfield(st, n), fieldnames(typeof(st))))
        extra = merge(extra, (; simulation_settings=SM.SimulationSettings(; merge(fields, settings)...)))
    end
    args = SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=effectors, control_rates=fill(10.0, length(effectors))), extra...)
    return run_simulation(args; return_results=true, isolate_state=false, kw...)
end

@testset "Power budget" begin
    @testset "records and mounting" begin
        @test_throws ArgumentError SolarArray(normal=SVector(0.0, 0.0, 0.0), area_m2=1.0, efficiency=0.3)
        @test_throws ArgumentError SolarArray(normal=SVector(1.0, 0.0, 0.0), area_m2=1.0, efficiency=1.5)
        @test_throws ArgumentError Battery(capacity_wh=0.0)
        @test_throws ArgumentError Battery(capacity_wh=1.0, initial_soc=1.2)
        @test_throws ArgumentError Battery(capacity_wh=1.0, recovery_soc=0.0)
        @test_throws ArgumentError PowerLoad(power_w=-1.0)
        @test SolarArray(normal=SVector(0.0, 2.0, 0.0), area_m2=1.0, efficiency=0.3).normal == SVector(0.0, 1.0, 0.0)
        root = SM.Link(root=true, m=10.0)
        panel = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0))
        sc = SM.SpacecraftModel(links=[root, panel], root=root, battery_energy_index=80)
        @test !SM.has_power_parts(sc) && !SM.has_battery(sc)
        @test_throws ArgumentError add_power!(sc, panel, Battery(capacity_wh=10.0))      # the battery mounts on the root
        b = add_power!(sc, root, Battery(capacity_wh=10.0))
        @test_throws ArgumentError add_power!(sc, root, Battery(capacity_wh=10.0))       # at most one
        @test_throws ArgumentError add_power!(sc, SM.Link(m=1.0), PowerLoad(power_w=1.0)) # not one of its links
        @test SM.has_battery(sc) && SM.battery(sc) === b
        @test SM.initial_battery_wh(sc) == 8.0                                           # battery_energy_index, 80 %
        add_power!(sc, panel, SolarArray(normal=SVector(1.0, 0.0, 0.0), area_m2=0.5, efficiency=0.2, max_power_w=100.0))
        add_power!(sc, root, PowerLoad(power_w=3.0)); add_power!(sc, panel, PowerLoad(power_w=1.5))
        @test SM.bus_load_w(sc) == 4.5
    end

    @testset "generation follows the array normal, the flux and the lit fraction" begin
        root = SM.Link(root=true, m=10.0)
        panel = SM.Link(m=1.0, r=MVector(0.0, 1.0, 0.0), q=MVector(0.0, 0.0, sin(π / 4), cos(π / 4)))   # turned 90 deg about z
        sc = SM.SpacecraftModel(links=[root, panel], root=root)
        add_power!(sc, panel, SolarArray(normal=SVector(1.0, 0.0, 0.0), area_m2=2.0, efficiency=0.25, max_power_w=500.0))
        flux = SM.solar_flux_w_m2(SM.AU_M)
        @test flux == 1361.0
        # The panel's +x faces body +y after its 90 degree turn.
        @test SM.solar_generation_w(sc, SVector(0.0, 1.0, 0.0), flux, 1.0) ≈ 500.0                       # capped
        @test SM.solar_generation_w(sc, SVector(0.0, 1.0, 0.0), 0.5 * flux, 1.0) ≈ 2.0 * 0.25 * 0.5 * flux
        @test SM.solar_generation_w(sc, normalize(SVector(0.0, 1.0, 1.0)), 0.5 * flux, 1.0) ≈ 2.0 * 0.25 * 0.5 * flux / sqrt(2)
        @test SM.solar_generation_w(sc, SVector(0.0, -1.0, 0.0), flux, 1.0) == 0.0                       # the back of the cells
        @test SM.solar_generation_w(sc, SVector(0.0, 1.0, 0.0), 0.5 * flux, 0.25) ≈ 2.0 * 0.25 * 0.5 * flux * 0.25   # penumbra
        @test SM.solar_generation_w(sc, SVector(0.0, 1.0, 0.0), flux, 0.0) == 0.0                        # umbra
        @test SM.mean_solar_distance_m(SM.Mars()) ≈ 1.5237 * SM.AU_M
    end

    @testset "battery rate, shedding tiers and the step-end update" begin
        b = Battery(capacity_wh=10.0, charge_efficiency=0.9, discharge_efficiency=0.8, max_charge_w=20.0)
        @test SM.battery_rate_wh_s(b, 10.0) ≈ 10.0 * 0.9 / 3600
        @test SM.battery_rate_wh_s(b, 50.0) ≈ 20.0 * 0.9 / 3600                 # charge limit
        @test SM.battery_rate_wh_s(b, -8.0) ≈ -8.0 / 0.8 / 3600
        b.full = true;  @test SM.battery_rate_wh_s(b, 10.0) == 0.0 && SM.battery_rate_wh_s(b, -8.0) < 0.0
        b.full = false; b.empty = true
        @test SM.battery_rate_wh_s(b, -8.0) == 0.0 && SM.battery_rate_wh_s(b, 10.0) > 0.0
        @test SM.shed_tier(100.0, 5.0, (10.0, 3.0, 2.0, 0.0)) == 0
        @test SM.shed_tier(12.0, 5.0, (10.0, 3.0, 2.0, 0.0)) == 1                # lasers off: 10 W
        @test SM.shed_tier(8.0, 5.0, (10.0, 3.0, 2.0, 0.0)) == 2                 # and rods: 7 W
        @test SM.shed_tier(6.0, 5.0, (10.0, 3.0, 2.0, 0.0)) == 3
        root = SM.Link(root=true, m=10.0)
        sc = SM.SpacecraftModel(links=[root], root=root)
        add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0)))
        add_actuator!(sc, root, Magnetorquer(axis=SVector(0.0, 1.0, 0.0), max_dipole_am2=1.0, power_at_max_dipole_w=3.0))
        add_actuator!(sc, root, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=2.0, efficiency=0.2))
        bat = add_power!(sc, root, Battery(capacity_wh=10.0))
        levels = [0.01, 1.0, 2.0]
        @test SM.shed_levels(sc, levels) === levels
        bat.shed_tier = 1; @test SM.shed_levels(sc, levels) == [0.01, 1.0, 0.0]
        bat.shed_tier = 2; @test SM.shed_levels(sc, levels) == [0.01, 0.0, 0.0]
        bat.shed_tier = 3; @test SM.shed_levels(sc, levels) == [0.0, 0.0, 0.0]
        bat.shed_tier = 0
        balance = SM.power_balance(sc, levels, 20.0)
        @test balance.load_w ≈ 3.0 + 10.0 && balance.net_w ≈ 20.0 - 13.0 && balance.by_tier == (10.0, 3.0, 0.0, 0.0)
        # Empty: shed to what the arrays alone carry; recharged to recovery_soc: back to normal.
        bat.empty = true
        @test SM.update_power_state!(sc, 0.0, 5.0, balance.by_tier) && bat.shed_tier == 1
        @test SM.power_balance(sc, levels, 5.0).load_w ≈ 3.0
        @test !SM.update_power_state!(sc, 0.4, 5.0, balance.by_tier) && bat.empty    # under 5 % of 10 Wh
        @test SM.update_power_state!(sc, 0.6, 5.0, balance.by_tier) && !bat.empty && bat.shed_tier == 0
        @test SM.battery_boundary!(sc, 0.0) && bat.empty && !SM.battery_boundary!(sc, 0.0)
    end

    @testset "a run charges in sunlight, sheds the laser when the battery is empty, and recovers" begin
        table = power_run(power_spacecraft()).table
        soc, gen, load, shed, laser = table.sc1_battery_soc, table.sc1_power_generated_w, table.sc1_power_load_w,
            table.sc1_power_shed_w, table.sc1_actuator_level_1
        lit = 24.5
        @test gen[1] ≈ 0.06 * 0.3 * SM.solar_flux_w_m2(SM.mean_solar_distance_m(PLANET) - ORBIT_R) rtol=1e-3
        @test all(g -> isapprox(g, gen[1]; rtol=1e-3) || g < 1e-9 || g < gen[1], gen)                        # lit or eclipsed
        # Charging at 9.5 W from half full: 2.5 Wh in about 947 s.
        k = findfirst(t -> t >= 450.0, table.time)
        @test soc[k] ≈ 0.5 + (gen[1] - 15.0) * table.time[k] / 3600 / 5.0 rtol=1e-3
        @test all(<=(1.0 + 1e-9), soc) && all(>=(-1e-9), soc)                                             # within its capacity
        full = findall(>=(1.0 - 1e-9), soc)
        @test !isempty(full) && all(i -> gen[i] > 15.0, full)                                             # full only in sunlight
        empty = findall(<=(1e-9), soc)
        @test !isempty(empty)
        @test all(i -> laser[i] == 0.0 && shed[i] ≈ 10.0 && load[i] ≈ 5.0, empty)                          # laser shed while empty
        @test all(i -> laser[i] == 2.0, findall(>(0.06), soc))                                             # on otherwise
        # After the eclipse the Sun carries the full load again and the battery recharges.
        last_empty = maximum(empty)
        @test any(i -> i > last_empty && soc[i] > 0.1 && laser[i] == 2.0, eachindex(soc))
    end

    @testset "an empty battery sheds the wheels' motors in the engine, not only in the columns" begin
        # Arrays facing away and an empty battery: everything sheddable is shed, so
        # the wheel commanded to 0.01 N m does not turn the body.
        shed = power_run(power_spacecraft(soc=0.0, wheel=true, facing=-1.0); mission_time=20.0, levels=[0.01, 2.0], return_solution=true)
        @test all(==(0.0), shed.table.sc1_actuator_level_1) && all(==(0.0), shed.table.sc1_actuator_level_2)
        @test norm(shed.solution.u[end].sc[1].ω) < 1e-12
        free = power_run(power_spacecraft(battery=false, wheel=true, facing=-1.0); mission_time=20.0, levels=[0.01, 2.0], return_solution=true)
        @test norm(free.solution.u[end].sc[1].ω) > 1e-3                         # without a battery nothing is shed
    end

    @testset "the load counts the actuators of every mounted control effector" begin
        # Two effectors each command the laser to 2 W, and the dynamics apply
        # both: the bus carries both draws (10 W each) besides the 5 W load.
        table = power_run(power_spacecraft(); mission_time=10.0, effectors=(HoldLevels([2.0]), HoldLevels([2.0]))).table
        @test table.sc1_power_load_w[1] ≈ 25.0
        @test table.sc1_power_net_w[1] ≈ table.sc1_power_generated_w[1] - 25.0
        @test table.sc1_actuator_power[1] ≈ 10.0                                 # the actuator columns follow the first effector
    end

    @testset "checkpoint segments keep an empty battery shedding until it recovers" begin
        # 12.25 W of arrays against a 15 W load: the battery starts empty, the
        # laser is shed while the arrays recharge it, and shedding lasts until
        # 20 % of its capacity (about 500 s), past the 50 s segments.
        sc() = power_spacecraft(soc=0.0, area=0.03, recovery_soc=0.2)
        dir = mktempdir()
        straight = power_run(sc(); mission_time=300.0).table
        split = power_run(sc(); mission_time=300.0, settings=(checkpoint_enabled=true, checkpoint_interval_s=50.0, checkpoint_directory=dir)).table
        @test all(==(0.0), straight.sc1_actuator_level_1)
        @test all(==(0.0), split.sc1_actuator_level_1)
        @test split.sc1_battery_wh[end] ≈ straight.sc1_battery_wh[end] rtol=1e-6
        @test 0.1 < straight.sc1_battery_soc[end] < 0.2                            # recharging, still below recovery
        # A run resumed from the file keeps shedding: the file holds the flags.
        resumed_dir = mktempdir()
        first_leg = power_run(sc(); mission_time=150.0, settings=(checkpoint_enabled=true, checkpoint_interval_s=50.0, checkpoint_directory=resumed_dir)).table
        resumed = power_run(sc(); mission_time=300.0, settings=(checkpoint_enabled=true, checkpoint_interval_s=50.0, checkpoint_directory=resumed_dir,
            resume_from_checkpoint=true, results=true, results_directory=resumed_dir)).table
        # It resumed at 150 s rather than starting again, and its table holds the
        # whole run: the first leg's rows, kept beside its checkpoint, then its own.
        @test TOML.parsefile(joinpath(resumed_dir, "simulation_results.manifest.toml"))["resumed_from_s"] == 150.0
        @test resumed.time[1] == 0.0
        @test isequal(resumed[resumed.time .< 150.0, :], first_leg[first_leg.time .< 150.0, :])
        after_resume = resumed.time .>= 150.0 - 1e-9
        @test any(after_resume) && all(==(0.0), resumed.sc1_actuator_level_1[after_resume])   # still shedding after the resume
        @test resumed.sc1_battery_wh[end] ≈ split.sc1_battery_wh[end] rtol=1e-6
        # Another battery configuration, or a file without the callbacks' state, is refused.
        @test_throws ArgumentError power_run(power_spacecraft(soc=0.0, area=0.03, recovery_soc=0.3); mission_time=300.0,
            settings=(resume_from_checkpoint=true, checkpoint_directory=resumed_dir))
        args = SM.SimConfig._with_configuration(make_example_config(planet=PLANET, spacecraft=sc(), mission_time=10.0, initial_time=EPOCH,
            dynamic_effectors=(), ephemerides_model=EPHEM, orientation_sim=true, keplerian=false, EI_km=100.0, verbose=false, results=false))
        flag = (shared_buffers=(callback_resume=Ref(false),),)
        @test_throws ArgumentError SpaceAGORA.SimulationEngine._restore_callback_checkpoint_state!(args, (t=50.0, runtime_state=nothing), flag)
        @test !flag.shared_buffers.callback_resume[]
    end

    @testset "power parts are opt-in; split solver modes refuse a battery" begin
        plain = SM.Link(root=true, m=10.0)
        sc = SM.SpacecraftModel(links=[plain], root=plain, inertia_tensor=plain.inertia,
            initial_condition=SM.CartesianInitialCondition(SVector(ORBIT_R, 0.0, 0.0), SVector(0.0, 7600.0, 0.0)), id=1)
        names_plain = names(power_run(sc; mission_time=10.0, levels=Float64[]).table)
        @test !any(n -> occursin("battery", n) || occursin("power_", n), names_plain)
        no_battery = names(power_run(power_spacecraft(battery=false); mission_time=10.0).table)
        @test "sc1_power_generated_w" in no_battery && "sc1_battery_soc" in no_battery
        @test_throws ArgumentError power_run(power_spacecraft(); mission_time=10.0, solver=:multirate)
    end
end

end # module PowerBudgetTests
