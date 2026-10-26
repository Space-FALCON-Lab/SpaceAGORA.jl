module InterlinkGroundStationTests

# Ground stations: the geodetic position, a station's state sample, the
# elevation mask, the station and its refusals, registration, pairs, the
# policies at a station end, a pass against the closed-form elevation, the
# planet's frame inside a step, an uplink, the setup checks, the viewer block,
# checkpoints, the stations' data counters and the stations with data stores,
# and a negative control.

using Test
using JSON
using Arrow
using DataFrames
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using OrdinaryDiffEq: DiscreteCallback

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const IM = SM.InterLinkModels
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const BUDGET = ShannonBudget(bandwidth_hz=1e6, frequency_hz=2.2e9)
const LINK = DataLink(transmit_power_w=1.0, budget=BUDGET, range_m=3000e3)
struct StationStoreInterrupted <: Exception end

# A station-compatible kind whose terminal-aware hooks cannot fall back silently.
struct StationWrenchProbe <: AbstractInterlink end
# A radio kind with different transmitter and receiver hardware requirements.
struct StationRoleProbe <: AbstractInterlink end
SpaceAGORA.link_range_m(::StationRoleProbe) = 1e8
SpaceAGORA.link_directed(::StationRoleProbe) = true
SpaceAGORA.terminal_supports(::RadioTerminal, ::StationRoleProbe) = true
SpaceAGORA.terminal_supports(terminal::RadioTerminal, ::StationRoleProbe, side::Symbol) =
    side === :transmitter ? terminal.gain_dbi == 0.0 : terminal.gain_dbi == 3.0

const STATION_HOOK_TERMINALS = Ref{Any}(nothing)
SpaceAGORA.link_range_m(::StationWrenchProbe) = 1e8
SpaceAGORA.link_directed(::StationWrenchProbe) = true
SpaceAGORA.link_applies_torque(::StationWrenchProbe) = true
SpaceAGORA.terminal_supports(::RadioTerminal, ::StationWrenchProbe) = true
SpaceAGORA.link_available(::StationWrenchProbe, a, b, ctx) = false
SpaceAGORA.link_available(::StationWrenchProbe, ta, tb, a, b, ctx) =
    (STATION_HOOK_TERMINALS[] = (ta, tb); ta isa RadioTerminal && tb isa RadioTerminal)
function SpaceAGORA.interlink_wrench(::StationWrenchProbe, ta, tb, a, b, ctx)
    @assert ta isa RadioTerminal && tb isa RadioTerminal
    STATION_HOOK_TERMINALS[] = (ta, tb)
    force = normalize(b.pos_ii - a.pos_ii)
    return force, SVector(1.0, 2.0, 3.0), -force, SVector(4.0, 5.0, 6.0)
end

radio(; kwargs...) = RadioTerminal(; boresight=[0.0, 0.0, 1.0], max_power_w=5.0, gain_dbi=0.0, efficiency=0.5, kwargs...)
station(; name="north", lat_deg=0.0, lon_deg=0.0, alt_m=0.0, min_elevation_deg=10.0, antennas=[radio()]) =
    GroundStation(; name, lat_deg, lon_deg, alt_m, min_elevation_deg, antennas)

# A 100 kg bus at `pos` with velocity `vel`, with a radio and, when asked, a battery.
function radio_sat(index; pos=SVector(RADIUS, 0.0, 0.0), vel=SVector(0.0, sqrt(PLANET.μ / RADIUS), 0.0),
                   terminal=radio(), n_terminal=0, battery=false)
    bus = Link(root=true, m=100.0)
    ic = CartesianInitialCondition(SVector{3, Float64}(pos), SVector{3, Float64}(vel); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=n_terminal)
    terminal === nothing || add_actuator!(sc, bus, terminal)
    battery && add_power!(sc, bus, Battery(capacity_wh=10.0, initial_soc=0.5))
    return sc
end

# A planet whose frame (J2000 to planet-fixed) is a turn of `θ` about the pole.
function turned_planet(θ)
    planet = make_no_gram_planet(:earth)
    planet.L_PI .= [cos(θ) sin(θ) 0.0; -sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
    return planet
end

# The inertial position `range_m` from station `st` at elevation `el` and azimuth
# `az` (rad, from north toward east), under the planet's frame.
function seen_from(st, planet, el, az, range_m)
    lat, lon = deg2rad(st.lat_deg), deg2rad(st.lon_deg)
    up = SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    east = SVector(-sin(lon), cos(lon), 0.0)
    north = cross(up, east)
    position = SM.Geodesy.geodetic_to_pp(lat, lon, st.alt_m, planet)
    direction = cos(el) * (cos(az) * north + sin(az) * east) + sin(el) * up
    return transpose(SMatrix{3, 3, Float64}(planet.L_PI)) * (position + range_m * direction)
end

state_at(positions...) = (sc=[(pos=p, vel=SVector(0.0, 7.5e3, 0.0), mass=100.0) for p in positions],)

# The closed-form Shannon rate of `BUDGET` at `d` between two 0 dBi radios at 500 K.
function shannon(d; power_w=1.0)
    λ = 299_792_458.0 / BUDGET.frequency_hz
    received = power_w * (λ / (4π * d))^2
    noise = 1.380649e-23 * 500.0 * BUDGET.bandwidth_hz
    return BUDGET.implementation_efficiency * BUDGET.bandwidth_hz * log2(1 + received / noise)
end

function simpson(f, a, b, n)
    h = (b - a) / n
    total = f(a) + f(b)
    for k in 1:n-1
        total += (isodd(k) ? 4 : 2) * f(a + k * h)
    end
    return total * h / 3
end

function run_args(spacecraft, model; duration=60.0, dt_max=10.0, policy=DataThroughputPolicy(),
                  settings=SimulationSettings(results=false, generate_plots=false))
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=false, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=dt_max),
        solver_config=SM.SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=policy)
end

run(spacecraft, model; kwargs...) = run_simulation(run_args(spacecraft, model; kwargs...);
    return_results=true, return_solution=true)
refusal(f) = try f(); nothing catch err err end

# The planet's frame at `t` s into a run, as the ephemerides give it.
const EPHEMERIDES = make_no_gram_environment(planet=PLANET).ephemerides_model
const ET0 = SM.EphemeridesModels.ephemerides_epoch_seconds(EPHEMERIDES, InitialTime())
frame(t) = SMatrix{3, 3, Float64}(SM.EphemeridesModels.ephemerides_frame_lpi(EPHEMERIDES, PLANET, ET0 + t))

# An equatorial station at longitude 0 and one spacecraft on a circular
# equatorial orbit, starting `lead_deg` west of it, with the link `link` (a
# downlink by default; `uplink=true` reverses it).
function pass_fixture(; lead_deg=25.0, st=station(), link=LINK, uplink=false, battery=false)
    L0 = frame(0.0)
    u = transpose(L0) * normalize(SM.Geodesy.geodetic_to_pp(0.0, 0.0, 0.0, PLANET))
    pole = transpose(L0) * SVector(0.0, 0.0, 1.0)
    φ = -deg2rad(lead_deg)
    r̂ = cos(φ) * u + sin(φ) * cross(pole, u)
    spacecraft = [radio_sat(1; pos=RADIUS * r̂, vel=sqrt(PLANET.μ / RADIUS) * cross(pole, r̂), battery)]
    model = InterlinkModel()
    add_ground_station!(model, st)
    uplink ? register_link!(model, spacecraft, (-1, 1), (1, 1), link) : register_link!(model, spacecraft, (1, 1), (-1, 1), link)
    return spacecraft, model
end

# The elevation of spacecraft 1 above the station of `pass_fixture` at `t`, from the solution.
const SITE = SM.Geodesy.geodetic_to_pp(0.0, 0.0, 0.0, PLANET)
pass_elevation(solution, t) = (d = frame(t) * SVector{3, Float64}(solution(t).sc[1].pos) - SITE; asin(d[1] / norm(d)))

# The data the link `key` of a `pass_fixture` run carried, by quadrature: from each
# history sample that selects it to the next one, at the closed-form rate between
# the spacecraft of the solution and the station under the planet's frame at each
# time, or, with `held`, under the frame at each step's start; up to `until`.
function carried(solution, key; held::Bool=false, until::Float64=Inf)
    history = solution.prob.p.args.interlink_model.history
    rate(t, L) = shannon(norm(SVector{3, Float64}(solution(t).sc[1].pos) - transpose(L) * SITE))
    total = 0.0
    for k in 1:length(history) - 1
        key in history[k].active || continue
        t0, t1 = history[k].time, min(history[k + 1].time, until)
        t0 < t1 || continue
        total += simpson(t -> rate(t, held ? frame(t0) : frame(t)), t0, t1, 64)
    end
    return total
end

@testset "Ground stations" begin
    @testset "station terminal-aware hooks retain frame and hardware" begin
        for uplink in (false, true)
            craft_radio = RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=5.0,
                gain_dbi=3.0, efficiency=0.5)
            station_radio = RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=5.0,
                gain_dbi=7.0, efficiency=0.5)
            spacecraft = [radio_sat(1; terminal=craft_radio)]
            model = InterlinkModel()
            st = station(antennas=[station_radio])
            add_ground_station!(model, st)
            ends = uplink ? ((-1, 1), (1, 1)) : ((1, 1), (-1, 1))
            register_link!(model, spacecraft, ends..., StationWrenchProbe())
            craft_terminal = only(mounted_link_terminals(only(spacecraft))).terminal
            expected = uplink ? (only(st.antennas), craft_terminal) : (craft_terminal, only(st.antennas))
            @test IM._edge_terminals(model, 1, spacecraft) === expected
            ctx = LinkContext(planet=turned_planet(0.0))
            state = state_at(SVector(RADIUS, 0.0, 0.0))
            @test update_availability!(model, spacecraft, state; context=ctx) == [true]
            @test STATION_HOOK_TERMINALS[] === expected
            @test SM.InterlinkTerminals.prepare_terminal_endpoints!(model, spacecraft)
            @test IM._edge_terminals(model, 1, spacecraft) === expected
            @test update_availability!(model, spacecraft, state; context=ctx) == [true]
            @test STATION_HOOK_TERMINALS[] === expected
            model.graph.active .= true
            explicit_frame = SMatrix{3, 3, Float64}(turned_planet(0.4).L_PI)
            station_pos = transpose(explicit_frame) * first(IM._station_site(st, ctx.planet))
            expected_force = normalize(station_pos - state.sc[1].pos)
            force, torque, turned = IM.interlink_wrench_on_spacecraft(model, state, 1;
                spacecraft, context=ctx, planet_frame=explicit_frame, is_active=[true])
            @test STATION_HOOK_TERMINALS[] === expected
            @test force ≈ expected_force rtol=1e-14
            @test torque == (uplink ? SVector(4.0, 5.0, 6.0) : SVector(1.0, 2.0, 3.0))
            @test turned
            @test interlink_force_on_spacecraft(model, state, 1; spacecraft, context=ctx,
                planet_frame=explicit_frame, is_active=[true]) == force
            @test IM.interlink_wrench_on_spacecraft(model, state, 1; spacecraft, context=ctx,
                planet_frame=explicit_frame, is_active=[false]) == (zeros(3), zeros(3), false)
            config = SE._interlink_checkpoint_configuration(run_args(spacecraft, model))
            @test hasproperty(config, :ground_stations) && hasproperty(config, :terminal_mounts)
            @test only(config.terminal_mounts).endpoint == (1, 1)
        end
    end

    @testset "geodetic position" begin
        for planet in (PLANET, make_no_gram_planet(:mars))
            @test SM.Geodesy.geodetic_to_pp(0.0, 0.0, 100.0, planet) ≈ [planet.Rp_e + 100.0, 0.0, 0.0] rtol=1e-15
            @test SM.Geodesy.geodetic_to_pp(0.0, π / 2, 0.0, planet) ≈ [0.0, planet.Rp_e, 0.0] atol=1e-6
            @test SM.Geodesy.geodetic_to_pp(π / 2, 0.0, 50.0, planet) ≈ [0.0, 0.0, planet.Rp_p + 50.0] atol=1e-6
            # The round trip through `rtolatlong`.
            for (lat, lon, alt) in ((42.3, -83.7, 270.0), (-60.0, 150.0, 3000.0), (89.0, 10.0, 0.0), (0.5, 179.0, -50.0))
                back = SM.FrameTransforms.rtolatlong(SM.Geodesy.geodetic_to_pp(deg2rad(lat), deg2rad(lon), alt, planet), planet)
                @test back[1] ≈ alt atol=1e-5
                @test back[2] ≈ deg2rad(lat) atol=1e-11
                @test back[3] ≈ deg2rad(lon) atol=1e-11
            end
        end
    end

    @testset "a station's state sample" begin
        planet = turned_planet(0.6)
        L = SMatrix{3, 3, Float64}(planet.L_PI)
        st = station(lon_deg=30.0)
        sample = IM._station_sample(st, 2, planet, nothing)
        site = SM.Geodesy.geodetic_to_pp(0.0, deg2rad(30.0), 0.0, planet)
        # Its inertial position is `L_PIᵀ r`, and it moves at `ω R` toward the planet-fixed east.
        @test L * sample.pos_ii ≈ site rtol=1e-14
        @test L * sample.vel_ii ≈ cross(planet.ω, site) rtol=1e-14
        @test norm(sample.vel_ii) ≈ planet.ω[3] * planet.Rp_e rtol=1e-12
        @test dot(L * sample.vel_ii, [-sin(deg2rad(30.0)), cos(deg2rad(30.0)), 0.0]) ≈ norm(sample.vel_ii) rtol=1e-12
        @test isnan(sample.mass_kg) && sample.q_ib === nothing && sample.ω_body === nothing
        @test sample.spacecraft === st && sample.spacecraft_index == -2
        # A frame given in place of the planet's current one.
        other = SMatrix{3, 3, Float64}(turned_planet(1.1).L_PI)
        @test other * IM._station_sample(st, 2, planet, other).pos_ii ≈ site rtol=1e-14
        err = refusal(() -> IM._station_sample(st, 1, nothing, nothing))
        @test err isa ArgumentError && occursin("needs the run's planet", err.msg)
    end

    @testset "the elevation mask" begin
        planet = turned_planet(0.3)
        for st in (station(lat_deg=60.0, lon_deg=-20.0, alt_m=500.0), station(lat_deg=-35.0, lon_deg=148.0))
            for (el, az) in ((0.2, 0.0), (deg2rad(10.05), 1.0), (deg2rad(9.95), 4.0), (1.4, 2.5), (-0.1, 3.0))
                # The geodetic up, not the geocentric radial: they differ by a few mrad at these latitudes.
                @test IM._station_elevation_rad(st, planet, nothing, seen_from(st, planet, el, az, 1500e3)) ≈ el atol=1e-12
            end
        end
        st = station(lat_deg=60.0, lon_deg=-20.0)
        spacecraft = [radio_sat(1)]
        model = InterlinkModel()
        add_ground_station!(model, st)
        register_link!(model, spacecraft, (1, 1), (-1, 1), LINK)
        available(el, az; range_m=1500e3, kwargs...) = only(update_availability!(model, spacecraft,
            state_at(seen_from(st, planet, el, az, range_m)); context=LinkContext(planet=planet), kwargs...))
        # Overhead and just above the mask: available, also through the planet's sphere of equatorial
        # radius, which a station at 60° latitude sits inside.
        @test available(π / 2, 0.0)
        @test available(deg2rad(10.05), 0.0) && available(deg2rad(10.05), π)
        # Just below the mask, north and south: not.
        @test !available(deg2rad(9.95), 0.0) && !available(deg2rad(9.95), π)
        # The link's own range applies too.
        @test !available(π / 3, 0.0; range_m=3100e3)
        # A station is always active; the spacecraft's own flag applies.
        @test available(π / 2, 0.0; is_active=[true]) && !available(π / 2, 0.0; is_active=[false])
        # A forbidden spacecraft-station pair.
        forbidden = InterlinkModel(forbidden_pairs=[(1, -1)])
        add_ground_station!(forbidden, st)
        register_link!(forbidden, spacecraft, (1, 1), (-1, 1), LINK)
        @test (-1, 1) in forbidden.forbidden_pairs
        @test !only(update_availability!(forbidden, spacecraft, state_at(seen_from(st, planet, π / 2, 0.0, 1500e3));
            context=LinkContext(planet=planet)))
    end

    @testset "the station" begin
        st = station(name="x", lat_deg=10, lon_deg=-200.0, alt_m=-30.0, min_elevation_deg=0.0)
        @test st.name == "x" && st.lat_deg === 10.0 && st.lon_deg == -200.0 && st.alt_m == -30.0 && st.min_elevation_deg == 0.0
        defaults = GroundStation(name="d", lat_deg=1.0, lon_deg=2.0, antennas=[radio()])
        @test defaults.alt_m == 0.0 && defaults.min_elevation_deg == 10.0
        for bad in ((name="",), (lat_deg=90.5,), (lat_deg=NaN,), (lon_deg=Inf,), (alt_m=NaN,), (min_elevation_deg=-1.0,),
                    (min_elevation_deg=90.0,), (antennas=RadioTerminal[],),
                    (antennas=[LaserTerminal(boresight=[1.0, 0.0, 0.0], max_power_w=1.0)],))
            @test_throws ArgumentError station(; bad...)
        end
        # The antennas are the station's own copy.
        list = [radio()]
        copied = station(antennas=list)
        push!(list, radio())
        @test length(copied.antennas) == 1
    end

    @testset "registration" begin
        spacecraft = [radio_sat(1), radio_sat(2), radio_sat(3; terminal=nothing, n_terminal=1)]
        model = InterlinkModel()
        # Two antennas with equal fields cannot be told apart, so these differ.
        dish = radio()
        @test add_ground_station!(model, station(name="a", antennas=[dish, radio(gain_dbi=3.0)])) == 1
        @test add_ground_station!(model, station(name="b")) == 2
        err = refusal(() -> add_ground_station!(model, station(name="a")))
        @test err isa ArgumentError && occursin("already has a ground station named \"a\"", err.msg)
        # A downlink keeps its order, an uplink takes the antenna itself.
        @test register_link!(model, spacecraft, (1, 1), (-1, 1), LINK) == ((1, 1), (-1, 1))
        @test register_link!(model, spacecraft, (-1, dish), (2, 1), LINK) == ((-1, 1), (2, 1))
        @test register_link!(model, spacecraft, (2, 1), (-2, 1), LINK) == ((2, 1), (-2, 1))
        @test model.graph.satellite_edges == [[1], [2, 3], Int[]]
        @test model.graph.vertices == [(1, 1), (-1, 1), (2, 1), (-2, 1)]
        for (first, second, link, needle) in (
                ((-1, 1), (-2, 1), LINK, "cannot join two ground stations"),
                ((1, 1), (-3, 1), LINK, "Invalid ground station 3"),
                ((1, 1), (-2, 2), LINK, "Invalid antenna 2 on ground station 2 (b)"),
                ((1, 1), (-1, radio(gain_dbi=7.0)), LINK, "not an antenna of ground station a"),
                ((3, 1), (-1, 2), LaserMomentumLink(), "Antenna 2 of ground station 1 (a), a RadioTerminal, cannot carry"),
                ((1, 1), (-1, 1), LINK, "already registered"))
            err = refusal(() -> register_link!(model, spacecraft, first, second, link))
            @test err isa ArgumentError && occursin(needle, err.msg)
        end
        # The station's transmitter at a level beyond its antenna's limits is refused at setup.
        loud = InterlinkModel()
        add_ground_station!(loud, station())
        register_link!(loud, spacecraft, (-1, 1), (1, 1), DataLink(transmit_power_w=10.0, budget=BUDGET, range_m=3e6))
        err = refusal(() -> SM.InterlinkTerminals.prepare_terminal_endpoints!(loud, spacecraft))
        @test err isa ArgumentError && occursin("antenna 1 of ground station 1 sets it to level 10.0", err.msg)
        # A station's links take no spacecraft's terminal endpoints, and its antennas are never shed.
        @test SM.InterlinkTerminals.prepare_terminal_endpoints!(model, spacecraft)
        @test model.terminal_endpoints[1] == [(1, 1, 1)] && model.terminal_endpoints[2] == [(2, 2, 1), (3, 1, 1)]
        @test SM.InterlinkTerminals.terminal_availability(model, spacecraft) == fill(true, 4)
    end

    @testset "pairs" begin
        spacecraft = [radio_sat(1), radio_sat(2)]
        model = InterlinkModel()
        add_ground_station!(model, station(name="a"))
        add_ground_station!(model, station(name="b", antennas=[radio(), radio()]))
        @test station_pairs(spacecraft, model) == [(1, -1), (1, -2), (2, -1), (2, -2)]
        @test station_pairs(spacecraft, model.ground_stations) == station_pairs(spacecraft, model)
        # Each station's links take its antennas in turn.
        @test register_links!(model, spacecraft, station_pairs(spacecraft, model), LINK) ==
            [((1, 1), (-1, 1)), ((1, 1), (-2, 1)), ((2, 1), (-1, 1)), ((2, 1), (-2, 2))]
        @test register_links!(model, spacecraft, reverse.(station_pairs(spacecraft, model)), LINK) ==
            [((-1, 1), (1, 1)), ((-2, 1), (1, 1)), ((-1, 1), (2, 1)), ((-2, 2), (2, 1))]
        err = refusal(() -> register_links!(model, spacecraft, [(1, -3)], LINK))
        @test err isa ArgumentError && occursin("Invalid ground station 3", err.msg)
        err = refusal(() -> register_links!(model, spacecraft, [(-1, -2)], LINK))
        @test err isa ArgumentError && occursin("cannot join two ground stations", err.msg)
    end

    @testset "automatic station hardware selection" begin
        # Shannon links skip a spacecraft's virtual terminal in both directions.
        spacecraft = [radio_sat(1; n_terminal=1)]
        model = InterlinkModel()
        add_ground_station!(model, station())
        @test register_links!(model, spacecraft, [(1, -1), (-1, 1)], LINK) ==
            [((1, 2), (-1, 1)), ((-1, 1), (1, 2))]
        # A station antenna must satisfy the same role-aware hook as a spacecraft terminal.
        for pair in ((1, -1), (-1, 1))
            @testset "direction $pair" begin
                spacecraft = [radio_sat(1)]
                add_actuator!(spacecraft[1], spacecraft[1].root, radio(gain_dbi=3.0))
                model = InterlinkModel()
                add_ground_station!(model, station(antennas=[radio(), radio(gain_dbi=3.0)]))
                @test register_links!(model, spacecraft, [pair], StationRoleProbe()) ==
                    [((pair[1], 1), (pair[2], 2))]
            end
        end
        # A receiver-only station cannot serve as a transmitter.
        model = InterlinkModel()
        add_ground_station!(model, station(antennas=[radio(gain_dbi=3.0)]))
        @test_throws ArgumentError register_links!(model, [radio_sat(1)], [(-1, 1)], StationRoleProbe())
    end

    @testset "the policies at a station end" begin
        planet = turned_planet(0.2)
        st = station(lat_deg=20.0, lon_deg=40.0)
        spacecraft = [radio_sat(1), radio_sat(2)]
        model = InterlinkModel()
        add_ground_station!(model, st)
        first_down = register_link!(model, spacecraft, (1, 1), (-1, 1), LINK)
        second_down = register_link!(model, spacecraft, (2, 1), (-1, 1), LINK)
        state = state_at(seen_from(st, planet, 1.2, 0.0, 800e3), seen_from(st, planet, 0.5, 2.0, 1200e3))
        context = LinkContext(time_s=0.0, planet=planet)
        scheduler = SchedulingPolicyModel(DataThroughputPolicy())
        # One antenna carries one link: the nearer spacecraft, with the higher rate, wins it.
        @test schedule_interlinks!(model, scheduler, spacecraft, state, planet.μ; context) == [first_down]
        @test scheduler.matching.score ≈ [shannon(800e3), shannon(1200e3)] rtol=1e-9
        # The candidate's station end: the station's sample and its antenna, with or without the spacecraft models.
        for models in (spacecraft, nothing)
            candidate = SM.SchedulingPolicies._link_candidate(model, 1,
                SchedulingContext(planet.μ, context, models, state, nothing))
            @test candidate.key == first_down && candidate.b.spacecraft === st && candidate.b.spacecraft_index == -1
            @test candidate.terminal_b === st.antennas[1]
            @test (candidate.terminal_a === nothing) == (models === nothing)
        end
        # The data the selected link carries counts at its spacecraft end (the run's hardware endpoints
        # give the spacecraft's radio).
        SM.InterlinkTerminals.prepare_terminal_endpoints!(model, spacecraft)
        sent, received = SM.data_link_rates_bps(model, spacecraft, state; context)
        @test sent ≈ [shannon(800e3), 0.0] rtol=1e-9
        @test received == [0.0, 0.0]
        # A station is never an orbital target and has no battery: these score a data link 0, without reading it.
        ctx = SchedulingContext(planet.μ, context, spacecraft, state, nothing)
        @test isnan(state_of_charge(ctx, -1))
        for policy in (SemimajorAxisRatePolicy(targets=:both), SemimajorAxisRatePolicy(targets=:lower_index),
                       SemimajorAxisRatePolicy(targets=[1, 2]), EccentricityRatePolicy(targets=:both),
                       EnergyBalancingPolicy(), PowerDeliveryPolicy())
            scored = SchedulingPolicyModel(policy)
            score_candidates!(model, scored, state, planet.μ; spacecraft, context)
            @test scored.matching.score == [0.0, 0.0]
        end
        mix = SchedulingPolicyModel(WeightedPolicy(PolicyTerm(DataThroughputPolicy(); scale=1e3),
            PolicyTerm(SemimajorAxisRatePolicy(targets=:both); scale=1.0)))
        score_candidates!(model, mix, state, planet.μ; spacecraft, context)
        @test mix.matching.score ≈ [shannon(800e3), shannon(1200e3)] ./ 1e3 rtol=1e-9
        # A caller's candidate force receives a station partner as its state sample.
        partners = Any[]
        force(target, partner, link) = (push!(partners, partner); SVector(0.0, 0.0, 0.0))
        score_candidates!(model, SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=:lower_index)), state, planet.μ;
            spacecraft, context, candidate_force=force)
        @test length(partners) == 2 && all(partner -> partner.spacecraft === st, partners)
    end

    @testset "a pass" begin
        spacecraft, model = pass_fixture()
        solution = run(spacecraft, model; duration=900.0, dt_max=30.0).solution
        history = solution.prob.p.args.interlink_model.history
        key = ((1, 1), (-1, 1))
        mask = deg2rad(10.0)
        # Available at exactly the accepted steps where the closed-form elevation clears the mask.
        elevations = [pass_elevation(solution, sample.time) for sample in history]
        flags = [key in sample.available for sample in history]
        clear = abs.(elevations .- mask) .> 1e-7
        @test flags[clear] == (elevations .>= mask)[clear]
        @test !flags[1] && any(flags) && !flags[end]
        window = only(contact_windows(history))
        rise, set = findfirst(flags), findlast(flags)
        @test window.key == key && window.start_s == history[rise].time && window.end_s == history[set + 1].time
        @test 100.0 < window.start_s < 170.0 && 720.0 < window.end_s < 790.0
        # The throughput policy selects it for the whole pass.
        @test only(contact_windows(history; active=true)).start_s == window.start_s
        # The link schedule shows the station end as -1.
        table = link_schedule(solution.prob.p.args.interlink_model)
        @test all(==(-1), table.b_spacecraft) && all(==(1), table.b_terminal) && all(==(1), table.a_spacecraft)
        # The downlink counts what it sent at the spacecraft end.
        @test solution.u[end].sc[1].data_sent_bits > 0 && solution.u[end].sc[1].data_received_bits == 0
    end

    @testset "the stations' counters" begin
        # A downlink over the pass: the station receives what the spacecraft sends, at the closed-form rate.
        spacecraft, model = pass_fixture()
        result = run(spacecraft, model; duration=900.0, dt_max=30.0)
        final = result.solution.u[end]
        @test propertynames(final.gs[1]) == (:received_bits, :sent_bits) && length(final.gs) == 1
        @test final.gs[1].received_bits ≈ final.sc[1].data_sent_bits rtol=1e-14
        @test final.gs[1].sent_bits == 0.0
        @test final.gs[1].received_bits ≈ carried(result.solution, ((1, 1), (-1, 1))) rtol=1e-7
        # The columns, from the integrated counters.
        table = result.table
        @test table.gs1_received_bits ≈ table.sc1_data_sent_bits rtol=1e-14
        @test all(iszero, table.gs1_sent_bits) && issorted(table.gs1_received_bits)
        @test table.gs1_received_bits[end] == final.gs[1].received_bits
        # An uplink: the station sends what the spacecraft receives.
        up = run(pass_fixture(uplink=true)...; duration=430.0, dt_max=30.0).solution.u[end]
        @test up.gs[1].sent_bits ≈ up.sc[1].data_received_bits rtol=1e-14
        @test up.gs[1].sent_bits > 0 && up.gs[1].received_bits == 0.0
        # Each station counts its own links: two stations, one spacecraft sending to the first only.
        two = pass_fixture()[2]
        add_ground_station!(two, station(name="south", lat_deg=-60.0))
        @test length(SE.build_initial_conditions(run_args(pass_fixture()[1], two)).gs) == 2
        # Without data links, or without stations, there are no counters.
        laser = InterlinkModel()
        add_ground_station!(laser, station())
        register_link!(laser, [radio_sat(1; terminal=nothing, n_terminal=1), radio_sat(2; terminal=nothing, n_terminal=1)],
            (1, 1), (2, 1), LaserMomentumLink())
        @test !IM.has_station_counters(laser) && IM.has_station_counters(model)
        plain = InterlinkModel()
        pair = [radio_sat(1), radio_sat(2)]
        register_link!(plain, pair, (1, 1), (2, 1), LINK)
        @test !hasproperty(SE.build_initial_conditions(run_args(pair, plain)), :gs)
    end

    @testset "stations with data stores" begin
        # The flows at one state: a station neither limits a link nor is limited.
        planet = turned_planet(0.2)
        st = station(lat_deg=20.0, lon_deg=40.0)
        fixed(rate) = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=rate), range_m=3e6)
        stored(capacity, held) = (sc = radio_sat(1; terminal=nothing, n_terminal=1);
            add_data!(sc, sc.root, DataStore(capacity_bits=capacity, initial_bits=held)); sc)
        spacecraft = [stored(1e9, 0.0), stored(1e9, 5e8), stored(2e9, 1e9)]
        add_data!(spacecraft[1], spacecraft[1].root, DataSource(rate_bps=1e3))
        SM.data_store(spacecraft[1]).empty = true
        SM.data_store(spacecraft[3]).full = true
        model = InterlinkModel()
        add_ground_station!(model, station(lat_deg=20.0, lon_deg=40.0, antennas=[radio(), radio(gain_dbi=1.0), radio(gain_dbi=2.0)]))
        register_link!(model, spacecraft, (1, 1), (-1, 1), fixed(4e3))   # out of an empty store: its production
        register_link!(model, spacecraft, (2, 1), (-1, 2), fixed(4e3))   # out of a store with data: the link's rate
        register_link!(model, spacecraft, (-1, 3), (3, 1), fixed(4e3))   # into a full store sending nothing: nothing
        model.graph.active .= true
        state = (sc=[(pos=seen_from(st, planet, 1.0, 0.0, 900e3), vel=SVector(0.0, 7.5e3, 0.0), mass=100.0,
                      data_stored_bits=SM.data_store(sc).initial_bits) for sc in spacecraft],)
        flows = SM.InterlinkTerminals._data_flows(model, spacecraft, state; context=LinkContext(time_s=0.0, planet=planet))
        @test flows.sent ≈ [1e3, 4e3, 0.0] && flows.received == [0.0, 0.0, 0.0]
        @test flows.station_received ≈ [5e3] && flows.station_sent == [0.0]
        # With room in the third store, the station's uplink carries its rate.
        SM.data_store(spacecraft[3]).full = false
        flows = SM.InterlinkTerminals._data_flows(model, spacecraft, state; context=LinkContext(time_s=0.0, planet=planet))
        @test flows.received ≈ [0.0, 0.0, 4e3] && flows.station_sent ≈ [4e3]
        # A station holds and has room without limit, for the backlog policy.
        ctx = SchedulingContext(planet.μ, LinkContext(planet=planet), spacecraft, state, nothing)
        @test stored_bits(ctx, -1) == Inf && store_room_bits(ctx, -1) == Inf
        @test stored_bits(SchedulingContext(planet.μ, LinkContext(planet=planet), nothing, state, nothing), -1) == Inf
        backlog = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=100.0))
        update_availability!(model, spacecraft, state; context=LinkContext(planet=planet))
        score_candidates!(model, backlog, state, planet.μ; spacecraft, context=LinkContext(planet=planet))
        # An empty producer scores its production; stored data reaches the rate cap; the station is unlimited.
        @test backlog.matching.score ≈ [1e3, 4e3, 4e3]

        # A run: a full store drains to the station under the backlog policy, which never limits it, and the
        # link drops once the store is empty.
        sc = radio_sat(1; pos=pass_fixture(lead_deg=10.0)[1][1].initial_condition.pos,
            vel=pass_fixture(lead_deg=10.0)[1][1].initial_condition.vel)
        add_data!(sc, sc.root, DataStore(capacity_bits=2e6, initial_bits=2e6))
        drain = InterlinkModel()
        add_ground_station!(drain, station())
        register_link!(drain, [sc], (1, 1), (-1, 1), LINK)
        result = run([sc], drain; duration=450.0, dt_max=30.0, policy=DataBacklogPolicy(horizon_s=60.0))
        final = result.solution.u[end]
        history = result.solution.prob.p.args.interlink_model.history
        @test final.sc[1].data_stored_bits < 1e-2
        @test final.gs[1].received_bits ≈ 2e6 atol=1e-2
        @test final.gs[1].received_bits ≈ final.sc[1].data_sent_bits rtol=1e-14
        @test history[1].active == [((1, 1), (-1, 1))] && isempty(history[end].active)
        # While the store held data, the link carried its full rate: the station did not limit it.
        before = findlast(sample -> result.solution(sample.time).sc[1].data_stored_bits > 1e5, history)
        t_before = history[before].time
        @test 60.0 < t_before < 440.0
        @test result.solution(t_before).sc[1].data_sent_bits ≈
            carried(result.solution, ((1, 1), (-1, 1)); until=t_before) rtol=1e-7
    end

    @testset "backlog counts a full receiver's current outflow" begin
        function relay_fixture(; draining=true)
            reference = pass_fixture(lead_deg=10.0)[1][1]
            sc = radio_sat(1; pos=reference.initial_condition.pos, vel=reference.initial_condition.vel)
            add_actuator!(sc, sc.root, radio())
            add_data!(sc, sc.root, DataStore(capacity_bits=1e7, initial_bits=1e7))
            model = InterlinkModel()
            add_ground_station!(model, station(antennas=[radio(), radio()]))
            link = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=1e4), range_m=3e6)
            register_link!(model, [sc], (-1, 1), (1, 1), link)
            draining && register_link!(model, [sc], (1, 2), (-1, 2), link)
            return [sc], model
        end
        uplink = ((-1, 1), (1, 1))
        for draining in (false, true)
            result = run(relay_fixture(; draining)...; duration=450.0, dt_max=5.0,
                policy=DataBacklogPolicy(horizon_s=60.0))
            table = result.table
            final = result.solution.u[end]
            chosen = [uplink in sample.active for sample in result.configuration.interlink_model.history]
            @test all(-1e-2 .<= table.sc1_data_stored_bits .<= 1e7 + 1e-2)
            @test all(iszero, table.sc1_data_lost_bits)
            @test final.sc[1].data_stored_bits ≈ 1e7 + final.gs[1].sent_bits - final.gs[1].received_bits atol=1e-2
            @test final.sc[1].data_received_bits ≈ final.gs[1].sent_bits atol=1e-2
            @test final.sc[1].data_sent_bits ≈ final.gs[1].received_bits atol=1e-2
            # Initialization cannot borrow outflow from its own first selection:
            # the saved initial row and final initial matching agree.
            @test !chosen[1]
            @test iszero(table.sc1_data_rate_rx_bps[1])
            if draining
                @test table.sc1_data_rate_tx_bps[1] == 1e4
                @test any(chosen)
                first_on = findfirst(chosen)
                @test all(chosen[first_on:end]) # no alternating selection after acquisition
                @test count(pair -> pair[1] != pair[2], zip(chosen[1:end-1], chosen[2:end])) <= 1
                @test final.gs[1].sent_bits > 4e6
            else
                @test !any(chosen)
                @test iszero(final.gs[1].sent_bits)
                @test final.sc[1].data_stored_bits ≈ 1e7 atol=1e-2
            end
        end
        # Score at an exactly full state: outflow is from the old matching,
        # with current availability and downstream flow caps, not the new choice.
        spacecraft, model = relay_fixture()
        args = run_args(spacecraft, model)
        state = SE.build_initial_conditions(args)
        SM.data_store(only(spacecraft)).full = true
        context = LinkContext(time_s=0.0, planet=PLANET)
        model.graph.available .= true
        model.graph.active .= [false, true]
        for policy in (DataBacklogPolicy(horizon_s=60.0),
                       WeightedPolicy(PolicyTerm(DataBacklogPolicy(horizon_s=60.0); scale=1.0)),
                       WeightedPolicy(PolicyTerm(WeightedPolicy(PolicyTerm(DataBacklogPolicy(horizon_s=60.0); scale=1.0)); scale=1.0)))
            scheduler = SchedulingPolicyModel(policy)
            score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context)
            @test scheduler.matching.score[1] == 1e4
            @test model.graph.active == [false, true] # scoring does not choose links
        end
        # A snapshot from before the first initialization pass stays authoritative.
        scheduler = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=60.0))
        score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context, previous_active=falses(2))
        @test scheduler.matching.score[1] == 0.0
        score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context, is_active=[false])
        @test scheduler.matching.score[1] == 0.0
        # Existing caller constructors remain usable without outflow credit.
        SP = SM.SchedulingPolicies
        ctx5 = SchedulingContext(PLANET.μ, context, spacecraft, state, nothing)
        ctx6 = SchedulingContext(PLANET.μ, context, spacecraft, state, nothing, ctx5._terminals)
        @test ctx5._outflows === nothing && ctx6._outflows === nothing
        model.graph.available[2] = false
        scheduler = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=60.0))
        score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context)
        @test scheduler.matching.score[1] == 0.0
        # Nominal rate is not room when the downstream store cannot drain.
        a = radio_sat(1)
        b = radio_sat(2)
        for sc in (a, b)
            add_actuator!(sc, sc.root, radio())
            add_data!(sc, sc.root, DataStore(capacity_bits=1e7, initial_bits=1e7))
        end
        chain = InterlinkModel()
        add_ground_station!(chain, station(antennas=[radio(), radio()]))
        fixed(rate) = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=rate), range_m=1e8)
        register_link!(chain, [a, b], (-1, 1), (1, 1), fixed(1e4))
        register_link!(chain, [a, b], (1, 2), (2, 1), fixed(1e4))
        register_link!(chain, [a, b], (2, 2), (-1, 2), fixed(2500.0))
        state2 = SE.build_initial_conditions(run_args([a, b], chain))
        foreach(sc -> (SM.data_store(sc).full = true), (a, b))
        chain.graph.available .= true
        chain.graph.active .= [false, true, false]
        scheduler2 = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=60.0))
        score_candidates!(chain, scheduler2, state2, PLANET.μ; spacecraft=[a, b], context)
        @test scheduler2.matching.score[1] == 0.0
        chain.graph.active[3] = true
        score_candidates!(chain, scheduler2, state2, PLANET.μ; spacecraft=[a, b], context)
        @test scheduler2.matching.score[1] == 2500.0
    end

    @testset "the planet's frame inside a step" begin
        # The rising half of the pass, before culmination: the station moves away from the spacecraft.
        spacecraft, model = pass_fixture()
        solution = run(spacecraft, model; duration=430.0, dt_max=30.0).solution
        exact, held = carried(solution, ((1, 1), (-1, 1))), carried(solution, ((1, 1), (-1, 1)); held=true)
        @test solution.u[end].sc[1].data_sent_bits ≈ exact rtol=1e-7
        # A frame held at each step's start, the negative control, misses by far more than that.
        @test abs(held - exact) > 1e-4 * exact
        # The engine reads the frame at the evaluation's time.
        p = solution.prob.p
        @test SM.SimulationCallbacks._station_planet_frame(p, 123.4) ≈ frame(123.4) rtol=1e-9
        @test SM.SimulationCallbacks._station_planet_frame(p, NaN) === nothing
    end

    @testset "an uplink" begin
        spacecraft, model = pass_fixture(uplink=true, battery=true)
        solution = run(spacecraft, model; duration=430.0, dt_max=30.0).solution
        final = solution.u[end].sc[1]
        # The spacecraft receives what the station sends, at the closed-form rate.
        @test final.data_received_bits > 0 && final.data_sent_bits == 0
        @test final.data_received_bits ≈ carried(solution, ((-1, 1), (1, 1))) rtol=1e-7
        # A receiving radio draws nothing, and a data link delivers no power.
        @test final.battery_wh == 5.0
    end

    @testset "setup checks" begin
        spacecraft, model = pass_fixture()
        empty!(model.ground_stations)
        err = refusal(() -> run_simulation(run_args(spacecraft, model)))
        @test err isa ArgumentError && occursin("references ground station 1", err.msg)
    end

    @testset "the viewer block" begin
        spacecraft = [radio_sat(1), radio_sat(2)]
        model = InterlinkModel()
        add_ground_station!(model, station())
        between = register_link!(model, spacecraft, (1, 1), (2, 1), LINK)
        down = register_link!(model, spacecraft, (2, 1), (-1, 1), LINK)
        push!(model.history, InterLinkHistory(0.0, [between, down], [between, down]))
        push!(model.history, InterLinkHistory(10.0, [between], [between]))
        table = link_schedule(model)
        @test sort(unique(table.b_spacecraft)) == [-1, 2]
        SV = SM.SceneVisualization
        scene = SV.build_visualization_scene(run_args(spacecraft, model))
        @test only(scene.ground_stations).position_m ≈ SITE
        @test only(scene.ground_stations).name == "north"
        values = ntuple(i -> getfield(scene, i), 10)
        legacy = SV.VisualizationScene(values..., scene.link_schedule, scene.ground_stations)
        @test legacy.ground_stations == scene.ground_stations && legacy.results_layout === :wide
        layout_only = SV.VisualizationScene(values..., scene.link_schedule, :long)
        @test isempty(layout_only.ground_stations) && layout_only.results_layout === :long
        @test table.b_name[findfirst(==(-1), table.b_spacecraft)] == "north"
        @test_throws ArgumentError SV.interlink_payload(table, 2)
        block = SV.interlink_payload(table, 2; ground_stations=scene.ground_stations)
        @test length(block["links"]) == 2 && block["links"][2]["b"] == -1
        @test block["active"]["count"] == 2 && block["contact"]["count"] == 2
        @test block["stations"] == [Dict("node" => -1, "name" => "north", "position_m" => collect(SITE))]
        mktempdir() do directory
            sidecar = joinpath(directory, "scene.json")
            SV.write_visualization_scene(sidecar, scene)
            @test SV.read_visualization_scene(sidecar) == scene
        end
    end

    @testset "saved station pass exports names, coordinates and windows" begin
        for uplink in (false, true), (layout, checkpoints) in ((:wide, false), (:long, false), (:long, true))
            mktempdir() do directory
                spacecraft, model = pass_fixture(; lead_deg=0.0, uplink)
                settings = SimulationSettings(results=true, results_directory=directory, generate_plots=false,
                    save_csv=true, save_visualization_scene=true, results_layout=layout,
                    checkpoint_enabled=checkpoints, checkpoint_interval_s=10.0,
                    checkpoint_directory=directory)
                result = run_simulation(run_args(spacecraft, model; duration=20.0, settings); return_results=true)
                prefix = joinpath(directory, "simulation_results")
                SV = SM.SceneVisualization
                scene = SV.read_visualization_scene(prefix * "_scene.json")
                @test scene.results_layout === layout
                @test SV._scene_from_dict(SV.scene_dict(scene)) == scene
                payload = JSON.parsefile(SV.export_viewer_payload(prefix; textures=false))
                table = link_schedule(result)
                disk = DataFrame(Arrow.Table(prefix * "_link_schedule.feather"))
                @test disk.a_name == table.a_name && disk.b_name == table.b_name
                @test payload["interlinks"] == JSON.parse(JSON.json(SV.interlink_payload(table, 1; ground_stations=scene.ground_stations)))
                link = only(payload["interlinks"]["links"])
                @test (link["a"], link["b"]) == (uplink ? (-1, 1) : (1, -1))
                @test only(payload["interlinks"]["stations"])["position_m"] ≈ SITE
                # The sidecar's rotation maps the saved fixed site to the same
                # inertial position used by the physical station at each sample.
                for (time, q) in zip(scene.planet.rotation_t_s, scene.planet.rotation_q_pi)
                    @test transpose(SM.rot(q)) * SITE ≈ transpose(frame(time)) * SITE atol=1e-6
                end
                if layout === :long
                    # Station geometry and a valid schedule must not conceal a
                    # mismatched long trajectory, including finalized streams.
                    trajectory = joinpath(directory, scene.results_feather)
                    original = read(trajectory)
                    rows = DataFrame(Arrow.Table(trajectory); copycols=true)
                    rows.pos_1 .+= 1000.0
                    Arrow.write(trajectory, rows)
                    @test_throws ArgumentError SV.export_viewer_payload(prefix; textures=false, out=prefix * "_stale.json")
                    @test !isfile(prefix * "_stale.json")
                    write(trajectory, original)
                end
            end
        end
    end

    @testset "resume initializes the station schedule at checkpoint time" begin
        # Choose a rising pass where the epoch-zero station frame sees the link
        # at the checkpoint, but the physically rotated station is still below
        # its elevation mask. The mission horizon is unchanged across the crash.
        checkpoint_time = 200.0
        mask = deg2rad(10.0)
        limit = acos(PLANET.Rp_e / RADIUS * cos(mask)) - mask
        mean_motion = sqrt(PLANET.μ / RADIUS^3)
        spin = norm(PLANET.ω)
        lead = rad2deg(limit + mean_motion * checkpoint_time - spin * checkpoint_time / 2)
        function args_at(directory; resume=false)
            settings = SimulationSettings(results=false, generate_plots=false, verbose=false,
                results_directory=directory, checkpoint_directory=directory,
                checkpoint_enabled=true, checkpoint_interval_s=100.0, resume_from_checkpoint=resume)
            return run_args(pass_fixture(lead_deg=lead)...; duration=400.0, dt_max=5.0, settings)
        end
        function observer(selections)
            DiscreteCallback((u, t, integrator) -> false, integrator -> nothing;
                initialize=(cb, u, t, integrator) -> begin
                    t == checkpoint_time && push!(selections,
                        copy(last(integrator.p.args.interlink_model.history).active))
                    nothing
                end, save_positions=(false, false))
        end
        uninterrupted = Any[]
        straight = run_simulation(args_at(mktempdir()); return_results=true,
            extra_callbacks=(observer(uninterrupted),))
        dir = mktempdir()
        crash = DiscreteCallback((u, t, integrator) -> t > 250.0,
            integrator -> error("station resume interruption"); save_positions=(false, false))
        @test_throws Exception run_simulation(args_at(dir); return_results=true, extra_callbacks=(crash,))
        selections = Any[]
        resumed = run_simulation(args_at(dir; resume=true); return_results=true,
            extra_callbacks=(observer(selections),))
        @test length(uninterrupted) == length(selections) == 1
        @test isempty(only(uninterrupted))
        @test selections == uninterrupted
        rows = findall(==(checkpoint_time), resumed.table.time)
        @test !isempty(rows)
        @test all(iszero, resumed.table.sc1_data_rate_tx_bps[rows])
        @test resumed.table.sc1_data_sent_bits[end] ≈ straight.table.sc1_data_sent_bits[end] rtol=1e-8
    end

    @testset "checkpoints" begin
        configuration = SE._interlink_checkpoint_configuration(run_args(pass_fixture()...))
        entry = only(configuration.ground_stations)
        @test keys(entry) == (:site, :antennas) && entry.site.min_elevation_deg.value == 10.0
        @test only(entry.antennas) == SE._interlink_terminal_identity(radio())
        # A changed site, mask or antenna is another configuration; a new name is not.
        for changed in (station(lat_deg=0.001), station(lon_deg=-0.001), station(alt_m=1.0), station(min_elevation_deg=10.5),
                        station(antennas=[radio(gain_dbi=1.0)]), station(antennas=[radio(), radio()]))
            @test !isequal(SE._interlink_checkpoint_configuration(run_args(pass_fixture(st=changed)...)), configuration)
        end
        @test isequal(SE._interlink_checkpoint_configuration(run_args(pass_fixture(st=station(name="renamed"))...)),
            configuration)
        # A checkpoint whose state lacks the station's counters is refused.
        args = run_args(pass_fixture()...)
        saved = (runtime_state=SE._interlink_checkpoint_runtime_state(args, [true]),
            u=(sc=[(pos=zeros(3), vel=zeros(3), mass=100.0)],))
        err = refusal(() -> SE._restore_interlink_checkpoint_state!(args, saved, [true]))
        @test err isa ArgumentError && occursin("one data counter pair per ground station (1 in this run)", err.msg)
        # The state's renaming for a checkpoint keeps the counters.
        u = SE.build_initial_conditions(args)
        u.gs[1].received_bits = 3.0
        @test SE._checkpoint_state(args, u).gs[1].received_bits == 3.0
        # A resumed run continues; one with a changed station is refused. The spacecraft is in view
        # from the start to the end, so the selection does not depend on the steps.
        overhead() = pass_fixture(lead_deg=10.0)
        straight = run(overhead()...; duration=450.0, dt_max=30.0).solution.u[end]
        mktempdir() do directory
            segments(resume) = SimulationSettings(results=resume, results_directory=directory, generate_plots=false,
                checkpoint_enabled=true, checkpoint_interval_s=100.0, checkpoint_directory=directory,
                resume_from_checkpoint=resume)
            run(overhead()...; duration=200.0, dt_max=30.0, settings=segments(false))
            err = refusal(() -> run(pass_fixture(lead_deg=10.0, st=station(min_elevation_deg=12.0))...; duration=450.0,
                dt_max=30.0, settings=segments(true)))
            @test err isa ArgumentError && occursin("ground stations", err.msg)
            resumed = run(overhead()...; duration=450.0, dt_max=30.0, settings=segments(true))
            @test all(sample -> ((1, 1), (-1, 1)) in sample.active, resumed.solution.prob.p.args.interlink_model.history)
            @test resumed.solution.u[end].sc[1].data_sent_bits ≈ straight.sc[1].data_sent_bits rtol=1e-8
            @test resumed.table.sc1_data_sent_bits[end] ≈ straight.sc[1].data_sent_bits rtol=1e-8
            # The station's counters continue too.
            @test resumed.solution.u[end].gs[1].received_bits ≈ straight.gs[1].received_bits rtol=1e-8
            @test resumed.table.gs1_received_bits[end] ≈ straight.gs[1].received_bits rtol=1e-8
        end
    end

    @testset "station and store checkpoints continue together" begin
        # The same checkpoint boundaries make the reference and resumed runs
        # comparable bit for bit. An empty producer stays selected on downlink;
        # the station remains an unlimited sender on uplink until the store fills.
        function stored_pass(uplink)
            spacecraft, model = pass_fixture(lead_deg=10.0; uplink)
            sc = only(spacecraft)
            add_data!(sc, sc.root, DataStore(capacity_bits=2e6,
                initial_bits=0.0))
            uplink || add_data!(sc, sc.root, DataSource(rate_bps=2e3))
            return spacecraft, model
        end
        settings(directory, resume=false) = SimulationSettings(results=false,
            generate_plots=false, checkpoint_enabled=true, checkpoint_interval_s=100.0,
            checkpoint_directory=directory, resume_from_checkpoint=resume)
        for uplink in (false, true)
            mktempdir() do reference_dir
                whole = run(stored_pass(uplink)...; duration=450.0, dt_max=30.0,
                    policy=DataBacklogPolicy(horizon_s=60.0), settings=settings(reference_dir))
                mktempdir() do split_dir
                    # Interrupt the same mission just after the 200 s checkpoint.
                    # Shortening its configured horizon to 200 s would test
                    # a changed mission as well as checkpoint restart.
                    interrupted_args = run_args(stored_pass(uplink)...; duration=450.0,
                        dt_max=30.0, policy=DataBacklogPolicy(horizon_s=60.0), settings=settings(split_dir))
                    stop = SM.SimulationCallbacks.PresetTimeCallback([201.0],
                        i -> throw(StationStoreInterrupted()); save_positions=(false, false))
                    @test_throws StationStoreInterrupted run_simulation(interrupted_args;
                        extra_callbacks=(stop,), return_results=true, return_solution=true)
                    saved = SE._load_checkpoint(interrupted_args)
                    @test saved.t == 200.0
                    first_counter = uplink ? saved.u.gs[1].sent_bits : saved.u.gs[1].received_bits
                    @test first_counter > 0.0
                    resumed = run(stored_pass(uplink)...; duration=450.0, dt_max=30.0,
                        policy=DataBacklogPolicy(horizon_s=60.0), settings=settings(split_dir, true))
                    @test resumed.solution.t == whole.solution.t
                    @test all(isequal(collect(a), collect(b)) for
                        (a, b) in zip(resumed.solution.u, whole.solution.u))
                    final = resumed.solution.u[end]
                    @test final.gs[1].sent_bits ≈ final.sc[1].data_received_bits rtol=1e-14
                    @test final.gs[1].received_bits ≈ final.sc[1].data_sent_bits rtol=1e-14
                    initial = 0.0
                    @test initial + final.sc[1].data_produced_bits + final.gs[1].sent_bits ≈
                        final.sc[1].data_stored_bits + final.sc[1].data_lost_bits +
                        final.gs[1].received_bits atol=1e-2
                    @test resumed.table.gs1_received_bits[end] == final.gs[1].received_bits
                    @test resumed.table.gs1_sent_bits[end] == final.gs[1].sent_bits
                    if !uplink
                        key = [((1, 1), (-1, 1))]
                        whole_history = whole.configuration.interlink_model.history
                        resumed_history = resumed.configuration.interlink_model.history
                        @test first(whole_history).active == key
                        @test all(sample -> sample.active == key, whole_history)
                        @test first(resumed_history).active == key
                        @test all(sample -> sample.active == key, resumed_history)
                        @test saved.u.sc[1].data_stored_bits ≈ 0.0 atol=1e-6
                        @test final.sc[1].data_stored_bits ≈ 0.0 atol=1e-6
                        @test final.sc[1].data_produced_bits ≈ 2e3 * 450.0 rtol=1e-9
                        @test final.sc[1].data_sent_bits ≈ final.sc[1].data_produced_bits rtol=1e-9
                        @test final.gs[1].received_bits > first_counter
                    end
                end
            end
        end
    end

    @testset "negative control" begin
        # Without stations: none on the model, no planet frame read, no checkpoint entry.
        model = InterlinkModel()
        @test isempty(model.ground_stations) && !IM.has_ground_stations(model) && !IM.has_ground_stations(nothing)
        @test SM.SimulationCallbacks._station_planet_frame((args=(interlink_model=model,),), 1.0) === nothing
        @test SM.SimulationCallbacks._station_planet_frame((args=(interlink_model=nothing,),), 1.0) === nothing
        spacecraft = [radio_sat(1), radio_sat(2)]
        register_link!(model, spacecraft, (1, 1), (2, 1), LINK)
        @test !haskey(SE._interlink_checkpoint_configuration(run_args(spacecraft, model)), :ground_stations)
        # A link between spacecraft keeps its occultation test.
        a = IM.link_endpoint((pos=[RADIUS, 0.0, 0.0],), nothing, 1)
        b = IM.link_endpoint((pos=[-RADIUS, 1.0, 0.0],), nothing, 2)
        @test !link_available(DataLink(transmit_power_w=1.0, budget=BUDGET, range_m=2e7), a, b, LinkContext(planet=PLANET))
        # A station does not see a spacecraft on the planet's far side.
        planet = turned_planet(0.0)
        st = station()
        far = InterlinkModel()
        add_ground_station!(far, st)
        register_link!(far, spacecraft[1:1], (1, 1), (-1, 1), DataLink(transmit_power_w=1.0, budget=BUDGET, range_m=2e7))
        @test !only(update_availability!(far, spacecraft[1:1], state_at(SVector(-RADIUS, 0.0, 0.0));
            context=LinkContext(planet=planet)))
    end
end

end # module InterlinkGroundStationTests
