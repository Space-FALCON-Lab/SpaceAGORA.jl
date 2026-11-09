include("run.jl")

module EquatorialGroundStations

using ..KuangLaterTasks
using SpaceAGORA
using SpaceAGORA.SimulationModel
using CSV, DataFrames, JSON, Base64, LinearAlgebra, Test, Dates, SHA

const K = KuangLaterTasks
const SM = SpaceAGORA.SimulationModel
const OUT = joinpath(@__DIR__, "results", "ground_stations")
const LONGITUDES = (0.0, 120.0, -120.0)
const NAMES = ("Equator_0", "Equator_120E", "Equator_120W")
const CAPACITY = 1e9
const INITIAL = 20e6
const PRODUCTION = 1e5
const RATE = 1e6
const HORIZON = 60.0
const RANGE = 5000e3
const MASK = 10.0
const BIT_TOL = 1.0
const KEY = K.LinkKey

function frame(args, t)
    ephemerides = args.environment_model.ephemerides_model
    epoch = SM.EphemeridesModels.ephemerides_epoch_seconds(ephemerides, args.initial_time)
    return SM.EphemeridesModels.ephemerides_frame_lpi(
        ephemerides, args.environment_model.planet, epoch + t)
end

site(s) = K.PLANET.Rp_e * [cosd(LONGITUDES[s]), sind(LONGITUDES[s]), 0.0]
key(i, s) = ((i, 1), (-s, 1))

function configuration()
    original = K.constellation_configuration(OUT)
    spacecraft = original.dynamics_model.spacecraft
    model = InterlinkModel(record_available=true)
    for (s, lon) in enumerate(LONGITUDES)
        antenna = RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=1.0,
            gain_dbi=0.0, efficiency=0.5, max_rate_bps=RATE)
        add_ground_station!(model, GroundStation(name=NAMES[s], lat_deg=0.0,
            lon_deg=lon, alt_m=0.0, min_elevation_deg=MASK, antennas=[antenna]))
    end
    for (i, sc) in enumerate(spacecraft)
        add_data!(sc, sc.root, DataStore(capacity_bits=CAPACITY, initial_bits=INITIAL))
        add_data!(sc, sc.root, DataSource(rate_bps=PRODUCTION))
        for s in 1:3
            register_link!(model, spacecraft, key(i, s)...,
                DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=RATE), range_m=RANGE))
        end
    end
    return K.configuration(spacecraft, model, DataBacklogPolicy(horizon_s=HORIZON);
        duration=original.mission_configuration.mission_time, dtmax=5.0,
        output=OUT, data_rate=10.0)
end

function candidates(args, state, time_s)
    rotation = frame(args, time_s)
    return DataFrame(map([(i, s) for i in 1:40 for s in 1:3]) do (i, s)
        fixed = rotation * state.sc[i].pos
        location = site(s)
        displacement = fixed - location
        distance = norm(displacement)
        elevation = asind(clamp(dot(displacement, location / norm(location)) / distance, -1.0, 1.0))
        available = distance <= RANGE && elevation >= MASK
        held = state.sc[i].data_stored_bits
        held <= CAPACITY * 1e-9 && (held = 0.0)
        score = available ? min(RATE, held / HORIZON + PRODUCTION) : 0.0
        (; time_s, spacecraft=i, station=s, station_name=NAMES[s], allowed=true,
            range_m=distance, elevation_deg=elevation, available, score_bps=score)
    end)
end

function geometry_probes(args)
    state = K.SE.build_initial_conditions(args)
    rotation = frame(args, 0.0)
    model = deepcopy(args.interlink_model)
    rows = NamedTuple[]
    for s in 1:3, elevation in (-5.0, 9.0, 11.0, 90.0)
        radius = norm(state.sc[1].pos)
        central_angle = acosd(K.PLANET.Rp_e / radius * cosd(elevation)) - elevation
        longitude = LONGITUDES[s] + central_angle
        state.sc[1].pos .= transpose(rotation) *
            (radius * [cosd(longitude), sind(longitude), 0.0])
        table = candidates(args, state, 0.0)
        row = only(eachrow(table[(table.spacecraft .== 1) .& (table.station .== s), :]))
        context = LinkContext(time_s=0.0, planet=deepcopy(K.PLANET))
        context.planet.L_PI .= rotation
        update_availability!(model, args.dynamics_model.spacecraft, state; context)
        actual = model.graph.available[model.graph.edge_index_lookup[key(1, s)]]
        @test isapprox(row.elevation_deg, elevation; atol=1e-10, rtol=0)
        @test actual == (elevation >= MASK && row.range_m <= RANGE)
        push!(rows, (; station=s, requested_elevation_deg=elevation,
            measured_elevation_deg=row.elevation_deg, row.range_m,
            expected_available=row.available, actual_available=actual))
    end
    CSV.write(joinpath(OUT, "geometry_probes.csv"), DataFrame(rows))
end

function verify_viewer_windows(block, schedule, history)
    decode(data) = collect(reinterpret(Float64, base64decode(data)))
    @test length(block["stations"]) == 3
    for s in 1:3
        @test block["stations"][s]["node"] == -s
        @test block["stations"][s]["position_m"] ≈ site(s) atol=1e-8
    end
    for kind in ("active", "contact")
        window = block[kind]
        starts, ends = decode(window["start_s"]), decode(window["end_s"])
        @test window["count"] == count(==(kind), schedule.window)
        for sample in history
            represented = Set{KEY}()
            for j in eachindex(starts)
                starts[j] <= sample.time < ends[j] || continue
                link = block["links"][window["index"][j] + 1]
                push!(represented, ((link["a"], link["a_terminal"]), (link["b"], link["b_terminal"])))
            end
            # A final sample has no following interval; the fixture viewer retains it explicitly.
            if sample.time < last(history).time
                @test represented == Set(kind == "active" ? sample.active : sample.available)
            end
        end
    end
end

function run()
    mkpath(OUT)
    status = Dict{String, Any}("status" => "running", "started_utc" => string(now(UTC)))
    K.write_toml(joinpath(OUT, "status.toml"), status)
    passed = false
    try
        @testset failfast=true "Three equatorial ground stations" begin
            args = configuration()
            initial = K.SE.build_initial_conditions(args)
            CSV.write(joinpath(OUT, "stations.csv"), DataFrame(station=1:3, name=collect(NAMES),
                latitude_deg=zeros(3), longitude_deg=collect(LONGITUDES), altitude_m=zeros(3),
                minimum_elevation_deg=fill(MASK, 3), antennas=ones(Int, 3)))
            CSV.write(joinpath(OUT, "candidates_before_run.csv"), candidates(args, initial, 0.0))
            for a in 1:3, b in (a + 1):3
                @test isapprox(acosd(dot(site(a), site(b)) / K.PLANET.Rp_e^2), 120; atol=1e-10)
            end
            geometry_probes(args)
            solution = K.run_checked(args)
            run_args = solution.prob.p.args
            model = run_args.interlink_model
            history = model.history
            states = K.history_states(solution)
            tables, accounting, network, audits, frames = DataFrame[], NamedTuple[], NamedTuple[], NamedTuple[], Any[]
            expected_ground = zeros(3)
            hand_uncertainty = zeros(3)
            max_balance_error = 0.0
            previous = Set{KEY}()
            for (k, (sample, state)) in enumerate(zip(history, states))
                table = candidates(run_args, state, sample.time)
                available = Set(key(row.spacecraft, row.station) for row in eachrow(table) if row.available)
                selected = Set(sample.active)
                @test Set(sample.available) == available
                @test issubset(selected, available)
                @test allunique(collect(Iterators.flatten(sample.active)))
                @test length(selected) <= 3
                # Visibility cones do not overlap: each component is one station's star.
                @test all(i -> count(row -> row.spacecraft == i && row.available, eachrow(table)) <= 1, 1:40)
                optimum = 0.0
                for s in 1:3
                    feasible = table[(table.station .== s) .& table.available, :]
                    best = maximum(feasible.score_bps; init=0.0)
                    optimum += best
                    chosen = [link for link in sample.active if link[2][1] == -s]
                    @test length(chosen) == (best > 0 ? 1 : 0)
                end
                @test isapprox(sum(sample.scores), optimum; atol=1e-5, rtol=1e-12)
                for (link, score) in zip(sample.active, sample.scores)
                    row = only(eachrow(table[(table.spacecraft .== link[1][1]) .&
                        (table.station .== -link[2][1]), :]))
                    @test isapprox(score, row.score_bps; atol=1e-5, rtol=1e-12)
                end
                table.selected = [key(row.spacecraft, row.station) in selected for row in eachrow(table)]
                push!(tables, table)
                if k > 1
                    dt = sample.time - history[k - 1].time
                    for link in history[k - 1].active
                        i, s = link[1][1], -link[2][1]
                        held = max(0.0, states[k - 1].sc[i].data_stored_bits)
                        upper = min(RATE * dt, held + PRODUCTION * dt)
                        lower = held <= CAPACITY * 1e-9 ? min(RATE, PRODUCTION) * dt : upper
                        expected_ground[s] += upper
                        # Tiny positive residuals may or may not have crossed the empty event.
                        hand_uncertainty[s] += upper - lower
                    end
                end
                for i in 1:40
                    sc = state.sc[i]
                    residual = INITIAL + sc.data_produced_bits + sc.data_received_bits -
                        sc.data_sent_bits - sc.data_lost_bits - sc.data_stored_bits
                    max_balance_error = max(max_balance_error, abs(residual))
                    @test abs(residual) <= BIT_TOL
                    @test -BIT_TOL <= sc.data_stored_bits <= CAPACITY + BIT_TOL
                    @test abs(sc.data_lost_bits) <= BIT_TOL
                    @test abs(sc.data_received_bits) <= BIT_TOL
                    @test isapprox(sc.data_produced_bits, PRODUCTION * sample.time; atol=BIT_TOL, rtol=0)
                    push!(accounting, (; time_s=sample.time, spacecraft=i,
                        stored_bits=sc.data_stored_bits, produced_bits=sc.data_produced_bits,
                        downlinked_bits=sc.data_sent_bits, lost_bits=sc.data_lost_bits,
                        balance_residual_bits=residual))
                end
                delivered = sum(station.received_bits for station in state.gs)
                @test isapprox(delivered, sum(sc.data_sent_bits for sc in state.sc); atol=BIT_TOL, rtol=0)
                residual = 40INITIAL + sum(sc.data_produced_bits - sc.data_stored_bits -
                    sc.data_lost_bits for sc in state.sc) - delivered
                @test abs(residual) <= BIT_TOL
                for s in 1:3
                    @test expected_ground[s] - hand_uncertainty[s] - BIT_TOL <=
                        state.gs[s].received_bits <= expected_ground[s] + BIT_TOL
                    @test state.gs[s].sent_bits == 0.0
                end
                push!(network, (; time_s=sample.time, delivered_bits=delivered,
                    stored_bits=sum(sc.data_stored_bits for sc in state.sc),
                    produced_bits=40PRODUCTION * sample.time,
                    lost_bits=sum(sc.data_lost_bits for sc in state.sc), residual_bits=residual,
                    gs1_received_bits=state.gs[1].received_bits,
                    gs2_received_bits=state.gs[2].received_bits,
                    gs3_received_bits=state.gs[3].received_bits))
                push!(audits, (; time_s=sample.time, available_count=length(available),
                    active_count=length(selected), selected_score_bps=sum(sample.scores),
                    optimal_score_bps=optimum, matching=K.matching_label(selected),
                    changed=sample.time > 0 && selected != previous))
                rotation = frame(run_args, sample.time)
                push!(frames, Dict("time_s" => sample.time,
                    "positions_m" => [collect(rotation * sc.pos) for sc in state.sc],
                    "active" => [K.key_label(link) for link in sample.active],
                    "available" => [K.key_label(link) for link in sample.available],
                    "delivered_bits" => [station.received_bits for station in state.gs]))
                previous = selected
            end
            CSV.write(joinpath(OUT, "candidate_scores.csv"), vcat(tables...))
            CSV.write(joinpath(OUT, "node_accounting.csv"), DataFrame(accounting))
            CSV.write(joinpath(OUT, "network_accounting.csv"), DataFrame(network))
            CSV.write(joinpath(OUT, "accepted_steps.csv"), DataFrame(audits))
            schedule = K.save_schedule(solution, OUT)
            summary = DataFrame(station=1:3, name=collect(NAMES),
                received_bits=[station.received_bits for station in solution.u[end].gs],
                hand_lower_bits=expected_ground - hand_uncertainty,
                hand_upper_bits=expected_ground, empty_threshold_uncertainty_bits=hand_uncertainty)
            @test all(>(0), summary.received_bits)
            @test K.schedule_changes(history) > 0
            @test any(row -> !row.available, eachrow(first(tables)))
            CSV.write(joinpath(OUT, "station_delivery.csv"), summary)
            SV = SM.SceneVisualization
            scene = SV.build_visualization_scene(run_args)
            SV.write_visualization_scene(joinpath(OUT, "scene.json"), scene)
            block = SV.interlink_payload(schedule, 40; ground_stations=scene.ground_stations)
            verify_viewer_windows(block, schedule, history)
            open(joinpath(OUT, "interlink_viewer_payload.json"), "w") do io
                JSON.print(io, block)
            end
            links = [Dict("key" => K.key_label(key(i, s)), "spacecraft" => i, "station" => s)
                for i in 1:40 for s in 1:3]
            payload = Dict("radius_m" => K.PLANET.Rp_e,
                "stations" => [Dict("name" => NAMES[s], "position_m" => site(s),
                    "longitude_deg" => LONGITUDES[s]) for s in 1:3],
                "links" => links, "frames" => frames)
            json = replace(JSON.json(payload), "</" => "<\\/")
            template = read(joinpath(@__DIR__, "ground_viewer.html"), String)
            write(joinpath(OUT, "viewer.html"), replace(template, "__GROUND_DATA__" => json))
            K.write_toml(joinpath(OUT, "settings.toml"), Dict(
                "duration_s" => solution.t[end], "spacecraft_count" => 40,
                "station_longitudes_deg" => collect(LONGITUDES), "station_latitude_deg" => 0.0,
                "elevation_mask_deg" => MASK, "range_m" => RANGE, "link_rate_bps" => RATE,
                "production_per_spacecraft_bps" => PRODUCTION, "initial_per_spacecraft_bits" => INITIAL,
                "capacity_per_spacecraft_bits" => CAPACITY, "backlog_horizon_s" => HORIZON,
                "dtmax_s" => 5.0, "reltol" => 1e-10, "abstol" => 1e-8,
                "accounting_atol_bits" => BIT_TOL, "virtual_spacecraft_terminals" => true,
                "ground_antennas_per_station" => 1, "record_available" => true,
                "source_sha256" => bytes2hex(sha256(read(@__FILE__))),
                "viewer_source_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "ground_viewer.html")))),
                "git_revision" => readchomp(`git -C $(K.ROOT) rev-parse HEAD`),
                "julia_version" => string(VERSION)))
            open(joinpath(@__DIR__, "GROUND_REPORT.md"), "w") do io
                println(io, "# Three equatorial ground stations\n\nGenerated: $(now(UTC)) UTC.\n")
                println(io, "This follow-up supplements, rather than replaces, [the original results](REPORT.md). ",
                    "Data now reaches actual ground-station endpoints, not a destination satellite.\n")
                println(io, "## Configuration\n\n",
                    "- Latitude 0 degrees; longitudes 0, 120 E and 120 W; altitude 0 m. ",
                    "Equal 120-degree longitude spacing was checked numerically.\n",
                    "- One receiving antenna per station, 10-degree minimum elevation, 5000 km range. ",
                    "Stations rotate with Earth using the simulation ephemerides frame. The mask ",
                    "enforces horizon clearance; station links do not use the spacecraft-to-spacecraft occultation test.\n",
                    "- Same 40 spacecraft orbital geometry as the original example: 30 at 1000 km and ",
                    "10 at 1100 km, circular equatorial orbits; 227 kg each. ",
                    "Each has one virtual data terminal, a 1 Gbit store initially holding 20 Mbit, ",
                    "and 0.1 Mbit/s continuous production.\n",
                    "- All 120 satellite-to-station downlinks are allowed; no intersatellite links ",
                    "or momentum forces in this separate data experiment. Fixed rate 1 Mbit/s; ",
                    "backlog horizon 60 s; no radio link-budget or battery limitation is asserted.\n",
                    "- Duration $(solution.t[end]) s (quarter target orbit), Tsit5, dtmax 5 s, ",
                    "reltol 1e-10, abstol 1e-8. No movie is generated.\n")
                println(io, "## Measured delivery\n\n| Station | Longitude | Received (Mbit) | Hand lower (Mbit) | Hand upper (Mbit) |\n|---|---:|---:|---:|---:|")
                for row in eachrow(summary)
                    println(io, "| $(row.name) | $(LONGITUDES[row.station]) | $(row.received_bits / 1e6) | $(row.hand_lower_bits / 1e6) | $(row.hand_upper_bits / 1e6) |")
                end
                println(io, "\nTotal Earth delivery: **$(sum(summary.received_bits) / 1e6) Mbit**. ",
                    "$(length(history)) scheduling samples; $(K.schedule_changes(history)) selected-set changes. ",
                    "All three stations receive data. No production loss occurred to the 1-bit tolerance.\n")
                println(io, "## Checks and interpretation\n\n",
                    "- Independently computed range/elevation for all 120 candidates at every accepted step. ",
                    "Separate tests at -5, 9, 11 and 90 degrees verify blocking and the 10-degree mask at every station.\n",
                    "- Actual availability and selected scores agree with those calculations. ",
                    "A spacecraft can see at most one station at once here, so feasible components ",
                    "are stars. The best matching is one maximum positive-score edge per station. ",
                    "This is checked without a Python cyclic matcher. Visible does not mean selected.\n",
                    "- Independent interval upper transfer: min(rate * dt, starting stored bits + production * dt), ",
                    "summed only over actually selected accepted-step intervals. When a store starts ",
                    "within its documented 1e-9-capacity empty threshold (1 bit here), the lower ",
                    "transfer is production * dt: a tiny positive residual may already be flagged empty. ",
                    "These interval uncertainties are accumulated explicitly, not hidden in solver tolerances. ",
                    "Every station counter lies between the resulting lower and upper integrals, ",
                    "with 1 bit numerical allowance. Final uncertainty widths: $(collect(hand_uncertainty)) bits.\n",
                    "- At every saved step: per-spacecraft balance, bounds, zero loss, and ",
                    "sum(spacecraft sent) = sum(ground received). Network balance is ",
                    "initial + produced = remaining spacecraft storage + Earth delivery + loss. ",
                    "Maximum node residual: $max_balance_error bit. Ground delivery is now ",
                    "outside the spacecraft storage total, unlike the earlier destination-satellite fixture.\n",
                    "- Contact boundaries and switches have accepted-step resolution, not exact ",
                    "continuously located visibility crossings. A 5-second step cap is not a claim ",
                    "that this scheduling behavior is timestep-converged. Busy-relay comparison ",
                    "and its refinement remain in the original separate fixture.\n")
                println(io, "## Viewer\n\nOpen [viewer.html](results/ground_stations/viewer.html). ",
                    "This is a self-contained, Earth-fixed north-polar **2D diagnostic viewer**, ",
                    "not the optional SpaceAGORAVisualizer 3D application. Green markers are stations; ",
                    "orange lines are active downlinks; dashed green lines are available but idle. ",
                    "The slider uses every accepted scheduling sample, so it cannot silently skip ",
                    "short recorded windows. Playback advances samples, not uniform physical time. ",
                    "At the final sample it shows the final recorded selection, without extrapolating beyond the run.\n")
                println(io, "The repository-native [viewer link payload](results/ground_stations/interlink_viewer_payload.json) ",
                    "and [scene](results/ground_stations/scene.json) are also saved. Their station coordinates ",
                    "and contact/active intervals were checked against history. ",
                    "The optional 3D renderer, terrain, weather, antenna pointing, atmospheric losses, ",
                    "power shedding and resumed histories are not validated by this fixture.\n")
                println(io, "The separate browser check writes [VIEWER_CHECK.md](VIEWER_CHECK.md) ",
                    "and a SHA-256-bound validation record after checking the actual rendered links, ",
                    "coordinates, delivery counters and controls at every recorded sample. ",
                    "Rerun that check after regenerating this viewer; see [README.md](README.md).\n")
                println(io, "## Evidence and rerun\n\n",
                    "- [Stations](results/ground_stations/stations.csv)\n",
                    "- [Candidates before propagation](results/ground_stations/candidates_before_run.csv)\n",
                    "- [Geometry probes](results/ground_stations/geometry_probes.csv)\n",
                    "- [All candidates and selections](results/ground_stations/candidate_scores.csv)\n",
                    "- [Recorded windows](results/ground_stations/link_schedule.csv)\n",
                    "- [Station delivery](results/ground_stations/station_delivery.csv)\n",
                    "- [Node accounting](results/ground_stations/node_accounting.csv)\n",
                    "- [Network accounting](results/ground_stations/network_accounting.csv)\n",
                    "- [Settings and hashes](results/ground_stations/settings.toml)\n\n",
                    "From the repository root:\n\n```sh\n",
                    "julia --startup-file=no --project=. kuang_later_tasks/ground_stations.jl\n```\n")
            end
        end
        passed = true
    finally
        status["status"] = passed ? "passed" : "failed"
        status["completed_utc"] = string(now(UTC))
        K.write_toml(joinpath(OUT, "status.toml"), status)
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    EquatorialGroundStations.run()
end
