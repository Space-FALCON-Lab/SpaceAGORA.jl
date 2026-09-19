using Test, SpaceAGORA, Arrow, DataFrames, JSON
module MissionDemoHelpersTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "common.jl"))
end
const MDH = MissionDemoHelpersTest

@testset "bounded mission searches" begin
    @test MDH._golden_min(t -> (t - 2)^2, 0.0, 5.0) ≈ 2.0 atol=0.001
    @test MDH._bisect(t -> 2 - t, 0.0, 5.0) ≈ 2.0 atol=0.0001
    @test MDH._bisect(t -> -t, 0.0, 5.0) == 0.0
    @test MDH._closest_approach_et(t -> (t - 2)^2, 0.0; window_s=5.0, step_s=1.0) ≈ 2.0 atol=0.001
    @test MDH._first_apoapsis_et_after(t -> 2 - t, 0.0; window_s=5.0, step_s=1.0) == 2.0
    @test first(MDH._search_grid(0.0, 5.0, 2.0)) == 0.0
    @test last(MDH._search_grid(0.0, 5.0, 2.0)) == 5.0
    @test maximum(diff(MDH._search_grid(0.0, 5.0, 2.0))) <= 2.0
    for bad in (0.0, -1.0, NaN, Inf)
        @test_throws ArgumentError MDH._golden_min(abs, -1.0, 1.0; tol=bad)
        @test_throws ArgumentError MDH._bisect(t -> -t, -1.0, 1.0; tol=bad)
        @test_throws ArgumentError MDH._closest_approach_et(abs, 0.0; window_s=bad)
        @test_throws ArgumentError MDH._search_grid(0.0, 5.0, bad)
    end
    @test_throws ArgumentError MDH._golden_min(t -> NaN, 0.0, 1.0)
    @test_throws ArgumentError MDH._golden_min(identity, 1.0, 0.0)
    @test_throws ArgumentError MDH._bisect(t -> 1.0, 0.0, 1.0)
    @test_throws ArgumentError MDH._bisect(t -> NaN, 0.0, 1.0)
    @test_throws ArgumentError MDH._closest_approach_et(identity, 0.0; window_s=5.0, step_s=1.0)
    @test_throws ArgumentError MDH._first_apoapsis_et_after(t -> 1.0, 0.0; window_s=5.0, step_s=1.0)
    @test_throws ArgumentError MDH._search_grid(0.0, 1e9, 1.0)
    @test_throws ArgumentError MDH._golden_min(t -> (t-1e16)^2, 1e16, 1e16+8; tol=eps())
end

@testset "fresh outputs and spacecraft-index references" begin
    mktempdir() do root
        out = joinpath(root, "results")
        options = MDH.viewer_demo_options("mission", 10.0; argv=["--output-dir",out,"--duration-s","5"])
        @test options.output_dir == out
        @test options.duration_s == 5.0
        mkpath(out); write(joinpath(out,"existing.txt"),"preserve me")
        @test_throws ArgumentError MDH.viewer_demo_options("mission",10.0;argv=["--output-dir",out])
        args=(simulation_settings=(results_directory=out,results=true,generate_filenames=false),)
        @test_throws ArgumentError MDH.run_or_reuse!(args,out)
        @test read(joinpath(out,"existing.txt"),String) == "preserve me"
        @test_throws ArgumentError MDH.run_or_reuse!(args,joinpath(root,"elsewhere"))
        prefix=joinpath(root,"simulation_results")
        df=DataFrame(time=[0.0,1.0],sc1_pos_1=[100.0,200.0],sc1_pos_2=zeros(2),sc1_pos_3=zeros(2),
                     sc2_pos_1=[1.0,2.0],sc2_pos_2=zeros(2),sc2_pos_3=zeros(2))
        Arrow.write(prefix*".feather",df)
        ref=(name="second",target=2,t_s=[0.0,1.0],pos_m=[1.0 2.0;0.0 0.0;0.0 0.0])
        @test MDH.reference_separation(prefix,ref) == [0.0,0.0]
        @test MDH.reference_separation(prefix,merge(ref,(target=1,))) == [99.0,198.0]
        @test_throws ArgumentError MDH.reference_separation(prefix,merge(ref,(target=3,)))
        @test_throws ArgumentError MDH.reference_separation(prefix,merge(ref,(t_s=[5.0,6.0],)))
        @test_throws ArgumentError MDH.reference_separation(prefix,merge(ref,(pos_m=fill(NaN,3,2),)))
        write(prefix*"_scene.json",JSON.json(Dict("epoch"=>Dict("et_start_s"=>0.0),"planet"=>Dict("name"=>"Earth"),"spacecraft"=>[Dict("id"=>7),Dict("id"=>91)])))
        mission=MDH.MissionSpice("fixture","-1","EARTH",String[])
        @test_throws ArgumentError MDH.spice_reference(mission,prefix;stride=0)
        @test_throws ArgumentError MDH.spice_reference(mission,prefix;target=91)
        @test_throws ArgumentError MDH.spice_reference(mission,prefix;opacity=NaN)
        @test_throws ArgumentError MDH.spice_reference(MDH.MissionSpice("bad center","-1","MARS",String[]),prefix)
        @test MDH.build_cdn_page("unused","unused","a","b","c","d","e") === nothing
        @test_throws ArgumentError MDH.ensure_mission_kernel("../outside.bsp","https://naif.jpl.nasa.gov/a")
        @test_throws ArgumentError MDH.ensure_mission_kernel("file.bsp","http://example.invalid/a")
    end
end

# Include each real driver in its own module. Their main guards must prevent
# network access, native loading, propagation, or output creation merely on load.
module ApolloMissionDriverTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "apollo11_lunar_orbit.jl"))
end
module MagellanMissionDriverTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "magellan_aerobraking.jl"))
end
module CassiniMissionDriverTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "cassini_titan_flyby.jl"))
end
module OdysseyMissionDriverTest
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "odyssey_aerobraking.jl"))
end
@testset "mission driver load boundaries" begin
    @test all(isdefined(m,:main) for m in (ApolloMissionDriverTest,MagellanMissionDriverTest,CassiniMissionDriverTest,OdysseyMissionDriverTest))
    @test !isdefined(MDH,:GRAMSuite)
    @test CassiniMissionDriverTest.FLYBYS["TA"].ca_guess == "2004-10-26T15:30:00"
    @test_throws ArgumentError CassiniMissionDriverTest.main(["all","--output-dir","not-used"])
    @test_throws ArgumentError CassiniMissionDriverTest.run_flyby("bad")
end

@testset "mission epoch agrees with the engine and saved SPICE reference" begin
    # Reuse the Earth starter-pack kernels already used by the regular suite.
    # The Moon is a moving SPK target without downloading a mission kernel or
    # constructing a native atmosphere. Missing files are an explicit skip;
    # malformed kernels, conversion failures and assertion failures are errors.
    spice_dir = MDH.SPICE_PATH
    planetary_override = strip(get(ENV, "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH", ""))
    planetary_candidates = isempty(planetary_override) ? (
        "spk/planets/de430.bsp", "spk/planets/de421.bsp",
        "spk/planets/de442s.bsp", "spk/planets/de442.bsp",
        "spk/planets/de440s.bsp", "spk/planets/de440_GRAM.bsp",
    ) : (planetary_override,)
    required_present = all(isfile(joinpath(spice_dir, name)) for name in
        ("lsk/naif0012.tls", "pck/pck00011.tpc")) &&
        any(isfile(joinpath(spice_dir, name)) for name in planetary_candidates)
    if !required_present
        @info "Skipping mission epoch SPICE regression: required starter-pack kernel files are absent" spice_dir
        @test_skip false
    else
        sm = MDH.SM
        planet = sm.Earth("", spice_dir)
        mission = MDH.MissionSpice("Moon SPICE clock reference", "MOON", "EARTH", String[])
        # Expected UTC strings are independent, explicit millisecond-rounded
        # answers. The four fractions reproduce the mission startup cases;
        # the final two exercise ordinary minute/day rollover, not leap seconds.
        cases = (
            ("1993-05-26T00:00:07.282648", "1993-05-26T00:00:07.283"),
            ("2001-11-06T10:05:12.693351", "2001-11-06T10:05:12.693"),
            ("2004-10-26T13:00:04.608162", "2004-10-26T13:00:04.608"),
            ("2005-04-16T16:41:45.772100", "2005-04-16T16:41:45.772"),
            ("2001-01-01T00:00:59.999600", "2001-01-01T00:01:00.000"),
            ("2001-01-01T23:59:59.999600", "2001-01-02T00:00:00.000"),
        )
        for (requested_utc, expected_utc) in cases
            @testset "$requested_utc" begin
                requested_et = MDH.et_of(requested_utc)
                initial_time = MDH.initial_time_of(requested_et)
                resolved_et = MDH.et_of(initial_time)
                expected_et = MDH.et_of(expected_utc)
                et_atol = 4 * eps(abs(expected_et))
                @test resolved_et ≈ expected_et rtol=0 atol=et_atol
                @test sm.ephemerides_time_seconds(initial_time, sm.SpiceEphemeridesModel()) ≈ expected_et rtol=0 atol=et_atol
                @test MDH.et_of(MDH.initial_time_of(resolved_et)) ≈ expected_et rtol=0 atol=et_atol
                @test abs(requested_et - expected_et) > 5e-5

                # Direct SPICE is the position/velocity oracle. A relative
                # tolerance at lunar distances could conceal the old metre or
                # sub-metre mismatch, so every comparison sets rtol=0.
                expected_state = lock(MDH.RuntimeServices.SPICE_LOCK) do
                    MDH.spkezr("MOON", expected_et, "J2000", "NONE", "EARTH")[1] .* 1e3
                end
                ic = MDH.cartesian_ic_at(mission, resolved_et)
                @test ic.pos ≈ expected_state[1:3] rtol=0 atol=1e-3
                @test ic.vel ≈ expected_state[4:6] rtol=0 atol=1e-8
                # Negative control: the original raw-search-ET initialization
                # must be distinguishable from the saved reference epoch.
                old_ic = MDH.cartesian_ic_at(mission, requested_et)
                @test MDH.norm(old_ic.pos - expected_state[1:3]) > 0.01
            end
        end

        mktempdir() do output
            requested_utc, expected_utc = cases[2]
            initial_time = MDH.initial_time_of(MDH.et_of(requested_utc))
            expected_et = MDH.et_of(expected_utc)
            et_atol = 4 * eps(abs(expected_et))
            ic = MDH.cartesian_ic_at(mission, MDH.et_of(initial_time))
            body = sm.Link(root=true, m=100.0, ref_area=1.0)
            craft = sm.SpacecraftModel(links=[body], root=body, initial_condition=ic, id=91)
            base = MDH.make_example_config(planet=planet, spacecraft=craft,
                mission_time=1.0, initial_time=initial_time, dynamic_effectors=(),
                density_model=sm.NoAtmosphereModel(), ephemerides_model=sm.SpiceEphemeridesModel(),
                orientation_sim=false, keplerian=false, verbose=false,
                results=true, results_directory=output)
            args = sm.SimConfig._with_configuration(base;
                mission_configuration=sm.MissionConfiguration(mission_type=sm.MissionTime,
                    mission_time=1.0, number_of_orbits=1, keplerian=false,
                    orientation_sim=false, num_steps_to_save=3, data_rate=0.5),
                solver_config=sm.SolverConfig(solver_mode=:tsit5))
            solution = SpaceAGORA.run_simulation(args; visualization=true, return_solution=true)
            @test solution.prob.p.shared_buffers.et_start[] ≈ expected_et rtol=0 atol=et_atol
            prefix = joinpath(output, "simulation_results")
            scene = JSON.parsefile(prefix * "_scene.json")
            saved = DataFrame(Arrow.Table(prefix * ".feather"))
            @test scene["epoch"]["et_start_s"] ≈ expected_et rtol=0 atol=et_atol
            @test scene["epoch"]["utc"] == expected_utc * "Z"
            @test saved.time[1] == 0.0
            @test saved.time[end] == 1.0
            saved_position = [saved[1, Symbol("sc1_pos_", axis)] for axis in 1:3]
            saved_velocity = [saved[1, Symbol("sc1_vel_", axis)] for axis in 1:3]
            @test saved_position ≈ ic.pos rtol=0 atol=1e-3
            @test saved_velocity ≈ ic.vel rtol=0 atol=1e-8
            reference = MDH.spice_reference(mission, prefix)
            @test reference.t_s == saved.time
            @test reference.pos_m[:, 1] ≈ saved_position rtol=0 atol=1e-3
            @test reference.vel_mps[:, 1] ≈ saved_velocity rtol=0 atol=1e-8
            @test first(MDH.reference_separation(prefix, reference)) <= 1e-3
            # This one-second force-free run checks the real clock/export path;
            # it does not assert that a free particle follows the Moon's orbit.
        end
    end
end
