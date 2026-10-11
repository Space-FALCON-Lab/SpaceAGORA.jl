module PassExponentialDensityTests

using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays
using LinearAlgebra
using Statistics: median

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels

const MU = 3.24858592e14           # Venus GM [m^3/s^2]
const RV = 6.0518e6

# fallback: a constant density, so fallback answers are recognisable
_fb() = SM.ConstantDensityModel(density_kg_m3=7.0e-9, temperature_k=222.0)

# Two profiles: flight passes 2951 (pdyn) and 2952 (dv)
_model(; offset=2949, ei=250.0e3) = EM.PassKeyedExponentialAtmosphereModel(
    _fb(), [2951, 2952], [0.01, NaN], [NaN, 0.5], [150.0e3, 140.0e3], [3500.0, 5000.0];
    counter_pass_offset=offset, entry_interface_m=ei)

# state at periapsis radius rp of an orbit with apoapsis radius ra
function _peri_state(rp, ra)
    a = 0.5 * (rp + ra)
    vp = sqrt(MU * (2 / rp - 1 / a))
    return SVector(rp, 0.0, 0.0), SVector(0.0, vp, 0.0), vp
end

@testset "per-pass exponential density" begin
    rp = RV + 150.0e3
    pos, vel, vp = _peri_state(rp, RV + 66000.0e3)
    m = _model()

    # profile formula: rho = 2 pdyn / Vp^2 * exp(-(h - hp)/H), Vp the osculating pericentre speed
    for h in (150.0e3, 140.0e3, 160.0e3)
        @test EM.pass_exponential_density(m, h, pos, vel, 2, MU) ≈ 2 * 0.01 / vp^2 * exp(-(h - 150.0e3) / 3500.0) rtol = 1e-12
    end
    # the same pericentre speed is recovered from a state on the inbound leg of the same orbit
    ra = RV + 66000.0e3; ecc = (ra - rp) / (ra + rp); pp = rp * (1 + ecc); nu = -deg2rad(3.0)
    rnu = pp / (1 + ecc * cos(nu)); vq = sqrt(MU / pp)
    posn = SVector(rnu * cos(nu), rnu * sin(nu), 0.0)
    veln = vq * SVector(-sin(nu) * 1.0 * cos(0.0) , (ecc + cos(nu)), 0.0)                      # perifocal velocity (x along periapsis)
    @test norm(cross(posn, veln)) ≈ norm(cross(pos, vel)) rtol = 1e-9
    @test EM.pass_exponential_density(m, 150.0e3, posn, veln, 2, MU) ≈ EM.pass_exponential_density(m, 150.0e3, pos, vel, 2, MU) rtol = 1e-9
    # a circular orbit at the entry state: Vp = circular speed
    pc = SVector(rp, 0.0, 0.0); vc = sqrt(MU / rp)
    @test EM.pass_exponential_density(m, 150.0e3, pc, SVector(0.0, vc, 0.0), 2, MU) ≈ 2 * 0.01 / vc^2 rtol = 1e-9
end

@testset "pass keying" begin
    pos, vel, vp = _peri_state(RV + 150.0e3, RV + 66000.0e3)
    m = _model()
    # counter 2 -> flight pass 2951 (pdyn 0.01, H 3.5 km), counter 3 -> 2952 (dv given, H 5 km)
    r2 = EM.pass_exponential_density(m, 145.0e3, pos, vel, 2, MU)
    r3 = EM.pass_exponential_density(m, 145.0e3, pos, vel, 3, MU)
    @test r2 ≈ 2 * 0.01 / vp^2 * exp(-(145.0e3 - 150.0e3) / 3500.0) rtol = 1e-12
    @test r3 != r2
    # dv mode inverts Damiani Eq. 4: dv = pdyn (Cd S / m) sqrt(2 pi H rp) / Vp
    rpp = norm(pos); H = 5000.0
    pdyn = 0.5 * 650.0 * vp / (2.2 * 10.4 * sqrt(2pi * H * rpp))
    @test r3 ≈ 2 * pdyn / vp^2 * exp(-(145.0e3 - 140.0e3) / H) rtol = 1e-12
    # the offset moves the key
    mo = _model(offset=2950)
    @test EM.pass_exponential_density(mo, 145.0e3, pos, vel, 1, MU) ≈ r2 rtol = 1e-12
end

@testset "fallback" begin
    pos, vel, _ = _peri_state(RV + 150.0e3, RV + 66000.0e3)
    m = _model()
    @test EM.pass_exponential_density(m, 145.0e3, pos, vel, 1, MU) === nothing     # no profile: counter 1 -> 2950
    @test EM.pass_exponential_density(m, 145.0e3, pos, vel, 4, MU) === nothing     # past the table (no nearest-pass substitution)
    @test EM.pass_exponential_density(m, 250.0e3, pos, vel, 2, MU) === nothing     # at the entry interface
    @test EM.pass_exponential_density(m, 400.0e3, pos, vel, 2, MU) === nothing     # above it
    @test EM.pass_exponential_density(m, 249.999e3, pos, vel, 2, MU) isa Float64   # just below it
    @test isnan(EM.pass_exponential_density(m, NaN, pos, vel, 2, MU))              # solver rejects the step
    # state-free queries answer with the fallback alone, temperature included
    rho, T, w = SM.getDensity(m, 145.0e3, 0.0, 0.0, 0.0, true, nothing)
    @test rho == 7.0e-9 && T == 222.0
    # an empty table is a pure pass-through (the logging control)
    me = EM.PassKeyedExponentialAtmosphereModel(_fb(), Int[], Float64[], Float64[], Float64[], Float64[]; counter_pass_offset=2949, entry_interface_m=250.0e3)
    @test EM.pass_exponential_density(me, 145.0e3, pos, vel, 2, MU) === nothing
    # construction guards
    @test_throws ArgumentError EM.PassKeyedExponentialAtmosphereModel(_fb(), [2952, 2951], [1.0, 1.0], [NaN, NaN], [1.0, 1.0], [3.5, 3.5]; counter_pass_offset=0, entry_interface_m=250.0e3)
    @test_throws ArgumentError EM.PassKeyedExponentialAtmosphereModel(_fb(), [1], [1.0], [1.0], [1.0], [3.5]; counter_pass_offset=0, entry_interface_m=250.0e3)
    @test_throws ArgumentError EM.PassKeyedExponentialAtmosphereModel(_fb(), [1], [NaN], [NaN], [1.0], [3.5]; counter_pass_offset=0, entry_interface_m=250.0e3)
end

const _MANIFEST_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "telemetry_validation", "manifests"))
_load(n) = only(SpaceAGORA.TelemetryVerification._load_scenarios_from_manifest(joinpath(_MANIFEST_DIR, "$n.toml")))
const _ORVM = ("vex_orvm_g115", "vex_orvm_gfit", "vex_orvm_f", "vex_orvm_vg", "vex_orvm_vgfit")

@testset "manifest keys" begin
    TV = SpaceAGORA.TelemetryVerification
    t(extra...) = Dict("atmosphere_truth" => Dict("atmosphere_model" => "gram_pass_exponential", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f", extra...))
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(t(), "ctx")                       # mode without file
    cfg = TV._parse_atmosphere_truth_config(t("pass_exponential_file" => "/x.csv", "pass_exponential_counter_offset" => 2949), "ctx")
    @test cfg.atmosphere_model == "gram_pass_exponential" && cfg.pass_exponential_counter_offset == 2949
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict(
        "atmosphere_model" => "GRAM", "atmosphere_dataset" => "d", "space_weather_model" => "s",
        "solar_flux_model" => "f", "pass_exponential_file" => "/x.csv")), "ctx")                     # file without mode
    # the shipped measured-drag manifest parses and points at an existing profile file
    sc = _load("vex_orvm_f")
    @test sc.atmosphere_truth.atmosphere_model == "gram_pass_exponential"
    @test isfile(sc.atmosphere_truth.pass_exponential_file)
    @test sc.atmosphere_truth.pass_exponential_counter_offset == 2948
end


@testset "periapsis pulse" begin
    SC = SM.SimulationCallbacks
    rp = RV + 150.0e3; ra = RV + 66000.0e3
    pos, vel, vp = _peri_state(rp, ra)
    P0, rp0 = SC.periapsis_pulse_osculating(pos, vel, MU)
    a = 0.5 * (rp + ra)
    @test P0 ≈ 2pi * sqrt(a^3 / MU) rtol = 1e-12
    @test rp0 ≈ rp rtol = 1e-12
    # a 0.04 m/s retrograde pulse at periapsis: period falls by 3 P dv Vp a / (mu) (first order), periapsis radius is unchanged to second order
    dv = 0.04
    P1, rp1 = SC.periapsis_pulse_osculating(pos, vel * (1 - dv / vp), MU)
    da = 2 * a^2 * vp * dv / MU
    @test (P1 - P0) ≈ -1.5 * P0 * da / a rtol = 1e-3
    @test abs(rp1 - rp0) < 1e-2
    # construction guards
    @test_throws ArgumentError SC.get_periapsis_pulse_callback(-0.04, 2951:3005)
    @test_throws ArgumentError SC.get_periapsis_pulse_callback(NaN, 2951:3005)
    @test SC.get_periapsis_pulse_callback(0.04, 2951:3005; counter_pass_offset=2949) isa Any
    # manifests: off unless the block says so; main's record manifests carry none of the new options
    for n in ("vex_venusgram", "odyssey_tolson", "odyssey_marsgram", "vex_orvm_g115", "vex_orvm_vg")
        sc = _load(n)
        @test sc.thruster_pulse_dv_mps == 0.0 && sc.calibration.cd_scale_min == 1.0 && !sc.calibration.fit_cd_scale
    end
    f = _load("vex_orvm_f")
    @test f.thruster_pulse_dv_mps == 0.04
    @test (f.thruster_pulse_first_pass, f.thruster_pulse_last_pass, f.thruster_pulse_counter_pass_offset) == (2950, 3005, 2948)
end


@testset "periapsis pulse in a propagation" begin
    # Regression: under DiffEqBase 7 the pulse callback and the orbit counter's apsis callback both locate
    # the periapsis root of r.v and re-detected it in turn forever after the first pulse. Over about
    # two orbits from just after apoapsis (orbit counting on), each pass gets exactly one pulse and the run completes.
    planet = Earth("", joinpath(@__DIR__, "..", "..", "..", "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE"))
    root = Link(root=true, m=140.0, ref_area=1.2)
    ic = InitialCondition(ra=planet.Rp_e + 900e3, rp=planet.Rp_e + 800e3, i=28.0, ω=15.0, Ω=20.0, ν=181.0)
    sc = SpacecraftModel(joints=Joint[], links=Link[root], root=root, instant_actuation=true, prop_mass=15.0,
        inertia_tensor=root.inertia, n_reaction_wheels=0, n_thrusters=0, initial_condition=ic, id=1)
    period = 2pi * sqrt((planet.Rp_e + 850e3)^3 / planet.μ)
    cfg = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionOrbits, keplerian=true, number_of_orbits=2,
            mission_time=4.0 * period, orientation_sim=false, num_steps_to_save=50),
        environment_model=EnvironmentModel(planet=planet, EI=120.0, density_model=ExponentialAtmosphereModel(planet),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0))
    text = withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        mktempdir() do tmp
            cd(tmp) do
                log_path = joinpath(tmp, "pulse_log.txt")
                open(log_path, "w") do f
                    redirect_stdout(f) do
                        run_simulation(cfg; return_solution=true,
                            extra_callbacks=(SM.SimulationCallbacks.get_periapsis_pulse_callback(0.04, 1:10),))
                    end
                end
                read(log_path, String)
            end
        end
    end
    pulses = [parse(Int, m.captures[1]) for m in eachmatch(r"periapsis_pulse sat=1 pass=(\d+)", text)]
    @test pulses == [1, 2]
end


@testset "density scale k" begin
    pos, vel, _ = _peri_state(RV + 150.0e3, RV + 66000.0e3)
    m1 = _model()
    mk(k) = EM.PassKeyedExponentialAtmosphereModel(_fb(), [2951, 2952], [0.01, NaN], [NaN, 0.5], [150.0e3, 140.0e3], [3500.0, 5000.0];
        counter_pass_offset=2949, entry_interface_m=250.0e3, rho_scale=k)
    base = EM.pass_exponential_density(m1, 145.0e3, pos, vel, 2, MU)
    @test m1.rho_scale === 1.0
    @test EM.pass_exponential_density(mk(1.0), 145.0e3, pos, vel, 2, MU) === base          # k = 1 is bit-identical
    @test EM.pass_exponential_density(mk(1.3), 145.0e3, pos, vel, 2, MU) ≈ 1.3 * base rtol = 1e-14
    @test EM.pass_exponential_density(mk(0.7), 145.0e3, pos, vel, 3, MU) ≈ 0.7 * EM.pass_exponential_density(m1, 145.0e3, pos, vel, 3, MU) rtol = 1e-14
    @test EM.pass_exponential_density(mk(1.3), 145.0e3, pos, vel, 1, MU) === nothing        # fallback passes are not scaled
    @test EM.pass_exponential_density(mk(1.3), 260.0e3, pos, vel, 2, MU) === nothing        # nor altitudes above the entry interface
    @test_throws ArgumentError mk(0.0)
    @test_throws ArgumentError mk(-1.0)
    @test_throws ArgumentError mk(NaN)
    TV = SpaceAGORA.TelemetryVerification
    t(extra...) = Dict("atmosphere_truth" => Dict("atmosphere_model" => "gram_pass_exponential", "atmosphere_dataset" => "d",
        "space_weather_model" => "s", "solar_flux_model" => "f", "pass_exponential_file" => "/x.csv", extra...))
    @test TV._parse_atmosphere_truth_config(t(), "ctx").pass_exponential_scale == 1.0
    @test TV._parse_atmosphere_truth_config(t("pass_exponential_scale" => 1.25), "ctx").pass_exponential_scale == 1.25
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(t("pass_exponential_scale" => 0.0), "ctx")
    @test_throws ArgumentError TV._parse_atmosphere_truth_config(Dict("atmosphere_truth" => Dict("atmosphere_model" => "GRAM",
        "atmosphere_dataset" => "d", "space_weather_model" => "s", "solar_flux_model" => "f", "pass_exponential_scale" => 1.2)), "ctx")
    @test _load("vex_orvm_f").atmosphere_truth.pass_exponential_scale == 1.0
end


@testset "pass-aligned comparison axis" begin
    TV = SpaceAGORA.TelemetryVerification
    tele = [0.006, 1.988, 3.024, 4.012, 5.047]
    ax = TV._pass_aligned_sim_axis(tele, 8)
    @test ax[1:5] == tele                                             # the compared events take the flight labels
    @test all(diff(ax) .> 0) && length(ax) == 8
    @test ax[6] ≈ tele[end] + median(diff(tele)) rtol = 1e-14         # events past the compared points extend at the median step
    @test TV._pass_aligned_sim_axis(tele, 3) == tele[1:3]             # fewer events than points
    # pairing: interpolating the simulation on this axis at the flight labels returns event i for flight point i, and the
    # error is simulated minus flight index for index; on the legacy axis (median step) a +1 orbit label offset pairs i with i + 1
    sim = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    flight = [11.0, 21.0, 31.0, 41.0, 51.0]
    _, rows = TV._compare_orbit_curve("t", "peri", tele, flight, sim; sim_axis=ax)
    @test rows.error_km ≈ fill(-1.0, 5) atol = 1e-9
    legacy_axis = tele[1] .+ median(diff(tele)) .* collect(0:7)
    _, rowsl = TV._compare_orbit_curve("t", "peri", tele, flight, sim; sim_axis=legacy_axis)
    @test !(rowsl.error_km ≈ fill(-1.0, 5))
    # manifest option: default legacy, refused when unknown or combined with epoch_orbit_offset
    @test TV._parse_comparison_axis(Dict{String, Any}(), "ctx") === :legacy
    @test TV._parse_comparison_axis(Dict("comparison_axis" => "pass_aligned"), "ctx") === :pass_aligned
    @test_throws ArgumentError TV._parse_comparison_axis(Dict("comparison_axis" => "bogus"), "ctx")
    @test_throws ArgumentError TV._parse_comparison_axis(Dict("comparison_axis" => "pass_aligned", "epoch_orbit_offset" => 19.0), "ctx")
    # main's record manifests keep the legacy axis and no skipped event
    for n in ("vex_venusgram", "odyssey_tolson", "odyssey_marsgram")
        sc = _load(n)
        @test sc.comparison_axis === :legacy && sc.pass_aligned_skipped_event == 0 && sc.pass_aligned_apoapsis_skipped_event == -1
    end
end


@testset "pass-aligned axis with a skipped orbit" begin
    TV = SpaceAGORA.TelemetryVerification
    tele = [0.006, 1.988, 3.024, 4.012, 5.047]
    ax = TV._pass_aligned_sim_axis(tele, 8; skipped_event=1)
    @test ax[1] == tele[1]                                  # event 0 with flight point 0
    @test ax[2] ≈ 0.5 * (tele[1] + tele[2])                 # the skipped event sits between points 0 and 1
    @test ax[3:6] == tele[2:5]                              # events i + 1 with flight points i >= 1
    @test all(diff(ax) .> 0) && length(ax) == 8
    @test TV._pass_aligned_sim_axis(tele, 8; skipped_event=0) == TV._pass_aligned_sim_axis(tele, 8)
    @test_throws ArgumentError TV._pass_aligned_sim_axis(tele, 8; skipped_event=5)
    # pairing: flight point 0 against event 0, point i >= 1 against event i + 1
    sim = [10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    flight = [11.0, 21.0, 31.0, 41.0, 51.0]
    _, rows = TV._compare_orbit_curve("t", "peri", tele, flight, sim; sim_axis=ax)
    @test rows.error_km ≈ fill(-1.0, 5) atol = 1e-9
    @test TV._parse_pass_aligned_skipped_event(Dict{String, Any}(), "ctx") == 0
    @test TV._parse_pass_aligned_skipped_event(Dict("comparison_axis" => "pass_aligned", "pass_aligned_skipped_event" => 1), "ctx") == 1
    @test_throws ArgumentError TV._parse_pass_aligned_skipped_event(Dict("pass_aligned_skipped_event" => 1), "ctx")
    @test_throws ArgumentError TV._parse_pass_aligned_skipped_event(Dict("comparison_axis" => "pass_aligned", "pass_aligned_skipped_event" => 0), "ctx")
    @test TV._parse_pass_aligned_apoapsis_skipped_event(Dict{String, Any}(), "ctx") == -1
    @test TV._parse_pass_aligned_apoapsis_skipped_event(Dict("comparison_axis" => "pass_aligned", "pass_aligned_apoapsis_skipped_event" => 0), "ctx") == 0
    @test_throws ArgumentError TV._parse_pass_aligned_apoapsis_skipped_event(Dict("pass_aligned_apoapsis_skipped_event" => 0), "ctx")
    @test_throws ArgumentError TV._parse_pass_aligned_apoapsis_skipped_event(Dict("comparison_axis" => "pass_aligned", "pass_aligned_apoapsis_skipped_event" => -2), "ctx")
end


@testset "Venus Express ORVM-referenced manifests" begin
    for n in _ORVM
        sc = _load(n)
        @test endswith(sc.telemetry_peri_path, "VEx/orvm/vex_orvm_periapsis.feather") && endswith(sc.telemetry_apo_path, "VEx/orvm/vex_orvm_apoapsis.feather")
        @test sc.comparison_axis === :pass_aligned && sc.pass_aligned_skipped_event == 0 && sc.pass_aligned_apoapsis_skipped_event == -1
        @test sc.epoch_orbit_offset === nothing
        orvm_ic = startswith(n, "vex_orvm_v")
        @test (sc.initial_state_j2000_m !== nothing) == orvm_ic
        @test sc.maneuver_orbit_numbers == (orvm_ic ? [6, 36, 41, 45] : [6, 37, 42, 46])
        @test sc.maneuver_delta_v_mps == [0.428, -0.177, -0.07, -0.05]
        @test sc.atmosphere_truth.atmosphere_model == (n == "vex_orvm_f" ? "gram_pass_exponential" : "GRAM")
    end
    cal(sc) = (sc.calibration.cr_min, sc.calibration.cr_max, sc.calibration.cr_steps, sc.calibration.fit_bias, sc.calibration.objective)
    @test cal(_load("vex_orvm_g115")) == cal(_load("vex_orvm_vg")) == cal(_load("vex_orvm_f")) == (1.15, 1.15, 1, true, "mean_nmae")
    @test cal(_load("vex_orvm_gfit")) == cal(_load("vex_orvm_vgfit")) == (1.15, 1.45, 3, true, "mean_nmae")
    # the measured-drag profile: counter c -> key c + 2948; counter 2 -> 2950 (point 0), counter 3 -> 2951 (absent:
    # VenusGRAM), counter 4 -> 2952, counter 57 -> 3005 (last)
    TV = SpaceAGORA.TelemetryVerification
    prof = TV.DataFrame(TV.CSV.File(_load("vex_orvm_f").atmosphere_truth.pass_exponential_file))
    keys = Set(prof.pass)
    @test (2 + 2948 in keys) && !(3 + 2948 in keys) && (4 + 2948 in keys) && (57 + 2948 in keys) && !(58 + 2948 in keys)
end

end # module
