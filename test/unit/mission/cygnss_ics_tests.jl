using Test
using JSON
using LinearAlgebra
using StaticArrays
using SPICE
using Arrow
using DataFrames

const CYGIC_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(CYGIC_REPO, "scripts", "dev", "viewer_demos", "cygnss_ics.jl"))
using .CygnssICs

const CYGIC_TLE_FILE = joinpath(CYGIC_REPO, "scripts", "dev", "viewer_demos", "cygnss_historical.tle")
const CYGIC_ICS_FILE = joinpath(CYGIC_REPO, "data", "telemetry", "CYGNSS", "constellation_ics_20250606.json")
const CYGIC_SPICE_PATH = joinpath(CYGIC_REPO, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const CYGIC_FM04_FEATHER = joinpath(CYGIC_REPO, "data", "telemetry", "CYGNSS", "cyg04_nasa_pvt_96hr.feather")

cygic_rz(θ) = SMatrix{3, 3, Float64}(cos(θ), -sin(θ), 0.0, sin(θ), cos(θ), 0.0, 0.0, 0.0, 1.0)
cygic_rx(θ) = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, cos(θ), -sin(θ), 0.0, sin(θ), cos(θ))
# d/dθ of cygic_rz, so that dR/dt = ω * cygic_drz(θ) for uniform spin.
cygic_drz(θ) = SMatrix{3, 3, Float64}(-sin(θ), -cos(θ), 0.0, cos(θ), -sin(θ), 0.0, 0.0, 0.0, 0.0)

@testset "CygnssICs" begin

    # =======================================================================
    # Frame conversions
    # =======================================================================
    @testset "frame conversion" begin
        ω = 7.2921158553e-5
        R = 6.8e6

        @testset "a point at rest on a rotating body moves at omega x r" begin
            for θ in (0.0, 0.7, 2.6, 5.9)
                l_pi = cygic_rz(θ)                     # inertial -> body-fixed
                r_f = SVector{3, Float64}(R, 0.0, 0.0)
                r_i, v_i = rotating_to_inertial(l_pi, SVector(0.0, 0.0, ω), r_f, zeros(3))
                @test norm(r_i) ≈ R rtol = 1e-12
                # Velocity is the rigid-body term: perpendicular to r and to the
                # spin axis, with magnitude ω|r|.
                @test norm(v_i) ≈ ω * R rtol = 1e-12
                @test abs(dot(v_i, r_i)) < 1e-6
                @test v_i ≈ cross(SVector(0.0, 0.0, ω), r_i) rtol = 1e-12
            end
        end

        @testset "the spin axis is body-fixed, not inertial" begin
            # With the body-fixed frame tilted away from the inertial z-axis,
            # forming the transport term as omega x r_inertial with the same
            # numbers treats a body-fixed vector as an inertial one. The two
            # answers must differ, and by about omega*R*sin(tilt).
            tilt = deg2rad(0.35)                        # J2000-to-2025 precession scale
            l_pi = cygic_rx(tilt) * cygic_rz(1.1)
            r_f = SVector{3, Float64}(R, 0.0, 0.0)
            w_f = SVector{3, Float64}(0.0, 0.0, ω)
            _, v_correct = rotating_to_inertial(l_pi, w_f, r_f, zeros(3))
            r_i = l_pi' * r_f
            v_wrong = cross(w_f, r_i)                   # the mistake, written out
            @test norm(v_correct - v_wrong) ≈ ω * R * sin(tilt) rtol = 0.05
            @test norm(v_correct - v_wrong) > 2.0       # m/s: not a rounding difference
        end

        @testset "exact derivative form matches uniform rotation" begin
            θ = 2.3
            l_pi = cygic_rz(θ)
            dl_pi = ω .* cygic_drz(θ)
            r_f = SVector{3, Float64}(4.0e6, -3.0e6, 2.5e6)
            v_f = SVector{3, Float64}(1200.0, 3400.0, -900.0)
            ra, va = rotating_to_inertial(l_pi, SVector(0.0, 0.0, ω), r_f, v_f)
            rb, vb = earth_fixed_to_inertial(l_pi, dl_pi, r_f, v_f)
            @test ra ≈ rb rtol = 1e-14
            @test va ≈ vb rtol = 1e-12
        end

        @testset "conversion is invertible" begin
            θ = 0.9
            l_pi = cygic_rz(θ)
            dl_pi = ω .* cygic_drz(θ)
            r_f = SVector{3, Float64}(5.4e6, -2.5e6, 3.3e6)
            v_f = SVector{3, Float64}(1657.0, 6635.0, 2375.0)
            r_i, v_i = earth_fixed_to_inertial(l_pi, dl_pi, r_f, v_f)
            # Inverse of [R' 0; dR' R'] is [R 0; -R dR' R  R].
            r_back = l_pi * r_i
            v_back = l_pi * v_i - l_pi * dl_pi' * r_back
            @test r_back ≈ r_f rtol = 1e-12
            @test v_back ≈ v_f rtol = 1e-10
        end

        @testset "TEME to J2000 is a pure rotation" begin
            # A quasi-inertial pair: one rotation carries r and v alike, so
            # every rotation invariant survives.
            Rm = cygic_rx(0.004) * cygic_rz(0.0011)
            r_t = SVector{3, Float64}(-3.89e6, -4.53e6, 3.29e6)
            v_t = SVector{3, Float64}(6277.0, -3670.0, 2359.0)
            r_j, v_j = teme_to_j2000(Rm, r_t, v_t)
            @test norm(r_j) ≈ norm(r_t) rtol = 1e-14
            @test norm(v_j) ≈ norm(v_t) rtol = 1e-14
            @test dot(r_j, v_j) ≈ dot(r_t, v_t) rtol = 1e-12
            @test norm(cross(r_j, v_j)) ≈ norm(cross(r_t, v_t)) rtol = 1e-13
            @test teme_to_j2000(Rm, r_t, v_t)[1] ≈ Rm * r_t rtol = 1e-14
        end

        @testset "an Earth-fixed velocity read as inertial is impossible" begin
            # This is the test that establishes the PVT columns' frame, written
            # as a check rather than an assumption. Take a real circular orbit
            # at CYGNSS radius and inclination, express it Earth-fixed, and show
            # that reading those numbers as inertial puts the semimajor axis
            # inside the Earth.
            μ = CYGNSS_MU_EARTH
            r_e = CYGNSS_R_EARTH_EQ
            rmag = 6.8185e6
            incl = deg2rad(34.94)
            r_i = SVector{3, Float64}(rmag, 0.0, 0.0)
            vc = sqrt(μ / rmag)
            v_i = SVector{3, Float64}(0.0, vc * cos(incl), vc * sin(incl))

            l_pi = cygic_rz(0.0)                        # frames aligned at this instant
            w_f = SVector{3, Float64}(0.0, 0.0, 7.2921158553e-5)
            r_f = l_pi * r_i
            v_f = l_pi * v_i - cross(w_f, r_f)          # inverse of rotating_to_inertial

            @test norm(v_f) < norm(v_i)                 # prograde: Earth-fixed is slower
            @test norm(v_i) - norm(v_f) > 300.0         # m/s, the scale actually seen

            a_if_inertial = 1.0 / (2.0 / rmag - norm(v_f)^2 / μ)
            @test a_if_inertial < r_e                   # below the surface: impossible
            a_correct = 1.0 / (2.0 / rmag - norm(v_i)^2 / μ)
            @test a_correct ≈ rmag rtol = 1e-10

            # And the round trip recovers the inertial state.
            r_back, v_back = rotating_to_inertial(l_pi, w_f, r_f, v_f)
            @test r_back ≈ r_i rtol = 1e-12
            @test v_back ≈ v_i rtol = 1e-12
        end
    end

    # =======================================================================
    # Classical elements
    # =======================================================================
    @testset "classical elements" begin
        μ = CYGNSS_MU_EARTH

        @testset "circular equatorial" begin
            rmag = 7.0e6
            r = SVector{3, Float64}(rmag, 0.0, 0.0)
            v = SVector{3, Float64}(0.0, sqrt(μ / rmag), 0.0)
            el = classical_elements(r, v)
            @test el.a_m ≈ rmag rtol = 1e-12
            @test el.e < 1e-12
            @test el.inclination_deg < 1e-10
            @test el.period_s ≈ 2π * sqrt(rmag^3 / μ) rtol = 1e-12
            @test el.altitude_m ≈ rmag - CYGNSS_R_EARTH_EQ rtol = 1e-12
        end

        @testset "inclined circular recovers RAAN and argument of latitude" begin
            rmag = 6.82e6
            for (incl_deg, raan_deg, u_deg) in ((34.94, 177.4, 57.3), (34.88, 201.6, 295.3), (35.0, 10.0, 350.0))
                i = deg2rad(incl_deg); Ω = deg2rad(raan_deg); u = deg2rad(u_deg)
                # Perifocal-to-inertial for a circular orbit, node at Ω.
                n̂ = SVector{3, Float64}(cos(Ω), sin(Ω), 0.0)
                ĥ = SVector{3, Float64}(sin(i) * sin(Ω), -sin(i) * cos(Ω), cos(i))
                ŵ = cross(ĥ, n̂)
                r = rmag .* (cos(u) .* n̂ .+ sin(u) .* ŵ)
                v = sqrt(μ / rmag) .* (-sin(u) .* n̂ .+ cos(u) .* ŵ)
                el = classical_elements(r, v)
                @test el.inclination_deg ≈ incl_deg rtol = 1e-9
                @test el.raan_deg ≈ raan_deg rtol = 1e-9
                @test el.arg_latitude_deg ≈ u_deg rtol = 1e-8
                @test el.e < 1e-12
            end
        end

        @testset "argument of latitude survives near-zero eccentricity" begin
            # arg_perigee + true_anomaly is ill-conditioned at e ~ 1e-4; the
            # node-referenced angle is not. Perturb a circular orbit and check
            # the angle barely moves.
            rmag = 6.82e6
            i = deg2rad(34.94); Ω = deg2rad(120.0); u = deg2rad(200.0)
            n̂ = SVector{3, Float64}(cos(Ω), sin(Ω), 0.0)
            ĥ = SVector{3, Float64}(sin(i) * sin(Ω), -sin(i) * cos(Ω), cos(i))
            ŵ = cross(ĥ, n̂)
            r = rmag .* (cos(u) .* n̂ .+ sin(u) .* ŵ)
            vdir = -sin(u) .* n̂ .+ cos(u) .* ŵ
            base = classical_elements(r, sqrt(μ / rmag) .* vdir)
            bumped = classical_elements(r, (sqrt(μ / rmag) * 1.0001) .* vdir)
            @test bumped.e > 1e-5
            @test bumped.arg_latitude_deg ≈ base.arg_latitude_deg atol = 1e-6
        end

        @testset "retrograde inclination exceeds 90 degrees" begin
            rmag = 7.0e6
            r = SVector{3, Float64}(rmag, 0.0, 0.0)
            v = SVector{3, Float64}(0.0, -sqrt(μ / rmag), 0.0)
            @test classical_elements(r, v).inclination_deg ≈ 180.0 rtol = 1e-10
        end
    end

    @testset "angle spread" begin
        @test angle_spread_deg([10.0, 20.0, 30.0]).span_deg ≈ 20.0
        # The wrap that a naive maximum-minus-minimum gets wrong.
        w = angle_spread_deg([359.0, 1.0])
        @test w.span_deg ≈ 2.0
        @test w.largest_gap_deg ≈ 358.0
        u = angle_spread_deg(collect(0.0:45.0:315.0))
        @test u.span_deg ≈ 315.0
        @test u.largest_gap_deg ≈ 45.0
        @test angle_spread_deg([12.0]).span_deg ≈ 0.0
        @test isnan(angle_spread_deg(Float64[]).span_deg)
        # Negative and over-360 inputs are wrapped, not rejected.
        @test angle_spread_deg([-1.0, 361.0]).span_deg ≈ 2.0
        @test wrap_360(-10.0) ≈ 350.0
    end

    # =======================================================================
    # Loader
    # =======================================================================
    @testset "loader" begin
        good = Dict(
            "epoch_utc" => "2025-06-06T00:00:00Z",
            "frame" => "J2000 Earth-centered inertial, meters and meters per second",
            "notes" => "n",
            "spacecraft" => [Dict(
                "name" => "CYGNSS FM01", "norad_id" => 41887,
                "r_ii_m" => [-4575510.123, 3629534.942, -3528444.131],
                "v_ii_m_s" => [-5440.1622, -5034.1701, 1863.5395],
                "provenance" => "telemetry", "source" => "s", "epoch_offset_s" => -1.0,
            )],
        )

        mktempdir() do dir
            write_doc(d, name) = (p = joinpath(dir, name); open(io -> JSON.print(io, d), p, "w"); p)

            @testset "round trip" begin
                ics = load_constellation_ics(write_doc(good, "good.json"))
                @test ics.epoch_utc == "2025-06-06T00:00:00Z"
                @test ics.notes == "n"
                @test length(ics.spacecraft) == 1
                s = ics.spacecraft[1]
                @test s.name == "CYGNSS FM01"
                @test s.norad_id == 41887
                @test s.provenance == "telemetry"
                @test s.epoch_offset_s ≈ -1.0
                @test s.r_ii_m ≈ SVector(-4575510.123, 3629534.942, -3528444.131)
                @test s.v_ii_m_s ≈ SVector(-5440.1622, -5034.1701, 1863.5395)
                @test s.r_ii_m isa SVector{3, Float64}
            end

            @testset "a missing norad id is allowed" begin
                d = deepcopy(good)
                delete!(d["spacecraft"][1], "norad_id")
                @test load_constellation_ics(write_doc(d, "nonorad.json")).spacecraft[1].norad_id === nothing
            end

            @testset "rejections" begin
                @test_throws ArgumentError load_constellation_ics(joinpath(dir, "absent.json"))

                for key in ("epoch_utc", "frame", "spacecraft")
                    d = deepcopy(good); delete!(d, key)
                    @test_throws ArgumentError load_constellation_ics(write_doc(d, "no_$key.json"))
                end
                for key in ("name", "r_ii_m", "v_ii_m_s", "provenance", "source", "epoch_offset_s")
                    d = deepcopy(good); delete!(d["spacecraft"][1], key)
                    @test_throws ArgumentError load_constellation_ics(write_doc(d, "sc_no_$key.json"))
                end

                d = deepcopy(good); d["spacecraft"] = []
                @test_throws ArgumentError load_constellation_ics(write_doc(d, "empty.json"))

                d = deepcopy(good); d["spacecraft"] = Dict("a" => 1)
                @test_throws ArgumentError load_constellation_ics(write_doc(d, "notarray.json"))

                # An unlabeled state is the failure the file exists to prevent.
                d = deepcopy(good); d["spacecraft"][1]["provenance"] = "guess"
                @test_throws ArgumentError load_constellation_ics(write_doc(d, "badprov.json"))

                d = deepcopy(good); d["spacecraft"][1]["r_ii_m"] = [1.0, 2.0]
                @test_throws ArgumentError load_constellation_ics(write_doc(d, "shortr.json"))
                d = deepcopy(good); d["spacecraft"][1]["v_ii_m_s"] = [1.0, 2.0, 3.0, 4.0]
                @test_throws ArgumentError load_constellation_ics(write_doc(d, "longv.json"))
            end
        end
    end

    # =======================================================================
    # The committed element-set snapshots
    # =======================================================================
    @testset "historical element sets" begin
        @test isfile(CYGIC_TLE_FILE)
        blocks = parse_tle_blocks(CYGIC_TLE_FILE)
        @test Set(keys(blocks)) == Set(["A", "B"])

        # NORAD ids CelesTrak's SATCAT gives for 2016-078; 41889 (FM06) is
        # absent from both captures because it had decayed by 2024-06-13.
        expected = Set([41884, 41885, 41886, 41887, 41888, 41890, 41891])
        for (name, recs) in blocks
            @test length(recs) == 7
            ids = Set{Int}()
            for rec in recs
                @test startswith(rec.name, "CYGFM")
                @test startswith(rec.line1, "1 ")
                @test startswith(rec.line2, "2 ")
                id1 = parse(Int, rec.line1[3:7])
                id2 = parse(Int, rec.line2[3:7])
                @test id1 == id2
                push!(ids, id1)
                # Modulo-10 checksum over each line, which catches a transcription
                # slip in the committed capture.
                for line in (rec.line1, rec.line2)
                    s = 0
                    for c in line[1:68]
                        isdigit(c) && (s += c - '0')
                        c == '-' && (s += 1)
                    end
                    @test s % 10 == parse(Int, line[69])
                end
                @test rec.block == name
            end
            @test ids == expected
            @test 41889 ∉ ids
        end

        @testset "malformed input is rejected" begin
            mktempdir() do dir
                p = joinpath(dir, "bad.tle")
                write(p, "BLOCK X cap\nCYGFM01\n1 41887U ...\nNOT A LINE TWO\n")
                @test_throws ArgumentError parse_tle_blocks(p)
                p2 = joinpath(dir, "short.tle")
                write(p2, "BLOCK X cap\nCYGFM01\n1 41887U ...\n")
                @test_throws ArgumentError parse_tle_blocks(p2)
                @test_throws ArgumentError parse_tle_blocks(joinpath(dir, "absent.tle"))
            end
        end
    end

    # =======================================================================
    # The generated artifact, when it is present.
    # data/telemetry/CYGNSS/ is gitignored, so this is skipped in a fresh
    # checkout and exercised on a machine that has run build_cygnss_ics.jl.
    # =======================================================================
    @testset "generated constellation file" begin
        if !isfile(CYGIC_ICS_FILE)
            @test_skip "constellation_ics_20250606.json absent (run scripts/dev/viewer_demos/build_cygnss_ics.jl)"
        else
            ics = load_constellation_ics(CYGIC_ICS_FILE)
            @test ics.epoch_utc == CYGNSS_EPOCH_UTC
            @test occursin("J2000", ics.frame)
            @test occursin("meters", ics.frame)

            # Seven, not eight: FM06 had decayed by this epoch.
            @test length(ics.spacecraft) == 7
            @test 41889 ∉ [s.norad_id for s in ics.spacecraft]
            @test !any(occursin("FM06", s.name) for s in ics.spacecraft)

            by_id = Dict(s.norad_id => s for s in ics.spacecraft)
            @test Set(keys(by_id)) == Set([41884, 41885, 41886, 41887, 41888, 41890, 41891])

            # FM01 and FM04 are the two with flight telemetry.
            @test by_id[41887].provenance == "telemetry"
            @test by_id[41885].provenance == "telemetry"
            @test all(by_id[i].provenance == "catalogue" for i in (41884, 41886, 41888, 41890, 41891))
            # Nothing may be unattributed, and nothing may be silently nominal.
            @test all(!isempty(s.source) for s in ics.spacecraft)

            for s in ics.spacecraft
                el = classical_elements(s.r_ii_m, s.v_ii_m_s)
                # The constellation as flown in June 2025, not the 510 km design
                # orbit: with no propulsion, drag had taken it to about 440 km.
                @test 400e3 < el.altitude_m < 500e3
                @test 34.5 < el.inclination_deg < 35.5
                @test 92 * 60 < el.period_s < 95 * 60
                @test el.e < 0.005
                @test norm(s.r_ii_m) > CYGNSS_R_EARTH_EQ
                # Near-circular, so speed is within a few m/s of circular speed.
                @test norm(s.v_ii_m_s) ≈ sqrt(CYGNSS_MU_EARTH / norm(s.r_ii_m)) rtol = 2e-3
            end

            # A catalog state was propagated across a real gap; a flight state
            # was not. This is what epoch_offset_s is for.
            @test abs(by_id[41887].epoch_offset_s) < 2.0
            @test abs(by_id[41885].epoch_offset_s) < 2.0
            @test all(abs(by_id[i].epoch_offset_s) > 3600.0 for i in (41884, 41886, 41888, 41890, 41891))

            # One plane, spread around it: that is the picture worth drawing.
            els = [classical_elements(s.r_ii_m, s.v_ii_m_s) for s in ics.spacecraft]
            @test angle_spread_deg([e.raan_deg for e in els]).span_deg < 60.0
            @test angle_spread_deg([e.arg_latitude_deg for e in els]).span_deg > 90.0
        end
    end

    # =======================================================================
    # The frame conversion against flight data, when both the telemetry mirror
    # and the SPICE kernels are present. The FM04 file carries an independently
    # produced inertial position, so the conversion can be checked rather than
    # trusted.
    # =======================================================================
    @testset "flight-data frame check" begin
        lsk = joinpath(CYGIC_SPICE_PATH, "lsk", "naif0012.tls")
        bpc = joinpath(CYGIC_SPICE_PATH, "pck", "earth_latest_high_prec.bpc")
        if !(isfile(CYGIC_FM04_FEATHER) && isfile(lsk) && isfile(bpc))
            @test_skip "FM04 telemetry or SPICE kernels absent"
        else
            let
                furnsh(lsk)
                furnsh(bpc)
                df = DataFrame(Arrow.Table(CYGIC_FM04_FEATHER))
                @test hasproperty(df, :pos_ii_1)
                worst = 0.0
                for k in 1:40_000:nrow(df)
                    M = sxform("ITRF93", "J2000", df.pvt_et_seconds[k])
                    l_pi = SMatrix{3, 3, Float64}(M[1:3, 1:3])'
                    dl_pi = SMatrix{3, 3, Float64}(M[4:6, 1:3])'
                    r_i, _ = earth_fixed_to_inertial(
                        l_pi, dl_pi,
                        SVector(df.sc_pos_x_pvt_m[k], df.sc_pos_y_pvt_m[k], df.sc_pos_z_pvt_m[k]),
                        SVector(0.0, 0.0, 0.0),
                    )
                    worst = max(worst, norm(r_i - SVector(df.pos_ii_1[k], df.pos_ii_2[k], df.pos_ii_3[k])))
                end
                # Sub-millimeter: the same ITRF93 -> J2000 transform that built
                # the file's own inertial columns.
                @test worst < 1e-3
            end
        end
    end
end
