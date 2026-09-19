using Test
using Arrow
using DataFrames
using JSON
using LinearAlgebra
using StaticArrays

const CYGTR_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

# `test/unit/runtests.jl` includes every unit test file into one shared scope, so
# `cygnss_ics_tests.jl` has already brought `CygnssICs` in by the time this file
# runs. Including it again would define a second module of the same name, and
# `using` both makes every shared export -- `classical_elements`,
# `CYGNSS_MU_EARTH` -- ambiguous in `Main` and undefined at the call site. The
# guard keeps this file runnable on its own as well as inside the suite.
isdefined(@__MODULE__, :CygnssICs) ||
    include(joinpath(CYGTR_REPO, "scripts", "dev", "viewer_demos", "cygnss_ics.jl"))
isdefined(@__MODULE__, :CygnssTracks) ||
    include(joinpath(CYGTR_REPO, "scripts", "dev", "viewer_demos", "cygnss_tracks.jl"))
using .CygnssICs
using .CygnssTracks

const CYGTR_TRACKS_FILE = joinpath(CYGTR_REPO, "data", "telemetry", "CYGNSS",
                                   "cygnss_constellation_tracks_20250606_96hr.feather")
const CYGTR_ICS_FILE = joinpath(CYGTR_REPO, "data", "telemetry", "CYGNSS",
                                "constellation_ics_20250606.json")
const CYGTR_COMPARISON_FILE = joinpath(CYGTR_REPO, "data", "telemetry", "CYGNSS",
                                       "catalog_vs_flown_20250606.json")

"A circular orbit in the x-y plane, tilted by `incl`, sampled at `t`."
function cygtr_circular(t::Real; radius = 6.82e6, mu = CYGNSS_MU_EARTH, incl = 0.0)
    n = sqrt(mu / radius^3)
    c, s = cos(n * t), sin(n * t)
    ci, si = cos(incl), sin(incl)
    r = SVector(radius * c, radius * s * ci, radius * s * si)
    v = SVector(-radius * n * s, radius * n * c * ci, radius * n * c * si)
    return (r, v)
end

"A small Feather track table, one circular orbit per spacecraft, each in its own plane."
function cygtr_write_table(path; names = ["A", "B"], times = 0.0:10.0:600.0)
    parts = DataFrame[]
    for (k, nm) in enumerate(names)
        rs = [cygtr_circular(t; incl = 0.2 * (k - 1)) for t in times]
        push!(parts, DataFrame(
            name = fill(nm, length(times)),
            norad_id = fill(Int32(40000 + k), length(times)),
            time = collect(float.(times)),
            pos_ii_1 = [p[1][1] for p in rs],
            pos_ii_2 = [p[1][2] for p in rs],
            pos_ii_3 = [p[1][3] for p in rs],
            vel_ii_1 = [p[2][1] for p in rs],
            vel_ii_2 = [p[2][2] for p in rs],
            vel_ii_3 = [p[2][3] for p in rs],
        ))
    end
    df = reduce(vcat, parts)
    Arrow.write(path, df)
    return df
end

@testset "CygnssTracks" begin
    @testset "reference decimation keeps coverage" begin
        sample = CygnssTracks.reference_sample_indices
        @test isempty(sample(Float64[], 10.0))
        @test isempty(sample([5.0, 6.0], 4.0))
        @test sample([5.0], 5.0) == [1]
        @test sample(0.0:10.0, 3.0; max_samples=2) == [1, 4]
        @test sample(0.0:10.0, 3.5; max_samples=3) == [1, 2, 4]
        @test_throws ArgumentError sample(0.0:10.0, 10.0; max_samples=1)
        @test_throws ArgumentError sample(0.0:10.0, Inf)
        for times in (collect(0.0:345599.0), [0.0, 1.0, 8.0, 100.0, 107.0])
            idx = sample(times, last(times); max_samples=3)
            @test length(idx) <= 3
            @test issorted(idx) && allunique(idx)
            @test times[first(idx)] == first(times)
            @test times[last(idx)] == last(times)
        end
    end


    # =======================================================================
    # Fitting a state out of a position arc
    # =======================================================================
    @testset "fit_state_from_arc" begin
        @testset "a polynomial inside the fit's own degree is recovered exactly" begin
            # Degree 3 in each axis, fitted at degree 5: the fit must reproduce
            # the generating polynomial's value and slope, not approximate them.
            coeffs = ([1.0e6, 2.0e3, -3.0, 0.4], [-5.0e5, 1.0e3, 7.0, -0.2], [2.0e5, -4.0e2, 1.5, 0.05])
            t = collect(-40.0:1.0:40.0)
            pos = [SVector(ntuple(c -> sum(coeffs[c][p + 1] * (ti - 7.0)^p for p in 0:3), 3)) for ti in t]
            r, v, rms = fit_state_from_arc(t, pos, 7.0; order = 5)
            @test rms < 1e-6
            for c in 1:3
                @test r[c] ≈ coeffs[c][1] rtol = 1e-9
                @test v[c] ≈ coeffs[c][2] rtol = 1e-7
            end
        end

        @testset "a circular orbit arc is recovered in position and velocity" begin
            t = collect(0.0:1.0:300.0)
            pos = [cygtr_circular(ti)[1] for ti in t]
            r0, v0 = cygtr_circular(0.0)
            r, v, rms = fit_state_from_arc(t, pos, 0.0; order = 5)
            # The epoch is the arc's left endpoint, exactly as in the real build.
            # A quintic over 300 s of a 93 minute orbit truncates at the
            # centimeter, two orders below the product's own 1 m quantization.
            @test norm(r - r0) < 0.1
            @test norm(v - v0) < 0.01
            @test rms < 0.05
        end

        @testset "a high-order fit over a long arc stays conditioned" begin
            # The regression test for scaling the fit variable. In raw seconds the
            # degree-9 Vandermonde over a 600 s arc spans nineteen decades and the
            # solve returns a state megameters from the data.
            t = collect(0.0:1.0:600.0)
            pos = [cygtr_circular(ti)[1] for ti in t]
            r0, v0 = cygtr_circular(0.0)
            for order in (5, 7, 9)
                r, v, rms = fit_state_from_arc(t, pos, 0.0; order = order)
                @test norm(r - r0) < 10.0      # meters, against a 6820 km radius
                @test norm(v - v0) < 0.2
                @test rms < 1.0
            end
        end

        @testset "a fit too poor to describe the arc says so in its residual" begin
            # A cubic cannot represent 38 degrees of orbit. The point is that the
            # residual reports the failure rather than the answer looking fine.
            t = collect(0.0:1.0:600.0)
            pos = [cygtr_circular(ti)[1] for ti in t]
            _, _, rms_cubic = fit_state_from_arc(t, pos, 0.0; order = 3)
            _, _, rms_quintic = fit_state_from_arc(t, pos, 0.0; order = 5)
            @test rms_cubic > 100.0
            @test rms_quintic < 1.0
        end

        @testset "malformed input is rejected" begin
            t = collect(0.0:1.0:10.0)
            pos = [cygtr_circular(ti)[1] for ti in t]
            @test_throws ArgumentError fit_state_from_arc(t[1:5], pos, 0.0)
            @test_throws ArgumentError fit_state_from_arc(t[1:4], pos[1:4], 0.0; order = 5)
            @test_throws ArgumentError fit_state_from_arc(zeros(9), pos[1:9], 0.0; order = 3)
        end
    end

    # =======================================================================
    # Comparing two states
    # =======================================================================
    @testset "rtn_offset" begin
        r = SVector(6.82e6, 0.0, 0.0)
        v = SVector(0.0, 7.6e3, 0.0)                 # h along +z

        @testset "each axis picks out its own offset" begin
            @test rtn_offset(r, v, SVector(500.0, 0.0, 0.0)).radial_m ≈ 500.0
            @test rtn_offset(r, v, SVector(0.0, 500.0, 0.0)).along_m ≈ 500.0
            @test rtn_offset(r, v, SVector(0.0, 0.0, 500.0)).cross_m ≈ 500.0
            @test abs(rtn_offset(r, v, SVector(500.0, 0.0, 0.0)).along_m) < 1e-9
            @test abs(rtn_offset(r, v, SVector(0.0, 500.0, 0.0)).cross_m) < 1e-9
        end

        @testset "the triad is orthonormal, so nothing leaks between axes" begin
            # An inclined, eccentric reference state and an arbitrary offset: the
            # three components must still add up in quadrature to the total.
            rr = SVector(4.1e6, -3.3e6, 2.7e6)
            vv = SVector(3.1e3, 6.2e3, -1.4e3)
            d = SVector(1234.5, -678.9, 246.8)
            o = rtn_offset(rr, vv, d)
            @test sqrt(o.radial_m^2 + o.along_m^2 + o.cross_m^2) ≈ norm(d) rtol = 1e-12
        end

        @testset "the sign convention is ahead-is-positive" begin
            # An offset along the direction of motion is a positive along-track.
            ahead = rtn_offset(r, v, normalize(v) * 1000.0)
            @test ahead.along_m > 0
        end
    end

    @testset "orbit_plane_angle_deg" begin
        r = SVector(6.82e6, 0.0, 0.0)
        v = SVector(0.0, 7.6e3, 0.0)
        # Identical states do not give exactly zero: acos of a dot product that
        # rounds to 1 returns the square root of the machine epsilon, about a
        # millionth of a degree. That is the function's resolution floor, and it
        # is four orders below the plane differences it is used to measure.
        @test orbit_plane_angle_deg(r, v, r, v) ≈ 0.0 atol = 1e-5
        # Tipping the velocity out of the plane by delta tips the plane by delta.
        # The angle is an arccosine of a near-unit dot product, so it carries
        # about half the available digits: at a thousandth of a degree the
        # absolute precision is around a millionth of a degree, which is what the
        # atol allows for. That is four orders below the angles being measured.
        for delta in (0.001, 0.05, 3.0)
            vd = SVector(0.0, 7.6e3 * cosd(delta), 7.6e3 * sind(delta))
            @test orbit_plane_angle_deg(r, v, r, vd) ≈ delta rtol = 1e-5 atol = 1e-6
        end
        # Scaling either state cannot change the angle between the planes.
        @test orbit_plane_angle_deg(r, v, 2r, 3v) ≈ 0.0 atol = 1e-5
    end

    @testset "along_track_time_s" begin
        # 7.6 km of along-track at 7.6 km/s is one second of timing error.
        @test along_track_time_s(7600.0, 7600.0) ≈ 1.0
        @test along_track_time_s(-7600.0, 7600.0) ≈ -1.0
    end

    # =======================================================================
    # The track table, against a synthetic file
    # =======================================================================
    @testset "track table" begin
        mktempdir() do dir
            path = joinpath(dir, "tracks.feather")
            cygtr_write_table(path)
            tracks = load_constellation_tracks(path)
            @test track_names(tracks) == ["A", "B"]

            @testset "a sample time returns that sample" begin
                for nm in ("A", "B"), t in (0.0, 120.0, 600.0)
                    r, v = track_state_at(tracks, nm, t)
                    row = findfirst(i -> tracks.table.name[i] == nm && tracks.table.time[i] == t,
                                    1:nrow(tracks.table))
                    @test r ≈ SVector(tracks.table.pos_ii_1[row], tracks.table.pos_ii_2[row],
                                      tracks.table.pos_ii_3[row]) rtol = 1e-12
                    @test v ≈ SVector(tracks.table.vel_ii_1[row], tracks.table.vel_ii_2[row],
                                      tracks.table.vel_ii_3[row]) rtol = 1e-12
                end
            end

            @testset "between samples it follows the orbit, not the chord" begin
                # Ten-second samples of a 93 minute orbit: the chord cuts about
                # 10 m off the arc at the midpoint, and Hermite must do far better.
                r, v = track_state_at(tracks, "A", 105.0)
                r_true, v_true = cygtr_circular(105.0)      # "A" is the untilted one
                a = findfirst(i -> tracks.table.name[i] == "A" && tracks.table.time[i] == 100.0,
                              1:nrow(tracks.table))
                chord = 0.5 .* (SVector(tracks.table.pos_ii_1[a], tracks.table.pos_ii_2[a],
                                        tracks.table.pos_ii_3[a]) .+
                                SVector(tracks.table.pos_ii_1[a + 1], tracks.table.pos_ii_2[a + 1],
                                        tracks.table.pos_ii_3[a + 1]))
                @test norm(chord - r_true) > 100.0
                @test norm(r - r_true) < 0.01
                @test norm(v - v_true) < 1e-4
            end

            @testset "a time outside the track is an error, not a clamp" begin
                @test_throws ArgumentError track_state_at(tracks, "A", -1.0)
                @test_throws ArgumentError track_state_at(tracks, "A", 601.0)
                @test_throws ArgumentError track_state_at(tracks, "Z", 10.0)
            end

            @testset "a malformed table is rejected" begin
                bad = joinpath(dir, "missing_column.feather")
                df = DataFrame(Arrow.Table(path))
                Arrow.write(bad, select(df, Not(:vel_ii_3)))
                @test_throws ArgumentError load_constellation_tracks(bad)

                unsorted = joinpath(dir, "unsorted.feather")
                shuffled = df[[2; 1; 3:nrow(df)], :]
                Arrow.write(unsorted, shuffled)
                @test_throws ArgumentError load_constellation_tracks(unsorted)

                @test_throws ArgumentError load_constellation_tracks(joinpath(dir, "absent.feather"))
            end
        end
    end

    # =======================================================================
    # The generated artifacts, when they are present. data/telemetry/CYGNSS/ is
    # gitignored, so this is skipped in a fresh checkout and exercised on a
    # machine that has run fetch_cygnss_l1_states.py and build_cygnss_ics.jl.
    # =======================================================================
    @testset "generated flight tracks" begin
        if !isfile(CYGTR_TRACKS_FILE)
            @test_skip "cygnss_constellation_tracks_20250606_96hr.feather absent " *
                       "(run scripts/dev/viewer_demos/build_cygnss_ics.jl)"
        else
            tracks = load_constellation_tracks(CYGTR_TRACKS_FILE)
            names = track_names(tracks)
            @test length(names) == 7
            @test !any(occursin("FM06", n) for n in names)

            for nm in names
                rows = tracks.table.name .== nm
                t = tracks.table.time[rows]
                @test first(t) ≈ 0.0 atol = 1e-6
                # The window is 96 hours; the last sample is at its far end.
                @test last(t) ≈ 96 * 3600.0 atol = 2.0
                # At 1 Hz over 96 hours, a few thousand dropped samples are the
                # product's own gaps; an order of magnitude more would mean the
                # fetch lost a granule.
                @test 330_000 < count(rows) <= 345_601

                r, v = track_state_at(tracks, nm, 0.0)
                el = classical_elements(r, v)
                @test 400e3 < el.altitude_m < 500e3
                @test 34.5 < el.inclination_deg < 35.5
                @test 92 * 60 < el.period_s < 95 * 60
            end

            if isfile(CYGTR_ICS_FILE)
                # The initial-conditions file and the track are two products of
                # one pipeline and must agree at the epoch. The state file is
                # fitted over a 300-second arc while the track carries the raw
                # 1 Hz samples, so they agree to the product's quantization and
                # not to the bit.
                ics = load_constellation_ics(CYGTR_ICS_FILE)
                @test length(ics.spacecraft) == 7
                @test all(s.provenance == "telemetry" for s in ics.spacecraft)
                for s in ics.spacecraft
                    s.name in names || continue
                    r, _ = track_state_at(tracks, s.name, 0.0)
                    @test norm(r - s.r_ii_m) < 25.0
                end
            end
        end
    end

    # =======================================================================
    # The fetcher's own decoder. It is Python, because the fetch is HTTP and
    # binary decoding, but it carries a network-free self-test and that test
    # belongs to this suite: the Level 1 states everything above rests on come
    # through it. Skipped where the interpreter or pyarrow is missing.
    # =======================================================================
    @testset "level 1 fetcher self-test" begin
        script = joinpath(CYGTR_REPO, "scripts", "dev", "viewer_demos",
                          "fetch_cygnss_l1_states.py")
        @test isfile(script)
        have_python = try
            success(pipeline(`python3 -c "import numpy, pyarrow"`; stdout = devnull, stderr = devnull))
        catch
            false
        end
        if !have_python
            @test_skip "python3 with numpy and pyarrow is not available"
        else
            out = IOBuffer()
            ok = success(pipeline(`python3 $script --self-test`; stdout = out, stderr = out))
            text = String(take!(out))
            ok || println(text)
            @test ok
            @test occursin("self-test: ok", text)
        end
    end

    @testset "catalog versus flown comparison" begin
        if !isfile(CYGTR_COMPARISON_FILE)
            @test_skip "catalog_vs_flown_20250606.json absent " *
                       "(run scripts/dev/viewer_demos/build_cygnss_ics.jl)"
        else
            doc = JSON.parsefile(CYGTR_COMPARISON_FILE)
            @test doc["epoch_utc"] == CYGNSS_EPOCH_UTC
            @test Set(doc["previously_catalog_fm"]) == Set([2, 3, 5, 7, 8])
            best = doc["captures"][doc["previous_file_capture"]]
            @test length(best["spacecraft"]) == 7
            for row in best["spacecraft"]
                # Each component must be consistent with the total it came from.
                total = sqrt(row["radial_m"]^2 + row["along_track_m"]^2 + row["cross_track_m"]^2)
                @test total ≈ row["position_error_m"] rtol = 1e-9
                # The finding this file records: the catalog error is a phase
                # error. The off-track part is a small fraction of the total and
                # the planes agree to well under a hundredth of a degree.
                off_track = sqrt(row["radial_m"]^2 + row["cross_track_m"]^2)
                @test off_track < 0.1 * row["position_error_m"]
                @test row["orbit_plane_deg"] < 0.01
                # Every element set in the chosen capture postdates the epoch and
                # was propagated backwards across days, not seconds.
                @test row["propagation_s"] < -3600.0
            end
        end
    end
end
