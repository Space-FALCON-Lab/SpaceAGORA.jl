using Test
using Dates
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const CYGSL_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(CYGSL_REPO, "scripts", "dev", "viewer_demos", "cygnss_slew_telemetry.jl"))
using .CygnssSlewTelemetry

const CYGSL_SM = SpaceAGORA.SimulationModel
const CYGSL_ADCS = joinpath(CYGSL_REPO, "data", "telemetry", "CYGNSS", "cyg01_slew_adcs.feather")
const CYGSL_PV = joinpath(CYGSL_REPO, "data", "telemetry", "CYGNSS", "cyg01_slew_pv_eci.feather")
const CYGSL_CONSTANTS = joinpath(CYGSL_REPO, "data", "telemetry", "CYGNSS", "cyg01_adcs_constants.toml")
# The leap-second kernel the epoch conversion needs: the starter pack, or the
# path the suite already honours. Absent means the SPICE-backed epoch tests skip.
const CYGSL_LSK = joinpath(get(ENV, "SPACEAGORA_SPICE_PATH",
    joinpath(CYGSL_REPO, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")), "lsk", "naif0012.tls")
const CYGSL_LSK_READY = isfile(CYGSL_LSK) && (lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
    CygnssSlewTelemetry.SPICE.furnsh(CYGSL_LSK)
    true
end)

"A unit quaternion from an axis and an angle, scalar-first, for the convention tests."
cygsl_axis_angle(axis, θ) = (n = normalize(SVector{3, Float64}(axis));
    SVector{4, Float64}(cos(θ / 2), (sin(θ / 2) .* n)...))

"""
The scalar-last quaternion of an inertial-to-body matrix, in SpaceAGORA's
convention: `A = (w^2 - |q|^2) I - 2 w [q x] + 2 q q'`, so the antisymmetric
part gives the vector components as `(A[2,3] - A[3,2]) / 4w` and its cyclic
partners. Only used on matrices whose trace is comfortably above -1.
"""
function cygsl_quaternion_of(a::AbstractMatrix)
    w = sqrt(max(0.0, 1 + a[1, 1] + a[2, 2] + a[3, 3])) / 2
    w > 1e-6 || error("cygsl_quaternion_of: the trace path is degenerate for this matrix")
    return SVector{4, Float64}((a[2, 3] - a[3, 2]) / (4w), (a[3, 1] - a[1, 3]) / (4w),
        (a[1, 2] - a[2, 1]) / (4w), w)
end

@testset "CygnssSlewTelemetry" begin

    # =======================================================================
    # The epoch
    # =======================================================================
    @testset "epoch conversion" begin
        @testset "the counter is ephemeris time past J2000" begin
            @test SLEW_TIME_ORIGIN == DateTime(2000, 1, 1, 12, 0, 0)
            @test slew_epoch_et(0.0) == 0.0
            @test slew_epoch_et(8.12811487477591e8) == 8.12811487477591e8
            @test slew_epoch_et(Float32(1.5)) === 1.5
        end

        if !CYGSL_LSK_READY
            @test_skip "leap-second kernel absent; the SPICE-backed epoch conversion is exercised where the starter pack is staged"
        else
            @testset "the FM01 export's first sample is 2025-10-04T00:56:58.295 UTC" begin
                # The export's own first time stamp, converted as ET through
                # the leap-second kernel. The GPS receiver's stamp in the same
                # row is 00:56:58.000, the age of the fix earlier.
                stamp = slew_epoch_utc(8.12811487477591e8)
                @test stamp == DateTime(2025, 10, 4, 0, 56, 58, 295)
                # The superseded reading, UTC seconds added to the calendar
                # origin, sat TT minus UTC later: 69.18 s in 2025.
                old = SLEW_TIME_ORIGIN + Millisecond(round(Int, 8.12811487477591e8 * 1000))
                @test Millisecond(69_180) <= old - stamp <= Millisecond(69_190)
                # And the reading from midnight of that day, twelve hours
                # earlier still, is what the SGP4 orbit-plane check rejects.
                midnight = DateTime(2000, 1, 1, 0, 0, 0) + Millisecond(round(Int, 8.12811487477591e8 * 1000))
                @test Date(midnight) == Date(2025, 10, 3)
            end

            @testset "the origin and sub-second resolution" begin
                # ET zero is 2000-01-01T12:00:00 TDB, 11:58:55.816 UTC.
                @test slew_epoch_utc(0.0) == DateTime(2000, 1, 1, 11, 58, 55, 816)
                @test slew_epoch_utc(1.5) == DateTime(2000, 1, 1, 11, 58, 57, 316)
                @test slew_epoch_utc(86_400.0) == DateTime(2000, 1, 2, 11, 58, 55, 816)
            end
        end
    end

    # =======================================================================
    # The attitude convention
    # =======================================================================
    @testset "quaternion convention" begin
        @testset "the mapping moves the scalar last and negates the vector" begin
            q = SVector{4, Float64}(0.5, 0.5, -0.5, 0.5)
            @test slew_quaternion_to_spaceagora(q) == SVector{4, Float64}(-0.5, 0.5, -0.5, 0.5)
        end

        @testset "rot of the mapped quaternion is the transpose of the scalar-first matrix" begin
            # This is the identity the mapping rests on: negating the vector
            # part conjugates the quaternion, and conjugating transposes the
            # attitude matrix. The nadir check in the loader measures which of
            # the two directions the export actually uses; this test pins that
            # the mapping delivers the other one, so the two statements together
            # fix the frame.
            for (axis, θ) in ((SVector(1.0, 0, 0), 0.4), (SVector(0.0, 1, 0), -1.1),
                              (SVector(1.0, -2.0, 0.5), 2.3), (SVector(0.3, 0.2, -0.9), 0.05))
                q = cygsl_axis_angle(axis, θ)
                mapped = slew_quaternion_to_spaceagora(q)
                @test norm(mapped) ≈ 1.0 atol = 1e-15
                @test CYGSL_SM.rot(mapped) ≈ slew_scalar_first_attitude_matrix(q)' atol = 1e-13
            end
        end

        @testset "the scalar-first matrix is a rotation about the stated axis" begin
            θ = 0.7
            q = cygsl_axis_angle(SVector(0.0, 0.0, 1.0), θ)
            a = slew_scalar_first_attitude_matrix(q)
            @test a ≈ SMatrix{3, 3, Float64}(cos(θ), -sin(θ), 0, sin(θ), cos(θ), 0, 0, 0, 1) atol = 1e-13
            @test det(a) ≈ 1.0 atol = 1e-13
            @test a * a' ≈ SMatrix{3, 3, Float64}(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0) atol = 1e-13
        end

        @testset "the mapping is its own inverse on the quaternion's own frame" begin
            q = cygsl_axis_angle(SVector(0.2, -0.7, 0.4), 1.9)
            mapped = slew_quaternion_to_spaceagora(q)
            # Mapping back: scalar first again, vector negated again.
            back = SVector{4, Float64}(mapped[4], -mapped[1], -mapped[2], -mapped[3])
            @test back ≈ q atol = 1e-15
        end
    end

    # =======================================================================
    # The body-rate and wheel-axis signs
    # =======================================================================
    @testset "body rate and wheel axes" begin
        @testset "the body rate is w_eci unchanged" begin
            w = SVector{3, Float64}(-1.0e-4, -1.2e-3, 7.0e-5)
            @test slew_body_rate(w) == w
        end

        @testset "the constants file's wheel-axis matrix enters negated" begin
            axes = [0.9 0.1 -0.2; -0.1 1.0 0.3; 0.2 -0.3 0.95]
            @test slew_wheel_axes(axes) == -axes
        end

        @testset "the two sign choices are not interchangeable in the torque law" begin
            # Negating both the body rate and the wheel momentum leaves the
            # algebraic conservation law `I w + H = const` untouched, which is
            # why the regression that produced the file's matrix could not see
            # which sign it had. The forward torque law is not invariant, and
            # this is the test that says so.
            inertia = SMatrix{3, 3, Float64}(1.4, 0, 0, 0, 0.82, 0, 0, 0, 1.95)
            ω = SVector{3, Float64}(1.0e-3, -5.0e-4, 2.0e-4)
            h = SVector{3, Float64}(1.2e-3, -4.0e-4, 8.0e-4)
            hdot = SVector{3, Float64}(2.0e-6, 1.0e-6, -3.0e-6)
            @test inertia * ω + h ≈ -(inertia * (-ω) + (-h)) atol = 1e-18
            τ_plus = CYGSL_SM.wheel_reaction_torque(h, hdot, ω)
            τ_minus = CYGSL_SM.wheel_reaction_torque(-h, -hdot, -ω)
            @test !isapprox(τ_plus, -τ_minus; atol=1e-12)
        end
    end

    # =======================================================================
    # The LVLH pointing angle
    # =======================================================================
    @testset "LVLH pointing angle" begin
        # A circular orbit in the equatorial plane: the spacecraft is on +x
        # moving toward +y, so nadir is -x and the orbit normal is +z. The
        # on-target body frame then has +z on -x, +y on -z and +x on +y.
        r = SVector{3, Float64}(7.0e6, 0.0, 0.0)
        v = SVector{3, Float64}(0.0, 7.5e3, 0.0)
        z_body = SVector{3, Float64}(-1.0, 0.0, 0.0)
        y_body = SVector{3, Float64}(0.0, 0.0, -1.0)
        x_body = cross(y_body, z_body)
        # Rows of the inertial-to-body matrix are the body axes in inertial.
        c_on_target = SMatrix{3, 3, Float64}(hcat(x_body, y_body, z_body)')
        @test x_body ≈ SVector{3, Float64}(0.0, 1.0, 0.0) atol = 1e-15

        q_on_target = cygsl_quaternion_of(c_on_target)

        @testset "the on-target attitude reads zero" begin
            @test CYGSL_SM.rot(q_on_target) ≈ c_on_target atol = 1e-12
            @test lvlh_pointing_angle_deg(q_on_target, r, v) < 1e-9
        end

        @testset "a known extra rotation reads its own angle" begin
            for θ_deg in (0.5, 10.0, 37.5)
                θ = deg2rad(θ_deg)
                # Rotate the body frame a further θ about its own x axis.
                extra = SMatrix{3, 3, Float64}(1.0, 0, 0, 0, cos(θ), -sin(θ), 0, sin(θ), cos(θ))
                q = cygsl_quaternion_of(extra * c_on_target)
                @test lvlh_pointing_angle_deg(q, r, v) ≈ θ_deg atol = 1e-8
            end
        end

        @testset "the sign of the quaternion does not matter" begin
            @test lvlh_pointing_angle_deg(-q_on_target, r, v) ≈ lvlh_pointing_angle_deg(q_on_target, r, v) atol = 1e-12
        end
    end

    @testset "attitude angle" begin
        q = cygsl_axis_angle(SVector(0.0, 0.0, 1.0), 0.0)
        for θ_deg in (0.0, 1.0, 45.0, 179.0)
            θ = deg2rad(θ_deg)
            p = SVector{4, Float64}(sin(θ / 2), 0.0, 0.0, cos(θ / 2))
            @test attitude_angle_deg(SVector{4, Float64}(0, 0, 0, 1), p) ≈ θ_deg atol = 1e-8
            @test attitude_angle_deg(SVector{4, Float64}(0, 0, 0, 1), -p) ≈ θ_deg atol = 1e-8
        end
        @test attitude_angle_deg(SVector{4, Float64}(0, 0, 0, 2.0), SVector{4, Float64}(0, 0, 0, 0.5)) < 1e-12
    end

    # =======================================================================
    # The loaders, when the gitignored telemetry is present
    # =======================================================================
    @testset "loaders" begin
        if !isfile(CYGSL_ADCS) || !isfile(CYGSL_PV)
            @test_skip "cyg01_slew_adcs.feather / cyg01_slew_pv_eci.feather absent " *
                "(gitignored); the FM01 export loader is exercised only where the telemetry is staged"
        else
            tel = load_slew_telemetry(CYGSL_ADCS, CYGSL_PV)
            @test length(tel.t_rel) == size(tel.q, 2) == size(tel.omega, 2)
            @test issorted(tel.t_rel)
            @test tel.nadir_spread < 0.15
            # The state channel is 1 Hz inside a 4 Hz file, so the distinct
            # fixes are about a quarter of the rows and never more than them.
            @test length(tel.pv_t_rel) == size(tel.pos_m, 2) == size(tel.vel_mps, 2)
            @test length(tel.pv_t_rel) < length(tel.t_rel)
            @test length(tel.pv_t_rel) > length(tel.t_rel) / 8
            @test all(k -> norm(tel.q[:, k]) ≈ 1.0, 1:100:size(tel.q, 2))
            # The quaternion is unwrapped: no sign flip between neighbors.
            @test all(k -> dot(tel.q[:, k], tel.q[:, k - 1]) > 0.0, 2:size(tel.q, 2))
            @test tel.t_abs[1] == 8.12811487477591e8
            if isfile(CYGSL_CONSTANTS)
                L = momentum_ledger(tel, load_slew_constants(CYGSL_CONSTANTS); window=(890.0, 1100.0))
                @test all(isfinite, L.drift_nm) && all(isfinite, L.gravity_gradient_nm)
                @test L.exchanged_nms > 0.0 && L.omitted_nms >= 0.0 && L.ratio == L.omitted_nms / L.exchanged_nms
                @test L.rod_duty === nothing && L.fixes >= 3
                @test_throws ArgumentError momentum_ledger(tel, load_slew_constants(CYGSL_CONSTANTS); window=(1100.0, 890.0))
            end
            CYGSL_LSK_READY && @test Date(slew_epoch_utc(tel.t_abs[1])) == Date(2025, 10, 4)
        end

        if !isfile(CYGSL_CONSTANTS)
            @test_skip "cyg01_adcs_constants.toml absent (gitignored); the ADCS constants loader is " *
                "exercised only where the file is staged"
        else
            c = load_slew_constants(CYGSL_CONSTANTS)
            @test c.inertia ≈ c.inertia' atol = 1e-12
            @test c.wheel_inertia > 0.0
            @test c.wheel_axes == -c.wheel_axes_from_file .* permutedims(c.effective_inertia_scale)
            @test all(isfinite, c.effective_inertia_scale) && all(>(0.0), c.effective_inertia_scale)
            # The file's inner lists are the matrix columns; if they were read
            # as rows the column norms would be nowhere near unity.
            @test all(0.8 .< [norm(c.wheel_axes_from_file[:, k]) for k in 1:3] .< 1.25)
            # inertia = momentum rating / speed rating, the file's own relation.
            @test c.wheel_inertia ≈ c.momentum_rating_nms / (c.speed_rating_rpm * 2pi / 60) rtol = 1e-3
        end
    end
end
