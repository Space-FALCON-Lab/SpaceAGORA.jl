include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyConfig.jl"))
using .StudyConfig
import Pkg
Pkg.activate(LIBRARY_ROOT; io=devnull)
using Arrow, DataFrames, LinearAlgebra, StaticArrays, Test, TOML

vec(table, sc, field, row) = SVector{3,Float64}(ntuple(
    k -> table[Symbol("sc$(sc)_$(field)_$k")][row], 3))

function interpolate(r0, v0, r1, v1, dt, s)
    r = (2s^3 - 3s^2 + 1) * r0 + (s^3 - 2s^2 + s) * dt * v0 +
        (-2s^3 + 3s^2) * r1 + (s^3 - s^2) * dt * v1
    v = ((6s^2 - 6s) * r0 + (-6s^2 + 6s) * r1) / dt +
        (3s^2 - 4s + 1) * v0 + (3s^2 - 2s) * v1
    r, v
end

function verify_saved_diagnostics(directory)
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    info["settings"]["schema_version"] == 1 || error(
        "This historical checker requires the dense 10 s output. For compact output, use production_pipeline_tests.jl verify-reference.")
    target = info["settings"]["target_index"]
    mass = info["settings"]["mass_kg"]
    force_size = info["settings"]["force_N"]
    laser = sort(DataFrame(Arrow.Table(joinpath(directory, "laser_firing.feather"))), :start_time_s)
    totals = zeros(5)
    work = zeros(target)
    previous_table = nothing
    previous_row = 0
    previous_time = 0.0
    cursor = 1
    max_rtn_error, max_signed_error, max_work_error = 0.0, 0.0, 0.0
    rows = 0
    for table in Arrow.Stream(joinpath(directory, "trajectory.feather"))
        work_columns = [table[Symbol("sc$(sc)_laser_work_J")] for sc in 1:target]
        for row in eachindex(table.time)
            time = table.time[row]
            if previous_table !== nothing
                dt = time - previous_time
                while cursor <= nrow(laser) && laser.end_time_s[cursor] <= previous_time
                    cursor += 1
                end
                event = cursor
                while event <= nrow(laser) && laser.start_time_s[event] < time
                    helper = laser.sc_a[event]
                    start = max(previous_time, laser.start_time_s[event])
                    stop = min(time, laser.end_time_s[event])
                    rt0, vt0 = vec(previous_table, target, "pos", previous_row), vec(previous_table, target, "vel", previous_row)
                    rt1, vt1 = vec(table, target, "pos", row), vec(table, target, "vel", row)
                    rh0, vh0 = vec(previous_table, helper, "pos", previous_row), vec(previous_table, helper, "vel", previous_row)
                    rh1, vh1 = vec(table, helper, "pos", row), vec(table, helper, "vel", row)
                    rate = function (t)
                        fraction = (t - previous_time) / dt
                        r, v = interpolate(rt0, vt0, rt1, vt1, dt, fraction)
                        rh, vh = interpolate(rh0, vh0, rh1, vh1, dt, fraction)
                        force = force_size * normalize(r - rh)
                        radial = normalize(r)
                        normal = normalize(cross(r, v))
                        transverse = cross(normal, radial)
                        R, T, N = dot(force, radial) / mass, dot(force, transverse) / mass, dot(force, normal) / mass
                        SVector(R, T, N, max(T, 0.0), min(T, 0.0), dot(force, v), -dot(force, vh))
                    end
                    # Split at T sign changes so the signed integrands are smooth.
                    knots = collect(range(start, stop; length=33))
                    boundaries = [start, stop]
                    for index in 1:length(knots)-1
                        lo, hi = knots[index], knots[index + 1]
                        a, b = rate(lo)[2], rate(hi)[2]
                        iszero(a) && push!(boundaries, lo)
                        if signbit(a) != signbit(b)
                            while hi - lo > 1e-7
                                middle = (lo + hi) / 2
                                value = rate(middle)[2]
                                if signbit(value) == signbit(a)
                                    lo, a = middle, value
                                else
                                    hi = middle
                                end
                            end
                            push!(boundaries, (lo + hi) / 2)
                        end
                    end
                    sort!(unique!(boundaries))
                    integral = zero(SVector{7,Float64})
                    pieces = 32
                    for index in 1:length(boundaries)-1
                        lo, hi = boundaries[index], boundaries[index + 1]
                        part = zero(integral)
                        for sample in 0:pieces
                            t = lo + (hi - lo) * sample / pieces
                            weight = sample in (0, pieces) ? 1 : isodd(sample) ? 4 : 2
                            part += weight * rate(t)
                        end
                        integral += part * ((hi - lo) / (3pieces))
                    end
                    totals .+= integral[1:5]
                    work[target] += integral[6]
                    work[helper] += integral[7]
                    event += 1
                end
            end
            measured = SVector(table.target_dv_R_mps[row], table.target_dv_T_mps[row], table.target_dv_N_mps[row])
            max_rtn_error = max(max_rtn_error, norm(measured - totals[1:3]))
            max_signed_error = max(max_signed_error,
                abs(table.target_dv_T_positive_mps[row] - totals[4]),
                abs(table.target_dv_T_negative_mps[row] - totals[5]))
            for sc in 1:target
                max_work_error = max(max_work_error, abs(work_columns[sc][row] - work[sc]))
            end
            previous_table, previous_row, previous_time = table, row, time
            rows += 1
        end
    end
    report = Dict{String,Any}("output_directory" => directory, "rows" => rows,
        "method" => "saved 10 s Cartesian Hermite interpolation, actual firing and T-zero boundaries, 32-subinterval Simpson quadrature on each smooth segment",
        "max_RTN_error_mps" => max_rtn_error, "max_signed_T_error_mps" => max_signed_error,
        "max_work_error_J" => max_work_error, "RTN_tolerance_mps" => 1e-8,
        "signed_T_tolerance_mps" => 1e-7, "work_tolerance_J" => 0.1 + 1e-8maximum(abs, work))
    report["passed"] = max_rtn_error <= report["RTN_tolerance_mps"] &&
        max_signed_error <= report["signed_T_tolerance_mps"] &&
        max_work_error <= report["work_tolerance_J"]
    root = joinpath(@__DIR__, "results", "production_pipeline", "saved_diagnostics")
    mkpath(root)
    atomic_toml(joinpath(root, info["run_id"] * ".toml"), report)
    TOML.print(stdout, report; sorted=true)
    @testset "Fresh-file independent RTN/work reconstruction" begin
        @test max_rtn_error <= report["RTN_tolerance_mps"]
        @test max_signed_error <= report["signed_T_tolerance_mps"]
        @test max_work_error <= report["work_tolerance_J"]
    end
end

length(ARGS) == 1 || error("Usage: julia verify_saved_production_diagnostics.jl CASE_DIRECTORY")
verify_saved_diagnostics(abspath(only(ARGS)))
