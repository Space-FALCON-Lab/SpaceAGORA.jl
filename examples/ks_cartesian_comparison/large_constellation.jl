# Large-constellation nonlinear KS versus Cartesian propagation (point mass + J2).
# Reuses comparison.jl's KS route. Errors are terminal (all spacecraft share the
# target clock, so the KS terminal event lands every spacecraft on t_f); dense
# sampling is skipped because it does not fit in memory at this size.
#   julia --project=examples/ks_cartesian_comparison examples/ks_cartesian_comparison/large_constellation.jl [planes slots orbits]
include(joinpath(@__DIR__, "comparison.jl"))
using .NativeCoordinateComparison
const M = NativeCoordinateComparison
using OrdinaryDiffEq, LinearAlgebra, CSV, DataFrames
using OrdinaryDiffEqLowOrderRK: DP5
const SVector = M.CS.SVector

"""Walker-delta constellation at the original shell-0 a, inclination and argument of periapsis."""
function walker_initial(planes, slots, e)
    base = M.walker27_constellation()
    a0, _, inc, _, arg = base.initial_elements[1, 5:9]
    a = a0 * 0.99 / (1 - e)  # same perigee-preserving rule as ten_initial
    n = planes * slots; x = zeros(6, n); planet = base.config.planet
    for P in 0:planes-1, S in 0:slots-1
        j = P * slots + S + 1
        node = 2pi * P / planes; Mean = 2pi * S / slots + 2pi * P / n  # phasing F = 1
        E = Mean
        for _ in 1:15; E -= (E - e * sin(E) - Mean) / (1 - e * cos(E)); end
        nu = 2atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
        r, v = M.FT.orbitalelemtorv(SVector(a, e, inc, node, arg, nu, 0.0), planet)
        x[:, j] .= vcat(r, v)
    end
    (; x, planet, a, e, period=2pi * sqrt(a^3 / planet.μ))
end

"""Same physics as M.cart_rhs!, without heap temporaries."""
function cart_rhs_static!(du, u, p, t)
    for j in 1:length(u)÷6
        k = 6j - 6; r = SVector(u[k+1], u[k+2], u[k+3])
        du[k+1] = u[k+4]; du[k+2] = u[k+5]; du[k+3] = u[k+6]
        acc = -p.μ * r / norm(r)^3 + M.KS.ks_j2_acceleration_si(r, p)
        du[k+4] = acc[1]; du[k+5] = acc[2]; du[k+6] = acc[3]
    end
end

"""Terminal Cartesian state (6 x n) and accepted step count."""
function propagate(c, route, tf, alg; tol=1e-10, step=0.0)
    common = (; reltol=tol, abstol=tol, adaptive=step == 0.0, maxiters=10^9, save_everystep=false, save_start=false)
    if route in (:cartesian, :cartesian_alloc)
        rhs = route == :cartesian ? cart_rhs_static! : M.cart_rhs!
        sol = solve(ODEProblem(rhs, vec(c.x), (0.0, tf), c.planet), alg; dt=step == 0.0 ? 1.0 : step, common...)
        X = reshape(sol.u[end], 6, :)
    else
        y = M.to_ks(c.x, c.planet)
        callback = ContinuousCallback((u, s, integrator) -> u[10] - tf, terminate!;
            abstol=1e-9, reltol=0.0, save_positions=(false, true))
        sol = solve(ODEProblem(M.ks_rhs!, y, (0.0, 2tf / (c.a * (1 - c.e))), c.planet), alg;
            dt=(step == 0.0 ? 1.0 : step) / c.a, callback, common...)
        X = M.to_cart(sol.u[end])
    end
    string(sol.retcode) in ("Success", "Terminated") || error("Propagation failed: $(sol.retcode)")
    clock = route == :ks_target_s ? abs(sol.u[end][10] - tf) : 0.0
    (; X, steps=sol.stats.naccept, rhs_calls=sol.stats.nf, clock)
end

function terminal_errors(X, R)
    dp = [norm(X[1:3, j] - R[1:3, j]) for j in axes(X, 2)]
    dv = [norm(X[4:6, j] - R[4:6, j]) for j in axes(X, 2)]
    (; max_position_m=maximum(dp), rms_position_m=sqrt(sum(abs2, dp) / length(dp)),
       max_velocity_m_s=maximum(dv), rms_velocity_m_s=sqrt(sum(abs2, dv) / length(dv)))
end

function large_runs(; planes=40, slots=25, orbits=50, steps=(5.0, 20.0, 100.0), tols=(1e-8, 1e-10, 1e-12),
        out=joinpath(M.OUTPUT, "large_constellation_$(planes*slots)sat_$(orbits)orb"))
    mkpath(out); rows = NamedTuple[]; truth_rows = NamedTuple[]
    for e in (0.01, 0.3)
        c = walker_initial(planes, slots, e); tf = orbits * c.period
        truth = propagate(c, :cartesian, tf, Vern9(); tol=2e-14).X
        refined = propagate(c, :cartesian, tf, Vern9(); tol=1e-14).X
        push!(truth_rows, (; e, nsat=planes * slots, orbits, terminal_errors(truth, refined)...))
        CSV.write(joinpath(out, "truth_refinement.csv"), DataFrame(truth_rows))
        for route in (:cartesian, :cartesian_alloc, :ks_target_s), (name, alg) in (("Tsit5", Tsit5()), ("DP5", DP5())),
            (mode, value) in vcat([("adaptive", t) for t in tols], [("fixed", s) for s in steps])
            tol = mode == "adaptive" ? value : 1e-12; step = mode == "fixed" ? value : 0.0
            propagate(c, route, min(tf, 100.0), alg; tol, step)  # warm compilation
            GC.gc(); stats = @timed propagate(c, route, tf, alg; tol, step)
            r = stats.value
            push!(rows, (; e, nsat=planes * slots, orbits, route=string(route), integrator=name, mode, setting=value,
                solve_s=stats.time, gc_s=stats.gctime, bytes=stats.bytes, steps=r.steps, rhs_calls=r.rhs_calls,
                clock_error_s=r.clock, terminal_errors(r.X, refined)...))
            CSV.write(joinpath(out, "nonlinear.csv"), DataFrame(rows))
            println("e=$e $route $name $mode $value: $(round(stats.time; digits=3)) s, ",
                "max terminal position error $(round(rows[end].max_position_m; sigdigits=3)) m"); flush(stdout)
        end
    end
    println("Results: $(abspath(out))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    M.BLAS.set_num_threads(1)
    args = parse.(Int, ARGS)
    isempty(args) ? large_runs() : large_runs(; planes=args[1], slots=args[2], orbits=args[3])
end
