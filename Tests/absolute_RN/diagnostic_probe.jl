function independent_rates(c, state, pairs)
    count = c.helpers + 1
    result = zeros(12 + count)
    isempty(pairs) && return result
    helper, target = only(pairs)
    target == count || error("Independent probe received the wrong target index")
    r, v = SVector{3}(state.sc[target].pos), SVector{3}(state.sc[target].vel)
    displacement = r - SVector{3}(state.sc[helper].pos)
    force = (c.power_W * B / LIGHT_SPEED) * displacement / norm(displacement)
    acceleration = force / state.sc[target].mass
    radial = r / norm(r)
    normal = normalize(cross(r, v))
    transverse = cross(normal, radial)
    result[1:3] .= (dot(acceleration, radial), dot(acceleration, transverse), dot(acceleration, normal))
    for axis in 1:3
        result[3 + axis] = max(result[axis], 0.0)
        result[6 + axis] = min(result[axis], 0.0)
        result[9 + axis] = abs(result[axis])
    end
    result[12 + target] = dot(force, v)
    result[12 + helper] = -dot(force, state.sc[helper].vel)
    result
end

function split_integral(rate, start, stop)
    grid = collect(range(start, stop; length=9))
    samples = rate.(grid)
    cuts = [start, stop]
    for axis in 1:3, index in 1:length(grid)-1
        left, right = grid[index], grid[index + 1]
        fl, fr = samples[index][axis], samples[index + 1][axis]
        iszero(fl) && push!(cuts, left)
        iszero(fr) && push!(cuts, right)
        fl * fr < 0 || continue
        for _ in 1:48
            middle = (left + right) / 2
            (middle == left || middle == right) && break
            fm = rate(middle)[axis]
            if iszero(fm)
                left = right = middle
                break
            elseif signbit(fm) == signbit(fl)
                left, fl = middle, fm
            else
                right = middle
            end
        end
        push!(cuts, (left + right) / 2)
    end
    sort!(unique!(cuts))
    fine, coarse = zeros(length(first(samples))), zeros(length(first(samples)))
    for index in 1:length(cuts)-1
        left, right = cuts[index], cuts[index + 1]
        width = right - left
        for k in 0:16
            values = rate(left + width * (k / 16))
            fine_weight = k in (0, 16) ? 1 : isodd(k) ? 4 : 2
            fine .+= (width / 48 * fine_weight) .* values
            if iseven(k)
                coarse_weight = k in (0, 16) ? 1 : isodd(k ÷ 2) ? 4 : 2
                coarse .+= (width / 24 * coarse_weight) .* values
            end
        end
    end
    (; fine, difference=abs.(fine - coarse))
end

function diagnostic_probe(c)
    count = c.helpers + 1
    totals = zeros(12 + count)
    net_error, branch_error, absolute_error = zeros(3), zeros(3), zeros(3)
    quadrature_difference = zeros(12)
    work_error = 0.0
    steps, next_sample = 0, 2
    times = trajectory_times(c)
    names = vcat([:time], [Symbol("target_dv_$(component)_mps") for component in ("R", "T", "N")],
        [Symbol("target_dv_$(component)_$(kind)_mps")
            for kind in ("positive", "negative", "absolute") for component in ("R", "T", "N")],
        [Symbol("sc$(sc)_laser_work_J") for sc in 1:count])
    reference = DataFrame([name => [0.0] for name in names])

    function observe(integrator, start, pairs)
        stop = Float64(integrator.t)
        function integral_to(endpoint)
            isempty(pairs) && return (fine=zeros(length(totals)), difference=zeros(length(totals)))
            split_integral(t -> independent_rates(c, integrator(t), pairs), start, endpoint)
        end
        while next_sample <= length(times) && times[next_sample] <= stop
            sample_time = times[next_sample]
            values = totals + integral_to(sample_time).fine
            push!(reference, (sample_time, values...))
            next_sample += 1
        end
        integral = integral_to(stop)
        totals .+= integral.fine
        quadrature_difference .+= integral.difference[1:12]
        target = integrator.u.sc[count]
        for (axis, component) in enumerate(("R", "T", "N"))
            positive = getproperty(target, Symbol("laser_dv_$(component)_positive"))
            negative = getproperty(target, Symbol("laser_dv_$(component)_negative"))
            net_error[axis] = max(net_error[axis], abs(target.laser_rtn[axis] - totals[axis]))
            branch_error[axis] = max(branch_error[axis], abs(positive - totals[3 + axis]),
                abs(negative - totals[6 + axis]))
            absolute_error[axis] = max(absolute_error[axis],
                abs(positive - negative - totals[9 + axis]))
        end
        for sc in 1:count
            work_error = max(work_error, abs(integrator.u.sc[sc].laser_work - totals[12 + sc]))
        end
        steps += 1
        nothing
    end

    function report()
        branch_tolerances = [1e-8, 1e-7, 1e-8]
        work_tolerance = 0.1 + 1e-8maximum(abs, totals[13:end])
        passed = steps > 0 && nrow(reference) == TRAJECTORY_ROWS &&
            all(net_error .<= 1e-8) && all(branch_error .<= branch_tolerances) &&
            all(absolute_error .<= branch_tolerances) &&
            maximum(quadrature_difference) <= 1e-9 && work_error <= work_tolerance
        Dict{String,Any}("passed" => passed, "steps" => steps, "reference_rows" => nrow(reference),
            "method" => "Held links and dense Cartesian states only; zero-crossing-split composite Simpson 16 versus 8 subintervals; abs(acceleration) integrated directly",
            "net_error_mps" => net_error, "branch_error_mps" => branch_error,
            "absolute_error_mps" => absolute_error, "component_order" => ["R", "T", "N"],
            "net_tolerance_mps" => 1e-8, "branch_tolerances_mps" => branch_tolerances,
            "summed_quadrature_difference_mps" => quadrature_difference,
            "quadrature_tolerance_mps" => 1e-9,
            "work_error_J" => work_error, "work_tolerance_J" => work_tolerance,
            "final_independent_net_mps" => totals[1:3],
            "final_independent_positive_mps" => totals[4:6],
            "final_independent_negative_mps" => totals[7:9],
            "final_independent_absolute_mps" => totals[10:12])
    end
    (; observe, report, reference)
end

function timing_probe(c)
    duration = mission_duration(c)
    warmup = min(duration / c.orbits, duration / 2)
    steady_start_wall, steady_start_simulation = 0.0, 0.0
    last_wall, last_simulation = 0.0, 0.0
    function observe(integrator, start, pairs)
        last_wall, last_simulation = time(), Float64(integrator.t)
        if steady_start_wall == 0 && last_simulation >= warmup
            steady_start_wall, steady_start_simulation = last_wall, last_simulation
        end
        nothing
    end
    function report()
        wall = last_wall - steady_start_wall
        simulated = last_simulation - steady_start_simulation
        steady_start_wall > 0 && simulated > 0 || error("No post-warmup timing window")
        Dict("warmup_simulated_s" => steady_start_simulation,
            "steady_wall_s" => wall, "steady_simulated_s" => simulated,
            "steady_wall_s_per_simulated_s" => wall / simulated)
    end
    (; observe, report)
end
