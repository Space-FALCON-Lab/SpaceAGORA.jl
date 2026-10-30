using Test
using LinearAlgebra

include(joinpath(@__DIR__, "../../../examples/periapsis_corridor_mean_elements.jl"))
using .PeriapsisCorridorMeanElements

const PCME = PeriapsisCorridorMeanElements

function _j2_acceleration(r, mu, radius, j2)
    r2 = sum(abs2, r)
    rmag = sqrt(r2)
    z2r2 = r[3]^2 / r2
    scale = 1.5 * j2 * mu * radius^2 / rmag^5
    return (-mu * r[1] / rmag^3 + scale * r[1] * (5z2r2 - 1),
            -mu * r[2] / rmag^3 + scale * r[2] * (5z2r2 - 1),
            -mu * r[3] / rmag^3 + scale * r[3] * (5z2r2 - 3))
end

function _j2_rk4_step(state, dt, mu, radius, j2)
    derivative(s) = (s[4], s[5], s[6], _j2_acceleration((s[1], s[2], s[3]), mu, radius, j2)...)
    add_scaled(s, k, scale) = ntuple(i -> s[i] + scale * k[i], 6)
    k1 = derivative(state)
    k2 = derivative(add_scaled(state, k1, dt / 2))
    k3 = derivative(add_scaled(state, k2, dt / 2))
    k4 = derivative(add_scaled(state, k3, dt))
    return ntuple(i -> state[i] + dt / 6 * (k1[i] + 2k2[i] + 2k3[i] + k4[i]), 6)
end

function _kepler_true_anomaly(mean_anomaly, eccentricity)
    E = mean_anomaly
    for _ in 1:12
        correction = (E - eccentricity * sin(E) - mean_anomaly) /
                     (1 - eccentricity * cos(E))
        E -= correction
        abs(correction) < 2e-15 && break
    end
    return 2atan(sqrt(1 + eccentricity) * sin(E / 2),
                 sqrt(1 - eccentricity) * cos(E / 2))
end

function _j2_map_truth_errors(mean, mu, radius, j2; steps=8192)
    n0 = sqrt(mu / mean.a^3)
    pbar = mean.a * (1 - mean.e^2)
    nbar = n0 * (1 + 0.75 * j2 * (radius / pbar)^2 * sqrt(1 - mean.e^2) *
                 (2 - 3sin(mean.i)^2))
    node_rate = -1.5 * nbar * j2 * (radius / pbar)^2 * cos(mean.i)
    periapsis_rate = 0.75 * nbar * j2 * (radius / pbar)^2 *
                     (4 - 5sin(mean.i)^2)
    period = 2π / nbar
    initial_osc = mean_to_osculating(mean; mu, radius, j2).elements
    r0, v0 = elements_to_cartesian(initial_osc; mu)
    state = (r0..., v0...)
    dt = period / steps
    max_position_error = max_velocity_error = 0.0
    periapsis_position_error = periapsis_velocity_error = 0.0
    apoapsis_position_error = apoapsis_velocity_error = 0.0
    initial_E = 2atan(sqrt(1 - mean.e) * sin(mean.nu / 2),
                      sqrt(1 + mean.e) * cos(mean.nu / 2))
    initial_M = initial_E - mean.e * sin(initial_E)
    periapsis_step = round(Int, mod(-initial_M, 2π) / nbar / dt)
    apoapsis_step = round(Int, mod(π - initial_M, 2π) / nbar / dt)
    sample_steps = Set(vcat(collect(0:(steps ÷ 32):steps),
                            [periapsis_step, apoapsis_step]))
    for step in 0:steps
        if step in sample_steps
            current_M = initial_M + nbar * step * dt
            current_nu = _kepler_true_anomaly(current_M, mean.e)
            secular_mean = merge(mean, (raan=mean.raan + node_rate * step * dt,
                omega=mean.omega + periapsis_rate * step * dt, nu=current_nu))
            predicted = mean_to_osculating(secular_mean; mu, radius, j2).elements
            rp, vp = elements_to_cartesian(predicted; mu)
            position_error = norm(collect(rp) - collect(state[1:3]))
            velocity_error = norm(collect(vp) - collect(state[4:6]))
            max_position_error = max(max_position_error, position_error)
            max_velocity_error = max(max_velocity_error, velocity_error)
            if step == periapsis_step
                periapsis_position_error = position_error
                periapsis_velocity_error = velocity_error
            elseif step == apoapsis_step
                apoapsis_position_error = position_error
                apoapsis_velocity_error = velocity_error
            end
        end
        step < steps && (state = _j2_rk4_step(state, dt, mu, radius, j2))
    end
    return (; max_position_error, max_velocity_error, periapsis_position_error,
            periapsis_velocity_error, apoapsis_position_error, apoapsis_velocity_error)
end

_radial_velocity_dot(state) = sum(state[1:3] .* state[4:6])

function _j2_radial_minimum_comparison(mean, mu, radius, j2;
                                      steps=65_536, map_steps=512)
    pbar = mean.a * (1 - mean.e^2)
    n0 = sqrt(mu / mean.a^3)
    nbar = n0 * (1 + 0.75 * j2 * (radius / pbar)^2 * sqrt(1 - mean.e^2) *
                 (2 - 3sin(mean.i)^2))
    node_rate = -1.5 * nbar * j2 * (radius / pbar)^2 * cos(mean.i)
    periapsis_rate = 0.75 * nbar * j2 * (radius / pbar)^2 *
                     (4 - 5sin(mean.i)^2)
    period = 2π / nbar
    E0 = 2atan(sqrt(1 - mean.e) * sin(mean.nu / 2),
               sqrt(1 + mean.e) * cos(mean.nu / 2))
    M0 = E0 - mean.e * sin(E0)
    initial_osc = mean_to_osculating(mean; mu, radius, j2).elements
    r0, v0 = elements_to_cartesian(initial_osc; mu)
    state = (r0..., v0...)
    dt = period / steps
    truth_radius = NaN
    truth_time = NaN
    for step in 1:steps
        next_state = _j2_rk4_step(state, dt, mu, radius, j2)
        if _radial_velocity_dot(state) ≤ 0 < _radial_velocity_dot(next_state)
            lo, hi = 0.0, dt
            for _ in 1:48
                mid = (lo + hi) / 2
                mid_state = _j2_rk4_step(state, mid, mu, radius, j2)
                if _radial_velocity_dot(mid_state) > 0
                    hi = mid
                else
                    lo = mid
                end
            end
            event_dt = (lo + hi) / 2
            event_state = _j2_rk4_step(state, event_dt, mu, radius, j2)
            truth_radius = norm(event_state[1:3])
            truth_time = (step - 1) * dt + event_dt
            break
        end
        state = next_state
    end
    isfinite(truth_radius) || error("failed to locate a Cartesian periapsis in one orbit")

    function map_radius(t)
        M = M0 + nbar * t
        nu = _kepler_true_anomaly(M, mean.e)
        mean_t = merge(mean, (raan=mean.raan + node_rate * t,
                              omega=mean.omega + periapsis_rate * t, nu=nu))
        osc = mean_to_osculating(mean_t; mu, radius, j2).elements
        return osc.a * (1 - osc.e^2) / (1 + osc.e * cos(osc.nu))
    end

    grid_t = range(0.0, period; length=map_steps + 1)
    grid_r = map(map_radius, grid_t)
    k = argmin(grid_r)
    left = grid_t[max(k - 1, 1)]
    right = grid_t[min(k + 1, length(grid_t))]
    ϕ = (sqrt(5) - 1) / 2
    x1, x2 = right - ϕ * (right - left), left + ϕ * (right - left)
    y1, y2 = map_radius(x1), map_radius(x2)
    for _ in 1:64
        if y1 < y2
            right, x2, y2 = x2, x1, y1
            x1 = right - ϕ * (right - left)
            y1 = map_radius(x1)
        else
            left, x1, y1 = x1, x2, y2
            x2 = left + ϕ * (right - left)
            y2 = map_radius(x2)
        end
    end
    map_time = (left + right) / 2
    map_radius_min = map_radius(map_time)
    return (; truth_radius, truth_time, map_radius=map_radius_min, map_time,
            radial_difference=map_radius_min - truth_radius)
end

@testset "J2 short-period mean/osculating conversion" begin
    μₑ = 3.986004418e14
    Rₑ = 6_378_136.3
    J₂ₑ = 1.08262668e-3

    @testset "high-precision first-order J2 reference case" begin
        # This independent reference was evaluated from the published
        # Vallado–Kwok short-period equations, with the parameters shown here.
        # (The upstream regression table is rounded and is therefore not used
        # as a numerical assertion.)
        mean = (a=7_190_982.0, e=0.001111, i=deg2rad(98.405),
                raan=deg2rad(90.0), omega=deg2rad(200.0), nu=deg2rad(45.0))
        mapped = mean_to_osculating(mean; mu=μₑ, radius=Rₑ, j2=J₂ₑ)

        @test mapped.elements.a ≈ 7_185_194.839432819 atol=1e-6
        @test mapped.elements.e ≈ 0.0006706278229540496 atol=1e-15
        @test mapped.elements.i ≈ 1.7175506481964766 atol=1e-14
        @test mapped.elements.raan ≈ 1.5707251346040267 atol=1e-14
        @test mapped.elements.omega ≈ 1.9014309490087546 atol=1e-14
        @test mapped.elements.nu ≈ 2.3745572583831631 atol=1e-14
        @test mapped.converged
        @test mapped.iterations == 0
    end

    @testset "high-eccentricity Earth aerobraking geometry round trips" begin
        rp = Rₑ + 97_500.0
        ra = Rₑ + 60_500_000.0
        mean = (a=(rp + ra) / 2, e=(ra - rp) / (ra + rp),
                i=deg2rad(50.0), raan=deg2rad(27.0),
                omega=deg2rad(90.0), nu=deg2rad(215.0))

        osculating = mean_to_osculating(mean; mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test !osculating.corridor_safety_certified
        recovered = osculating_to_mean(osculating.elements;
            mu=μₑ, radius=Rₑ, j2=J₂ₑ)

        @test osculating.converged
        @test recovered.converged
        @test recovered.iterations > 0
        @test !recovered.corridor_safety_certified
        @test recovered.residual_norm ≤ 1e-10
        @test recovered.elements.a ≈ mean.a rtol=2e-11
        @test recovered.elements.e ≈ mean.e atol=2e-12
        @test recovered.elements.i ≈ mean.i atol=2e-12
        @test PCME.angle_difference(recovered.elements.raan, mean.raan) ≈ 0 atol=2e-12
        @test PCME.angle_difference(recovered.elements.omega, mean.omega) ≈ 0 atol=2e-12
        @test PCME.angle_difference(recovered.elements.nu, mean.nu) ≈ 0 atol=2e-12

        # This map only accounts for J2 short-period terms. The example's J3
        # coefficient belongs to the truth model and is intentionally separate.
        @test isfinite(osculating.elements.a * (1 - osculating.elements.e))

        # Sweep phase across periapsis and apoapsis for both this high-e orbit
        # and a low-e reference orbit; inversion must remain well conditioned.
        low_e = (a=7_000_000.0, e=0.005, i=deg2rad(50.0), raan=0.3,
                 omega=1.4, nu=0.0)
        for base in (mean, low_e), phase in range(0.0, 2π; length=17)[1:end-1]
            phased = merge(base, (nu=phase,))
            mapped_phase = mean_to_osculating(phased; mu=μₑ, radius=Rₑ, j2=J₂ₑ)
            recovered_phase = osculating_to_mean(mapped_phase.elements;
                mu=μₑ, radius=Rₑ, j2=J₂ₑ)
            @test recovered_phase.converged
            @test !recovered_phase.corridor_safety_certified
            @test recovered_phase.residual_norm ≤ 1e-12
            @test recovered_phase.elements.e ≈ phased.e atol=2e-11
            @test PCME.angle_difference(recovered_phase.elements.nu, phase) ≈ 0 atol=2e-10
        end

        just_after_wrap = mean_to_osculating(merge(low_e, (nu=1e-10,));
            mu=μₑ, radius=Rₑ, j2=J₂ₑ).elements
        just_before_wrap = mean_to_osculating(merge(low_e, (nu=2π - 1e-10,));
            mu=μₑ, radius=Rₑ, j2=J₂ₑ).elements
        @test just_after_wrap.a ≈ just_before_wrap.a atol=1e-5
        @test just_after_wrap.e ≈ just_before_wrap.e atol=1e-12
        @test just_after_wrap.i ≈ just_before_wrap.i atol=1e-12
        @test PCME.angle_difference(just_after_wrap.raan, just_before_wrap.raan) ≈ 0 atol=1e-12
        @test PCME.angle_difference(just_after_wrap.omega, just_before_wrap.omega) ≈ 0 atol=2e-10
        @test PCME.angle_difference(just_after_wrap.nu, just_before_wrap.nu) ≈ 0 atol=3e-10

        # Independent Cartesian RK4 truth under central + J2 acceleration on
        # the high-e aerobraking seed. Secular rates and anomalies are evaluated
        # separately from first-order Kozai rates documented in SatelliteToolbox
        # j2.jl; the map supplies only the phase-local short-period correction.
        truth_mean = mean
        n0 = sqrt(μₑ / truth_mean.a^3)
        pbar = truth_mean.a * (1 - truth_mean.e^2)
        nbar = n0 * (1 + 0.75 * J₂ₑ * (Rₑ / pbar)^2 * sqrt(1 - truth_mean.e^2) *
                     (2 - 3sin(truth_mean.i)^2))
        node_rate = -1.5 * nbar * J₂ₑ * (Rₑ / pbar)^2 * cos(truth_mean.i)
        periapsis_rate = 0.75 * nbar * J₂ₑ * (Rₑ / pbar)^2 * (4 - 5sin(truth_mean.i)^2)
        period = 2π / nbar
        truth_osculating = mean_to_osculating(truth_mean;
            mu=μₑ, radius=Rₑ, j2=J₂ₑ).elements
        r0, v0 = elements_to_cartesian(truth_osculating; mu=μₑ)
        state = (r0..., v0...)
        steps = 8192
        dt = period / steps
        max_position_error = 0.0
        max_velocity_error = 0.0
        periapsis_position_error = 0.0
        periapsis_velocity_error = 0.0
        apoapsis_position_error = 0.0
        apoapsis_velocity_error = 0.0
        initial_E = 2atan(sqrt(1 - truth_mean.e) * sin(truth_mean.nu / 2),
                          sqrt(1 + truth_mean.e) * cos(truth_mean.nu / 2))
        initial_M = initial_E - truth_mean.e * sin(initial_E)
        periapsis_step = round(Int, mod(-initial_M, 2π) / nbar / dt)
        apoapsis_step = round(Int, mod(π - initial_M, 2π) / nbar / dt)
        sample_steps = Set(vcat(collect(0:(steps ÷ 32):steps),
                                [periapsis_step, apoapsis_step]))
        for step in 0:steps
            t = step * dt
            if step in sample_steps
                current_M = initial_M + nbar * t
                current_nu = _kepler_true_anomaly(current_M, truth_mean.e)
                secular_mean = merge(truth_mean, (raan=truth_mean.raan + node_rate * t,
                                             omega=truth_mean.omega + periapsis_rate * t,
                                             nu=current_nu))
                predicted = mean_to_osculating(secular_mean;
                    mu=μₑ, radius=Rₑ, j2=J₂ₑ).elements
                rp, vp = elements_to_cartesian(predicted; mu=μₑ)
                position_error = norm(collect(rp) - collect(state[1:3]))
                velocity_error = norm(collect(vp) - collect(state[4:6]))
                max_position_error = max(max_position_error, position_error)
                max_velocity_error = max(max_velocity_error, velocity_error)
                if step == periapsis_step
                    periapsis_position_error = position_error
                    periapsis_velocity_error = velocity_error
                elseif step == apoapsis_step
                    apoapsis_position_error = position_error
                    apoapsis_velocity_error = velocity_error
                end
            end
            step < steps && (state = _j2_rk4_step(state, dt, μₑ, Rₑ, J₂ₑ))
        end
        @info "high-e J2 map phase errors" max_position_error max_velocity_error periapsis_position_error periapsis_velocity_error apoapsis_position_error apoapsis_velocity_error
        # The vector residual above includes along-track phase error and is not
        # itself a periapsis-radius error. Compare minimum radii separately.
        radial = _j2_radial_minimum_comparison(mean, μₑ, Rₑ, J₂ₑ)
        @info "high-e J2 map physical radial minimum comparison" radial.truth_radius radial.map_radius radial.radial_difference radial.truth_time radial.map_time
        @test isfinite(radial.radial_difference)
        @test abs(radial.radial_difference) < 2_000.0
        coarser_radial = _j2_radial_minimum_comparison(mean, μₑ, Rₑ, J₂ₑ;
                                                       steps=32_768)
        @test abs(coarser_radial.truth_radius - radial.truth_radius) < 10.0
    end

    @testset "periapsis-radius accuracy for the vacuum seed" begin
        # Compare the actual Cartesian minimum radius with the minimum radius
        # of the mapped secular trajectory. This is the corridor observable;
        # full Cartesian phase error is reported separately above.
        mean = (a=Rₑ + 700_000.0, e=0.001, i=deg2rad(50.0),
                raan=deg2rad(27.0), omega=deg2rad(90.0), nu=deg2rad(215.0))
        mapped = mean_to_osculating(mean; mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test !mapped.corridor_safety_certified
        radial = _j2_radial_minimum_comparison(mean, μₑ, Rₑ, J₂ₑ)
        @info "vacuum-seed physical radial minimum comparison" radial.truth_radius radial.map_radius radial.radial_difference radial.truth_time radial.map_time
        @test abs(radial.radial_difference) < 2_000.0
    end

    @testset "e=0.15 periapsis-radius diagnostic" begin
        mean = (a=8_200_000.0, e=0.15, i=deg2rad(50.0), raan=deg2rad(27.0),
                omega=deg2rad(90.0), nu=deg2rad(215.0))
        radial = _j2_radial_minimum_comparison(mean, μₑ, Rₑ, J₂ₑ)
        @info "e=0.15 radial minimum diagnostic" radial.truth_radius radial.map_radius radial.radial_difference radial.truth_time radial.map_time
        @test isfinite(radial.radial_difference)
        @test abs(radial.radial_difference) < 2_000.0
        mapped = mean_to_osculating(mean; mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test !mapped.corridor_safety_certified
    end

    @testset "e=0.15 full-state phase residual diagnostic" begin
        # This is not a periapsis-radius error or a corridor safety assertion.
        mean = (a=8_200_000.0, e=0.15, i=deg2rad(50.0), raan=deg2rad(27.0),
                omega=deg2rad(90.0), nu=deg2rad(215.0))
        mapped = mean_to_osculating(mean; mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test !mapped.corridor_safety_certified
        errors = _j2_map_truth_errors(mean, μₑ, Rₑ, J₂ₑ)
        @info "e=0.15 J2 map phase residuals (not periapsis-radius errors)" errors.max_position_error errors.max_velocity_error errors.periapsis_position_error errors.periapsis_velocity_error
        @test errors.max_position_error > 2_000.0
        @test errors.max_velocity_error > 2.0
        @test errors.periapsis_position_error < 2_000.0
        @test errors.periapsis_velocity_error < 2.0
    end

    @testset "independent first-order radial-rate fixtures" begin
        # e=0.15 rows from SatelliteToolboxPropagators v1.2.0's pinned
        # Osculating Radial Rate fixture. It derives mean elements from a
        # six-orbit Hann-windowed Cartesian RK4 propagation; see source notes.
        reference = [
            (7181478.5102035254, 0.1496306659, 1.7175806949, 3.4879659757, 0.0019215108,  0.0000000140),
            (7200571.5026989058, 0.1514567721, 1.7174260880, 3.4954699797, 0.5839455326, 798.7165162242),
            (7198483.9469006751, 0.1498757265, 1.7174115022, 3.5011015424, 1.2624022277, 1129.5557297251),
            (7185794.5462470576, 0.1497891434, 1.7175398742, 3.4898013120, 2.1280406405, 798.7165162118),
            (7185378.2525120089, 0.1500913695, 1.7175510083, 3.4918752282, 3.1398949178, -0.0000000076),
            (7194644.1613620706, 0.1490904054, 1.7174331680, 3.4859319271, 4.1607463286, -798.7165162234),
            (7198483.2383694621, 0.1504233653, 1.7174239326, 3.4812668302, 5.0209647203, -1129.5557297248),
            (7183802.8619856136, 0.1501249115, 1.7175679770, 3.4934212469, 5.6916574363, -798.7165162062),
        ]
        μₑ_egm2008 = 3.986004415e14
        Rₑ_egm2008 = 6_378_137.0
        J₂ₑ_egm2008 = 0.0010826261738522227
        for (a, e, inc, omega, M, radial_rate_reference) in reference
            mean_anomaly = mod(M, 2π)
            nu = _kepler_true_anomaly(mean_anomaly, e)
            mean = (a=a, e=e, i=inc, raan=deg2rad(90.0), omega=omega, nu=nu)
            osculating = mean_to_osculating(mean; mu=μₑ_egm2008,
                radius=Rₑ_egm2008, j2=J₂ₑ_egm2008).elements
            r, v = elements_to_cartesian(osculating; mu=μₑ_egm2008)
            radial_rate = sum(r .* v) / norm(r)
            @test radial_rate ≈ radial_rate_reference atol=0.35
        end
    end

    @testset "two-body limit and input domain" begin
        elements = (a=7_000_000.0, e=0.1, i=0.8, raan=1.1, omega=2.0, nu=3.4)
        mapped = mean_to_osculating(elements; mu=μₑ, radius=Rₑ, j2=0.0)
        @test mapped.elements.a ≈ elements.a
        @test mapped.elements.e ≈ elements.e
        @test mapped.elements.i ≈ elements.i
        @test PCME.angle_difference(mapped.elements.raan, elements.raan) ≈ 0 atol=1e-14
        @test PCME.angle_difference(mapped.elements.omega, elements.omega) ≈ 0 atol=1e-14
        @test PCME.angle_difference(mapped.elements.nu, elements.nu) ≈ 0 atol=1e-14

        @test_throws ArgumentError mean_to_osculating(
            merge(elements, (e=1.0,)); mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test_throws ArgumentError mean_to_osculating(
            merge(elements, (a=Rₑ,)); mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test_throws ArgumentError mean_to_osculating(
            merge(elements, (e=0.0,)); mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test_throws ArgumentError mean_to_osculating(
            merge(elements, (i=0.0,)); mu=μₑ, radius=Rₑ, j2=J₂ₑ)

        # Critical inclination remains mathematically supported by the local
        # conversion; the first corridor case avoids it because secular apsidal
        # behavior there needs separate study.
        critical_i = acos(sqrt(1 / 5))
        critical = mean_to_osculating(merge(elements, (i=critical_i,));
            mu=μₑ, radius=Rₑ, j2=J₂ₑ)
        @test critical.converged
        @test !critical.corridor_safety_certified

        @test_throws ErrorException osculating_to_mean(
            merge(elements, (a=elements.a + 100_000.0,)); mu=μₑ,
            radius=Rₑ, j2=J₂ₑ, max_iterations=1, tolerance=1e-30)
    end
end
