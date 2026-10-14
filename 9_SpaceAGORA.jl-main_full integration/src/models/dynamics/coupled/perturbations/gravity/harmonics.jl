"""
    calcForceTorque(model::GravitationalHarmonicsModel, x, param, i)

Compute spherical-harmonics gravity using coefficients stored in the model's internal fully
normalized convention. By default this includes the inverse-square central term, so it should
not be paired with a separate `InverseSquaredGravityModel` for the same primary body.
"""
# Inner kernel: compute harmonics force/torque given a pre-computed inertial→planet-fixed
# rotation matrix. Called directly by the flat-batch parallel region (which already computed
# L_PI once serially) to avoid per-satellite cache-key allocation inside the parallel loop.
# The one harmonics kernel. `_harmonics_calcforcetorque_with_lpi` (legacy RHS
# path), `wrench` (typed path, gravity backbone) and the flat-route pre-pass all
# call this one compiled body, so every route and thread count agrees bit for
# bit. Never inlined and no fast-math: the compiler must not reassociate,
# contract or fuse differently per call site or instruction set.
# Compiled once and never inlined: StaticArrays products use `muladd`, whose
# fusion LLVM decides per compilation context, so the same expression written in
# two functions can round differently. One body, one codegen, one answer.
# Per-satellite frame terms, as ONE compiled body shared by the scalar kernel and
# the batched one. `L_PI * pos_ii` is a StaticArrays product and StaticArrays uses
# `muladd`, which LLVM may contract to an FMA or not depending on the surrounding
# context. Computing it at two call sites therefore yields `u` values differing in
# the last bit, and `A[2,1] = u * sqrt_3` propagates that through the whole
# Legendre recurrence. Measured before this was shared: scalar A[2,1]
# 0.92297461883779863 against batched 0.92297461883779841.
@noinline function _harmonics_frame_terms(
    L_PI::SMatrix{3, 3, Float64, 9},
    pos_ii::SVector{3, Float64},
)::Tuple{SVector{3, Float64}, Float64, Float64, Float64, Float64}
    rVec_cart = L_PI * pos_ii
    inv_r = 1.0 / norm(rVec_cart)
    return rVec_cart, inv_r, rVec_cart[1] * inv_r, rVec_cart[2] * inv_r, rVec_cart[3] * inv_r
end

# The harmonics back-transform, as ONE compiled body shared by the scalar kernel
# and the batched one. It is `@noinline` for the same reason the kernels are:
# these are StaticArrays products, StaticArrays uses `muladd`, and `muladd` is
# free to contract to an FMA or not. Two call sites that each got their own
# inlined copy could therefore round differently, which is exactly the class of
# difference the flat-route parity probe exists to catch.
@noinline function _harmonics_back_transform(
    a1::Float64, a2::Float64, a3::Float64, a4::Float64,
    s::Float64, t::Float64, u::Float64, inv_r::Float64,
    mass::Float64, rVec_cart::SVector{3, Float64},
    L_PI::SMatrix{3, 3, Float64, 9}, gm::Float64, include_central::Bool,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    g_pp_generic = SVector{3, Float64}(-a1 - s * a4, -a2 - t * a4, -a3 - u * a4)
    g_pp = if include_central
        g_pp_generic - gm * inv_r^3 * rVec_cart
    else
        g_pp_generic
    end
    force_ii = mass * L_PI' * g_pp
    return force_ii, g_pp_generic
end

# No `@noinline` here any more: every contraction-ambiguous operation this body
# used to contain (the StaticArrays products) now lives in _harmonics_frame_terms
# and _harmonics_back_transform, which are themselves single shared bodies. What
# remains is plain `*` and `+`, which LLVM will not fuse without fast-math, so
# inlining this at its several call sites cannot change the answer. Verified by
# the bit-exactness harness, which still reports zero mismatches.
function _harmonics_scalar_force_ii(
    model::GravitationalHarmonicsModel,
    workspace::HarmonicsScratchWorkspace,
    pos_ii::SVector{3, Float64},
    mass::Float64,
    L_PI::SMatrix{3, 3, Float64, 9};
    include_central::Bool=model.include_central,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    rVec_cart, inv_r, s, t, u = _harmonics_frame_terms(L_PI, pos_ii)
    A = workspace.A
    R = workspace.R
    I = workspace.I

    RE = model.reference_radius_m
    L = model.L
    M = model.M

    # NOTE: A[1,1] is already initialised to 1.0 in _make_harmonics_scratch_workspace.
    # We only need to refresh the sub-diagonal entries that depend on the current
    # position (u = z/r changes every call).
    A[2, 1] = u * sqrt_3
    @inbounds for n = 1:L+1
        row = n + 1
        A[row + 1, row] = u * model.sqrt_2n_plus_3[n] * A[row, row]
    end

    # Longitude trig recurrence: R[j] = Re[(s+it)^(j-1)], I[j] = Im[(s+it)^(j-1)]
    R[1] = 1.0
    I[1] = 0.0
    Rn = 1.0
    In = 0.0
    @inbounds for j = 2:(M + 2)
        Rn, In = s * Rn - t * In, s * In + t * Rn
        R[j] = Rn
        I[j] = In
    end

    ρ = RE * inv_r
    a1 = a2 = a3 = a4 = 0.0

    max_recur_row = 2
    ρ_np1 = -model.gm_m3s2 * inv_r * ρ
    @inbounds for l = 1:L
        row = l + 1

        if row > max_recur_row
            jmax = min(max(M, 1) + 1, l - 1)
            @inbounds for j = 1:jmax
                A[row, j] = u * model.N1[row, j] * A[row - 1, j] - model.N2[row, j] * A[row - 2, j]
            end
            max_recur_row = row
        end

        next_row = row + 1
        if next_row > max_recur_row
            jmax_next = min(max(M, 1) + 1, l)
            @inbounds for j = 1:jmax_next
                A[next_row, j] = u * model.N1[next_row, j] * A[next_row - 1, j] - model.N2[next_row, j] * A[next_row - 2, j]
            end
            max_recur_row = next_row
        end

        ρ_np1 *= ρ
        rr = ρ_np1 / RE
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0

        C0 = model.C[row, 1]
        # I[1] == 0.0 always; R[1] == 1.0 always — exploit this to save two multiplies.
        D0 = C0 * sqrt_2
        sum3 += model.VR01[row, 1] * A[row, 2] * D0
        sum4 += model.VR11[row, 1] * A[row + 1, 2] * D0

        active_orders = model.active_orders_by_degree[row]
        @inbounds for idx in eachindex(active_orders)
            m = active_orders[idx]
            j = m + 1
            C = model.C[row, j]
            S = model.S[row, j]
            R_term = R[j - 1]
            I_term = I[j - 1]
            Rj  = R[j]
            Ij  = I[j]
            D = (C * Rj  + S * Ij)  * sqrt_2
            E = (C * R_term + S * I_term) * sqrt_2
            F = (S * R_term - C * I_term) * sqrt_2

            mA = m * A[row, j]
            sum1 += mA * E
            sum2 += mA * F
            sum3 += model.VR01[row, j] * A[row, j + 1] * D
            sum4 += model.VR11[row, j] * A[row + 1, j + 1] * D
        end

        a1 += rr * sum1
        a2 += rr * sum2
        a3 += rr * sum3
        a4 -= rr * sum4
    end

    return _harmonics_back_transform(
        a1, a2, a3, a4, s, t, u, inv_r, mass, rVec_cart,
        L_PI, model.gm_m3s2, include_central,
    )
end

@inline function _harmonics_calcforcetorque_with_lpi(
    model::GravitationalHarmonicsModel,
    x::AbstractVector{Float64},
    param::ODEParams,
    i::Int64,
    L_PI::SMatrix{3, 3, Float64, 9},
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    mass = Float64(x[7])

    workspace = _harmonics_workspace_for_sat!(model, param, i)
    force_ii, g_pp_generic = _harmonics_scalar_force_ii(model, workspace, pos_ii, mass, L_PI)
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    if _DEBUG_COMPARE_J2[] && model.L == 2 && model.M == 0
        C20 = model.C[3, 1]
        if isfinite(C20) && C20 != 0.0
            rVec_cart = L_PI * pos_ii
            x_pp, y_pp, z_pp = rVec_cart
            r = norm(rVec_cart)
            RE = model.reference_radius_m
            r2 = r * r
            r4 = r2 * r2
            z2 = z_pp * z_pp
            J2 = -sqrt(5.0) * C20

            common = 1.5 * model.gm_m3s2 * J2 * RE^2 / r4
            g_pp_analytic = SVector{3, Float64}(
                common * (x_pp / r) * (5.0 * z2 / r2 - 1.0),
                common * (y_pp / r) * (5.0 * z2 / r2 - 1.0),
                common * (z_pp / r) * (5.0 * z2 / r2 - 3.0)
            )

            Δg_pp = g_pp_generic - g_pp_analytic
            rel = norm(Δg_pp) / max(norm(g_pp_analytic), 1e-30)
            prev = Base.Threads.atomic_add!(_J2_COMPARE_LOG_COUNT, 1)
            if prev < 20
                @info "J2 generic-vs-analytic comparison" sat=i et=et r_m=r c20=C20 j2=J2 g_generic_pp=norm(g_pp_generic) g_analytic_pp=norm(g_pp_analytic) delta_pp=norm(Δg_pp) rel_pp=rel
            end
        end
    end
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function calcForceTorque(model::GravitationalHarmonicsModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    L_PI = _harmonics_lpi_at!(model, param, et)
    return _harmonics_calcforcetorque_with_lpi(model, x, param, i, L_PI)
end

@inline environment_requirements(::GravitationalHarmonicsModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::GravitationalHarmonicsModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("GravitationalHarmonicsModel wrench requires env.planet_frame."))

    workspace = _make_harmonics_scratch_workspace(model)
    force_ii, _ = _harmonics_scalar_force_ii(model, workspace, x.pos_ii, x.mass_kg, planet_frame.l_pi)
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline gravity_backbone_structure(::GravitationalHarmonicsModel) = :position_only_static_gravity

@inline function gravity_backbone_acceleration_ii(
    model::GravitationalHarmonicsModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    force_ii, _ = wrench(model, x, env, t)
    return force_ii / x.mass_kg
end

# Relative gravity: the central term by Encke, the non-central part (J2 and higher harmonics) differenced
# plainly. Without `include_central` the model has no point-mass term and the whole difference is plain.
@inline function gravity_backbone_relative_acceleration_ii(
    model::GravitationalHarmonicsModel,
    x_base::StateSample, ρ, x_far::StateSample, env_base::EnvironmentSample, env_far::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    pb = env_base.planet_frame
    pf = env_far.planet_frame
    (pb === nothing || pf === nothing) && throw(ArgumentError("GravitationalHarmonicsModel relative gravity requires env.planet_frame."))
    workspace = _make_harmonics_scratch_workspace(model)
    f_far, _ = _harmonics_scalar_force_ii(model, workspace, x_far.pos_ii, 1.0, pf.l_pi; include_central=false)
    f_base, _ = _harmonics_scalar_force_ii(model, workspace, x_base.pos_ii, 1.0, pb.l_pi; include_central=false)
    noncentral = f_far - f_base
    model.include_central || return noncentral
    return encke_point_mass_difference(model.gm_m3s2, x_base.pos_ii, ρ) + noncentral
end
