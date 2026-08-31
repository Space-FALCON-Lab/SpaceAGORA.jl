#=
"""
    Condensed prediction form.

    The QP is written in terms of area increments over the horizon. The output
    perturbation is deltaY = H*deltaU + Mx*deltaX0.
"""
=#
function condensed_form(A_hist, B_hist, C_hist, D_hist)
    N  = length(C_hist)
    nx = size(A_hist[1], 1)
    ny = size(C_hist[1], 1)
    nu = size(B_hist[1], 2)

    Φ_hist = Vector{Matrix{Float64}}(undef, N)
    Φ_hist[1] = Matrix(I, nx, nx)
    for k = 2:N
        Φ_hist[k] = A_hist[k-1] * Φ_hist[k-1]
    end

    Mx = zeros(N*ny, nx)
    for k = 1:N
        rk = (k-1)*ny + 1 : k*ny
        Mx[rk, :] .= C_hist[k] * Φ_hist[k]
    end

    H = zeros(N*ny, N*nu)
    for k = 1:N
        rk = (k-1)*ny + 1 : k*ny
        ck = (k-1)*nu + 1 : k*nu
        H[rk, ck] .= C_hist[k] * B_hist[k] + D_hist[k]
    end

    for j = 1:N-1
        F = copy(B_hist[j])
        for k = j+1:N
            if k > j+1
                F = A_hist[k-1] * F
            end
            rk = (k-1)*ny + 1 : k*ny
            cj = (j-1)*nu + 1 : j*nu
            H[rk, cj] .= C_hist[k] * F
        end
    end

    return Φ_hist, Mx, H
end

function finite_difference_matrix(N::Int)
    D = zeros(N-1, N)
    for k = 1:N-1
        D[k,k]   = -1.0
        D[k,k+1] =  1.0
    end
    return D
end

function build_output_selector(N::Int, ny::Int; yidx::Int)
    S = zeros(N, N*ny)
    for k = 1:N
        S[k, (k-1)*ny + yidx] = 1.0
    end
    return S
end

function predict_outputs(H, Mx, δX0, Ybar, ny::Int, ΔU)
    N = size(Ybar, 1)
    ΔY = H * ΔU + Mx * δX0
    ΔYmat = reshape(ΔY, ny, N)'
    return Ybar .+ ΔYmat
end
