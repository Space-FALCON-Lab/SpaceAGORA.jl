module TwoBodyLinearization

using LinearAlgebra

"""
    two_body_f(x::AbstractVector{<:Real}, u::AbstractVector{<:Real}, μ::Real) -> Vector{Float64}

Two-body dynamics function f(x, u) = [v; a_grav + u].

# Arguments
- `x::AbstractVector{<:Real}`: State vector [r; v] (6 elements)
- `u::AbstractVector{<:Real}`: Control input (3 elements)
- `μ::Real`: Gravitational parameter

# Returns
- `Vector{Float64}`: Time derivative of state
"""
function two_body_f(x::AbstractVector{<:Real}, u::AbstractVector{<:Real}, μ::Real)
    @assert length(x) == 6
    @assert length(u) == 3
    r = @view x[1:3]
    v = @view x[4:6]
    rnorm = norm(r)
    a_grav = -μ * r / (rnorm^3)
    xdot = similar(x, Float64)
    xdot[1:3] .= v
    xdot[4:6] .= a_grav .+ u
    return xdot
end

"""
    two_body_fd_rk4(x::AbstractVector{<:Real}, u::AbstractVector{<:Real}, μ::Real, h::Real) -> Vector{Float64}

RK4 integration step for two-body dynamics.

# Arguments
- `x::AbstractVector{<:Real}`: State vector [r; v] (6 elements)
- `u::AbstractVector{<:Real}`: Control input (3 elements)
- `μ::Real`: Gravitational parameter
- `h::Real`: Time step

# Returns
- `Vector{Float64}`: State at time t + h
"""
function two_body_fd_rk4(x::AbstractVector{<:Real}, u::AbstractVector{<:Real}, μ::Real, h::Real)
    x0 = Float64.(x)
    u0 = Float64.(u)
    k1 = two_body_f(x0,               u0, μ)
    k2 = two_body_f(x0 .+ 0.5h .* k1, u0, μ)
    k3 = two_body_f(x0 .+ 0.5h .* k2, u0, μ)
    k4 = two_body_f(x0 .+ h   .* k3,  u0, μ)
    return x0 .+ (h / 6.0) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
end

"""
    discrete_taylor_linearize_fd(fd::Function, xbar::AbstractVector{<:Real}, ubar::AbstractVector{<:Real}; epsx::Real=1e-6, epsu::Real=1e-6) -> Tuple{Matrix{Float64}, Matrix{Float64}}

Linearize a discrete-time function using finite differences.

# Arguments
- `fd::Function`: Discrete-time dynamics function fd(x, u)
- `xbar::AbstractVector{<:Real}`: Reference state
- `ubar::AbstractVector{<:Real}`: Reference control
- `epsx::Real=1e-6`: Finite difference step for state
- `epsu::Real=1e-6`: Finite difference step for control

# Returns
- `Tuple{Matrix{Float64}, Matrix{Float64}}`: (A, B) matrices where x_{k+1} ≈ A x_k + B u_k
"""
function discrete_taylor_linearize_fd(
    fd::Function,
    xbar::AbstractVector{<:Real},
    ubar::AbstractVector{<:Real};
    epsx::Real = 1e-6,
    epsu::Real = 1e-6,
)
    n = length(xbar)
    m = length(ubar)
    x0 = Float64.(xbar)
    u0 = Float64.(ubar)

    f0 = fd(x0, u0)
    @assert length(f0) == n

    A = zeros(Float64, n, n)
    B = zeros(Float64, n, m)

    for i in 1:n
        dx = zeros(Float64, n)
        dx[i] = epsx
        fp = fd(x0 .+ dx, u0)
        fm = fd(x0 .- dx, u0)
        A[:, i] .= (fp .- fm) ./ (2epsx)
    end

    for j in 1:m
        du = zeros(Float64, m)
        du[j] = epsu
        fp = fd(x0, u0 .+ du)
        fm = fd(x0, u0 .- du)
        B[:, j] .= (fp .- fm) ./ (2epsu)
    end

    return A, B
end

"""
    van_loan_discretize(Ac::AbstractMatrix{<:Real}, Bc::AbstractMatrix{<:Real}, h::Real) -> Tuple{Matrix{Float64}, Matrix{Float64}}

Discretize continuous-time linear system using Van Loan's method.

For continuous-time system ẋ = Ac x + Bc u, the discrete-time system is:
x_{k+1} = Ad x_k + Bd u_k

# Arguments
- `Ac::AbstractMatrix{<:Real}`: Continuous-time state matrix
- `Bc::AbstractMatrix{<:Real}`: Continuous-time input matrix
- `h::Real`: Time step

# Returns
- `Tuple{Matrix{Float64}, Matrix{Float64}}`: (Ad, Bd) discrete-time matrices
"""
function van_loan_discretize(Ac::AbstractMatrix{<:Real}, Bc::AbstractMatrix{<:Real}, h::Real)
    n = size(Ac, 1)
    m = size(Bc, 2)
    @assert size(Ac, 2) == n
    @assert size(Bc, 1) == n

    M = zeros(Float64, n + m, n + m)
    M[1:n, 1:n] .= Ac
    M[1:n, n+1:n+m] .= Bc
    M .*= h

    E = exp(M)
    Ad = E[1:n, 1:n]
    Bd = E[1:n, n+1:n+m]
    return Ad, Bd
end

"""
    two_body_Ac(r::AbstractVector{<:Real}, μ::Real) -> Matrix{Float64}

Compute the continuous-time state matrix for two-body dynamics at position r.

Ac = [0 I; ∂a_grav/∂r 0]

# Arguments
- `r::AbstractVector{<:Real}`: Position vector (3 elements)
- `μ::Real`: Gravitational parameter

# Returns
- `Matrix{Float64}`: 6×6 continuous-time state matrix
"""
function two_body_Ac(r::AbstractVector{<:Real}, μ::Real)
    @assert length(r) == 3
    r2 = dot(r, r)
    r1 = sqrt(r2)
    r3 = r1^3
    r5 = r1^5

    I3 = Matrix{Float64}(I, 3, 3)
    rrT = Float64.(r) * Float64.(r)'
    dadr = -μ * (I3 / r3 - 3.0 * rrT / r5)

    Ac = zeros(6, 6)
    Ac[1:3, 4:6] .= I3
    Ac[4:6, 1:3] .= dadr
    return Ac
end

"""
    two_body_AdBd(xbar_or_rbar::AbstractVector{<:Real}, μ::Real, h::Real; method::Symbol=:vanloan, ubar::AbstractVector{<:Real}=zeros(3), epsx::Real=1e-6, epsu::Real=1e-6) -> Tuple{Matrix{Float64}, Matrix{Float64}}

Compute discrete-time state transition matrices for two-body dynamics.

# Arguments
- `xbar_or_rbar::AbstractVector{<:Real}`: Reference state (6 elements) or position (3 elements)
- `μ::Real`: Gravitational parameter
- `h::Real`: Time step
- `method::Symbol=:vanloan`: Discretization method (:vanloan or :discrete_taylor_fd)
- `ubar::AbstractVector{<:Real}=zeros(3)`: Reference control (for discrete_taylor_fd)
- `epsx::Real=1e-6`: Finite difference step for state (for discrete_taylor_fd)
- `epsu::Real=1e-6`: Finite difference step for control (for discrete_taylor_fd)

# Returns
- `Tuple{Matrix{Float64}, Matrix{Float64}}`: (Ad, Bd) discrete-time matrices
"""
function two_body_AdBd(
    xbar_or_rbar::AbstractVector{<:Real},
    μ::Real,
    h::Real;
    method::Symbol = :vanloan,
    ubar::AbstractVector{<:Real} = zeros(3),
    epsx::Real = 1e-6,
    epsu::Real = 1e-6,
)
    if method == :vanloan
        rbar = length(xbar_or_rbar) == 3 ? xbar_or_rbar : @view xbar_or_rbar[1:3]
        Ac = two_body_Ac(rbar, μ)
        Bc = zeros(Float64, 6, 3)
        Bc[4:6, 1:3] .= Matrix{Float64}(I, 3, 3)
        return van_loan_discretize(Ac, Bc, h)
    elseif method == :discrete_taylor_fd
        @assert length(xbar_or_rbar) == 6 "paper-style discrete Taylor needs full xbar (6-state)"
        xbar = xbar_or_rbar
        fd = (x, u) -> two_body_fd_rk4(x, u, μ, h)
        return discrete_taylor_linearize_fd(fd, xbar, ubar; epsx=epsx, epsu=epsu)
    else
        error("Unknown method=$method. Use :vanloan or :discrete_taylor_fd")
    end
end

"""
    cw_coefficients(a::Real, μ::Real) -> NamedTuple

Compute Clohessy-Wiltshire coefficients for circular orbit.

# Arguments
- `a::Real`: Semi-major axis [m]
- `μ::Real`: Gravitational parameter [m^3/s^2]

# Returns
- `NamedTuple`: (n, ω0, τ) where n is mean motion, ω0 is orbital frequency, τ is orbital period
"""
function cw_coefficients(a::Real, μ::Real)
    n = sqrt(μ / a^3)  # Mean motion
    ω0 = n              # Orbital frequency
    τ = 2π / n          # Orbital period
    return (n=n, ω0=ω0, τ=τ)
end

export two_body_f, two_body_fd_rk4, discrete_taylor_linearize_fd,
       van_loan_discretize, two_body_Ac, two_body_AdBd, cw_coefficients

end # module TwoBodyLinearization
