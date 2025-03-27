using SPICE
using LoopVectorization
using AssociatedLegendrePolynomials
using LinearAlgebra
# import .config

# Define delta function
δ(x,y) = ==(x,y)
δ(x) = δ(x,0)
function gravity_n_bodies(et, pos_ii, p, n_body)

    primary_body_name = p.name
    n_body_name = n_body.name

    if cmp(lowercase(primary_body_name), "mars") == 0 || cmp(lowercase(primary_body_name), "jupiter") == 0 || cmp(lowercase(primary_body_name), "saturn") == 0 || cmp(lowercase(primary_body_name), "uranus") == 0 || cmp(lowercase(primary_body_name), "neptune") == 0 || cmp(lowercase(primary_body_name), "earth") == 0
        primary_body_name *= "_barycenter"
    end

    if cmp(lowercase(n_body_name), "mars") == 0 || cmp(lowercase(n_body_name), "jupiter") == 0 || cmp(lowercase(n_body_name), "saturn") == 0 || cmp(lowercase(n_body_name), "uranus") == 0 || cmp(lowercase(n_body_name), "neptune") == 0 || cmp(lowercase(n_body_name), "earth") == 0
        n_body_name *= "_barycenter"
    end

    pos_primary_k = p.J2000_to_pci*spkpos(n_body_name, et, "J2000", "none", primary_body_name)[1] * 1e3
    pos_spacecraft_k = pos_primary_k - pos_ii
    pos_spacecraft_k_mag = norm(pos_spacecraft_k)

    g = n_body.μ * ((pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3))
    return g
end

function eclipse_area_calc(r_sat::Vector{Float64}, r_sun::Vector{Float64}, A::Float64, rp::Float64)
    """
    Calculate the exposed area of the satellite. Translated from Python to Julia.

    Parameters 
    ----------
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    r_sun : Vectro{Float64}
        Position vector of the Sun relative to the planet.
    A : Float64
        Area of the satellite.
    rp : Float64
        Radius of the planet.


    Returns
    -------
    A_exp : Float64
        Exposed area of the satellite.
    
    """
    shadow = "none"
    rs = 6.9634e8 # Radius of the Sun in meters 
    RP = norm(r_sun) # Distance from Sun to the planet 
    alpha_umb = asin((rs - rp) / RP)  # Umbra angle
    alpha_pen = asin((rs + rp) / RP)  # Penumbra angle

    if dot(r_sun, r_sat) < 0 # if the angle is greater than 90 degrees, satellite is potentially in an eclipse
        # Compute the satellite's horizontal and vertical distances
        sigma = acos(dot(-r_sat, r_sun) / (norm(r_sat) * norm(r_sun)))
        sat_horiz = norm(r_sat) * cos(sigma)
        sat_vert = norm(r_sat) * sin(sigma)
        # Determine the eclipse conditions
        x = rp / sin(alpha_pen)
        pen_vert = tan(alpha_pen) * (x + sat_horiz)

        if sat_vert <= pen_vert # if true, the satellite is in partial shadow(penumbra)
            shadow = "penumbra"
            y = rp / sin(alpha_umb)
            umb_vert = tan(alpha_umb) * (y - sat_horiz)
            if sat_vert <= umb_vert
                shadow = "umbra"
            end
        end
    end

    if shadow == "none"
        A_exp = A
    elseif shadow == "penumbra"
        A_exp = A * (1 - (1 - (sat_vert / pen_vert))^2)
    else 
        A_exp = 0
    end
    
    return A_exp
end

function srp(p, p_srp_unscaled::Float64, cR::Float64, A_sat::Float64, m::Float64, r_sat::Vector{Float64}, time_et::Float64)
    """
    Calculate acceleration due to solar radiation pressure. 

    Parameters 
    ----------
    p : Struct
        Contains planetary parameters, including equatorial radius.
    p_srp_unscaled : Float64
        Solar radiation pressure at 1 AU.
    cR : Float64 
        Coefficient of reflectivity.
    A_sat : Float64
        Area of the satellite
    m : Float64
        Mass of the satellite
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    time_et : String
        Ephemeris time of the simulation, used for SPICE calculations.
    
    Returns
    -------
    srp_accel : Vector{Float64}
        Acceleration due to solar radiation pressure.
    """
    rp = p.Rp_e # Equatorial radius of the planet
    r_sun, ltime = spkpos("SUN", time_et, "ECLIPJ2000", "NONE", uppercase(p.name))
    r_sun = r_sun .* 1e3 # Convert from km to m

    A_exp = eclipse_area_calc(r_sat, r_sun, A_sat, rp)
    p_srp = p_srp_unscaled * norm(r_sun - r_sat) / (1.496e11) # Scale SRP for distance from Sun
    r_sat_to_sun = r_sun - r_sat
    r_sat_to_sun_mag = norm(r_sat_to_sun)

    srp_accel = -p_srp * cR * A_exp / m * (r_sat_to_sun / r_sat_to_sun_mag)
    return srp_accel
end

"""
    alf_IDR(x::Real, N::Integer, M::Integer)

Fully normalized Associate Legendre Functions (fnALFs) from degree 0 to N and order 0 to M,
evaluate at x ∈ [0, 1]. The ALFs are computed through Increasing Degree Recursion (IDR)
formulas that keep either the degree or the order fixed. All ALFs for M > N are zero by
definition.
"""
function fnALF_IDR!(A, x, N, M)
    # Allocate the output array.
    # TODO: this should be done in-place to avoid allocating a new array at each run.
    # A = zeros(N+1,M+1)

    # Precompute sqrt
    ξ = √(1 - x^2)
    # Initialize top left element
    A[1,1] = 1;

    # Sectorials
    if M > 0
        @inbounds for i=2:min(N+1,M+1)
            # For ease of notation
            n = i - 1

            # TODO: preallocate fn's for speed.
            fn = √( ((1 + δ(1,n))*(2n + 1)) / 2n )
            A[i,i] = fn * ξ * A[i-1,i-1]

        end
    end

    # Zonals and tesserals
    for j = 1:M+1
        @inbounds for i = j+1:N+1
            # For ease of notation
            m = j - 1
            n = i - 1
            #TODO: preallocate gnm, hnm for speed
            gnm = √(((2n + 1)*(2n-1)) / ((n+m)*(n-m)))

            if i == j + 1
                A[i,j] = gnm * x * A[i-1,j]
            else
                hnm = √(((2n + 1)*(n-m-1)*(n+m-1))/((2n-3)*(n+m)*(n-m)))
                A[i,j] = gnm * x * A[i-1,j] - hnm * A[i-2,j]
            end
            
        end      
    end

    return A
end

"""
    gradU_sph(rVec_sph, μ, RE, Clm, Slm, L, M)

Gradient of the geopotential in spherical coordinates. The inputs are the state vector in an
Earth-centred, Earth-fixed frame (ECEF), `rVec_sph`, the gravitational parameter `μ`, the
arrays of normalized spherical harmonic coefficients `Clm, Slm`, and the desired truncation
degree `L` and truncation order `M`.

Use stable, increasing-degree recursions to compute fully normalized associated Legendre
functions, and recursions for the trigonometric functions of the longitude.

The realization of the ECEF frame depends on the geopotential model.
"""
function gradU_sph!(rVec_sph, μ, RE, Clm, Slm, L, M, A)
    # Unpack spherical coordinates
    r, φ, λ = rVec_sph[1:3]

    # fnALFs until degree L and order M+1 (one more order is needed for ∂U/∂φ) 
    sinφ = sin(φ)
    P_LM = fnALF_IDR!(A, sinφ, L, M+1)

    # Recursions for trig functions of λ, φ
    if M == 0
        sinmλ = [0.0]
        cosmλ = [1.0]
    elseif M == 1
        sinmλ = [0.0; sin(λ)]
        cosmλ = [1.0; cos(λ)]
    elseif M > 1
        sinmλ = [0.0; sin(λ); zeros(M-1)]
        cosmλ = [1.0; cos(λ); zeros(M-1)]
    end
    @inbounds for j = 3:M+1
        sinmλ[j] = 2cosmλ[2] * sinmλ[j-1] - sinmλ[j-2]
        cosmλ[j] = 2cosmλ[2] * cosmλ[j-1] - cosmλ[j-2]
    end
    tanφ = tan(φ)

    # ∇U
    # TODO: reformulate with lumped coefficients
    ∂U∂r, ∂U∂φ, ∂U∂λ = (0.0, 0.0, 0.0)
    REoverR = RE/r; REoverRl = REoverR^(L+1)
    # ∂U∂r_list = zeros(L)
    # ∂U∂φ_list = zeros(L)
    # ∂U∂λ_list = zeros(L)
    @inbounds for l = L:-1:2
        i = l + 1
        ∂U∂r_l, ∂U∂φ_l, ∂U∂λ_l = (0.0, 0.0, 0.0)
        M_max = min(M, l)
        # ∂U∂r_l_list = zeros(M_max+1)
        # ∂U∂φ_l_list = zeros(M_max+1)
        # ∂U∂λ_l_list = zeros(M_max+1)
        @inbounds @simd for m = M_max:-1:0
            j = m + 1
            ∂U∂r_l += P_LM[i,j] * (Clm[i,j] * cosmλ[j] + Slm[i,j] * sinmλ[j])
            
            ∂U∂λ_l += m * P_LM[i,j] * ( -Clm[i,j] * sinmλ[j] + Slm[i,j] * cosmλ[j] )
            
            Π_ratio = √( ((l+m+1)*(l-m))/(1 + δ(m,0)) )
            ∂U∂φ_l += 
            ( Π_ratio * P_LM[i,j+1] - m*tanφ*P_LM[i,j] ) * 
            (Clm[i,j] * cosmλ[j] + Slm[i,j] * sinmλ[j])

            # ∂U∂r_l_list[m+1] = P_LM[i,j] * (Clm[i,j] * cosmλ[j] + Slm[i,j] * sinmλ[j])
            
            # ∂U∂λ_l_list[m+1] = m * P_LM[i,j] * ( -Clm[i,j] * sinmλ[j] + Slm[i,j] * cosmλ[j] )
            
            # Π_ratio = √( ((l+m+1)*(l-m))/(1 + δ(m,0)) )
            # ∂U∂φ_l_list[m+1] = 
            # ( Π_ratio * P_LM[i,j+1] - m*tanφ*P_LM[i,j] ) * 
            # (Clm[i,j] * cosmλ[j] + Slm[i,j] * sinmλ[j])
        end
        # ∂U∂r_l = sum(∂U∂r_l_list)
        # ∂U∂φ_l = sum(∂U∂φ_l_list)
        # ∂U∂λ_l = sum(∂U∂λ_l_list)
        REoverRl /= REoverR
        ∂U∂r += REoverRl * (l+1) * ∂U∂r_l
        ∂U∂φ += REoverRl * ∂U∂φ_l
        ∂U∂λ += REoverRl * ∂U∂λ_l
        # ∂U∂r_list[l-1] = REoverRl * (l+1) * ∂U∂r_l
        # ∂U∂φ_list[l-1] = REoverRl * ∂U∂φ_l
        # ∂U∂λ_list[l-1] = REoverRl * ∂U∂λ_l
    end
    # ∂U∂r = sum(∂U∂r_list)
    # ∂U∂φ = sum(∂U∂φ_list)
    # ∂U∂λ = sum(∂U∂λ_list)
    # Multiply by outer powers of 1/r
    ∂U∂r *= -μ/r^2
    ∂U∂φ *= μ/r
    ∂U∂λ *= μ/r

    return [∂U∂r; ∂U∂φ; ∂U∂λ]
end

"""
    acc_NSG(rVec_cart::Vector, ∇U_sph::Vector)

Acceleration in Cartesian coordinates in the ECEF frame from the geopotential gradient in
spherical coordinates, `∇U_sph`, and the state vector in Cartesian coordinates in the ECEF
frame, `rVec_cart`.
"""
function acc_NSG(rVec_cart, ∇U_sph)
    # Unpack arguments
    x, y, z = rVec_cart
    ∂U∂r, ∂U∂φ, ∂U∂λ = ∇U_sph

    # Preliminary quantities
    r = √(x^2 + y^2 + z^2)
    r_eq = √(x^2 + y^2)
    # println("r_eq: ", r_eq)
    # println("r: ", r)
    # Acceleration - equatorial components
    a = zeros(3)
    α = 1.0/r * ∂U∂r - z/(r^2 * r_eq) * ∂U∂φ
    β = 1.0/r_eq^2 * ∂U∂λ
    # println("α: ", α)
    # println("β: ", β)
    a[1] = α*x - β*y
    a[2] = β*x + α*y

    # Acceleration - vertical component
    a[3] = z/r * ∂U∂r + r_eq/r^2 * ∂U∂φ
    # println("a: ", a)
    # println("r: ", r)
    # println("norm(a): ", norm(a))
    return a
end

function acc_gravity_pines!(rVec_cart::SVector{3, Float64}, Clm::Matrix{Float64}, Slm::Matrix{Float64}, L::Int64, M::Int64, μ::Float64, RE::Float64, planet)
    """
        Calculate the acceleration due to gravity using the Pines method.
        (S. Pines, “Uniform Representation of the Gravitational Potential and its Derivatives,” AIAA Journal,
        vol. 11, no. 11, pp. 1508–1511, 1973.))

        Parameters
        ----------
        rVec_cart : Vector{Float64}
            Position vector of the satellite in the ECEF frame.
        Clm : Matrix{Float64}
            Array with the cosine harmonics.
        Slm : Matrix{Float64}
            Array with the sine harmonics.
        L : Int64
            Maximum degree of the spherical harmonics.
        M : Int64
            Maximum order of the spherical harmonics.
        μ : Float64
            Gravitational parameter of the planet.
        RE : Float64    
            Equatorial radius of the planet.
        planet : config.Planet
            Planet object with the spherical harmonics and other parameters.
    """
    x, y, z = rVec_cart
    r = norm(rVec_cart)
    s = x/r
    t = y/r
    u = z/r
    

    VR01 = planet.VR01
    VR11 = planet.VR11
    N1 = planet.N1
    N2 = planet.N2
    A = planet.A
    R = planet.Re
    I = planet.Im
    sqrt_2 = sqrt(2)

    A[2, 1] = u*sqrt(3)
    # Fill the off diagonal elements of A
    @simd for n = 1:L+1
        i = n + 1
        A[i+1, i] = u*sqrt(2*n+3)*A[i, i]
    end
    # Fill the rest of A
    for m = 0:M+1
        j = m + 1
        for l = m+2:L+1
            i = l + 1
            A[i, j] = u*N1[i, j]*A[i-1, j] - N2[i, j]*A[i-2, j]
        end
        R[j] = m == 0 ? 1 : s*R[j-1] - t*I[j-1]
        I[j] = m == 0 ? 0 : s*I[j-1] + t*R[j-1]
    end

    ρ = RE/r
    ρ_np1 = -μ/r * ρ
    a1, a2, a3, a4 = (0.0, 0.0, 0.0, 0.0)
    sum1, sum2, sum3, sum4 = (0.0, 0.0, 0.0, 0.0)
    for l = 1:L
        i = l + 1
        ρ_np1 *= ρ
        sum1 = 0
        sum2 = 0
        sum3 = 0
        sum4 = 0
        @simd for m = 0:min(l, M)
            j = m + 1
            C = Clm[i, j]
            S = Slm[i, j]
            D =              (C*R[j]   + S*I[j])   * sqrt_2
            E = m == 0 ? 0 : (C*R[j-1] + S*I[j-1]) * sqrt_2
            F = m == 0 ? 0 : (S*R[j-1] - C*I[j-1]) * sqrt_2

            Avv00 = A[i, j]
            Avv01 = VR01[i, j]*A[i, j+1]
            Avv11 = VR11[i, j]*A[i+1, j+1]

            sum1 += m * Avv00 * E
            sum2 += m * Avv00 * F
            sum3 +=     Avv01 * D
            sum4 +=     Avv11 * D
        end
        rr = ρ_np1/RE
        a1 += rr * sum1
        a2 += rr * sum2
        a3 += rr * sum3
        a4 -= rr * sum4
    end
    g = -[a1 + s*a4; a2 + t*a4; a3 + u*a4]
    return g
end
# function acc_gravity_pines(rVec_cart, Clm, Slm, L, M, μ, RE)
#     """
#         Calculate the acceleration due to gravity using the Pines method.
#         (S. Pines, “Uniform Representation of the Gravitational Potential and its Derivatives,” AIAA Journal,
#         vol. 11, no. 11, pp. 1508–1511, 1973.))

#         Parameters
#         ----------
#         rVec_cart : Vector{Float64}
#             Position vector of the satellite in the ECEF frame.
#         latitude : Float64
#             Latitude of the satellite in radians.
#     """
#     # Define the dimensionless coordinates
#     x, y, z = rVec_cart
#     r = norm(rVec_cart)
#     s = x/r
#     t = y/r
#     u = z/r
#     println("s: ", s)
#     println("t: ", t)
#     println("u: ", u)
#     # Precompute R, I, and A using recurrence relations
#     R = zeros(L+1)
#     I = zeros(L+1)
#     A = zeros(L+2, L+2)
#     ρ = zeros(L+2)

#     # Compute ρ TODO: precompute this
#     ρ[1] = μ/r
#     Re_over_r = RE/r
#     for l = 2:L+2
#         ρ[l] = ρ[l-1]*Re_over_r
#     end

#     # Compute R and I
#     R[1] = 1
#     A[1,1] = 1

#     for m = 1:M-1
#         R[m+1] = s * R[m] - t * I[m]
#         I[m+1] = s * I[m] + t * R[m]
#     end


#     # Compute A
#     for l = 1:L+1
#         i = l + 1
#         for m = l:-1:0
#             j = m + 1
#             # TODO: precompute these
#             if m == l
#                 A[i, j] = 2*l-1
#                 # A[i, j] = √((2*l+1)*(2-δ(l))/(2*l)/(2-δ(l-1)))*A[i-1, i-1]
#             elseif m == l-1
#                 A[i, j] = u*A[i, i]
#             #     A[i, j] = u*√((2*l)*(2-δ(l-1))/(2-δ(l)))*A[i,i]
#             else
#                 # N1 = √((2*l+1)*(2*l-1)/(l+m)/(l-m))
#                 # N2 = √((l+m-1)*(2*l+1)*(l-m-1)/(2*l-3)/(l+m)/(l-m))
#                 # A[i, j] = N1*u*A[i-1, j] - N2*A[i-2, j]
#                 A[i, j] = 1/(l-m)*(u*A[i, j+1] - A[i-1, j+1])
#             end
#         end
#     end
#     A = LowerTriangular(A)
#     println("A: ", A)

#     # Compute a1-a4
#     a1, a2, a3, a4 = (0.0, 0.0, 0.0, 0.0)
#     for l = 0:L
#         i = l + 1
#         for m = 1:min(l-1, M)
#             j = m + 1
#             Nlm_lm1 = √((l-m)*(2-δ(m))*(l+m+1)/(2-δ(m+1)))
#             Nlm_l1m1 = √((l+m+2)*(l+m+1)*(2*l+1)*(2-δ(m))/(2*l+3)/(2-δ(m+1)))
#             D = Clm[i,j]*R[j] + Slm[i,j]*I[j]
#             a1 += ρ[i+1]/RE*m*A[i,j]*(Clm[i,j]*R[j-1] + Slm[i,j]*I[j-1])
#             a2 += ρ[i+1]/RE*m*A[i,j]*(Slm[i,j]*R[j-1] - Clm[i,j]*I[j-1])
#             a3 += ρ[i+1]/RE*Nlm_lm1*A[i,j+1]*D
#             a4 += ρ[i+1]/RE*Nlm_l1m1*A[i+1,j+1]*D
#         end
#     end

#     # Compute the acceleration
#     g = [a1 + s*a4; a2 + t*a4; a3 + u*a4]
#     return g
# end

