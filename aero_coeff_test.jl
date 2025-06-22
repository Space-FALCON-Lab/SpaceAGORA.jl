using LinearAlgebra
using Plots
using SpecialFunctions
plotly()

function CA(S, α, β, σ, T, Tw, lx, ly, lz)
    """
    Calculate the aerodynamic coefficient CA for a blunt body.
    # Arguments
    - 'S': Surface area, m²
    - 'α': Angle of attack, rad
    - 'β': Angle of sideslip, rad
    - 'σ': Reflection coefficient
    - 'T': Temperature, K
    - 'Tw': Wall temperature, K
    - 'lx': Length in x-direction, m
    - 'ly': Length in y-direction, m
    - 'lz': Length in z-direction, m
    # Returns
    - 'CA': Aerodynamic coefficient CA
    """
    σN = σ
    σT = σ

    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)
    CA = ((2-σN)/(S√(π))*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
            (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
            σT*cosα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
            lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    return CA
end

function CS(S, α, β, σ, T, Tw, lx, ly, lz)
    """
    Calculate the aerodynamic coefficient CN for a blunt body.
    # Arguments
    - 'S': Surface area, m²
    - 'α': Angle of attack, rad
    - 'β': Angle of sideslip, rad
    - 'σ': Reflection coefficient
    - 'T': Temperature, K
    - 'Tw': Wall temperature, K
    - 'lx': Length in x-direction, m
    - 'ly': Length in y-direction, m
    - 'lz': Length in z-direction, m
    # Returns
    - 'CN': Aerodynamic coefficient CN
    """
    σN = σ
    σT = σ

    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)
    
    CS = lx/ly*(((2-σN)/(S*√(π))*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
            (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
            (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
            σT*sinβ*(1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    return CS
end

function CN(S, α, β, σ, T, Tw, lx, ly, lz)
    """
    Calculate the aerodynamic coefficient CN for a blunt body.
    # Arguments
    - 'S': Surface area, m²
    - 'α': Angle of attack, rad
    - 'β': Angle of sideslip, rad
    - 'σ': Reflection coefficient
    - 'T': Temperature, K
    - 'Tw': Wall temperature, K
    - 'lx': Length in x-direction, m
    - 'ly': Length in y-direction, m
    - 'lz': Length in z-direction, m
    # Returns
    - 'CN': Aerodynamic coefficient CN
    """
    σN = σ
    σT = σ

    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)

    CN = lx/lz*((((2-σN)/(S*√(π))*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
            (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
            (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
            σT*sinα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
            (1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))
    return CN
end

function normalcoefficient(S, aoa, sigma)
        CN = 1 / (S^2) * ((((2 - sigma) / sqrt(pi)) * S * sin(aoa) + sigma / 2) * exp(-(S * sin(aoa))^2) + ((2 - sigma) * ((S * sin(aoa))^2 + 0.5) + sigma / 2 * sqrt(pi) * (S * sin(aoa))) * (1 + erf(S * sin(aoa))))
        return CN
end

function axialcoefficient(S, aoa, sigma)
    CA = ((sigma * cos(aoa)) / (sqrt(pi) * S)) * (exp(-(S * sin(aoa))^2) + sqrt(pi) * (S * sin(aoa)) * (1 + erf(S * sin(aoa))))
    return CA
end

T = 973
V = 7500
Tw = 300
α = deg2rad(0.0)
β = deg2rad(0.0)
# σ_list = [0, 0.2, 0.4, 0.6, 0.8, 1]
σ_list = [0.9]
lx = 0.01
ly = 3.76/2
lz = 1.93
gas_constant = 287.05  # Specific gas constant for air, J/(kg·K)
S = V/√(2*gas_constant*T)

min_angle = -180
max_angle = 360
# min_angle = deg2rad(min_angle)
# max_angle = deg2rad(max_angle)
p = plot()
for σ in σ_list
    ca_list = Float64[]
    cn_list = Float64[]
    cs_list = Float64[]
    for β = min_angle:max_angle
        β_rad = deg2rad(β)
        ca = CA(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        cn = CN(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        cs = CS(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        push!(ca_list, ca)
        push!(cn_list, cn)
        push!(cs_list, cs)
    end
    plot!(p, [min_angle:max_angle], ca_list, label="σ = $σ, CA", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
    plot!(p, [min_angle:max_angle], cn_list, label="σ = $σ, CN", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
    plot!(p, [min_angle:max_angle], cs_list, label="σ = $σ, CS", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
end
display(p)

p = plot()
for σ in σ_list
    ca_list = Float64[]
    cn_list = Float64[]
    cs_list = Float64[]
    ca_old_list = Float64[]
    cn_old_list = Float64[]
    for α = min_angle:max_angle
        α_rad = deg2rad(α)
        ca = CA(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        cn = CN(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        cs = CS(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        ca_old = axialcoefficient(S, α_rad, σ)
        cn_old = normalcoefficient(S, α_rad, σ)
        push!(ca_list, ca)
        push!(cn_list, cn)
        push!(cs_list, cs)
        push!(ca_old_list, ca_old)
        push!(cn_old_list, cn_old)
    end
    plot!(p, [min_angle:max_angle], ca_old_list, linestyle=:dash, label="σ = $σ, CA (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    plot!(p, [min_angle:max_angle], cn_old_list, linestyle=:dash, label="σ = $σ, CN (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    plot!(p, [min_angle:max_angle], ca_list, label="σ = $σ, CA", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    plot!(p, [min_angle:max_angle], cn_list, label="σ = $σ, CN", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    # plot!(p, [min_angle:max_angle], cs_list, label="σ = $σ, CS", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
end
display(p)

p = plot()
for σ in σ_list
    cd_list = Float64[]
    cl_list = Float64[]
    cs_list = Float64[]
    for α in min_angle:max_angle
        α_rad = -deg2rad(α)
        ca = CA(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        cn = CN(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        cs = CS(S, α_rad, β, σ, T, Tw, lx, ly, lz)
        cd = cos(α_rad)*cos(β)*ca + sin(β)*cs - sin(α_rad)*cos(β)*cn
        cl = sin(α_rad)*ca + cos(α_rad)*cn
        cs = cos(α_rad)*sin(β)*ca - cos(β)*cs - sin(α_rad)*sin(β)*cn
        push!(cd_list, cd)
        push!(cl_list, cl)
        push!(cs_list, cs)
    end
    plot!(p, [min_angle:max_angle], cd_list, label="σ = $σ, CD", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    plot!(p, [min_angle:max_angle], cl_list, label="σ = $σ, CL", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
    plot!(p, [min_angle:max_angle], cs_list, label="σ = $σ, CS", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
end
display(p)

p = plot()
for σ in σ_list
    cd_list = Float64[]
    cl_list = Float64[]
    cs_list = Float64[]
    for β in min_angle:max_angle
        β_rad = deg2rad(β)
        ca = CA(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        cn = CN(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        cs = CS(S, α, β_rad, σ, T, Tw, lx, ly, lz)
        cd = cos(α)*cos(β_rad)*ca + sin(β_rad)*cs - sin(α)*cos(β_rad)*cn
        cl = sin(α)*ca + cos(α)*cn
        cs = cos(α)*sin(β_rad)*ca + cos(β_rad)*cs - sin(α)*sin(β_rad)*cn
        push!(cd_list, cd)
        push!(cl_list, cl)
        push!(cs_list, cs)
    end
    plot!(p, [min_angle:max_angle], cd_list, label="σ = $σ, CD", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
    plot!(p, [min_angle:max_angle], cl_list, label="σ = $σ, CL", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
    plot!(p, [min_angle:max_angle], cs_list, label="σ = $σ, CS", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
end
display(p)