struct SolarRadiationPressureModel <: AbstractForceTorqueModel
    Cr::Float64 # Reflectivity coefficient
    A::Float64  # Cross-sectional area in m^2
    AU_m::Float64 # Astronomical unit in meters
    direct::Bool
    albedo::Bool
    ir::Bool
    planet_albedo::Float64
    planet_ir_flux_w_m2::Float64
end

function SolarRadiationPressureModel(
    Cr::Real,
    A::Real,
    AU_m::Real=149_597_870_700.0;
    direct::Bool=true,
    albedo::Bool=false,
    ir::Bool=false,
    planet_albedo::Real=0.3,
    planet_ir_flux_w_m2::Real=237.0,
)
    Cr_f = Float64(Cr)
    A_f = Float64(A)
    AU_f = Float64(AU_m)
    planet_albedo_f = Float64(planet_albedo)
    planet_ir_flux_f = Float64(planet_ir_flux_w_m2)
    Cr_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.Cr must be >= 0, got $Cr_f."))
    A_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.A must be >= 0 m^2, got $A_f."))
    AU_f > 0.0 || throw(ArgumentError("SolarRadiationPressureModel.AU_m must be > 0 m, got $AU_f."))
    0.0 <= planet_albedo_f <= 1.0 || throw(ArgumentError("SolarRadiationPressureModel.planet_albedo must be in [0, 1], got $planet_albedo_f."))
    planet_ir_flux_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.planet_ir_flux_w_m2 must be >= 0, got $planet_ir_flux_f."))
    return SolarRadiationPressureModel(Cr_f, A_f, AU_f, direct, albedo, ir, planet_albedo_f, planet_ir_flux_f)
end
