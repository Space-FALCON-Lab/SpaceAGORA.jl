module LegacyModelCodes

export LegacyGravityModelCode, LegacyDensityModelCode, LegacyAerodynamicModelCode, LegacyThermalModelCode, LegacyThrustControlCode
export LegacyGravityConstant, LegacyGravityInverseSquared, LegacyGravityInverseSquaredJ2, LegacyGravityGRAM
export LegacyDensityConstant, LegacyDensityExponential, LegacyDensityNoDensity, LegacyDensityGRAM, LegacyDensityNRLMSISE
export LegacyAeroCdClConstant, LegacyAeroDiffusive, LegacyAeroNoBallisticAxial
export LegacyThermalConvectiveRadiative, LegacyThermalMaxwellian
export LegacyThrustNone, LegacyThrustAerobrakingManeuver, LegacyThrustDragPassageFiring

@enum LegacyGravityModelCode::Int8 begin
    LegacyGravityConstant = 0
    LegacyGravityInverseSquared = 1
    LegacyGravityInverseSquaredJ2 = 2
    LegacyGravityGRAM = 3
end

@enum LegacyDensityModelCode::Int8 begin
    LegacyDensityConstant = 0
    LegacyDensityExponential = 1
    LegacyDensityNoDensity = 2
    LegacyDensityGRAM = 3
    LegacyDensityNRLMSISE = 4
end

@enum LegacyAerodynamicModelCode::Int8 begin
    LegacyAeroCdClConstant = 0
    LegacyAeroDiffusive = 1
    LegacyAeroNoBallisticAxial = 2
end

@enum LegacyThermalModelCode::Int8 begin
    LegacyThermalConvectiveRadiative = 1
    LegacyThermalMaxwellian = 2
end

@enum LegacyThrustControlCode::Int8 begin
    LegacyThrustNone = 0
    LegacyThrustAerobrakingManeuver = 1
    LegacyThrustDragPassageFiring = 2
end

@inline function _compat_enum_parse(::Type{T}, x::T) where {T <: Enum}
    return x
end
@inline function _compat_enum_parse(::Type{LegacyGravityModelCode}, x::Integer)
    x == 0 && return LegacyGravityConstant
    x == 1 && return LegacyGravityInverseSquared
    x == 2 && return LegacyGravityInverseSquaredJ2
    x == 3 && return LegacyGravityGRAM
    throw(ArgumentError("Invalid gravity model code $x."))
end
@inline function _compat_enum_parse(::Type{LegacyDensityModelCode}, x::Integer)
    x == 0 && return LegacyDensityConstant
    x == 1 && return LegacyDensityExponential
    x == 2 && return LegacyDensityNoDensity
    x == 3 && return LegacyDensityGRAM
    x == 4 && return LegacyDensityNRLMSISE
    throw(ArgumentError("Invalid density model code $x."))
end
@inline function _compat_enum_parse(::Type{LegacyAerodynamicModelCode}, x::Integer)
    x == 0 && return LegacyAeroCdClConstant
    x == 1 && return LegacyAeroDiffusive
    x == 2 && return LegacyAeroNoBallisticAxial
    throw(ArgumentError("Invalid aerodynamic model code $x."))
end
@inline function _compat_enum_parse(::Type{LegacyThermalModelCode}, x::Integer)
    x == 1 && return LegacyThermalConvectiveRadiative
    x == 2 && return LegacyThermalMaxwellian
    throw(ArgumentError("Invalid thermal model code $x."))
end
@inline function _compat_enum_parse(::Type{LegacyThrustControlCode}, x::Integer)
    x == 0 && return LegacyThrustNone
    x == 1 && return LegacyThrustAerobrakingManeuver
    x == 2 && return LegacyThrustDragPassageFiring
    throw(ArgumentError("Invalid thrust control code $x."))
end

# Keep legacy numeric comparisons (`ip.gm == 1`) working while fields are typed enums.
@inline Base.:(==)(lhs::LegacyGravityModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
@inline Base.:(==)(lhs::Integer, rhs::LegacyGravityModelCode) = rhs == lhs
@inline Base.:(==)(lhs::LegacyDensityModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
@inline Base.:(==)(lhs::Integer, rhs::LegacyDensityModelCode) = rhs == lhs
@inline Base.:(==)(lhs::LegacyAerodynamicModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
@inline Base.:(==)(lhs::Integer, rhs::LegacyAerodynamicModelCode) = rhs == lhs
@inline Base.:(==)(lhs::LegacyThermalModelCode, rhs::Integer) = Int(lhs) == Int(rhs)
@inline Base.:(==)(lhs::Integer, rhs::LegacyThermalModelCode) = rhs == lhs
@inline Base.:(==)(lhs::LegacyThrustControlCode, rhs::Integer) = Int(lhs) == Int(rhs)
@inline Base.:(==)(lhs::Integer, rhs::LegacyThrustControlCode) = rhs == lhs

end # module LegacyModelCodes
