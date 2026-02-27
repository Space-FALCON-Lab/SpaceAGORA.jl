const BODY_NO = 0
const BODY_VENUS = 1
const BODY_EARTH = 2
const BODY_MARS = 3
const BODY_JUPITER = 4
const BODY_SATURN = 5
const BODY_URANUS = 6
const BODY_NEPTUNE = 7
const BODY_TITAN = 8

const GAS_ARGON = 0
const GAS_HELIUM = 1
const GAS_HYDROGEN = 2
const GAS_DIHYDROGEN = 3
const GAS_NITROGEN = 4
const GAS_DINITROGEN = 5
const GAS_OXYGEN = 6
const GAS_DIOXYGEN = 7
const GAS_METHANE = 8
const GAS_CARBON_MONOXIDE = 9
const GAS_CARBON_DIOXIDE = 10
const GAS_OZONE = 11
const GAS_NITROUS_OXIDE = 12
const GAS_WATER = 13

struct GramTimeC
    year::Cint
    month::Cint
    day::Cint
    hour::Cint
    minute::Cint
    seconds::Cdouble
    scale::Cint
    frame::Cint
end

struct PositionC
    height::Cdouble
    latitude::Cdouble
    longitude::Cdouble
    elapsedTime::Cdouble
    isPlanetoCentric::Cint
    totalRadius::Cdouble
    latitudeRadius::Cdouble
    surfaceHeight::Cdouble
    gravity::Cdouble
end

struct DynamicsStateC
    temperature::Cdouble
    pressure::Cdouble
    density::Cdouble
    pressureScaleHeight::Cdouble
    densityScaleHeight::Cdouble
    speedOfSound::Cdouble
    referenceDensity::Cdouble
    referenceTemperature::Cdouble
    referencePressure::Cdouble
    sigmaLevel::Cdouble
    pressureAtSurface::Cdouble
    pressureAltitude::Cdouble
end

struct DensityStateC
    density::Cdouble
    lowDensity::Cdouble
    highDensity::Cdouble
    densityPerturbation::Cdouble
    perturbedDensity::Cdouble
    densityStandardDeviation::Cdouble
    perturbedSpeedOfSound::Cdouble
    relativeStepSize::Cdouble
    referenceDensity::Cdouble
    densityDeviation::Cdouble
    lowDensityDeviation::Cdouble
    highDensityDeviation::Cdouble
    perturbedDensityDeviation::Cdouble
    updateStatus::Cint
end

struct WindsStateC
    ewWind::Cdouble
    nsWind::Cdouble
    verticalWind::Cdouble
    ewWindPerturbation::Cdouble
    nsWindPerturbation::Cdouble
    verticalWindPerturbation::Cdouble
    perturbedEWWind::Cdouble
    perturbedNSWind::Cdouble
    perturbedVerticalWind::Cdouble
    ewStandardDeviation::Cdouble
    nsStandardDeviation::Cdouble
    verticalStandardDeviation::Cdouble
end

struct GasesStateC
    totalNumberDensity::Cdouble
    averageMolecularWeight::Cdouble
    compressibilityFactor::Cdouble
    specificGasConstant::Cdouble
    specificHeatRatio::Cdouble
end

struct ConstituentGasC
    moleFraction::Cdouble
    massFraction::Cdouble
    numberDensity::Cdouble
    averageMolecularWeight::Cdouble
    specificHeatCapacity::Cdouble
end

struct PerturbationStateC
    densityRandomNumber::Cdouble
    ewWindRandomNumber::Cdouble
    nsWindRandomNumber::Cdouble
    verticalWindRandomNumber::Cdouble
end

struct EphemerisStateC
    solarTime::Cdouble
    longitudeSun::Cdouble
    subsolarLatitude::Cdouble
    subsolarLongitude::Cdouble
    orbitalRadius::Cdouble
    oneWayLightTime::Cdouble
    solarZenithAngle::Cdouble
    secondsPerSol::Cdouble
end

struct TrajectoryPointC
    position::PositionC
    dynamics::DynamicsStateC
    winds::WindsStateC
    ephemeris::EphemerisStateC
end

struct MarsStateC
    planetoGraphicHeight::Cdouble
    planetoGraphicLatitude::Cdouble
    referenceHeight::Cdouble
    referenceRadius::Cdouble
    groundTemperature::Cdouble
    thermosphereBaseHeight::Cdouble
    thermosphereBaseTemperature::Cdouble
    exosphericTemperature::Cdouble
    f1PeakHeight::Cdouble
    albedo::Cdouble
    heightOffset::Cdouble
    localHeightOffset::Cdouble
    dustOpticalDepth::Cdouble
    dustColumnArealDensity::Cdouble
    dustMixingRatio::Cdouble
    dustMassDensity::Cdouble
    dustNumberDensity::Cdouble
    iceIsPresent::Cint
end

struct DailyDynamicsStateC
    temperatureDaily::Cdouble
    pressureDaily::Cdouble
    densityDaily::Cdouble
    ewWindDaily::Cdouble
    nsWindDaily::Cdouble
    densityMin::Cdouble
    densityMax::Cdouble
    temperatureMin::Cdouble
    temperatureMax::Cdouble
end

struct EarthStateC
    perturbedTemperature::Cdouble
    temperaturePerturbation::Cdouble
    temperatureStandardDeviation::Cdouble
    perturbedPressure::Cdouble
    pressurePerturbation::Cdouble
    pressureStandardDeviation::Cdouble
    vaporPressure::Cdouble
    vaporPressureSD::Cdouble
    vaporDensity::Cdouble
    vaporDensitySD::Cdouble
    dewPoint::Cdouble
    dewPointSD::Cdouble
    relativeHumidity::Cdouble
    relativeHumiditySD::Cdouble
    geodeticLatitude::Cdouble
    rraWeight::Cdouble
    rraSiteName::NTuple{6, Cchar}
    windSpeed::Cdouble
    windSpeedStandardDeviation::Cdouble
    windCorrelation::Cdouble
    severityLevel::Cint
end

struct EarthPertsC
    presPertSmall::Cdouble
    densPertSmall::Cdouble
    tempPertSmall::Cdouble
    ewWindPertSmall::Cdouble
    nsWindPertSmall::Cdouble
    presSDSmall::Cdouble
    densSDSmall::Cdouble
    tempSDSmall::Cdouble
    ewWindSDSmall::Cdouble
    nsWindSDSmall::Cdouble
    presPertLarge::Cdouble
    densPertLarge::Cdouble
    tempPertLarge::Cdouble
    ewWindPertLarge::Cdouble
    nsWindPertLarge::Cdouble
    presSDLarge::Cdouble
    densSDLarge::Cdouble
    tempSDLarge::Cdouble
    ewWindSDLarge::Cdouble
    nsWindSDLarge::Cdouble
end

struct EarthSurfaceC
    windSpeedAtSurface::Cdouble
    windSpeedSDAtSurface::Cdouble
    temperatureAtSurface::Cdouble
    temperatureSDAtSurface::Cdouble
    pressureSDAtSurface::Cdouble
    densitySDAtSurface::Cdouble
    densityAtSurface::Cdouble
    ewWindAtSurface::Cdouble
    nsWindAtSurface::Cdouble
    ewWindSDAtSurface::Cdouble
    nsWindSDAtSurface::Cdouble
    windCorrelationAtSurface::Cdouble
end

struct EarthBoundaryLayerC
    landCode::Cint
    surfaceRoughness::Cdouble
    netRadiationIndex::Cdouble
    stability::Cdouble
    inverseLength::Cdouble
    frictionVelocity::Cdouble
    BVFrequencySquare::Cdouble
    boundaryLayerDepth::Cdouble
    neutralBoundaryLayerDepth::Cdouble
    unstableBLFactor::Cdouble
    sigmaRatio::Cdouble
    sigmaW::Cdouble
    metersAboveSurface::Cdouble
    perturbedWindSpeedAtSurface::Cdouble
end

mutable struct Atmosphere
    ptr::Ptr{Cvoid}
    body::Cint
    function Atmosphere(ptr::Ptr{Cvoid}, body::Integer = BODY_NO)
        ptr == C_NULL && error("createAtmosphere_C returned NULL.")
        obj = new(ptr, Cint(body))
        finalizer(obj) do x
            if x.ptr != C_NULL
                ccall((:deleteAtmosphere_C, _libgram()), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        return obj
    end
end

const _LIBGRAM = Ref{String}("")
const _SPICE_PATH = Ref{String}("")
const _REGISTERED_BODIES = Set{Int}()
