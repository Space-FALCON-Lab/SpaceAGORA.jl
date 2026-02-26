!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!! \ingroup F_Mars
MODULE GramStructs

  USE ISO_C_BINDING
  TYPE, BIND(C) :: GramTime_F
    INTEGER(C_INT) :: year      = 0      !< A four digit year
    INTEGER(C_INT) :: month     = 0      !< 1 through 12
    INTEGER(C_INT) :: day       = 0      !< 1 through 31
    INTEGER(C_INT) :: hour      = 0      !< 0 through 23
    INTEGER(C_INT) :: minute    = 0      !< 0 through 59
    REAL(C_DOUBLE) :: seconds   = 0.0D0  !< 0 through 60, may contain fractions
    INTEGER(C_INT) :: scale     = 0      !< 0 = TDT, 1 = UTC, 2 = TDB
    INTEGER(C_INT) :: frame     = 0      !< 0 = PET, 1 = ERT
  END TYPE

  TYPE, BIND(C) :: Position_F
    REAL(C_DOUBLE) :: height         = 0.0D0  !< \copydoc GRAM::Position::height
    REAL(C_DOUBLE) :: latitude       = 0.0D0  !< \copydoc GRAM::Position::latitude
    REAL(C_DOUBLE) :: longitude      = 0.0D0  !< \copydoc GRAM::Position::longitude
    REAL(C_DOUBLE) :: elapsedTime    = 0.0D0  !< \copydoc GRAM::Position::elapsedTime
    INTEGER(C_INT) :: isPlanetoCentric = 1    !< 1 if planeto-centric, 0 if planeto-graphic

    REAL(C_DOUBLE) :: totalRadius    = 0.0D0 !< \copydoc GRAM::Position::totalRadius
    REAL(C_DOUBLE) :: latitudeRadius = 0.0D0 !< \copydoc GRAM::Position::latitudeRadius
    REAL(C_DOUBLE) :: surfaceHeight  = 0.0D0 !< \copydoc GRAM::Position::surfaceHeight
    REAL(C_DOUBLE) :: gravity        = 0.0D0 !< \copydoc GRAM::Position::gravity
  END TYPE

  TYPE, BIND(C) :: EphemerisState_F
    REAL(C_DOUBLE) :: solarTime         = 0.0D0 !< \copydoc GRAM::EphemerisState::solarTime
    REAL(C_DOUBLE) :: longitudeSun      = 0.0D0 !< \copydoc GRAM::EphemerisState::longitudeSun
    REAL(C_DOUBLE) :: subsolarLatitude  = 0.0D0 !< \copydoc GRAM::EphemerisState::subsolarLatitude
    REAL(C_DOUBLE) :: subsolarLongitude = 0.0D0 !< \copydoc GRAM::EphemerisState::subsolarLongitude
    REAL(C_DOUBLE) :: orbitalRadius     = 0.0D0 !< \copydoc GRAM::EphemerisState::orbitalRadius
    REAL(C_DOUBLE) :: oneWayLightTime   = 0.0D0 !< \copydoc GRAM::EphemerisState::oneWayLightTime
    REAL(C_DOUBLE) :: solarZenithAngle  = 0.0D0 !< \copydoc GRAM::EphemerisState::solarZenithAngle
    REAL(C_DOUBLE) :: secondsPerSol     = 0.0D0 !< \copydoc GRAM::EphemerisState::secondsPerSol
  END TYPE

  TYPE, BIND(C) :: DynamicsState_F
    REAL(C_DOUBLE) :: temperature           = 0.0D0 !< \copydoc GRAM::AtmosphereState::temperature
    REAL(C_DOUBLE) :: pressure              = 0.0D0 !< \copydoc GRAM::AtmosphereState::pressure
    REAL(C_DOUBLE) :: density               = 0.0D0 !< \copydoc GRAM::AtmosphereState::density
    REAL(C_DOUBLE) :: pressureScaleHeight   = 0.0D0 !< \copydoc GRAM::AtmosphereState::pressureScaleHeight
    REAL(C_DOUBLE) :: densityScaleHeight    = 0.0D0 !< \copydoc GRAM::AtmosphereState::densityScaleHeight
    REAL(C_DOUBLE) :: speedOfSound          = 0.0D0 !< \copydoc GRAM::AtmosphereState::speedOfSound
    REAL(C_DOUBLE) :: referenceDensity      = 0.0D0 !< \copydoc GRAM::AtmosphereState::referenceDensity
    REAL(C_DOUBLE) :: referenceTemperature  = 0.0D0 !< \copydoc GRAM::AtmosphereState::referenceTemperature
    REAL(C_DOUBLE) :: referencePressure     = 0.0D0 !< \copydoc GRAM::AtmosphereState::referencePressure
    REAL(C_DOUBLE) :: sigmaLevel            = 0.0D0 !< \copydoc GRAM::AtmosphereState::sigmaLevel
    REAL(C_DOUBLE) :: pressureAtSurface     = 0.0D0 !< \copydoc GRAM::AtmosphereState::pressureAtSurface
    REAL(C_DOUBLE) :: pressureAltitude      = 0.0D0 !< \copydoc GRAM::AtmosphereState::pressureAltitude
  END TYPE

  TYPE, BIND(C) :: DensityState_F
    REAL(C_DOUBLE) :: density                   = 0.0D0 !< \copydoc GRAM::AtmosphereState::density
    REAL(C_DOUBLE) :: lowDensity                = 0.0D0 !< \copydoc GRAM::AtmosphereState::lowDensity
    REAL(C_DOUBLE) :: highDensity               = 0.0D0 !< \copydoc GRAM::AtmosphereState::highDensity
    REAL(C_DOUBLE) :: densityPerturbation       = 0.0D0 !< \copydoc GRAM::AtmosphereState::densityPerturbation
    REAL(C_DOUBLE) :: perturbedDensity          = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedDensity
    REAL(C_DOUBLE) :: densityStandardDeviation  = 0.0D0 !< \copydoc GRAM::AtmosphereState::densityStandardDeviation
    REAL(C_DOUBLE) :: perturbedSpeedOfSound     = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedSpeedOfSound
    REAL(C_DOUBLE) :: relativeStepSize          = 0.0D0 !< \copydoc GRAM::AtmosphereState::relativeStepSize
    REAL(C_DOUBLE) :: referenceDensity          = 0.0D0 !< \copydoc GRAM::AtmosphereState::referenceDensity
    REAL(C_DOUBLE) :: densityDeviation          = 0.0D0 !< \copydoc GRAM::AtmosphereState::densityDeviation
    REAL(C_DOUBLE) :: lowDensityDeviation       = 0.0D0 !< \copydoc GRAM::AtmosphereState::lowDensityDeviation
    REAL(C_DOUBLE) :: highDensityDeviation      = 0.0D0 !< \copydoc GRAM::AtmosphereState::highDensityDeviation
    REAL(C_DOUBLE) :: perturbedDensityDeviation = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedDensityDeviation
    INTEGER(C_INT) :: updateStatus              = 0     !< \copydoc GRAM::AtmosphereState::updateStatus
  END TYPE

  TYPE, BIND(C) :: WindsState_F
    REAL(C_DOUBLE) :: ewWind                    = 0.0D0 !< \copydoc GRAM::AtmosphereState::ewWind
    REAL(C_DOUBLE) :: nsWind                    = 0.0D0 !< \copydoc GRAM::AtmosphereState::nsWind
    REAL(C_DOUBLE) :: verticalWind              = 0.0D0 !< \copydoc GRAM::AtmosphereState::verticalWind
    REAL(C_DOUBLE) :: ewWindPerturbation        = 0.0D0 !< \copydoc GRAM::AtmosphereState::ewWindPerturbation
    REAL(C_DOUBLE) :: nsWindPerturbation        = 0.0D0 !< \copydoc GRAM::AtmosphereState::nsWindPerturbation
    REAL(C_DOUBLE) :: verticalWindPerturbation  = 0.0D0 !< \copydoc GRAM::AtmosphereState::verticalWindPerturbation
    REAL(C_DOUBLE) :: perturbedEWWind           = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedEWWind
    REAL(C_DOUBLE) :: perturbedNSWind           = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedNSWind
    REAL(C_DOUBLE) :: perturbedVerticalWind     = 0.0D0 !< \copydoc GRAM::AtmosphereState::perturbedVerticalWind
    REAL(C_DOUBLE) :: ewStandardDeviation       = 0.0D0 !< \copydoc GRAM::AtmosphereState::ewStandardDeviation
    REAL(C_DOUBLE) :: nsStandardDeviation       = 0.0D0 !< \copydoc GRAM::AtmosphereState::nsStandardDeviation
    REAL(C_DOUBLE) :: verticalStandardDeviation = 0.0D0 !< \copydoc GRAM::AtmosphereState::verticalStandardDeviation
  END TYPE

  TYPE, BIND(C) :: GasesState_F
    REAL(C_DOUBLE) :: totalNumberDensity      = 0.0D0 !< \copydoc GRAM::AtmosphereState::totalNumberDensity
    REAL(C_DOUBLE) :: averageMolecularWeight  = 0.0D0 !< \copydoc GRAM::AtmosphereState::averageMolecularWeight
    REAL(C_DOUBLE) :: compressibilityFactor   = 0.0D0 !< \copydoc GRAM::AtmosphereState::compressibilityFactor
    REAL(C_DOUBLE) :: specificGasConstant     = 0.0D0 !< \copydoc GRAM::AtmosphereState::specificGasConstant
    REAL(C_DOUBLE) :: specificHeatRatio       = 0.0D0 !< \copydoc GRAM::AtmosphereState::specificHeatRatio
  END TYPE

  TYPE, BIND(C) :: ConstituentGas_F
    REAL(C_DOUBLE) :: moleFraction           = 0.0D0 !< \copydoc GRAM::ConstituentGas::moleFraction
    REAL(C_DOUBLE) :: massFraction           = 0.0D0 !< \copydoc GRAM::ConstituentGas::massFraction
    REAL(C_DOUBLE) :: numberDensity          = 0.0D0 !< \copydoc GRAM::ConstituentGas::numberDensity
    REAL(C_DOUBLE) :: averageMolecularWeight = 0.0D0 !< \copydoc GRAM::ConstituentGas::averageMolecularWeight
    REAL(C_DOUBLE) :: specificHeatCapacity   = 0.0D0 !< \copydoc GRAM::ConstituentGas::specificHeatCapacity
  END TYPE

  TYPE, BIND(C) :: PerturbationState_F
    REAL(C_DOUBLE) :: densityRandomNumber      = 0.0D0 !< \copydoc GRAM::PerturbationState::densityRandomNumber
    REAL(C_DOUBLE) :: ewWindRandomNumber       = 0.0D0 !< \copydoc GRAM::PerturbationState::ewWindRandomNumber
    REAL(C_DOUBLE) :: nsWindRandomNumber       = 0.0D0 !< \copydoc GRAM::PerturbationState::nsWindRandomNumber
    REAL(C_DOUBLE) :: verticalWindRandomNumber = 0.0D0 !< \copydoc GRAM::PerturbationState::verticalWindRandomNumber
  END TYPE

END MODULE
