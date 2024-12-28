!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Mars The FORTRAN Interface for Mars
!! @{
!! \brief A listing of the FORTRAN interface functions for MarsGRAM.
!!
!! The FORTRAN interface for Mars is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE MarsGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Mars/examples/Mars_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: mars
!! mars = createAtmosphere_M(dataPath)
!! call setPosition_M(mars, &pos)
!! call update_M(mars)
!! 
!! // In C++
!! MarsAtmosphere mars(dataPath);
!! mars.setPosition(pos);
!! mars.update(); 
!! \endcode
!! \example Mars/examples/Mars_F.f90

MODULE MarsStructs

  USE ISO_C_BINDING

  TYPE, BIND(C) :: DailyDynamicsState_F
    REAL(C_DOUBLE) :: temperatureDaily         !< \copydoc GRAM::MarsAtmosphereState::temperatureDaily
    REAL(C_DOUBLE) :: pressureDaily            !< \copydoc GRAM::MarsAtmosphereState::pressureDaily
    REAL(C_DOUBLE) :: densityDaily             !< \copydoc GRAM::MarsAtmosphereState::densityDaily
    REAL(C_DOUBLE) :: ewWindDaily              !< \copydoc GRAM::MarsAtmosphereState::ewWindDaily
    REAL(C_DOUBLE) :: nsWindDaily              !< \copydoc GRAM::MarsAtmosphereState::nsWindDaily
    REAL(C_DOUBLE) :: densityMin               !< \copydoc GRAM::MarsAtmosphereState::densityMin
    REAL(C_DOUBLE) :: densityMax               !< \copydoc GRAM::MarsAtmosphereState::densityMax
    REAL(C_DOUBLE) :: temperatureMin           !< \copydoc GRAM::MarsAtmosphereState::temperatureMin
    REAL(C_DOUBLE) :: temperatureMax           !< \copydoc GRAM::MarsAtmosphereState::temperatureMax
  END TYPE

  TYPE, BIND(C) :: MarsState_F
    REAL(C_DOUBLE) :: planetoGraphicHeight         !< \copydoc GRAM::MarsAtmosphereState::planetoGraphicHeight
    REAL(C_DOUBLE) :: planetoGraphicLatitude       !< \copydoc GRAM::MarsAtmosphereState::planetoGraphicLatitude
    REAL(C_DOUBLE) :: referenceHeight              !< \copydoc GRAM::MarsAtmosphereState::referenceHeight
    REAL(C_DOUBLE) :: referenceRadius              !< \copydoc GRAM::MarsAtmosphereState::referenceRadius
    REAL(C_DOUBLE) :: groundTemperature            !< \copydoc GRAM::MarsAtmosphereState::groundTemperature
    REAL(C_DOUBLE) :: thermosphereBaseHeight       !< \copydoc GRAM::MarsAtmosphereState::thermosphereBaseHeight
    REAL(C_DOUBLE) :: thermosphereBaseTemperature  !< \copydoc GRAM::MarsAtmosphereState::thermosphereBaseTemperature
    REAL(C_DOUBLE) :: exosphericTemperature        !< \copydoc GRAM::MarsAtmosphereState::exosphericTemperature
    REAL(C_DOUBLE) :: f1PeakHeight                 !< \copydoc GRAM::MarsAtmosphereState::f1PeakHeight
    REAL(C_DOUBLE) :: albedo                       !< \copydoc GRAM::MarsAtmosphereState::albedo
    REAL(C_DOUBLE) :: heightOffset                 !< \copydoc GRAM::MarsAtmosphereState::heightOffset
    REAL(C_DOUBLE) :: localHeightOffset            !< \copydoc GRAM::MarsAtmosphereState::localHeightOffset
    REAL(C_DOUBLE) :: dustOpticalDepth             !< \copydoc GRAM::MarsAtmosphereState::dustOpticalDepth
    REAL(C_DOUBLE) :: dustColumnArealDensity       !< \copydoc GRAM::MarsAtmosphereState::dustColumnArealDensity
    REAL(C_DOUBLE) :: dustMixingRatio              !< \copydoc GRAM::MarsAtmosphereState::dustMixingRatio
    REAL(C_DOUBLE) :: dustMassDensity              !< \copydoc GRAM::MarsAtmosphereState::dustMassDensity
    REAL(C_DOUBLE) :: dustNumberDensity            !< \copydoc GRAM::MarsAtmosphereState::dustNumberDensity
    INTEGER(C_INT) :: iceIsPresent                 !< \copydoc GRAM::MarsAtmosphereState::iceIsPresent
  END TYPE
END MODULE

MODULE MarsGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_M()
    SUBROUTINE tryGetSpicePath_M(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_M")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc tryGetDataPaths_M()
    SUBROUTINE tryGetDataPaths_M(spicePath, dataPath, bufferSize) BIND(C, NAME="tryGetDataPaths_M")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: dataPath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_M()
    SUBROUTINE setSpiceLsk_M(lsk) BIND(C, NAME="setSpiceLsk_M")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_M()
    SUBROUTINE setSpicePck_M(pck) BIND(C, NAME="setSpicePck_M")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_M()
    SUBROUTINE setSpiceKernel_M(bsp) BIND(C, NAME="setSpiceKernel_M")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_M()
    SUBROUTINE initialize_M(spicePath) BIND(C, NAME="initialize_M")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_M()
    SUBROUTINE loadSpiceFile_M(fileName) BIND(C, NAME="loadSpiceFile_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_M()
    FUNCTION createAtmosphere_M(dataPath) BIND(C, NAME="createAtmosphere_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: dataPath(*)
      TYPE(C_PTR) :: createAtmosphere_M
    END FUNCTION

    !> \copydoc copyAtmosphere_M()
    FUNCTION copyAtmosphere_M(atmos) BIND(C, NAME="copyAtmosphere_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_M
    END FUNCTION

    !> \copydoc deleteAtmosphere_M()
    SUBROUTINE deleteAtmosphere_M(atmos) BIND(C, NAME="deleteAtmosphere_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setPlanetaryRadii_M()
    SUBROUTINE setPlanetaryRadii_M(atmos, eqrad, prad) BIND(C, NAME="setPlanetaryRadii_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: eqrad
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: prad
    END SUBROUTINE

    !> \copydoc setMapYear_M()
    SUBROUTINE setMapYear_M(atmos, year) BIND(C, NAME="setMapYear_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: year
    END SUBROUTINE

    !> \copydoc setHeightOffsetModel_M()
    SUBROUTINE setHeightOffsetModel_M(atmos, model, hgtOffset) BIND(C, NAME="setHeightOffsetModel_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: model
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: hgtOffset
    END SUBROUTINE

    !> \copydoc setHeightAboveSurface_M()
    SUBROUTINE setHeightAboveSurface_M(atmos, heightAboveSurface) BIND(C, NAME="setHeightAboveSurface_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: heightAboveSurface
    END SUBROUTINE

    !> \copydoc setMGCMDustLevels_M()
    SUBROUTINE setMGCMDustLevels_M(atmos, constantDustLevel, minDustLevel, maxDustLevel) BIND(C, NAME="setMGCMDustLevels_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: constantDustLevel
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minDustLevel
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: maxDustLevel
    END SUBROUTINE

    !> \copydoc setF107_M()
    SUBROUTINE setF107_M(atmos, f107) BIND(C, NAME="setF107_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: f107
    END SUBROUTINE

    !> \copydoc setPerturbationWaveLengthScale_M()
    SUBROUTINE setPerturbationWaveLengthScale_M(atmos, wscale) BIND(C, NAME="setPerturbationWaveLengthScale_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: wscale
    END SUBROUTINE

    !> \copydoc setMOLAHeights_M()
    SUBROUTINE setMOLAHeights_M(atmos, isMola) BIND(C, NAME="setMOLAHeights_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: isMola
    END SUBROUTINE

    !> \copydoc setMinMax_M()
    SUBROUTINE setMinMax_M(atmos, minMax) BIND(C, NAME="setMinMax_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: minMax
    END SUBROUTINE

    !> \copydoc setDustStorm_M()
    SUBROUTINE setDustStorm_M(atmos, longitudeSun, duration, intensity, maxRadius, lat, lon) BIND(C, NAME="setDustStorm_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: longitudeSun
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: duration
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: intensity
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: maxRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: lat
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: lon
    END SUBROUTINE

    !> \copydoc setDustDensity_M()
    SUBROUTINE setDustDensity_M(atmos, nu, diameter, density) BIND(C, NAME="setDustDensity_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nu
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: diameter
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: density
    END SUBROUTINE

    !> \copydoc setExosphericTemperature_M()
    SUBROUTINE setExosphericTemperature_M(atmos, hoffset, factor) BIND(C, NAME="setExosphericTemperature_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: hoffset
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: factor
    END SUBROUTINE

    !> \copydoc setWaveDefaults_M()
    SUBROUTINE setWaveDefaults_M(atmos, date, scale, mean, a1, p1, r1, a2, p2, r2, a3, p3, r3) BIND(C, NAME="setWaveDefaults_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: date
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: scale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: mean
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: a1
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: p1
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: r1
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: a2
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: p2
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: r2
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: a3
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: p3
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: r3
    END SUBROUTINE

    !> \copydoc setWaveFile_M()
    SUBROUTINE setWaveFile_M(waveFile) BIND(C, NAME="setWaveFile_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: waveFile(*)
    END SUBROUTINE

    !> \copydoc setWindScales_M()
    SUBROUTINE setWindScales_M(atmos, meanWinds, boundaryLayerWinds) BIND(C, NAME="setWindScales_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: meanWinds
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: boundaryLayerWinds
    END SUBROUTINE


  !> \copydoc setStartTime_M()
    SUBROUTINE setStartTime_M(atmos, time) BIND(C, NAME="setStartTime_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setSeed_M()
    SUBROUTINE setSeed_M(atmos, seed) BIND(C, NAME="setSeed_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_M()
    SUBROUTINE setMinRelativeStepSize_M(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_M()
    SUBROUTINE setPerturbationScales_M(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_M()
    SUBROUTINE addAuxiliaryAtmosphere_M(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_M()
    SUBROUTINE setAuxiliaryValues_M(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_M()
    SUBROUTINE setPosition_M(atmos, position) BIND(C, NAME="setPosition_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_M()
    SUBROUTINE setDelta_M(atmos, position) BIND(C, NAME="setDelta_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_M()
    SUBROUTINE setPerturbationAction_M(atmos, action) BIND(C, NAME="setPerturbationFactors_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_M()
    SUBROUTINE setEphemerisState_M(atmos, state) BIND(C, NAME="setEphemerisState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_M()
    SUBROUTINE setEphemerisFastModeOn_M(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_M()
    SUBROUTINE setSubsolarUpdateTime_M(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_M()
    INTEGER(C_INT) FUNCTION update_M(atmos) BIND(C, NAME="update_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_M()
    FUNCTION getErrorMessage_M(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_M
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_M()
    SUBROUTINE getStartTime_M(atmos, time) BIND(C, NAME="getStartTime_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getVersionString_M()
    FUNCTION getVersionString_M(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_M")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_M
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_M()
    SUBROUTINE getPosition_M(atmos, position) BIND(C, NAME="getPosition_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_M()
    SUBROUTINE getDynamicsState_M(atmos, state) BIND(C, NAME="getDynamicsState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_M()
    SUBROUTINE getDensityState_M(atmos, state) BIND(C, NAME="getDensityState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_M()
    SUBROUTINE getWindsState_M(atmos, state) BIND(C, NAME="getWindsState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_M()
    SUBROUTINE getGasesState_M(atmos, state, argon, carbonDioxide, carbonMonoxide, dihydrogen, dinitrogen, dioxygen, helium, hydrogen, oxygen, water) BIND(C, NAME="getGasesState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: argon
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonDioxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonMonoxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: dihydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: dinitrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: dioxygen
      TYPE(ConstituentGas_F), INTENT(OUT) :: helium
      TYPE(ConstituentGas_F), INTENT(OUT) :: hydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: oxygen
      TYPE(ConstituentGas_F), INTENT(OUT) :: water
    END SUBROUTINE

    !> \copydoc getEphemerisState_M()
    SUBROUTINE getEphemerisState_M(atmos, state) BIND(C, NAME="getEphemerisState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_M()
    SUBROUTINE getPerturbationState_M(atmos, state) BIND(C, NAME="getPerturbationState_M")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDailyDynamicsState_M()
    SUBROUTINE getDailyDynamicsState_M(atmos, state) BIND(C, NAME="getDailyDynamicsState_M")
      USE GramStructs
      USE MarsStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DailyDynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getMarsState_M()
    SUBROUTINE getMarsState_M(atmos, state) BIND(C, NAME="getMarsState_M")
      USE GramStructs
      USE MarsStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(MarsState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


