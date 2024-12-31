!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Earth The FORTRAN Interface for Earth
!! @{
!! \brief A listing of the FORTRAN interface functions for EarthGRAM.
!!
!! The FORTRAN interface for Earth is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE EarthGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Earth/examples/Earth_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: earth
!! earth = createAtmosphere_E(dataPath)
!! call setPosition_E(earth, &pos)
!! call update_E(earth)
!! 
!! // In C++
!! EarthAtmosphere earth(dataPath);
!! earth.setPosition(pos);
!! earth.update(); 
!! \endcode
!! \example Earth/examples/Earth_F.f90

MODULE EarthStructs

  USE ISO_C_BINDING

  TYPE, BIND(C) :: EarthState_F
    REAL(C_DOUBLE) :: perturbedTemperature;       !< \copydoc GRAM::EarthAtmosphereState::perturbedTemperature
    REAL(C_DOUBLE) :: temperaturePerturbation;    !< \copydoc GRAM::EarthAtmosphereState::temperaturePerturbation
    REAL(C_DOUBLE) :: temperatureStandardDeviation;!< \copydoc GRAM::AtmosphereState::temperatureStandardDeviation
    REAL(C_DOUBLE) :: perturbedPressure;          !< \copydoc GRAM::EarthAtmosphereState::perturbedPressure
    REAL(C_DOUBLE) :: pressurePerturbation;       !< \copydoc GRAM::EarthAtmosphereState::pressurePerturbation
    REAL(C_DOUBLE) :: pressureStandardDeviation;  !< \copydoc GRAM::AtmosphereState::pressureStandardDeviation
    
    REAL(C_DOUBLE) :: vaporPressure;              !< \copydoc GRAM::EarthAtmosphereState::vaporPressure
    REAL(C_DOUBLE) :: vaporPressureSD;            !< \copydoc GRAM::EarthAtmosphereState::vaporPressureSD
    REAL(C_DOUBLE) :: vaporDensity;               !< \copydoc GRAM::EarthAtmosphereState::vaporDensity
    REAL(C_DOUBLE) :: vaporDensitySD;             !< \copydoc GRAM::EarthAtmosphereState::vaporDensitySD
    REAL(C_DOUBLE) :: dewPoint;                   !< \copydoc GRAM::EarthAtmosphereState::dewPoint
    REAL(C_DOUBLE) :: dewPointSD;                 !< \copydoc GRAM::EarthAtmosphereState::dewPointSD
    REAL(C_DOUBLE) :: relativeHumidity;           !< \copydoc GRAM::EarthAtmosphereState::relativeHumidity
    REAL(C_DOUBLE) :: relativeHumiditySD;         !< \copydoc GRAM::EarthAtmosphereState::relativeHumiditySD
    
    REAL(C_DOUBLE) :: geodeticLatitude;           !< \copydoc GRAM::EarthAtmosphereState::geodeticLatitude
    REAL(C_DOUBLE) :: rraWeight;                  !< \copydoc GRAM::EarthAtmosphereState::rraWeight
    CHARACTER(C_CHAR) :: rraSiteName(6);          !< \copydoc GRAM::EarthAtmosphereState::rraSiteName
    
    REAL(C_DOUBLE) :: windSpeed;                  !< \copydoc GRAM::EarthAtmosphereState::windSpeed
    REAL(C_DOUBLE) :: windSpeedStandardDeviation; !< \copydoc GRAM::EarthAtmosphereState::windSpeedStandardDeviation
    REAL(C_DOUBLE) :: windCorrelation;            !< \copydoc GRAM::EarthAtmosphereState::windCorrelation
    INTEGER(C_INT) :: severityLevel;              !< \copydoc GRAM::EarthAtmosphereState::severityLevel
  END TYPE
  
  TYPE, BIND(C) :: EarthPerts_F
    REAL(C_DOUBLE) :: presPertSmall;              !< \copydoc GRAM::EarthAtmosphereState::presPertSmall
    REAL(C_DOUBLE) :: densPertSmall;              !< \copydoc GRAM::EarthAtmosphereState::densPertSmall
    REAL(C_DOUBLE) :: tempPertSmall;              !< \copydoc GRAM::EarthAtmosphereState::tempPertSmall
    REAL(C_DOUBLE) :: ewWindPertSmall;            !< \copydoc GRAM::EarthAtmosphereState::ewWindPertSmall
    REAL(C_DOUBLE) :: nsWindPertSmall;            !< \copydoc GRAM::EarthAtmosphereState::nsWindPertSmall

    REAL(C_DOUBLE) :: presSDSmall;                !< \copydoc GRAM::EarthAtmosphereState::presSDSmall
    REAL(C_DOUBLE) :: densSDSmall;                !< \copydoc GRAM::EarthAtmosphereState::densSDSmall
    REAL(C_DOUBLE) :: tempSDSmall;                !< \copydoc GRAM::EarthAtmosphereState::tempSDSmall
    REAL(C_DOUBLE) :: ewWindSDSmall;              !< \copydoc GRAM::EarthAtmosphereState::ewWindSDSmall
    REAL(C_DOUBLE) :: nsWindSDSmall;              !< \copydoc GRAM::EarthAtmosphereState::nsWindSDSmall

    REAL(C_DOUBLE) :: presPertLarge;              !< \copydoc GRAM::EarthAtmosphereState::presPertLarge
    REAL(C_DOUBLE) :: densPertLarge;              !< \copydoc GRAM::EarthAtmosphereState::densPertLarge
    REAL(C_DOUBLE) :: tempPertLarge;              !< \copydoc GRAM::EarthAtmosphereState::tempPertLarge
    REAL(C_DOUBLE) :: ewWindPertLarge;            !< \copydoc GRAM::EarthAtmosphereState::ewWindPertLarge
    REAL(C_DOUBLE) :: nsWindPertLarge;            !< \copydoc GRAM::EarthAtmosphereState::nsWindPertLarge

    REAL(C_DOUBLE) :: presSDLarge;                !< \copydoc GRAM::EarthAtmosphereState::presSDLarge
    REAL(C_DOUBLE) :: densSDLarge;                !< \copydoc GRAM::EarthAtmosphereState::densSDLarge
    REAL(C_DOUBLE) :: tempSDLarge;                !< \copydoc GRAM::EarthAtmosphereState::tempSDLarge
    REAL(C_DOUBLE) :: ewWindSDLarge;              !< \copydoc GRAM::EarthAtmosphereState::ewWindSDLarge
    REAL(C_DOUBLE) :: nsWindSDLarge;              !< \copydoc GRAM::EarthAtmosphereState::nsWindSDLarge
  END TYPE

  TYPE, BIND(C) :: EarthSurface_F
    REAL(C_DOUBLE) :: windSpeedAtSurface;         !< \copydoc GRAM::EarthAtmosphereState::windSpeedAtSurface
    REAL(C_DOUBLE) :: windSpeedSDAtSurface;       !< \copydoc GRAM::EarthAtmosphereState::windSpeedSDAtSurface
    REAL(C_DOUBLE) :: temperatureAtSurface;       !< \copydoc GRAM::EarthAtmosphereState::temperatureAtSurface
    REAL(C_DOUBLE) :: temperatureSDAtSurface;     !< \copydoc GRAM::EarthAtmosphereState::temperatureSDAtSurface
    REAL(C_DOUBLE) :: pressureSDAtSurface;        !< \copydoc GRAM::EarthAtmosphereState::pressureSDAtSurface
    REAL(C_DOUBLE) :: densitySDAtSurface;         !< \copydoc GRAM::EarthAtmosphereState::densitySDAtSurface
    REAL(C_DOUBLE) :: densityAtSurface;           !< \copydoc GRAM::EarthAtmosphereState::densityAtSurface
    REAL(C_DOUBLE) :: ewWindAtSurface;            !< \copydoc GRAM::EarthAtmosphereState::ewWindAtSurface
    REAL(C_DOUBLE) :: nsWindAtSurface;            !< \copydoc GRAM::EarthAtmosphereState::nsWindAtSurface
    REAL(C_DOUBLE) :: ewWindSDAtSurface;          !< \copydoc GRAM::EarthAtmosphereState::ewWindSDAtSurface
    REAL(C_DOUBLE) :: nsWindSDAtSurface;          !< \copydoc GRAM::EarthAtmosphereState::nsWindSDAtSurface
    REAL(C_DOUBLE) :: windCorrelationAtSurface;   !< \copydoc GRAM::EarthAtmosphereState::windCorrelationAtSurface
  END TYPE

  TYPE, BIND(C) :: EarthBoundaryLayer_F
    INTEGER(C_INT) :: landCode;                      !< \copydoc GRAM::EarthAtmosphereState::landCode
    REAL(C_DOUBLE) :: surfaceRoughness;              !< \copydoc GRAM::EarthAtmosphereState::surfaceRoughness
    REAL(C_DOUBLE) :: netRadiationIndex;             !< \copydoc GRAM::EarthAtmosphereState::netRadiationIndex
    REAL(C_DOUBLE) :: stability;                     !< \copydoc GRAM::EarthAtmosphereState::stability
    REAL(C_DOUBLE) :: inverseLength;                 !< \copydoc GRAM::EarthAtmosphereState::inverseLength
    REAL(C_DOUBLE) :: frictionVelocity;              !< \copydoc GRAM::EarthAtmosphereState::frictionVelocity
    REAL(C_DOUBLE) :: BVFrequencySquare;             !< \copydoc GRAM::EarthAtmosphereState::BVFrequencySquare
    REAL(C_DOUBLE) :: boundaryLayerDepth;            !< \copydoc GRAM::EarthAtmosphereState::boundaryLayerDepth
    REAL(C_DOUBLE) :: neutralBoundaryLayerDepth;     !< \copydoc GRAM::EarthAtmosphereState::neutralBoundaryLayerDepth
    REAL(C_DOUBLE) :: unstableBLFactor;              !< \copydoc GRAM::EarthAtmosphereState::unstableBLFactor
    REAL(C_DOUBLE) :: sigmaRatio;                    !< \copydoc GRAM::EarthAtmosphereState::sigmaRatio
    REAL(C_DOUBLE) :: sigmaW;                        !< \copydoc GRAM::EarthAtmosphereState::sigmaW
    REAL(C_DOUBLE) :: metersAboveSurface;            !< \copydoc GRAM::EarthAtmosphereState::metersAboveSurface
    REAL(C_DOUBLE) :: perturbedWindSpeedAtSurface;   !< \copydoc GRAM::EarthAtmosphereState::perturbedWindSpeedAtSurface
  END TYPE
END MODULE

MODULE EarthGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_E()
    SUBROUTINE tryGetSpicePath_E(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_E")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc tryGetDataPaths_E()
    SUBROUTINE tryGetDataPaths_E(spicePath, dataPath, bufferSize) BIND(C, NAME="tryGetDataPaths_E")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: dataPath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_E()
    SUBROUTINE setSpiceLsk_E(lsk) BIND(C, NAME="setSpiceLsk_E")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_E()
    SUBROUTINE setSpicePck_E(pck) BIND(C, NAME="setSpicePck_E")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_E()
    SUBROUTINE setSpiceKernel_E(bsp) BIND(C, NAME="setSpiceKernel_E")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_E()
    SUBROUTINE initialize_E(spicePath) BIND(C, NAME="initialize_E")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_E()
    SUBROUTINE loadSpiceFile_E(fileName) BIND(C, NAME="loadSpiceFile_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_E()
    FUNCTION createAtmosphere_E(dataPath) BIND(C, NAME="createAtmosphere_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: dataPath(*)
      TYPE(C_PTR) :: createAtmosphere_E
    END FUNCTION

    !> \copydoc copyAtmosphere_E()
    FUNCTION copyAtmosphere_E(atmos) BIND(C, NAME="copyAtmosphere_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_E
    END FUNCTION

    !> \copydoc deleteAtmosphere_E()
    SUBROUTINE deleteAtmosphere_E(atmos) BIND(C, NAME="deleteAtmosphere_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setThermosphereModel_E()
    SUBROUTINE setThermosphereModel_E(atmos, model) BIND(C, NAME="setThermosphereModel_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: model
    END SUBROUTINE

    !> \copydoc setSurfaceRoughness_E()
    SUBROUTINE setSurfaceRoughness_E(atmos, z0) BIND(C, NAME="setSurfaceRoughness_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: z0
    END SUBROUTINE

    !> \copydoc setUseNCEP_E()
    SUBROUTINE setUseNCEP_E(atmos, useNCEP) BIND(C, NAME="setUseNCEP_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: useNCEP
    END SUBROUTINE

    !> \copydoc setNCEPPath_E()
    SUBROUTINE setNCEPPath_E(atmos, path) BIND(C, NAME="setNCEPPath_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: path(*)
    END SUBROUTINE

    !> \copydoc setNCEPParameters_E()
    SUBROUTINE setNCEPParameters_E(atmos, NCEPYear, NCEPHour) BIND(C, NAME="setNCEPParameters_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: NCEPYear
      INTEGER(C_INT), VALUE, INTENT(IN) :: NCEPHour
    END SUBROUTINE

    !> \copydoc setMERRA2Path_E()
    SUBROUTINE setMERRA2Path_E(atmos, path) BIND(C, NAME="setMERRA2Path_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: path(*)
    END SUBROUTINE

    !> \copydoc setMERRA2Parameters_E()
    SUBROUTINE setMERRA2Parameters_E(atmos, M2Hour, latMin, latMax, lonMin, lonMax) BIND(C, NAME="setMERRA2Parameters_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: M2Hour
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: latMin
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: latMax
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: lonMin
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: lonMax
    END SUBROUTINE

    !> \copydoc setRRAParameters_E()
    SUBROUTINE setRRAParameters_E(atmos, year, innerRadius, outerRadius) BIND(C, NAME="setRRAParameters_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: year
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
    END SUBROUTINE

    !> \copydoc setRRAPath_E()
    SUBROUTINE setRRAPath_E(atmos, path) BIND(C, NAME="setRRAPath_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: path(*)
    END SUBROUTINE

    !> \copydoc setRRASiteList_E()
    SUBROUTINE setRRASiteList_E(atmos, fileName) BIND(C, NAME="setRRASiteList_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc setUseRRA_E()
    SUBROUTINE setUseRRA_E(atmos, useFlag) BIND(C, NAME="setUseRRA_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: useFlag
    END SUBROUTINE

    !> \copydoc setSolarParameters_E()
    SUBROUTINE setSolarParameters_E(atmos, dailyF10, meanF10, ap) BIND(C, NAME="setSolarParameters_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dailyF10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: meanF10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ap
    END SUBROUTINE

    !> \copydoc setJB2008Parameters_E()
    SUBROUTINE setJB2008Parameters_E(atmos, dailyS10, meanS10, dailyXM10, meanXM10, dailyY10, meanY10, dstdtc) BIND(C, NAME="setJB2008Parameters_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dailyS10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: meanS10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dailyXM10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: meanXM10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dailyY10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: meanY10
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dstdtc
    END SUBROUTINE

    !> \copydoc setInitialPerturbations_E()
    SUBROUTINE setInitialPerturbations_E(atmos, densPert, tempPert, ewPert, nsPert, verticalPert) BIND(C, NAME="setInitialPerturbations_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densPert
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: tempPert
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewPert
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsPert
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: verticalPert
    END SUBROUTINE

    !> \copydoc setPatchiness_E()
    SUBROUTINE setPatchiness_E(atmos, usePatchiness) BIND(C, NAME="setPatchiness_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: usePatchiness
    END SUBROUTINE

    !> \copydoc setStartTime_E()
    SUBROUTINE setStartTime_E(atmos, time) BIND(C, NAME="setStartTime_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setSeed_E()
    SUBROUTINE setSeed_E(atmos, seed) BIND(C, NAME="setSeed_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setPerturbationScales_E()
    SUBROUTINE setPerturbationScales_E(atmos, randomScale, horizontalWindScale, verticalWindScale) BIND(C, NAME="setPerturbationScales_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: randomScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: horizontalWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: verticalWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_E()
    SUBROUTINE addAuxiliaryAtmosphere_E(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_E()
    SUBROUTINE setAuxiliaryValues_E(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_E()
    SUBROUTINE setPosition_E(atmos, position) BIND(C, NAME="setPosition_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_E()
    SUBROUTINE setDelta_E(atmos, position) BIND(C, NAME="setDelta_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_E()
    SUBROUTINE setPerturbationAction_E(atmos, action) BIND(C, NAME="setPerturbationFactors_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_E()
    SUBROUTINE setEphemerisState_E(atmos, state) BIND(C, NAME="setEphemerisState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_E()
    SUBROUTINE setEphemerisFastModeOn_E(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_E()
    SUBROUTINE setSubsolarUpdateTime_E(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_E()
    INTEGER(C_INT) FUNCTION update_E(atmos) BIND(C, NAME="update_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_E()
    FUNCTION getErrorMessage_E(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_E
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_E()
    SUBROUTINE getStartTime_E(atmos, time) BIND(C, NAME="getStartTime_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getVersionString_E()
    FUNCTION getVersionString_E(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_E")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_E
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_E()
    SUBROUTINE getPosition_E(atmos, position) BIND(C, NAME="getPosition_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_E()
    SUBROUTINE getDynamicsState_E(atmos, state) BIND(C, NAME="getDynamicsState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_E()
    SUBROUTINE getDensityState_E(atmos, state) BIND(C, NAME="getDensityState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_E()
    SUBROUTINE getWindsState_E(atmos, state) BIND(C, NAME="getWindsState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_E()
    SUBROUTINE getGasesState_E(atmos, state, argon, carbonDioxide, carbonMonoxide, dinitrogen, dioxygen, helium, hydrogen, methane, nitrogen, nitrousOxide, oxygen, ozone, water) BIND(C, NAME="getGasesState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: argon
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonDioxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonMonoxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: dinitrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: dioxygen
      TYPE(ConstituentGas_F), INTENT(OUT) :: helium
      TYPE(ConstituentGas_F), INTENT(OUT) :: hydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: methane
      TYPE(ConstituentGas_F), INTENT(OUT) :: nitrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: nitrousOxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: oxygen
      TYPE(ConstituentGas_F), INTENT(OUT) :: ozone
      TYPE(ConstituentGas_F), INTENT(OUT) :: water
    END SUBROUTINE

    !> \copydoc getEphemerisState_E()
    SUBROUTINE getEphemerisState_E(atmos, state) BIND(C, NAME="getEphemerisState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_E()
    SUBROUTINE getPerturbationState_E(atmos, state) BIND(C, NAME="getPerturbationState_E")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getEarthState_E()
    SUBROUTINE getEarthState_E(atmos, state) BIND(C, NAME="getEarthState_E")
      USE GramStructs
      USE EarthStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EarthState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getEarthPerts_E()
    SUBROUTINE getEarthPerts_E(atmos, state) BIND(C, NAME="getEarthPerts_E")
      USE GramStructs
      USE EarthStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EarthPerts_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getEarthSurface_E()
    SUBROUTINE getEarthSurface_E(atmos, state) BIND(C, NAME="getEarthSurface_E")
      USE GramStructs
      USE EarthStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EarthSurface_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getEarthBoundaryLayer_E()
    SUBROUTINE getEarthBoundaryLayer_E(atmos, state) BIND(C, NAME="getEarthBoundaryLayer_E")
      USE GramStructs
      USE EarthStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EarthBoundaryLayer_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


