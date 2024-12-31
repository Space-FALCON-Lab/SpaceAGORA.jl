!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Venus The FORTRAN Interface for Venus
!! @{
!! \brief A listing of the FORTRAN interface functions for VenusGRAM.
!!
!! The FORTRAN interface for Venus is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE VenusGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Venus/examples/Venus_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: venus
!! venus = createAtmosphere_V(dataPath)
!! call setPosition_V(venus, &pos)
!! call update_V(venus)
!! 
!! // In C++
!! VenusAtmosphere venus(dataPath);
!! venus.setPosition(pos);
!! venus.update(); 
!! \endcode
!! \example Venus/examples/Venus_F.f90

MODULE VenusGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_V()
    SUBROUTINE tryGetSpicePath_V(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_V")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_V()
    SUBROUTINE setSpiceLsk_V(lsk) BIND(C, NAME="setSpiceLsk_V")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_V()
    SUBROUTINE setSpicePck_V(pck) BIND(C, NAME="setSpicePck_V")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_V()
    SUBROUTINE setSpiceKernel_V(bsp) BIND(C, NAME="setSpiceKernel_V")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_V()
    SUBROUTINE initialize_V(spicePath) BIND(C, NAME="initialize_V")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_V()
    SUBROUTINE loadSpiceFile_V(fileName) BIND(C, NAME="loadSpiceFile_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_V()
    FUNCTION createAtmosphere_V() BIND(C, NAME="createAtmosphere_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR) :: createAtmosphere_V
    END FUNCTION

    !> \copydoc copyAtmosphere_V()
    FUNCTION copyAtmosphere_V(atmos) BIND(C, NAME="copyAtmosphere_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_V
    END FUNCTION

    !> \copydoc deleteAtmosphere_V()
    SUBROUTINE deleteAtmosphere_V(atmos) BIND(C, NAME="deleteAtmosphere_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setStartTime_V()
    SUBROUTINE setStartTime_V(atmos, time) BIND(C, NAME="setStartTime_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setSeed_V()
    SUBROUTINE setSeed_V(atmos, seed) BIND(C, NAME="setSeed_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_V()
    SUBROUTINE setMinRelativeStepSize_V(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_V()
    SUBROUTINE setPerturbationScales_V(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_V()
    SUBROUTINE addAuxiliaryAtmosphere_V(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_V()
    SUBROUTINE setAuxiliaryValues_V(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_V()
    SUBROUTINE setPosition_V(atmos, position) BIND(C, NAME="setPosition_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_V()
    SUBROUTINE setDelta_V(atmos, position) BIND(C, NAME="setDelta_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_V()
    SUBROUTINE setPerturbationAction_V(atmos, action) BIND(C, NAME="setPerturbationFactors_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_V()
    SUBROUTINE setEphemerisState_V(atmos, state) BIND(C, NAME="setEphemerisState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_V()
    SUBROUTINE setEphemerisFastModeOn_V(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_V()
    SUBROUTINE setSubsolarUpdateTime_V(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_V()
    INTEGER(C_INT) FUNCTION update_V(atmos) BIND(C, NAME="update_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_V()
    FUNCTION getErrorMessage_V(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_V
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_V()
    SUBROUTINE getStartTime_V(atmos, time) BIND(C, NAME="getStartTime_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getVersionString_V()
    FUNCTION getVersionString_V(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_V")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_V
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_V()
    SUBROUTINE getPosition_V(atmos, position) BIND(C, NAME="getPosition_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_V()
    SUBROUTINE getDynamicsState_V(atmos, state) BIND(C, NAME="getDynamicsState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_V()
    SUBROUTINE getDensityState_V(atmos, state) BIND(C, NAME="getDensityState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_V()
    SUBROUTINE getWindsState_V(atmos, state) BIND(C, NAME="getWindsState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_V()
    SUBROUTINE getGasesState_V(atmos, state, hydrogen, helium, oxygen, nitrogen, dinitrogen, carbonMonoxide, carbonDioxide) BIND(C, NAME="getGasesState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: hydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: helium
      TYPE(ConstituentGas_F), INTENT(OUT) :: oxygen
      TYPE(ConstituentGas_F), INTENT(OUT) :: nitrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: dinitrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonMonoxide
      TYPE(ConstituentGas_F), INTENT(OUT) :: carbonDioxide
    END SUBROUTINE

    !> \copydoc getEphemerisState_V()
    SUBROUTINE getEphemerisState_V(atmos, state) BIND(C, NAME="getEphemerisState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_V()
    SUBROUTINE getPerturbationState_V(atmos, state) BIND(C, NAME="getPerturbationState_V")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


