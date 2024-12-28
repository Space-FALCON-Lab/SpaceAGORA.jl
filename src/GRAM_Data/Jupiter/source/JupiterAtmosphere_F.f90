! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Global Reference Atmospheric Model (GRAM) Framework
!
! No recipient of this code should forward copies outside of the United 
! States without explicit approval by NASA Marshall Space Flight Center.
! 
! Module: Jupiter-GRAM
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> \defgroup F_Jupiter The FORTRAN Interface for Jupiter
!! @{
!! \brief A listing of the FORTRAN interface functions for JupiterGRAM.
!!
!! The FORTRAN interface for Jupiter is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE JupiterGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Jupiter/examples/Jupiter_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: jupiter
!! jupiter = createAtmosphere_J(dataPath)
!! call setPosition_J(jupiter, &pos)
!! call update_J(jupiter)
!! 
!! // In C++
!! JupiterAtmosphere uranus(dataPath);
!! jupiter.setPosition(pos);
!! jupiter.update(); 
!! \endcode
!! \example Jupiter/examples/Jupiter_F.f90

MODULE JupiterGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_J()
    SUBROUTINE tryGetSpicePath_J(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_J")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_J()
    SUBROUTINE setSpiceLsk_J(lsk) BIND(C, NAME="setSpiceLsk_J")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_J()
    SUBROUTINE setSpicePck_J(pck) BIND(C, NAME="setSpicePck_J")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_J()
    SUBROUTINE setSpiceKernel_J(bsp) BIND(C, NAME="setSpiceKernel_J")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_J()
    SUBROUTINE initialize_J(spicePath) BIND(C, NAME="initialize_J")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_J()
    SUBROUTINE loadSpiceFile_J(fileName) BIND(C, NAME="loadSpiceFile_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_J()
    FUNCTION createAtmosphere_J() BIND(C, NAME="createAtmosphere_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR) :: createAtmosphere_J
    END FUNCTION

    !> \copydoc copyAtmosphere_J()
    FUNCTION copyAtmosphere_J(atmos) BIND(C, NAME="copyAtmosphere_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_J
    END FUNCTION

    !> \copydoc deleteAtmosphere_J()
    SUBROUTINE deleteAtmosphere_J(atmos) BIND(C, NAME="deleteAtmosphere_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setStartTime_J()
    SUBROUTINE setStartTime_J(atmos, time) BIND(C, NAME="setStartTime_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setSeed_J()
    SUBROUTINE setSeed_J(atmos, seed) BIND(C, NAME="setSeed_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_J()
    SUBROUTINE setMinRelativeStepSize_J(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_J()
    SUBROUTINE setPerturbationScales_J(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_J()
    SUBROUTINE addAuxiliaryAtmosphere_J(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_J()
    SUBROUTINE setAuxiliaryValues_J(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_J()
    SUBROUTINE setPosition_J(atmos, position) BIND(C, NAME="setPosition_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_J()
    SUBROUTINE setDelta_J(atmos, position) BIND(C, NAME="setDelta_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_J()
    SUBROUTINE setPerturbationAction_J(atmos, action) BIND(C, NAME="setPerturbationFactors_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_J()
    SUBROUTINE setEphemerisState_J(atmos, state) BIND(C, NAME="setEphemerisState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_J()
    SUBROUTINE setEphemerisFastModeOn_J(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_J()
    SUBROUTINE setSubsolarUpdateTime_J(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_J()
    INTEGER(C_INT) FUNCTION update_J(atmos) BIND(C, NAME="update_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_J()
    FUNCTION getErrorMessage_J(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_J
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_J()
    SUBROUTINE getStartTime_J(atmos, time) BIND(C, NAME="getStartTime_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getVersionString_J()
    FUNCTION getVersionString_J(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_J")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_J
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_J()
    SUBROUTINE getPosition_J(atmos, position) BIND(C, NAME="getPosition_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_J()
    SUBROUTINE getDynamicsState_J(atmos, state) BIND(C, NAME="getDynamicsState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_J()
    SUBROUTINE getDensityState_J(atmos, state) BIND(C, NAME="getDensityState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_J()
    SUBROUTINE getWindsState_J(atmos, state) BIND(C, NAME="getWindsState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_J()
    SUBROUTINE getGasesState_J(atmos, state) BIND(C, NAME="getGasesState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getEphemerisState_J()
    SUBROUTINE getEphemerisState_J(atmos, state) BIND(C, NAME="getEphemerisState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_J()
    SUBROUTINE getPerturbationState_J(atmos, state) BIND(C, NAME="getPerturbationState_J")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


