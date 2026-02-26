!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Uranus The FORTRAN Interface for Uranus
!! @{
!! \brief A listing of the FORTRAN interface functions for UranusGRAM.
!!
!! The FORTRAN interface for Uranus is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE UranusGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Uranus/examples/Uranus_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: uranus
!! uranus = createAtmosphere_U(dataPath)
!! call setPosition_U(uranus, &pos)
!! call update_U(uranus)
!! 
!! // In C++
!! UranusAtmosphere uranus(dataPath);
!! uranus.setPosition(pos);
!! uranus.update(); 
!! \endcode
!! \example Uranus/examples/Uranus_F.f90

MODULE UranusGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_U()
    SUBROUTINE tryGetSpicePath_U(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_U")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_U()
    SUBROUTINE setSpiceLsk_U(lsk) BIND(C, NAME="setSpiceLsk_U")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_U()
    SUBROUTINE setSpicePck_U(pck) BIND(C, NAME="setSpicePck_U")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_U()
    SUBROUTINE setSpiceKernel_U(bsp) BIND(C, NAME="setSpiceKernel_U")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_U()
    SUBROUTINE initialize_U(spicePath) BIND(C, NAME="initialize_U")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_U()
    SUBROUTINE loadSpiceFile_U(fileName) BIND(C, NAME="loadSpiceFile_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_U()
    FUNCTION createAtmosphere_U() BIND(C, NAME="createAtmosphere_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR) :: createAtmosphere_U
    END FUNCTION

    !> \copydoc copyAtmosphere_U()
    FUNCTION copyAtmosphere_U(atmos) BIND(C, NAME="copyAtmosphere_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_U
    END FUNCTION

    !> \copydoc deleteAtmosphere_U()
    SUBROUTINE deleteAtmosphere_U(atmos) BIND(C, NAME="deleteAtmosphere_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setStartTime_U()
    SUBROUTINE setStartTime_U(atmos, time) BIND(C, NAME="setStartTime_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setSeed_U()
    SUBROUTINE setSeed_U(atmos, seed) BIND(C, NAME="setSeed_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_U()
    SUBROUTINE setMinRelativeStepSize_U(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_U()
    SUBROUTINE setPerturbationScales_U(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_U()
    SUBROUTINE addAuxiliaryAtmosphere_U(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_U()
    SUBROUTINE setAuxiliaryValues_U(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_U()
    SUBROUTINE setPosition_U(atmos, position) BIND(C, NAME="setPosition_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_U()
    SUBROUTINE setDelta_U(atmos, position) BIND(C, NAME="setDelta_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_U()
    SUBROUTINE setPerturbationAction_U(atmos, action) BIND(C, NAME="setPerturbationFactors_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_U()
    SUBROUTINE setEphemerisState_U(atmos, state) BIND(C, NAME="setEphemerisState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_U()
    SUBROUTINE setEphemerisFastModeOn_U(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_U()
    SUBROUTINE setSubsolarUpdateTime_U(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_U()
    INTEGER(C_INT) FUNCTION update_U(atmos) BIND(C, NAME="update_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_U()
    FUNCTION getErrorMessage_U(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_U
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_U()
    SUBROUTINE getStartTime_U(atmos, time) BIND(C, NAME="getStartTime_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getVersionString_U()
    FUNCTION getVersionString_U(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_U")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_U
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_U()
    SUBROUTINE getPosition_U(atmos, position) BIND(C, NAME="getPosition_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_U()
    SUBROUTINE getDynamicsState_U(atmos, state) BIND(C, NAME="getDynamicsState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_U()
    SUBROUTINE getDensityState_U(atmos, state) BIND(C, NAME="getDensityState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_U()
    SUBROUTINE getWindsState_U(atmos, state) BIND(C, NAME="getWindsState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_U()
    SUBROUTINE getGasesState_U(atmos, state, dihydrogen, methane, helium) BIND(C, NAME="getGasesState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: dihydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: methane
      TYPE(ConstituentGas_F), INTENT(OUT) :: helium
    END SUBROUTINE

    !> \copydoc getEphemerisState_U()
    SUBROUTINE getEphemerisState_U(atmos, state) BIND(C, NAME="getEphemerisState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_U()
    SUBROUTINE getPerturbationState_U(atmos, state) BIND(C, NAME="getPerturbationState_U")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


