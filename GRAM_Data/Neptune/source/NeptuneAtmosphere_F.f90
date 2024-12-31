!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Neptune The FORTRAN Interface for Neptune
!! @{
!! \brief A listing of the FORTRAN interface functions for NeptuneGRAM.
!!
!! The FORTRAN interface for Neptune is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE NeptuneGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Neptune/examples/Neptune_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: neptune
!! neptune = createAtmosphere_N(dataPath)
!! call setPosition_N(neptune, &pos)
!! call update_N(neptune)
!! 
!! // In C++
!! NeptuneAtmosphere neptune(dataPath);
!! neptune.setPosition(pos);
!! neptune.update(); 
!! \endcode
!! \example Neptune/examples/Neptune_F.f90

MODULE NeptuneGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_N()
    SUBROUTINE tryGetSpicePath_N(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_N")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_N()
    SUBROUTINE setSpiceLsk_N(lsk) BIND(C, NAME="setSpiceLsk_N")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_N()
    SUBROUTINE setSpicePck_N(pck) BIND(C, NAME="setSpicePck_N")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_N()
    SUBROUTINE setSpiceKernel_N(bsp) BIND(C, NAME="setSpiceKernel_N")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_N()
    SUBROUTINE initialize_N(spicePath) BIND(C, NAME="initialize_N")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_N()
    SUBROUTINE loadSpiceFile_N(fileName) BIND(C, NAME="loadSpiceFile_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_N()
    FUNCTION createAtmosphere_N() BIND(C, NAME="createAtmosphere_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR) :: createAtmosphere_N
    END FUNCTION

    !> \copydoc copyAtmosphere_N()
    FUNCTION copyAtmosphere_N(atmos) BIND(C, NAME="copyAtmosphere_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_N
    END FUNCTION

    !> \copydoc deleteAtmosphere_N()
    SUBROUTINE deleteAtmosphere_N(atmos) BIND(C, NAME="deleteAtmosphere_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setStartTime_N()
    SUBROUTINE setStartTime_N(atmos, time) BIND(C, NAME="setStartTime_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setDinitrogenMoleFraction_N()
    SUBROUTINE setDinitrogenMoleFraction_N(atmos, n2mf) BIND(C, NAME="setDinitrogenMoleFraction_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: n2mf
    END SUBROUTINE

    !> \copydoc setMinMaxFactor_N()
    SUBROUTINE setMinMaxFactor_N(atmos, factor, computeFlag) BIND(C, NAME="setMinMaxFactor_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: factor
      INTEGER(C_INT), VALUE, INTENT(IN) :: computeFlag
    END SUBROUTINE

    !> \copydoc setSeed_N()
    SUBROUTINE setSeed_N(atmos, seed) BIND(C, NAME="setSeed_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_N()
    SUBROUTINE setMinRelativeStepSize_N(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_N()
    SUBROUTINE setPerturbationScales_N(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_N()
    SUBROUTINE addAuxiliaryAtmosphere_N(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_N()
    SUBROUTINE setAuxiliaryValues_N(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_N()
    SUBROUTINE setPosition_N(atmos, position) BIND(C, NAME="setPosition_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_N()
    SUBROUTINE setDelta_N(atmos, position) BIND(C, NAME="setDelta_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_N()
    SUBROUTINE setPerturbationAction_N(atmos, action) BIND(C, NAME="setPerturbationFactors_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_N()
    SUBROUTINE setEphemerisState_N(atmos, state) BIND(C, NAME="setEphemerisState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_N()
    SUBROUTINE setEphemerisFastModeOn_N(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_N()
    SUBROUTINE setSubsolarUpdateTime_N(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_N()
    INTEGER(C_INT) FUNCTION update_N(atmos) BIND(C, NAME="update_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_N()
    FUNCTION getErrorMessage_N(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_N
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_N()
    SUBROUTINE getStartTime_N(atmos, time) BIND(C, NAME="getStartTime_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getMinMaxFactor_N()
    SUBROUTINE getMinMaxFactor_N(atmos, factor, computeFlag) BIND(C, NAME="getMinMaxFactor_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(OUT) :: factor
      INTEGER(C_INT), INTENT(OUT) :: computeFlag
    END SUBROUTINE

    !> \copydoc getVersionString_N()
    FUNCTION getVersionString_N(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_N")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_N
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_N()
    SUBROUTINE getPosition_N(atmos, position) BIND(C, NAME="getPosition_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_N()
    SUBROUTINE getDynamicsState_N(atmos, state) BIND(C, NAME="getDynamicsState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_N()
    SUBROUTINE getDensityState_N(atmos, state) BIND(C, NAME="getDensityState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_N()
    SUBROUTINE getWindsState_N(atmos, state) BIND(C, NAME="getWindsState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_N()
    SUBROUTINE getGasesState_N(atmos, state, dihydrogen, methane, helium, dinitrogen) BIND(C, NAME="getGasesState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: dihydrogen
      TYPE(ConstituentGas_F), INTENT(OUT) :: methane
      TYPE(ConstituentGas_F), INTENT(OUT) :: helium
      TYPE(ConstituentGas_F), INTENT(OUT) :: dinitrogen
    END SUBROUTINE

    !> \copydoc getEphemerisState_N()
    SUBROUTINE getEphemerisState_N(atmos, state) BIND(C, NAME="getEphemerisState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_N()
    SUBROUTINE getPerturbationState_N(atmos, state) BIND(C, NAME="getPerturbationState_N")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


