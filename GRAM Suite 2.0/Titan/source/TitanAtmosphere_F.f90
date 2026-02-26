!
!To change this license header, choose License Headers in Project Properties.
!To change this template file, choose Tools | Templates
!and open the template in the editor.
!

!> \defgroup F_Titan The FORTRAN Interface for Titan
!! @{
!! \brief A listing of the FORTRAN interface functions for TitanGRAM.
!!
!! The FORTRAN interface for Titan is declared in the following:
!! \code{cpp}
!!   USE, INTRINSIC :: ISO_C_BINDING
!!   USE GramStructs
!!   USE TitanGRAM
!! \endcode 
!! An example using the FORTRAN interface can be found \ref Titan/examples/Titan_F.f90 "here".
!! <br><br>The FORTRAN interface is a wrapper around the C++ GRAM library. The design of the interface
!! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
!! FORTRAN and C++ code below.
!! \code{cpp}
!! ! In FORTRAN
!! TYPE(C_PTR) :: titan
!! titan = createAtmosphere_T(dataPath)
!! call setPosition_T(titan, &pos)
!! call update_T(titan)
!! 
!! // In C++
!! TitanAtmosphere titan(dataPath);
!! titan.setPosition(pos);
!! titan.update(); 
!! \endcode
!! \example Titan/examples/Titan_F.f90


MODULE TitanGRAM
  INTERFACE

    !> \copydoc tryGetSpicePath_T()
    SUBROUTINE tryGetSpicePath_T(spicePath, bufferSize) BIND(C, NAME="tryGetSpicePath_T")
      USE ISO_C_BINDING, ONLY : C_CHAR, C_INT
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: spicePath(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: bufferSize
    END SUBROUTINE

    !> \copydoc setSpiceLsk_T()
    SUBROUTINE setSpiceLsk_T(lsk) BIND(C, NAME="setSpiceLsk_T")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: lsk(*)
    END SUBROUTINE

    !> \copydoc setSpicePck_T()
    SUBROUTINE setSpicePck_T(pck) BIND(C, NAME="setSpicePck_T")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: pck(*)
    END SUBROUTINE

    !> \copydoc setSpiceKernel_T()
    SUBROUTINE setSpiceKernel_T(bsp) BIND(C, NAME="setSpiceKernel_T")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: bsp(*)
    END SUBROUTINE

    !> \copydoc initialize_T()
    SUBROUTINE initialize_T(spicePath) BIND(C, NAME="initialize_T")
      USE ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: spicePath(*)
    END SUBROUTINE

    !> \copydoc loadSpiceFile_T()
    SUBROUTINE loadSpiceFile_T(fileName) BIND(C, NAME="loadSpiceFile_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
    END SUBROUTINE

    !> \copydoc createAtmosphere_T()
    FUNCTION createAtmosphere_T() BIND(C, NAME="createAtmosphere_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR) :: createAtmosphere_T
    END FUNCTION

    !> \copydoc copyAtmosphere_T()
    FUNCTION copyAtmosphere_T(atmos) BIND(C, NAME="copyAtmosphere_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(C_PTR) :: copyAtmosphere_T
    END FUNCTION

    !> \copydoc deleteAtmosphere_T()
    SUBROUTINE deleteAtmosphere_T(atmos) BIND(C, NAME="deleteAtmosphere_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END SUBROUTINE

    !> \copydoc setStartTime_T()
    SUBROUTINE setStartTime_T(atmos, time) BIND(C, NAME="setStartTime_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(IN) :: time
    END SUBROUTINE

    !> \copydoc setMethaneMoleFraction_T()
    SUBROUTINE setMethaneMoleFraction_T(atmos, mmf) BIND(C, NAME="setMethaneMoleFraction_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: mmf
    END SUBROUTINE

    !> \copydoc setModelType_T()
    SUBROUTINE setModelType_T(atmos, type) BIND(C, NAME="setModelType_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: type
    END SUBROUTINE

    !> \copydoc setMinMaxFactor_T()
    SUBROUTINE setMinMaxFactor_T(atmos, factor, computeFlag) BIND(C, NAME="setMinMaxFactor_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: factor
      INTEGER(C_INT), VALUE, INTENT(IN) :: computeFlag
    END SUBROUTINE

    !> \copydoc setSeed_T()
    SUBROUTINE setSeed_T(atmos, seed) BIND(C, NAME="setSeed_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: seed
    END SUBROUTINE

    !> \copydoc setMinRelativeStepSize_T()
    SUBROUTINE setMinRelativeStepSize_T(atmos, minRelativeStepSize) BIND(C, NAME="setMinRelativeStepSize_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: minRelativeStepSize
    END SUBROUTINE

    !> \copydoc setPerturbationScales_T()
    SUBROUTINE setPerturbationScales_T(atmos, densityScale, ewWindScale, nsWindScale) BIND(C, NAME="setPerturbationScales_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: densityScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ewWindScale
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nsWindScale
    END SUBROUTINE

    !> \copydoc addAuxiliaryAtmosphere_T()
    SUBROUTINE addAuxiliaryAtmosphere_T(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive) BIND(C, NAME="addAuxiliaryAtmosphere_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: fileName(*)
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: innerRadius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: outerRadius
      INTEGER(C_INT), VALUE, INTENT(IN) :: isEastLongitudePositive
    END SUBROUTINE

    !> \copydoc setAuxiliaryValues_T()
    SUBROUTINE setAuxiliaryValues_T(atmos, dens, pres, temp, ew, ns) BIND(C, NAME="setAuxiliaryValues_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT, C_CHAR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: dens
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: pres
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: temp
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ew
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: ns
    END SUBROUTINE

    !> \copydoc setPosition_T()
    SUBROUTINE setPosition_T(atmos, position) BIND(C, NAME="setPosition_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setDelta_T()
    SUBROUTINE setDelta_T(atmos, position) BIND(C, NAME="setDelta_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(IN) :: position
    END SUBROUTINE

    !> \copydoc setPerturbationAction_T()
    SUBROUTINE setPerturbationAction_T(atmos, action) BIND(C, NAME="setPerturbationFactors_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: action
    END SUBROUTINE

    !> \copydoc setEphemerisState_T()
    SUBROUTINE setEphemerisState_T(atmos, state) BIND(C, NAME="setEphemerisState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(IN) :: state
    END SUBROUTINE

    !> \copydoc setEphemerisFastModeOn_T()
    SUBROUTINE setEphemerisFastModeOn_T(atmos, flag) BIND(C, NAME="setEphemerisFastModeOn_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      INTEGER(C_INT), VALUE, INTENT(IN) :: flag
    END SUBROUTINE

    !> \copydoc setSubsolarUpdateTime_T()
    SUBROUTINE setSubsolarUpdateTime_T(atmos, utime) BIND(C, NAME="setSubsolarUpdateTime_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(IN) :: utime
    END SUBROUTINE

    !> \copydoc update_T()
    INTEGER(C_INT) FUNCTION update_T(atmos) BIND(C, NAME="update_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
    END FUNCTION

    !> \copydoc getErrorMessage_T()
    FUNCTION getErrorMessage_T(atmos, message, bufferSize) BIND(C, NAME="getErrorMessage_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getErrorMessage_T
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: message(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getStartTime_T()
    SUBROUTINE getStartTime_T(atmos, time) BIND(C, NAME="getStartTime_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GramTime_F), INTENT(OUT) :: time
    END SUBROUTINE

    !> \copydoc getMinMaxFactor_T()
    SUBROUTINE getMinMaxFactor_T(atmos, factor, computeFlag) BIND(C, NAME="getMinMaxFactor_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_DOUBLE, C_INT
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      REAL(C_DOUBLE), INTENT(OUT) :: factor
      INTEGER(C_INT), INTENT(OUT) :: computeFlag
    END SUBROUTINE

    !> \copydoc getVersionString_T()
    FUNCTION getVersionString_T(atmos, versionString, bufferSize) BIND(C, NAME="getVersionString_T")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_SIZE_T, C_PTR, C_CHAR
      INTEGER(C_SIZE_T) :: getVersionString_T
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      CHARACTER(KIND=C_CHAR), INTENT(OUT) :: versionString(*)
      INTEGER(C_SIZE_T), VALUE, INTENT(IN) :: bufferSize
    END FUNCTION

    !> \copydoc getPosition_T()
    SUBROUTINE getPosition_T(atmos, position) BIND(C, NAME="getPosition_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(Position_F), INTENT(OUT) :: position
    END SUBROUTINE

    !> \copydoc getDynamicsState_T()
    SUBROUTINE getDynamicsState_T(atmos, state) BIND(C, NAME="getDynamicsState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DynamicsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getDensityState_T()
    SUBROUTINE getDensityState_T(atmos, state) BIND(C, NAME="getDensityState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(DensityState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getWindsState_T()
    SUBROUTINE getWindsState_T(atmos, state) BIND(C, NAME="getWindsState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(WindsState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getGasesState_T()
    SUBROUTINE getGasesState_T(atmos, state, argon, methane, dinitrogen) BIND(C, NAME="getGasesState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(GasesState_F), INTENT(OUT) :: state
      TYPE(ConstituentGas_F), INTENT(OUT) :: argon
      TYPE(ConstituentGas_F), INTENT(OUT) :: methane
      TYPE(ConstituentGas_F), INTENT(OUT) :: dinitrogen
    END SUBROUTINE

    !> \copydoc getEphemerisState_T()
    SUBROUTINE getEphemerisState_T(atmos, state) BIND(C, NAME="getEphemerisState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(EphemerisState_F), INTENT(OUT) :: state
    END SUBROUTINE

    !> \copydoc getPerturbationState_T()
    SUBROUTINE getPerturbationState_T(atmos, state) BIND(C, NAME="getPerturbationState_T")
      USE GramStructs
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE, INTENT(IN) :: atmos
      TYPE(PerturbationState_F), INTENT(OUT) :: state
    END SUBROUTINE

  END INTERFACE

END MODULE

!> @}


