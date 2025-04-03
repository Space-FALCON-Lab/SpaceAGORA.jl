!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.

PROGRAM Titan_F
  USE, INTRINSIC :: ISO_C_BINDING
  USE GramStructs
  USE TitanGRAM

  TYPE(C_PTR) :: titan, titan2
  TYPE(GramTime_F) :: time
  TYPE(Position_F) :: pos, pos2
  TYPE(DynamicsState_F) :: dstate, dstate2
  CHARACTER(KIND=C_CHAR, LEN = 500) :: spicePath
  CHARACTER(KIND=C_CHAR, LEN = 100) :: version
  INTEGER(C_SIZE_T) :: bufferSize
  INTEGER(C_SIZE_T) :: stringLength
  INTEGER(C_INT) :: errorFlag;
  TYPE(GasesState_F) :: gstate, gstate2
  TYPE(ConstituentGas_F) :: argon, argon2
  TYPE(ConstituentGas_F) :: methane, methane2
  TYPE(ConstituentGas_F) :: dinitrogen, dinitrogen2
  TYPE(DensityState_F) ::  rstate, rstate2
  TYPE(WindsState_F) :: wstate, wstate2
  TYPE(EphemerisState_F) :: estate, estate2
  CHARACTER(LEN=*), PARAMETER  :: FMT = "(A, T27, ES13.6, 1X, ES13.6)"

  write(*,*) "================="
  write(*,*) "A FORTRAN Example"
  write(*,*) "================="
  write(*,*)

  ! The first call must be to initialize the Spice path.
  call tryGetSpicePath_T(spicePath, 500);
  call initialize_T(spicePath)

  ! Create a Titan atmosphere.  A pointer is returned that is our
  ! handle on the titan model.  More than one Titan atmosphere
  ! can be created.  (Also see deleteAtmosphere_T() below.)
  titan = createAtmosphere_T()
  titan2 = createAtmosphere_T()

  ! Get and print the version string.
  ! On input, bufferSize is the length of the character array.
  ! On output, bufferSize is the length of the version string.
  bufferSize = 100
  stringLength = getVersionString_T(titan, version, bufferSize)
  write(*,*) version(1:stringLength)
  write(*,*)

  ! Compare two different models
  call setModelType_T(titan, 1);   ! Yelle97
  call setMinMaxFactor_T(titan, 0.0d0, 1);
  call setModelType_T(titan2, 2);  ! GCM95

  ! Set the initial seed, the relative step size, and the perturbation scale factors.
  call setSeed_T(titan, 1001)
  call setMinRelativeStepSize_T(titan, 0.5D0)
  call setPerturbationScales_T(titan, 1.5D0, 1.5D0, 1.5D0)
  call setSeed_T(titan2, 1001)
  call setMinRelativeStepSize_T(titan2, 0.5D0)
  call setPerturbationScales_T(titan2, 1.5D0, 1.5D0, 1.5D0)

  ! Set the start time.
  time%year = 2020
  time%month = 3
  time%day = 15
  time%hour = 0
  time%minute = 0
  time%seconds = 0.0D0
  time%frame = 1
  time%scale = 1
  call setStartTime_T(titan, time)
  call setStartTime_T(titan2, time)

  ! Set the position.
  pos%height = 50.0
  pos%latitude = 22.0
  pos%longitude = 48.0
  pos%elapsedTime = 100.0
  pos%isPlanetoCentric = 1
  call setPosition_T(titan, pos)
  call setPosition_T(titan2, pos)

  ! Update the model.
  errorFlag = update_T(titan)
  errorFlag = update_T(titan2)

  ! The position and other output is returned in the Position_F structure.
  call getPosition_T(titan, pos)
  call getPosition_T(titan2, pos2)
  write(*,*) "                            Titan #1      Titan #2"
  write(*,FMT) "Total Radius:", pos%totalRadius, pos2%totalRadius
  write(*,FMT) "Gravity:", pos%gravity, pos2%gravity
  write(*,*)

  ! The ephemeris values are returned in the EphemerisState_F structure
  call getEphemerisState_T(titan, estate)
  call getEphemerisState_T(titan2, estate2)
  write(*,FMT) "Solar Time:", estate%solarTime, estate2%solarTime
  write(*,FMT) "Longitude of the Sun:", estate%longitudeSun, estate2%longitudeSun
  write(*,*)

  ! The atmosphere state is returned in multiple structures:
  ! DynamicsState_F, DensityState_F, WindsState_F, GasesState_F
  call getDynamicsState_T(titan, dstate)
  call getDynamicsState_T(titan2, dstate2)
  write(*,FMT) "Temperature:", dstate%temperature, dstate2%temperature
  write(*,FMT) "Pressure:", dstate%pressure, dstate2%pressure
  write(*,FMT) "Density:", dstate%density, dstate2%density
  write(*,FMT) "Pressure Scale Height:", dstate%pressureScaleHeight, dstate2%pressureScaleHeight
  write(*,FMT) "Density Scale Height:", dstate%densityScaleHeight, dstate2%densityScaleHeight
  write(*,*)

  ! Get and print perturbed density
  call getDensityState_T(titan, rstate)
  call getDensityState_T(titan2, rstate2)
  write(*,FMT) "Mean Density:", rstate%density, rstate2%density
  write(*,FMT) "Perturbed Density:", rstate%perturbedDensity, rstate2%perturbedDensity
  write(*,FMT) "Perturbation Percent:", rstate%densityPerturbation, rstate2%densityPerturbation
  write(*,*)

  ! Get and print gases.
  call getGasesState_T(titan, gstate, argon, methane, dinitrogen)
  call getGasesState_T(titan2, gstate2, argon2, methane2, dinitrogen2)
  write(*,FMT) "Average Molecular Weight:", gstate%averageMolecularWeight, gstate2%averageMolecularWeight
  write(*,FMT) "Argon Mole Fraction:", argon%moleFraction, argon2%moleFraction
  write(*,FMT) "Methane Mole Fraction:", methane%moleFraction, methane2%moleFraction
  write(*,FMT) "Dinitrogen Mole Fraction:", dinitrogen%moleFraction, dinitrogen2%moleFraction
  write(*,*)

  ! Get and print wind data
  call getWindsState_T(titan, wstate)
  call getWindsState_T(titan2, wstate2)
  write(*,FMT) "EW Wind:", wstate%ewWind, wstate2%ewWind
  write(*,FMT) "Perturbed EW Wind:", wstate%perturbedEWWind, wstate2%perturbedEWWind

  ! Memory is freed up in the deleteAtmosphere call.
  call deleteAtmosphere_T(titan)
  call deleteAtmosphere_T(titan2)

END PROGRAM
