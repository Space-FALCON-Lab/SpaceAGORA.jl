!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.

PROGRAM Venus_F
  USE, INTRINSIC :: ISO_C_BINDING
  USE GramStructs
  USE VenusGRAM

  TYPE(C_PTR) :: venus, venus2
  TYPE(GramTime_F) :: time
  TYPE(Position_F) :: pos, pos2
  TYPE(DynamicsState_F) :: dstate, dstate2
  CHARACTER(KIND=C_CHAR, LEN = 500) :: spicePath
  CHARACTER(KIND=C_CHAR, LEN = 100) :: version
  INTEGER(C_SIZE_T) :: bufferSize
  INTEGER(C_SIZE_T) :: stringLength
  INTEGER(C_INT) :: errorFlag;
  TYPE(GasesState_F) :: gstate, gstate2
  TYPE(ConstituentGas_F) :: hydrogen, hydrogen2
  TYPE(ConstituentGas_F) :: helium, helium2
  TYPE(ConstituentGas_F) :: oxygen, oxygen2
  TYPE(ConstituentGas_F) :: nitrogen, nitrogen2
  TYPE(ConstituentGas_F) :: dinitrogen, dinitrogen2
  TYPE(ConstituentGas_F) :: carbonMonoxide, carbonMonoxide2
  TYPE(ConstituentGas_F) :: carbonDioxide, carbonDioxide2
  TYPE(DensityState_F) ::  rstate, rstate2
  TYPE(WindsState_F) :: wstate, wstate2
  TYPE(EphemerisState_F) :: estate, estate2
  CHARACTER(LEN=*), PARAMETER  :: FMT = "(A, T32, ES13.6, 1X, ES13.6)"

  write(*,*) "================="
  write(*,*) "A FORTRAN Example"
  write(*,*) "================="
  write(*,*)

  ! The first call must be to initialize the Spice path.
  call tryGetSpicePath_V(spicePath, 500);
  call initialize_V(spicePath)

  ! Create a Venus atmosphere.  A pointer is returned that is our
  ! handle on the venus model.  More than one Venus atmosphere
  ! can be created.  (Also see deleteAtmosphere_V() below.)
  venus = createAtmosphere_V()
  venus2 = createAtmosphere_V()

  ! Get and print the version string.
  ! On input, bufferSize is the length of the character array.
  ! On output, bufferSize is the length of the version string.
  bufferSize = 100
  stringLength = getVersionString_V(venus, version, bufferSize)
  write(*,*) version(1:stringLength)
  write(*,*)

  ! Set the initial seed, the relative step size, and the perturbation scale factors.
  call setSeed_V(venus, 1001)
  call setMinRelativeStepSize_V(venus, 0.5D0)
  call setPerturbationScales_V(venus, 1.5D0, 1.5D0, 1.5D0)
  call setSeed_V(venus2, 1001)
  call setMinRelativeStepSize_V(venus2, 0.5D0)
  call setPerturbationScales_V(venus2, 1.5D0, 1.5D0, 1.5D0)

  ! Set the start time.
  time%year = 2020
  time%month = 3
  time%day = 15
  time%hour = 0
  time%minute = 0
  time%seconds = 0.0D0
  time%frame = 1
  time%scale = 1
  call setStartTime_V(venus, time)
  call setStartTime_V(venus2, time)

  ! Set the position.
  pos%height = 50.0
  pos%latitude = 22.0
  pos%longitude = 48.0
  pos%elapsedTime = 100.0
  pos%isPlanetoCentric = 1
  call setPosition_V(venus, pos)
  pos%height = 200.0
  call setPosition_V(venus2, pos)

  ! Update the model.
  errorFlag = update_V(venus)
  errorFlag = update_V(venus2)

  ! The position and other output is returned in the Position_F structure.
  call getPosition_V(venus, pos)
  call getPosition_V(venus2, pos2)
  write(*,*) "                            Venus #1      Venus #2"
  write(*,FMT) "Total Radius:", pos%totalRadius, pos2%totalRadius
  write(*,FMT) "Gravity:", pos%gravity, pos2%gravity
  write(*,*)

  ! The ephemeris values are returned in the EphemerisState_F structure
  call getEphemerisState_V(venus, estate)
  call getEphemerisState_V(venus2, estate2)
  write(*,FMT) "Solar Time:", estate%solarTime, estate2%solarTime
  write(*,FMT) "Longitude of the Sun:", estate%longitudeSun, estate2%longitudeSun
  write(*,*)

  ! The atmosphere state is returned in multiple structures:
  ! DynamicsState_F, DensityState_F, WindsState_F, GasesState_F
  call getDynamicsState_V(venus, dstate)
  call getDynamicsState_V(venus2, dstate2)
  write(*,FMT) "Temperature:", dstate%temperature, dstate2%temperature
  write(*,FMT) "Pressure:", dstate%pressure, dstate2%pressure
  write(*,FMT) "Density:", dstate%density, dstate2%density
  write(*,FMT) "Pressure Scale Height:", dstate%pressureScaleHeight, dstate2%pressureScaleHeight
  write(*,FMT) "Density Scale Height:", dstate%densityScaleHeight, dstate2%densityScaleHeight
  write(*,*)

  ! Get and print perturbed density
  call getDensityState_V(venus, rstate)
  call getDensityState_V(venus2, rstate2)
  write(*,FMT) "Mean Density:", rstate%density, rstate2%density
  write(*,FMT) "Perturbed Density:", rstate%perturbedDensity, rstate2%perturbedDensity
  write(*,FMT) "Perturbation Percent:", rstate%densityPerturbation, rstate2%densityPerturbation
  write(*,*)

  ! Get and print gases.
  call getGasesState_V(venus, gstate, hydrogen, helium, oxygen, nitrogen, dinitrogen, carbonMonoxide, carbonDioxide)
  call getGasesState_V(venus2, gstate2, hydrogen2, helium2, oxygen2, nitrogen2, dinitrogen2, carbonMonoxide2, carbonDioxide2)
  write(*,FMT) "Average Molecular Weight:", gstate%averageMolecularWeight, gstate2%averageMolecularWeight
  write(*,FMT) "Hydrogen Mole Fraction:", hydrogen%moleFraction, hydrogen2%moleFraction
  write(*,FMT) "Helium Mole Fraction:", helium%moleFraction, helium2%moleFraction
  write(*,FMT) "Oxygen Mole Fraction:", oxygen%moleFraction, oxygen2%moleFraction
  write(*,FMT) "Nitrogen Mole Fraction:", nitrogen%moleFraction, nitrogen2%moleFraction
  write(*,FMT) "Dinitrogen Mole Fraction:", dinitrogen%moleFraction, dinitrogen2%moleFraction
  write(*,FMT) "Carbon Monoxide Mole Fraction:", carbonMonoxide%moleFraction, carbonMonoxide2%moleFraction
  write(*,FMT) "Carbon Dioxide Mole Fraction:", carbonDioxide%moleFraction, carbonDioxide2%moleFraction
  write(*,*)

  ! Get and print wind data
  call getWindsState_V(venus, wstate)
  call getWindsState_V(venus2, wstate2)
  write(*,FMT) "EW Wind:", wstate%ewWind, wstate2%ewWind
  write(*,FMT) "NS Wind:", wstate%nsWind, wstate2%nsWind
  write(*,FMT) "Perturbed EW Wind:", wstate%perturbedEWWind, wstate2%perturbedEWWind
  write(*,FMT) "Perturbed NS Wind:", wstate%perturbedNSWind, wstate2%perturbedNSWind

  ! Memory is freed up in the deleteAtmosphere call.
  call deleteAtmosphere_V(venus)
  call deleteAtmosphere_V(venus2)

END PROGRAM
