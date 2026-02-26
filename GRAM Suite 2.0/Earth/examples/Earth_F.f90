!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.

PROGRAM Earth_F
  USE, INTRINSIC :: ISO_C_BINDING
  USE GramStructs
  USE EarthStructs
  USE EarthGRAM

  TYPE(C_PTR) :: earth, earth2
  TYPE(GramTime_F) :: time
  TYPE(Position_F) :: pos, pos2
  TYPE(DynamicsState_F) :: dstate, dstate2
  CHARACTER(KIND=C_CHAR, LEN = 500) :: spicePath
  CHARACTER(KIND=C_CHAR, LEN = 500) :: dataPath
  CHARACTER(KIND=C_CHAR, LEN = 100) :: version
  INTEGER(C_SIZE_T) :: bufferSize
  INTEGER(C_SIZE_T) :: stringLength
  INTEGER(C_INT) :: errorFlag;
  TYPE(GasesState_F) :: gstate, gstate2
  TYPE(ConstituentGas_F) :: argon, argon2
  TYPE(ConstituentGas_F) :: carbonDioxide, carbonDioxide2
  TYPE(ConstituentGas_F) :: carbonMonoxide, carbonMonoxide2
  TYPE(ConstituentGas_F) :: dinitrogen, dinitrogen2
  TYPE(ConstituentGas_F) :: dioxygen, dioxygen2
  TYPE(ConstituentGas_F) :: helium, helium2
  TYPE(ConstituentGas_F) :: hydrogen, hydrogen2
  TYPE(ConstituentGas_F) :: methane, methane2
  TYPE(ConstituentGas_F) :: nitrogen, nitrogen2
  TYPE(ConstituentGas_F) :: nitrousOxide, nitrousOxide2
  TYPE(ConstituentGas_F) :: oxygen, oxygen2
  TYPE(ConstituentGas_F) :: ozone, ozone2
  TYPE(ConstituentGas_F) :: water, water2
  TYPE(DensityState_F) ::  rstate, rstate2
  TYPE(WindsState_F) :: wstate, wstate2
  TYPE(EphemerisState_F) :: estate, estate2
  TYPE(EarthState_F) :: eatmos, eatmos2
  CHARACTER(LEN=*), PARAMETER  :: FMT = "(A, T27, ES13.6, 1X, ES13.6)"

  write(*,*) "================="
  write(*,*) "A FORTRAN Example"
  write(*,*) "================="
  write(*,*)

  ! The first call must be to initialize the Spice path.
  call tryGetDataPaths_E(spicePath, dataPath, 500);
  call initialize_E(spicePath)

  ! Create a Earth atmosphere.  A pointer is returned that is our
  ! handle on the earth model.  More than one Earth atmosphere
  ! can be created.  (Also see deleteAtmosphere_E() below.)
  earth = createAtmosphere_E(dataPath)
  earth2 = createAtmosphere_E(dataPath)

  ! Get and print the version string.
  ! On input, bufferSize is the length of the character array.
  ! On output, bufferSize is the length of the version string.
  bufferSize = 100
  stringLength = getVersionString_E(earth, version, bufferSize)
  write(*,*) version(1:stringLength)
  write(*,*)

  ! Set the MERRA-2 hour and data extents
  call setMERRA2Parameters_E(earth, 1, 0.0D0, 40.0D0, 20.0D0, 60.0D0);
  call setMERRA2Parameters_E(earth2, 9, 0.0D0, 40.0D0, 20.0D0, 60.0D0);

  ! For NCEP model, uncomment lines below
  !call setUseNCEP_E(earth, true);
  !call setUseNCEP_E(earth2, true);
  !call setNCEPParameters_E(earth, 9715, 1);
  !call setNCEPParameters_E(earth2, 9715, 5);

  ! Set the initial seed, the relative step size, and the perturbation scale factors.
  call setSeed_E(earth, 1001)
  call setPerturbationScales_E(earth, 1.5D0, 1.5D0, 1.5D0)
  call setSeed_E(earth2, 1001)
  call setPerturbationScales_E(earth2, 1.5D0, 1.5D0, 1.5D0)

  ! Set the start time.
  time%year = 2020
  time%month = 12
  time%day = 15
  time%hour = 0
  time%minute = 0
  time%seconds = 0.0D0
  time%frame = 1
  time%scale = 1
  call setStartTime_E(earth, time)
  call setStartTime_E(earth2, time)

  ! Set the position.
  pos%height = 5.0
  pos%latitude = 22.0
  pos%longitude = 48.0
  pos%elapsedTime = 100.0
  pos%isPlanetoCentric = 1
  call setPosition_E(earth, pos)
  call setPosition_E(earth2, pos)

  ! Update the model.
  errorFlag = update_E(earth)
  errorFlag = update_E(earth2)

  ! The position and other output is returned in the Position_F structure.
  call getPosition_E(earth, pos)
  call getPosition_E(earth2, pos2)
  write(*,*) "                            Earth #1      Earth #2"
  write(*,FMT) "Total Radius:", pos%totalRadius, pos2%totalRadius
  write(*,FMT) "Gravity:", pos%gravity, pos2%gravity
  write(*,*)

  ! The ephemeris values are returned in the EphemerisState_F structure
  call getEphemerisState_E(earth, estate)
  call getEphemerisState_E(earth2, estate2)
  write(*,FMT) "Solar Time:", estate%solarTime, estate2%solarTime
  write(*,FMT) "Longitude of the Sun:", estate%longitudeSun, estate2%longitudeSun
  write(*,*)

  ! The atmosphere state is returned in multiple structures:
  ! DynamicsState_F, DensityState_F, WindsState_F, GasesState_F
  call getDynamicsState_E(earth, dstate)
  call getDynamicsState_E(earth2, dstate2)
  write(*,FMT) "Temperature:", dstate%temperature, dstate2%temperature
  write(*,FMT) "Pressure:", dstate%pressure, dstate2%pressure
  write(*,FMT) "Density:", dstate%density, dstate2%density
  write(*,FMT) "Pressure Scale Height:", dstate%pressureScaleHeight, dstate2%pressureScaleHeight
  write(*,FMT) "Density Scale Height:", dstate%densityScaleHeight, dstate2%densityScaleHeight
  write(*,*)

  ! Get and print perturbed density
  call getDensityState_E(earth, rstate)
  call getDensityState_E(earth2, rstate2)
  call getEarthState_E(earth, eatmos);
  call getEarthState_E(earth2, eatmos2);
  write(*,FMT) "Perturbed Density:", rstate%perturbedDensity, rstate2%perturbedDensity
  write(*,FMT) "Perturbed Pressure:", eatmos%perturbedPressure, eatmos2%perturbedPressure
  write(*,FMT) "Perturbed Temperature:", eatmos%perturbedTemperature, eatmos2%perturbedTemperature
  write(*,FMT) "Density Perturbation:", rstate%densityPerturbation, rstate2%densityPerturbation
  write(*,FMT) "Pressure Perturbation:", eatmos%pressurePerturbation, eatmos2%pressurePerturbation
  write(*,FMT) "Temperature Perturbation:", eatmos%temperaturePerturbation, eatmos2%temperaturePerturbation
  write(*,*)

  ! Get and print gases.
  call getGasesState_E(earth, gstate, argon, carbonDioxide, carbonMonoxide, dinitrogen, dioxygen, helium, hydrogen, methane, nitrogen, nitrousOxide, oxygen, ozone, water)
  call getGasesState_E(earth2, gstate2, argon2, carbonDioxide2, carbonMonoxide2, dinitrogen2, dioxygen2, helium2, hydrogen2, methane2, nitrogen2, nitrousOxide2, oxygen2, ozone2, water2)
  write(*,FMT) "Average Molecular Weight:", gstate%averageMolecularWeight, gstate2%averageMolecularWeight
  write(*,FMT) "O2 Mole Fraction:", dioxygen%moleFraction * 100.0, dioxygen2%moleFraction * 100.0
  write(*,FMT) "N2 Mole Fraction:", dinitrogen%moleFraction * 100.0, dinitrogen2%moleFraction * 100.0
  write(*,FMT) "CO2 Mole Fraction:", carbonDioxide%moleFraction * 100.0, carbonDioxide2%moleFraction * 100.0
  write(*,FMT) "Helium Mole Fraction:", helium%moleFraction * 100.0, helium2%moleFraction * 100.0
  write(*,FMT) "Argon Mole Fraction:", argon%moleFraction * 100.0, argon2%moleFraction * 100.0
  write(*,*)

  ! Get and print wind data
  call getWindsState_E(earth, wstate)
  call getWindsState_E(earth2, wstate2)
  write(*,FMT) "EW Wind:", wstate%ewWind, wstate2%ewWind
  write(*,FMT) "NS Wind:", wstate%nsWind, wstate2%nsWind
  write(*,FMT) "Vertical Wind:", wstate%verticalWind, wstate2%verticalWind
  write(*,FMT) "Perturbed EW Wind:", wstate%perturbedEWWind, wstate2%perturbedEWWind
  write(*,FMT) "Perturbed NS Wind:", wstate%perturbedNSWind, wstate2%perturbedNSWind
  write(*,FMT) "Perturbed Vertical Wind:", wstate%perturbedVerticalWind, wstate2%perturbedVerticalWind
  write(*,*)

  ! Print Earth specific data
  write(*,FMT) "Vapor Pressure:", eatmos%vaporPressure, eatmos2%vaporPressure
  write(*,FMT) "Vapor Density:", eatmos%vaporDensity, eatmos2%vaporDensity
  write(*,FMT) "Dew Point:", eatmos%dewPoint, eatmos2%dewPoint
  write(*,FMT) "Relative Humidity:", eatmos%relativeHumidity, eatmos2%relativeHumidity

  ! Memory is freed up in the deleteAtmosphere call.
  call deleteAtmosphere_E(earth)
  call deleteAtmosphere_E(earth2)

END PROGRAM
