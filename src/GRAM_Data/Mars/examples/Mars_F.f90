!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.

PROGRAM Mars_F
  USE, INTRINSIC :: ISO_C_BINDING
  USE GramStructs
  USE MarsStructs
  USE MarsGRAM

  TYPE(C_PTR) :: mars, mars2
  TYPE(GramTime_F) :: time
  TYPE(Position_F) :: pos, pos2
  TYPE(DynamicsState_F) :: dstate, dstate2
  TYPE(DailyDynamicsState_F) :: ddstate, ddstate2
  TYPE(MarsState_F) ::  mstate, mstate2

  CHARACTER(KIND=C_CHAR, LEN = 500) :: spicePath
  CHARACTER(KIND=C_CHAR, LEN = 100) :: version
  CHARACTER(KIND=C_CHAR, LEN = 100) :: dataPath
  INTEGER(C_SIZE_T) :: bufferSize
  INTEGER(C_SIZE_T) :: stringLength
  INTEGER(C_INT) :: errorFlag;
  TYPE(GasesState_F) :: gstate, gstate2
  TYPE(ConstituentGas_F) :: carbonDioxide, carbonDioxide2
  TYPE(ConstituentGas_F) :: dinitrogen, dinitrogen2
  TYPE(ConstituentGas_F) :: argon, argon2
  TYPE(ConstituentGas_F) :: carbonMonoxide, carbonMonoxide2
  TYPE(ConstituentGas_F) :: dioxygen, dioxygen2
  TYPE(ConstituentGas_F) :: dihydrogen, dihydrogen2
  TYPE(ConstituentGas_F) :: hydrogen, hydrogen2
  TYPE(ConstituentGas_F) :: oxygen, oxygen2
  TYPE(ConstituentGas_F) :: helium, helium2
  TYPE(ConstituentGas_F) :: water, water2
  TYPE(DensityState_F) ::  rstate, rstate2
  TYPE(WindsState_F) :: wstate, wstate2
  TYPE(EphemerisState_F) :: estate, estate2
  CHARACTER(LEN=*), PARAMETER  :: FMT = "(A, T33, ES13.6, 1X, ES13.6)"

  write(*,*) "================="
  write(*,*) "A FORTRAN Example"
  write(*,*) "================="
  write(*,*)

  ! The first call must be to initialize the Spice path.
  call tryGetDataPaths_M(spicePath, dataPath, 500);
  call initialize_M(spicePath)

  ! Create a Mars atmosphere.  A pointer is returned that is our
  ! handle on the mars model.  More than one Mars atmosphere
  ! can be created.  (Also see deleteAtmosphere_M() below.)
  mars = createAtmosphere_M(dataPath)
  mars2 = createAtmosphere_M(dataPath)

  ! Get and print the version string.
  ! On input, bufferSize is the length of the character array.
  ! On output, bufferSize is the length of the version string.
  bufferSize = 100
  stringLength = getVersionString_M(mars, version, bufferSize)
  write(*,*) version(1:stringLength)
  write(*,*)

  ! Mars specific settings.
  call setMapYear_M(mars, 0)
  call setMGCMDustLevels_M(mars, 2.0D0, 0.0D0, 0.0D0)
  call setMapYear_M(mars2, 2)

  ! Set the initial seed, the relative step size, and the perturbation scale factors.
  call setSeed_M(mars, 1001)
  call setPerturbationScales_M(mars, 1.5D0, 1.5D0, 1.5D0)
  call setMinRelativeStepSize_M(mars, 0.5D0)
  call setSeed_M(mars2, 1001)
  call setPerturbationScales_M(mars2, 1.5D0, 1.5D0, 1.5D0)
  call setMinRelativeStepSize_M(mars2, 0.5D0)

  ! Set the start time.
  time%year = 2020
  time%month = 3
  time%day = 15
  time%hour = 0
  time%minute = 0
  time%seconds = 0.0D0
  time%frame = 1
  time%scale = 1
  call setStartTime_M(mars, time)
  call setStartTime_M(mars2, time)

  ! Set the position.
  pos%height = 2.0D0
  pos%latitude = 22.0D0
  pos%longitude = 48.0D0
  pos%elapsedTime = 100.0D0
  pos%isPlanetoCentric = 1
  call setPosition_M(mars, pos)
  call setPosition_M(mars2, pos)

  ! Update the model.
  errorFlag = update_M(mars)
  errorFlag = update_M(mars2)

  ! The position and other output is returned in the Position_F structure.
  call getPosition_M(mars, pos)
  call getPosition_M(mars2, pos2)
  write(*,*) "                                   Mars #1       Mars #2"
  write(*,FMT) "Total Radius:", pos%totalRadius, pos2%totalRadius
  write(*,FMT) "Gravity:", pos%gravity, pos2%gravity
  write(*,*)

  ! The ephemeris values are returned in the EphemerisState_F structure
  call getEphemerisState_M(mars, estate)
  call getEphemerisState_M(mars2, estate2)
  write(*,FMT) "Solar Time:", estate%solarTime, estate2%solarTime
  write(*,FMT) "Longitude of the Sun:", estate%longitudeSun, estate2%longitudeSun
  write(*,*)

  ! The atmosphere state is returned in multiple structures:
  ! DynamicsState_F, DensityState_F, WindsState_F, GasesState_F
  call getDynamicsState_M(mars, dstate)
  call getDynamicsState_M(mars2, dstate2)
  write(*,FMT) "Temperature:", dstate%temperature, dstate2%temperature
  write(*,FMT) "Pressure:", dstate%pressure, dstate2%pressure
  write(*,FMT) "Density:", dstate%density, dstate2%density
  write(*,FMT) "Pressure Scale Height:", dstate%pressureScaleHeight, dstate2%pressureScaleHeight
  write(*,FMT) "Density Scale Height:", dstate%densityScaleHeight, dstate2%densityScaleHeight
  write(*,*)

  ! Get and print perturbed density
  call getDensityState_M(mars, rstate)
  call getDensityState_M(mars2, rstate2)
  write(*,FMT) "Mean Density:", rstate%density, rstate2%density
  write(*,FMT) "Perturbed Density:", rstate%perturbedDensity, rstate2%perturbedDensity
  write(*,FMT) "Perturbation Percent:", rstate%densityPerturbation, rstate2%densityPerturbation
  write(*,*)

  ! Get and print gases.
  call getGasesState_M(mars, gstate, argon, carbonDioxide, carbonMonoxide, dihydrogen, dinitrogen, dioxygen, helium, hydrogen, oxygen, water);
  call getGasesState_M(mars2, gstate2, argon2, carbonDioxide2, carbonMonoxide2, dihydrogen2, dinitrogen2, dioxygen2, helium2, hydrogen2, oxygen2, water2);
  write(*,FMT) "Average Molecular Weight:", gstate%averageMolecularWeight, gstate2%averageMolecularWeight
  write(*,FMT) "Carbon Dioxide Mole Fraction:", carbonDioxide%moleFraction * 100.0D0, carbonDioxide2%moleFraction * 100.0D0
  write(*,FMT) "Dinitrogen Mole Fraction:", dinitrogen%moleFraction * 100.0D0, dinitrogen2%moleFraction * 100.0D0
  write(*,FMT) "Argon Mole Fraction:", argon%moleFraction * 100.0D0, argon2%moleFraction * 100.0D0
  write(*,FMT) "Carbon Monoxide Mole Fraction:", carbonMonoxide%moleFraction * 100.0D0, carbonMonoxide2%moleFraction * 100.0D0
  write(*,FMT) "Dioxygen Mole Fraction:", dioxygen%moleFraction * 100.0D0, dioxygen2%moleFraction * 100.0D0
  write(*,FMT) "Dihydrogen Mole Fraction:", dihydrogen%moleFraction * 100.0D0, dihydrogen2%moleFraction * 100.0D0
  write(*,FMT) "Hydrogen Mole Fraction:", hydrogen%moleFraction * 100.0D0, hydrogen2%moleFraction * 100.0D0
  write(*,FMT) "Oxygen Mole Fraction:", oxygen%moleFraction * 100.0D0, oxygen2%moleFraction * 100.0D0
  write(*,FMT) "Helium Mole Fraction:", helium%moleFraction * 100.0D0, helium2%moleFraction * 100.0D0
  write(*,FMT) "Water Vapor Mole Fraction:", water%moleFraction * 100.0D0, water2%moleFraction * 100.0D0
  write(*,*)

  ! Get and print wind data
  call getWindsState_M(mars, wstate)
  call getWindsState_M(mars2, wstate2)
  write(*,FMT) "EW Wind:", wstate%ewWind, wstate2%ewWind
  write(*,FMT) "NS Wind:", wstate%nsWind, wstate2%nsWind
  write(*,FMT) "Vertical Wind:", wstate%verticalWind, wstate2%verticalWind
  write(*,FMT) "Perturbed EW Wind:", wstate%perturbedEWWind, wstate2%perturbedEWWind
  write(*,FMT) "Perturbed NS Wind:", wstate%perturbedNSWind, wstate2%perturbedNSWind
  write(*,*)

  ! Get and print some Mars specific data
  call getMarsState_M(mars, mstate)
  call getMarsState_M(mars2, mstate2)
  write(*,FMT) "Dust Optical Depth:", mstate%dustOpticalDepth, mstate2%dustOpticalDepth
  call getDailyDynamicsState_M(mars, ddstate)
  call getDailyDynamicsState_M(mars2, ddstate2)
  write(*,FMT) "Daily Mean Temperature:", ddstate%temperatureDaily, ddstate2%temperatureDaily
  write(*,FMT) "Daily Maximum Temperature:", ddstate%temperatureMax, ddstate2%temperatureMax
  write(*,FMT) "Daily Minimum Temperature:", ddstate%temperatureMin, ddstate2%temperatureMin
  write(*,*)

  ! Memory is freed up in the deleteAtmosphere call.
  call deleteAtmosphere_M(mars)
  call deleteAtmosphere_M(mars2)

END PROGRAM
