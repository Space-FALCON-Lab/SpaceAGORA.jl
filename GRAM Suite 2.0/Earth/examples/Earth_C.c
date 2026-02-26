//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "EarthAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  char dataPath[500];
  tryGetDataPaths_E(spicePath, dataPath, 500);
  initialize_E(spicePath);

  // Create a Earth atmosphere.  A pointer is returned that is our
  // handle on the Earth model.  More than one Earth atmosphere
  // can be created.  (Also see deleteAtmosphere_E() below.)
  EarthAtmosphere_C* earth = createAtmosphere_E(dataPath);
  EarthAtmosphere_C* earth2 = createAtmosphere_E(dataPath);

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_E(earth, version, bufferSize);
  printf("%s\n\n", version);

  // Set the MERRA-2 hour and data extents
  setMERRA2Parameters_E(earth, 1, 0.0, 40.0, 20.0, 60.0);
  setMERRA2Parameters_E(earth2, 9, 0.0, 40.0, 20.0, 60.0);

  // For NCEP model, uncomment lines below
  //setUseNCEP_E(earth, true);
  //setUseNCEP_E(earth2, true);
  //setNCEPParameters_E(earth, 9715, 1);
  //setNCEPParameters_E(earth2, 9715, 5);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_E(earth, 1001);
  setPerturbationScales_E(earth, 1.5, 1.5, 1.5);
  setSeed_E(earth2, 1001);
  setPerturbationScales_E(earth2, 1.5, 1.5, 1.5);

  // Set the start time.
  GramTime_C time;
  time.year = 2020;
  time.month = 12;
  time.day = 15;
  time.hour = 0;
  time.minute = 0;
  time.seconds = 0.0;
  time.frame = 0;
  time.scale = 1;
  setStartTime_E(earth, &time);
  setStartTime_E(earth2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 5.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_E(earth, &pos);
  setPosition_E(earth2, &pos);


  // Update the model.
  update_E(earth); 
  update_E(earth2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_E(earth, &pos);
  getPosition_E(earth2, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                              Earth #1        Earth #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_E(earth, &estate);
  getEphemerisState_E(earth2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_E(earth, &dstate);
  getDynamicsState_E(earth2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_E(earth, &rstate);
  getDensityState_E(earth2, &rstate2);
  EarthState_C eatmos, eatmos2;
  getEarthState_E(earth, &eatmos);
  getEarthState_E(earth2, &eatmos2);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbed Pressure:", eatmos.perturbedPressure, eatmos2.perturbedPressure);
  printf(fmt, "Perturbed Temperature:", eatmos.perturbedTemperature, eatmos2.perturbedTemperature);
  printf(fmt, "Density Perturbation:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf(fmt, "Pressure Perturbation:", eatmos.pressurePerturbation * 100.0, eatmos2.pressurePerturbation * 100.0);
  printf(fmt, "Temperature Perturbation:", eatmos.temperaturePerturbation * 100.0, eatmos2.temperaturePerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C argon, carbonDioxide, carbonMonoxide, dinitrogen, dioxygen, 
    helium, hydrogen, methane, nitrogen, nitrousOxide, oxygen, ozone, water;
  ConstituentGas_C argon2, carbonDioxide2, carbonMonoxide2, dinitrogen2, dioxygen2, 
    helium2, hydrogen2, methane2, nitrogen2, nitrousOxide2, oxygen2, ozone2, water2;
  getGasesState_E(earth, &gstate, &argon, &carbonDioxide, &carbonMonoxide, &dinitrogen, &dioxygen, 
    &helium, &hydrogen, &methane, &nitrogen, &nitrousOxide, &oxygen, &ozone, &water);
  getGasesState_E(earth2, &gstate2, &argon2, &carbonDioxide2, &carbonMonoxide2, &dinitrogen2, &dioxygen2, 
    &helium2, &hydrogen2, &methane2, &nitrogen2, &nitrousOxide2, &oxygen2, &ozone2, &water2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "O2 Mole Fraction:", dioxygen.moleFraction * 100.0, dioxygen2.moleFraction * 100.0);
  printf(fmt, "N2 Mole Fraction:", dinitrogen.moleFraction * 100.0, dinitrogen2.moleFraction * 100.0);
  printf(fmt, "CO2 Mole Fraction:", carbonDioxide.moleFraction * 100.0, carbonDioxide2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf(fmt, "Argon Mole Fraction:", argon.moleFraction * 100.0, argon2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_E(earth, &wstate);
  getWindsState_E(earth2, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "NS Wind:", wstate.nsWind, wstate2.nsWind);
  printf(fmt, "Vertical Wind:", wstate.verticalWind, wstate2.verticalWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);
  printf(fmt, "Perturbed NS Wind:", wstate.perturbedNSWind, wstate2.perturbedNSWind);
  printf(fmt, "Perturbed Vertical Wind:", wstate.perturbedVerticalWind, wstate2.perturbedVerticalWind);
  printf("\n");

  // Print Earth specific data
  printf(fmt, "Vapor Pressure:", eatmos.vaporPressure, eatmos2.vaporPressure);
  printf(fmt, "Vapor Density:", eatmos.vaporDensity, eatmos2.vaporDensity);
  printf(fmt, "Dew Point:", eatmos.dewPoint, eatmos2.dewPoint);
  printf(fmt, "Relative Humidity:", eatmos.relativeHumidity, eatmos2.relativeHumidity);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_E(earth);
  deleteAtmosphere_E(earth2);

  return 0;
}