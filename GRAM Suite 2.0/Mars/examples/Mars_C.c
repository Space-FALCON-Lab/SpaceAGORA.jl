//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "MarsAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500] = "";
  char dataPath[500] = "";
  tryGetDataPaths_M(spicePath, dataPath, 500);
  initialize_M(spicePath);

  // Create a Mars atmosphere.  A pointer is returned that is our
  // handle on the Mars model.  More than one Mars atmosphere
  // can be created.  (Also see deleteAtmosphere_M() below.)
  MarsAtmosphere_C* mars = createAtmosphere_M(dataPath);
  // Since the data path is critical, make sure creation worked.
  if (mars == NULL) {
    size_t bufferSize = 1000;
    char msg[1000];
    getErrorMessage_M(msg, bufferSize);
    printf("%s\n\n", msg);
    exit(0);
  }
  MarsAtmosphere_C* mars2 = createAtmosphere_M(dataPath);
  if (mars2 == NULL) {
    size_t bufferSize = 1000;
    char msg[1000];
    getErrorMessage_M(msg, bufferSize);
    printf("%s\n\n", msg);
    exit(0);
  }

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_M(mars, version, bufferSize);
  printf("%s\n\n", version);

  // Mars specific settings.
  setMapYear_M(mars, 0);
  setMGCMDustLevels_M(mars, 2.0, 0.0, 0.0);
  setMapYear_M(mars2, 2);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_M(mars, 1001);
  setPerturbationScales_M(mars, 1.5, 1.5, 1.5, 1.5);
  setMinRelativeStepSize_M(mars, 0.5);
  setSeed_M(mars2, 1001);
  setPerturbationScales_M(mars2, 1.5, 1.5, 1.5, 1.5);
  setMinRelativeStepSize_M(mars2, 0.5);

  // Set the start time.
  GramTime_C time;
  time.year = 2020;
  time.month = 3;
  time.day = 15;
  time.hour = 0;
  time.minute = 0;
  time.seconds = 0.0;
  time.frame = 1;
  time.scale = 1;
  setStartTime_M(mars, &time);
  setStartTime_M(mars2, &time);

  // Set the position.
  Position_C pos = { 0 };
  pos.height = 2.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_M(mars, &pos);
  setPosition_M(mars2, &pos);

  // Update the model.
  update_M(mars); 
  update_M(mars2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_M(mars, &pos);
  getPosition_M(mars2, &pos2);
  const char* fmt = "%-33s %14.6e  %14.6e\n";
  printf("                                      Mars #1         Mars #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_M(mars, &estate);
  getEphemerisState_M(mars2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_M(mars, &dstate);
  getDynamicsState_M(mars2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_M(mars, &rstate);
  getDensityState_M(mars2, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C carbonDioxide, carbonDioxide2;
  ConstituentGas_C dinitrogen, dinitrogen2;
  ConstituentGas_C argon, argon2;
  ConstituentGas_C carbonMonoxide, carbonMonoxide2;
  ConstituentGas_C dioxygen, dioxygen2;
  ConstituentGas_C dihydrogen, dihydrogen2;
  ConstituentGas_C hydrogen, hydrogen2;
  ConstituentGas_C oxygen, oxygen2;
  ConstituentGas_C helium, helium2;
  ConstituentGas_C water, water2;
  getGasesState_M(mars, &gstate, &argon, &carbonDioxide, &carbonMonoxide, &dihydrogen, &dinitrogen,
                  &dioxygen, &helium, &hydrogen, &oxygen, &water);
  getGasesState_M(mars2, &gstate2, &argon2, &carbonDioxide2, &carbonMonoxide2, &dihydrogen2, &dinitrogen2,
                  &dioxygen2, &helium2, &hydrogen2, &oxygen2, &water2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "Carbon Dioxide Mole Fraction:", carbonDioxide.moleFraction * 100.0, carbonDioxide2.moleFraction * 100.0);
  printf(fmt, "Dinitrogen Mole Fraction:", dinitrogen.moleFraction * 100.0, dinitrogen2.moleFraction * 100.0);
  printf(fmt, "Argon Mole Fraction:", argon.moleFraction * 100.0, argon2.moleFraction * 100.0);
  printf(fmt, "Carbon Monoxide Mole Fraction:", carbonMonoxide.moleFraction * 100.0, carbonMonoxide2.moleFraction * 100.0);
  printf(fmt, "Dioxygen Mole Fraction:", dioxygen.moleFraction * 100.0, dioxygen2.moleFraction * 100.0);
  printf(fmt, "Dihydrogen Mole Fraction:", dihydrogen.moleFraction * 100.0, dihydrogen2.moleFraction * 100.0);
  printf(fmt, "Hydrogen Mole Fraction:", hydrogen.moleFraction * 100.0, hydrogen2.moleFraction * 100.0);
  printf(fmt, "Oxygen Mole Fraction:", oxygen.moleFraction * 100.0, oxygen2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf(fmt, "Water Vapor Mole Fraction:", water.moleFraction * 100.0, water2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_M(mars, &wstate);
  getWindsState_M(mars2, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "NS Wind:", wstate.nsWind, wstate2.nsWind);
  printf(fmt, "Vertical Wind:", wstate.verticalWind, wstate2.verticalWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);
  printf(fmt, "Perturbed NS Wind:", wstate.perturbedNSWind, wstate2.perturbedNSWind);
  printf("\n");

  // Get and print some Mars specific data
  MarsState_C mstate, mstate2;
  getMarsState_M(mars, &mstate);
  getMarsState_M(mars2, &mstate2);
  printf(fmt, "Dust Optical Depth:", mstate.dustOpticalDepth, mstate2.dustOpticalDepth);
  DailyDynamicsState_C ddstate, ddstate2;
  getDailyDynamicsState_M(mars, &ddstate);
  getDailyDynamicsState_M(mars2, &ddstate2);
  printf(fmt, "Daily Mean Temperature:", ddstate.temperatureDaily, ddstate2.temperatureDaily);
  printf(fmt, "Daily Maximum Temperature:", ddstate.temperatureMax, ddstate2.temperatureMax);
  printf(fmt, "Daily Minimum Temperature:", ddstate.temperatureMin, ddstate2.temperatureMin);
  printf("\n");

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_M(mars);
  deleteAtmosphere_M(mars2);
}