//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "JupiterAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  tryGetSpicePath_J(spicePath, 500);
  initialize_J(spicePath);

  // Create a Jupiter atmosphere.  A pointer is returned that is our
  // handle on the Jupiter model.  More than one Jupiter atmosphere
  // can be created.  (Also see deleteAtmosphere_J() below.)
  JupiterAtmosphere_C* jupiter = createAtmosphere_J();
  JupiterAtmosphere_C* jupiter2 = createAtmosphere_J();

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_J(jupiter, version, bufferSize);
  printf("%s\n\n", version);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_J(jupiter, 1001);
  setMinRelativeStepSize_J(jupiter, 0.5);
  setPerturbationScales_J(jupiter, 1.5, 1.5, 1.5);
  setSeed_J(jupiter2, 1001);
  setMinRelativeStepSize_J(jupiter2, 0.5);
  setPerturbationScales_J(jupiter2, 1.5, 1.5, 1.5);

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
  setStartTime_J(jupiter, &time);
  setStartTime_J(jupiter2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_J(jupiter, &pos);
  pos.height = 1000.0;
  setPosition_J(jupiter2, &pos);


  // Update the model.
  update_J(jupiter); 
  update_J(jupiter2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_J(jupiter, &pos);
  getPosition_J(jupiter2, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                             Jupiter #1      Jupiter #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_J(jupiter, &estate);
  getEphemerisState_J(jupiter2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_J(jupiter, &dstate);
  getDynamicsState_J(jupiter2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  //DensityState_C rstate, rstate2;
  //getDensityState_J(jupiter, &rstate);
  //getDensityState_J(jupiter2, &rstate2);
  //printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  //printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  //printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  //printf("\n");

  // Get and print gases.
  // JupiterGRAM has no species model at this time.  
  //GasesState_C gstate, gstate2;
  //getGasesState_J(jupiter, &gstate);
  //getGasesState_J(jupiter2, &gstate2);
  //printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);

  // Get and print wind data
  //WindsState_C wstate, wstate2;
  //getWindsState_J(jupiter, &wstate);
  //getWindsState_J(jupiter2, &wstate2);
  //printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  //printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_J(jupiter);
  deleteAtmosphere_J(jupiter2);

  return 0;
}