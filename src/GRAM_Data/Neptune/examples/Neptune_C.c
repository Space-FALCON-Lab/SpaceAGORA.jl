//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "NeptuneAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  tryGetSpicePath_N(spicePath, 500);
  initialize_N(spicePath);

  // Create a Neptune atmosphere.  A pointer is returned that is our
  // handle on the neptune model.  More than one Neptune atmosphere
  // can be created.  (Also see deleteAtmosphere_N() below.)
  NeptuneAtmosphere_C* neptune = createAtmosphere_N();
  NeptuneAtmosphere_C* neptune2 = createAtmosphere_N();

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_N(neptune, version, bufferSize);
  printf("%s\n\n", version);

  // Set the minMaxFactor for the models
  setMinMaxFactor_N(neptune, 0.5, 1);
  setMinMaxFactor_N(neptune2, -0.5, 1);

  // Set dinitrogen mole fraction for neptune2 only
  setDinitrogenMoleFraction_N(neptune2, 0.002);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_N(neptune, 1001);
  setMinRelativeStepSize_N(neptune, 0.5);
  setPerturbationScales_N(neptune, 1.5, 1.5, 1.5);
  setSeed_N(neptune2, 3333);
  setMinRelativeStepSize_N(neptune2, 0.5);
  setPerturbationScales_N(neptune2, 0.5, 0.5, 0.5);

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
  setStartTime_N(neptune, &time);
  setStartTime_N(neptune2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_N(neptune, &pos);
  pos.latitude = -22.0;
  setPosition_N(neptune2, &pos);


  // Update the model.
  update_N(neptune); 
  update_N(neptune2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_N(neptune, &pos);
  getPosition_N(neptune2, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                             Neptune #1      Neptune #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_N(neptune, &estate);
  getEphemerisState_N(neptune2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_N(neptune, &dstate);
  getDynamicsState_N(neptune2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_N(neptune, &rstate);
  getDensityState_N(neptune2, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C dihydrogen, dihydrogen2;
  ConstituentGas_C methane, methane2;
  ConstituentGas_C helium, helium2;
  ConstituentGas_C dinitrogen, dinitrogen2;
  getGasesState_N(neptune, &gstate, &dihydrogen, &methane, &helium, &dinitrogen);
  getGasesState_N(neptune2, &gstate2, &dihydrogen2, &methane2, &helium2, &dinitrogen2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "Dihydrogen Mole Fraction:", dihydrogen.moleFraction * 100.0, dihydrogen2.moleFraction * 100.0);
  printf(fmt, "Methane Mole Fraction:", methane.moleFraction * 100.0, methane2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf(fmt, "Dinitrogen Mole Fraction:", dinitrogen.moleFraction * 100.0, dinitrogen2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_N(neptune, &wstate);
  getWindsState_N(neptune2, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_N(neptune);
  deleteAtmosphere_N(neptune2);

  return 0;
}