//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "TitanAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  tryGetSpicePath_T(spicePath, 500);
  initialize_T(spicePath);

  // Create a Titan atmosphere.  A pointer is returned that is our
  // handle on the Titan model.  More than one Titan atmosphere
  // can be created.  (Also see deleteAtmosphere_T() below.)
  TitanAtmosphere_C* titan = createAtmosphere_T();
  TitanAtmosphere_C* titan2 = createAtmosphere_T();

  size_t bufferSize = 100;
  char version[100];
  getVersionString_T(titan, version, bufferSize);
  printf("%s\n\n", version);
  
  // Compare two different models
  setModelType_T(titan, 1);  // Yelle97
  setMinMaxFactor_T(titan, 0.0, 1);
  setModelType_T(titan2, 2); // GCM95

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_T(titan, 1001);
  setMinRelativeStepSize_T(titan, 0.5);
  setPerturbationScales_T(titan, 1.5, 1.5, 1.5);
  setSeed_T(titan2, 1001);
  setMinRelativeStepSize_T(titan2, 0.5);
  setPerturbationScales_T(titan2, 1.5, 1.5, 1.5);

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
  setStartTime_T(titan, &time);
  setStartTime_T(titan2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_T(titan, &pos);
  setPosition_T(titan2, &pos);


  // Update the model.
  update_T(titan); 
  update_T(titan2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_T(titan, &pos);
  getPosition_T(titan2, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                              Titan #1        Titan #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_T(titan, &estate);
  getEphemerisState_T(titan2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_T(titan, &dstate);
  getDynamicsState_T(titan2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_T(titan, &rstate);
  getDensityState_T(titan2, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C argon, argon2;
  ConstituentGas_C methane, methane2;
  ConstituentGas_C dinitrogen, dinitrogen2;
  getGasesState_T(titan, &gstate, &argon, &methane, &dinitrogen);
  getGasesState_T(titan2, &gstate2, &argon2, &methane2, &dinitrogen2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "Argon Mole Fraction:", argon.moleFraction * 100.0, argon2.moleFraction * 100.0);
  printf(fmt, "Methane Mole Fraction:", methane.moleFraction * 100.0, methane2.moleFraction * 100.0);
  printf(fmt, "Dinitrogen Mole Fraction:", dinitrogen.moleFraction * 100.0, dinitrogen2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_T(titan, &wstate);
  getWindsState_T(titan2, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_T(titan);
  deleteAtmosphere_T(titan2);

  return 0;
}