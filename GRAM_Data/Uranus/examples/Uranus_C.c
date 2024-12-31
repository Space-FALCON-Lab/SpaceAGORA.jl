//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "UranusAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  tryGetSpicePath_U(spicePath, 500);
  initialize_U(spicePath);

  // Create a Uranus atmosphere.  A pointer is returned that is our
  // handle on the Uranus model.  More than one Uranus atmosphere
  // can be created.  (Also see deleteAtmosphere_U() below.)
  UranusAtmosphere_C* uranus = createAtmosphere_U();
  UranusAtmosphere_C* uranus2 = createAtmosphere_U();

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_U(uranus, version, bufferSize);
  printf("%s\n\n", version);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_U(uranus, 1001);
  setMinRelativeStepSize_U(uranus, 0.5);
  setPerturbationScales_U(uranus, 1.5, 1.5, 1.5);
  setSeed_U(uranus2, 1001);
  setMinRelativeStepSize_U(uranus2, 0.5);
  setPerturbationScales_U(uranus2, 1.5, 1.5, 1.5);

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
  setStartTime_U(uranus, &time);
  setStartTime_U(uranus2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_U(uranus, &pos);
  pos.height = 1000.0;
  setPosition_U(uranus2, &pos);


  // Update the model.
  update_U(uranus); 
  update_U(uranus2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_U(uranus, &pos);
  getPosition_U(uranus2, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                             Uranus #1       Uranus #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_U(uranus, &estate);
  getEphemerisState_U(uranus2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_U(uranus, &dstate);
  getDynamicsState_U(uranus2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_U(uranus, &rstate);
  getDensityState_U(uranus2, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C dihydrogen, dihydrogen2;
  ConstituentGas_C methane, methane2;
  ConstituentGas_C helium, helium2;
  getGasesState_U(uranus, &gstate, &dihydrogen, &methane, &helium);
  getGasesState_U(uranus2, &gstate2, &dihydrogen2, &methane2, &helium2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "Dihydrogen Mole Fraction:", dihydrogen.moleFraction * 100.0, dihydrogen2.moleFraction * 100.0);
  printf(fmt, "Methane Mole Fraction:", methane.moleFraction * 100.0, methane2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  //WindsState_C wstate, wstate2;
  //getWindsState_U(uranus, &wstate);
  //getWindsState_U(uranus2, &wstate2);
  //printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  //printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_U(uranus);
  deleteAtmosphere_U(uranus2);

  return 0;
}