//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "NeptuneAtmosphere_C.h"
#include "UranusAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  // One call must be made for each planet to be created.
  char spicePath[500];
  tryGetSpicePath_C(spicePath, 500);
  initialize_N(spicePath);
  initialize_U(spicePath);

  // Create atmospheres.  A pointer is returned that is our handle
  // on the model.  More than one atmosphere can be created.  
  // Allocated memory is freed below using deleteAtmosphere_C().
  // The handles could also be stored as void pointers. In that case, 
  // any call to a _N or _U function would require casting.
  NeptuneAtmosphere_C* neptune = (NeptuneAtmosphere_C*)createAtmosphere_C(NEPTUNE_C);
  UranusAtmosphere_C* uranus = (UranusAtmosphere_C*)createAtmosphere_C(URANUS_C);

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_N(neptune, version, bufferSize);
  printf("%s\n", version);
  getVersionString_U(uranus, version, bufferSize);
  printf("%s\n\n", version);

  // Model specific settings must use the planet specific interface
  setMinMaxFactor_N(neptune, 0.5, 1);


  // Set the perturbation scale factor, relative step size, and the initial seed.
  setSeed_C(neptune, 1423);
  setPerturbationScales_C(neptune, 1.5, 0.5, 0.5, 0.0);
  setPerturbationScales_C(uranus, 1.5, 0.5, 0.5, 0.0);

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
  setStartTime_C(neptune, &time);
  setStartTime_C(uranus, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_C(neptune, &pos);
  pos.latitude = -22.0;
  setPosition_C(uranus, &pos);


  // Update the model.
  update_C(neptune); 
  update_C(uranus);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_C(neptune, &pos);
  getPosition_C(uranus, &pos2);
  const char* fmt = "%-25s %14.6e  %14.6e\n";
  printf("                               Neptune        Uranus\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_C(neptune, &estate);
  getEphemerisState_C(uranus, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_C(neptune, &dstate);
  getDynamicsState_C(uranus, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_C(neptune, &rstate);
  getDensityState_C(uranus, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  // The generic interface can only retrieve the state structure.  
  GasesState_C gstate, gstate2;
  getGasesState_C(neptune, &gstate);
  getGasesState_C(uranus, &gstate2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  
  // If species data is required, use the planet specific interface.
  ConstituentGas_C dihydrogen, dihydrogen2;
  ConstituentGas_C methane, methane2;
  ConstituentGas_C helium, helium2;
  ConstituentGas_C dinitrogen;
  getGasesState_N(neptune, &gstate, &dihydrogen, &methane, &helium, &dinitrogen);
  getGasesState_U(uranus, &gstate2, &dihydrogen2, &methane2, &helium2);
  printf(fmt, "Dihydrogen Mole Fraction:", dihydrogen.moleFraction * 100.0, dihydrogen2.moleFraction * 100.0);
  printf(fmt, "Methane Mole Fraction:", methane.moleFraction * 100.0, methane2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf(fmt, "Dinitrogen Mole Fraction:", dinitrogen.moleFraction * 100.0, 0.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_C(neptune, &wstate);
  getWindsState_C(uranus, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_C(neptune);
  deleteAtmosphere_C(uranus);

  return 0;
}