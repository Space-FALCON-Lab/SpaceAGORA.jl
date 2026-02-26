//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "VenusAtmosphere_C.h"

int main(int argc, char** argv)
{
  printf("===========\n");
  printf("A C Example\n");
  printf("===========\n\n");

  // The first calls must be to initialize the Spice path.
  char spicePath[500];
  tryGetSpicePath_V(spicePath, 500);
  initialize_V(spicePath);

  // Create a Venus atmosphere.  A pointer is returned that is our
  // handle on the Venus model.  More than one Venus atmosphere
  // can be created.  (Also see deleteAtmosphere_V() below.)
  VenusAtmosphere_C* venus = createAtmosphere_V();
  VenusAtmosphere_C* venus2 = createAtmosphere_V();

  // Get and print the version string.
  size_t bufferSize = 100;
  char version[100];
  getVersionString_V(venus, version, bufferSize);
  printf("%s\n\n", version);

  // Set the initial seed, the relative step size, and the perturbation scale factors.
  setSeed_V(venus, 1001);
  setMinRelativeStepSize_V(venus, 0.5);
  setPerturbationScales_V(venus, 1.5, 1.5, 1.5);
  setSeed_V(venus2, 1001);
  setMinRelativeStepSize_V(venus2, 0.5);
  setPerturbationScales_V(venus2, 1.5, 1.5, 1.5);

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
  setStartTime_V(venus, &time);
  setStartTime_V(venus2, &time);

  // Set the position.
  Position_C pos;
  pos.height = 50.0;
  pos.latitude = 22.0;
  pos.longitude = 48.0;
  pos.elapsedTime = 100.0;
  pos.isPlanetoCentric = 1;
  setPosition_V(venus, &pos);
  pos.height = 200.0;
  setPosition_V(venus2, &pos);


  // Update the model.
  update_V(venus); 
  update_V(venus2);

  // The position and other output is returned in the Position_C structure.
  Position_C pos2;
  getPosition_V(venus, &pos);
  getPosition_V(venus2, &pos2);
  const char* fmt = "%-30s %14.6e  %14.6e\n";
  printf("                                   Venus #1        Venus #2\n");
  printf(fmt, "Total Radius:", pos.totalRadius, pos2.totalRadius);
  printf(fmt, "Gravity:", pos.gravity, pos2.gravity);
  printf("\n");

  // The ephemeris values are returned in the EphemerisState_C structure
  EphemerisState_C estate, estate2;
  getEphemerisState_V(venus, &estate);
  getEphemerisState_V(venus2, &estate2);
  printf(fmt, "Solar Time:", estate.solarTime, estate2.solarTime);
  printf(fmt, "Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  printf("\n");

  // The atmosphere state is returned in multiple structures:
  // DynamicsState_C, DensityState_C, WindsState_C, GasesState_C
  DynamicsState_C dstate, dstate2;
  getDynamicsState_V(venus, &dstate);
  getDynamicsState_V(venus2, &dstate2);
  printf(fmt, "Temperature:", dstate.temperature, dstate2.temperature);
  printf(fmt, "Pressure:", dstate.pressure, dstate2.pressure);
  printf(fmt, "Density:", dstate.density, dstate2.density);
  printf(fmt, "Pressure Scale Height:", dstate.pressureScaleHeight, dstate2.pressureScaleHeight);
  printf(fmt, "Density Scale Height:", dstate.densityScaleHeight, dstate2.densityScaleHeight);
  printf("\n");

  // Get and print perturbed density
  DensityState_C rstate, rstate2;
  getDensityState_V(venus, &rstate);
  getDensityState_V(venus2, &rstate2);
  printf(fmt, "Mean Density:", rstate.density, rstate2.density);
  printf(fmt, "Perturbed Density:", rstate.perturbedDensity, rstate2.perturbedDensity);
  printf(fmt, "Perturbation Percent:", rstate.densityPerturbation * 100.0, rstate2.densityPerturbation * 100.0);
  printf("\n");

  // Get and print gases.
  GasesState_C gstate, gstate2;
  ConstituentGas_C hydrogen, hydrogen2;
  ConstituentGas_C helium, helium2;
  ConstituentGas_C oxygen, oxygen2;
  ConstituentGas_C nitrogen, nitrogen2;
  ConstituentGas_C dinitrogen, dinitrogen2;
  ConstituentGas_C carbonMonoxide, carbonMonoxide2;
  ConstituentGas_C carbonDioxide, carbonDioxide2;
  getGasesState_V(venus, &gstate, &hydrogen, &helium, &oxygen, &nitrogen, &dinitrogen, &carbonMonoxide, &carbonDioxide);
  getGasesState_V(venus2, &gstate2, &hydrogen2, &helium2, &oxygen2, &nitrogen2, &dinitrogen2, &carbonMonoxide2, &carbonDioxide2);
  printf(fmt, "Average Molecular Weight:", gstate.averageMolecularWeight, gstate2.averageMolecularWeight);
  printf(fmt, "Hydrogen Mole Fraction:", hydrogen.moleFraction * 100.0, hydrogen2.moleFraction * 100.0);
  printf(fmt, "Helium Mole Fraction:", helium.moleFraction * 100.0, helium2.moleFraction * 100.0);
  printf(fmt, "Oxygen Mole Fraction:", oxygen.moleFraction * 100.0, oxygen2.moleFraction * 100.0);
  printf(fmt, "Nitrogen Mole Fraction:", nitrogen.moleFraction * 100.0, nitrogen2.moleFraction * 100.0);
  printf(fmt, "Dinitrogen Mole Fraction:", dinitrogen.moleFraction * 100.0, dinitrogen2.moleFraction * 100.0);
  printf(fmt, "Carbon Monoxide Mole Fraction:", carbonMonoxide.moleFraction * 100.0, carbonMonoxide2.moleFraction * 100.0);
  printf(fmt, "Carbon Dioxide Mole Fraction:", carbonDioxide.moleFraction * 100.0, carbonDioxide2.moleFraction * 100.0);
  printf("\n");

  // Get and print wind data
  WindsState_C wstate, wstate2;
  getWindsState_V(venus, &wstate);
  getWindsState_V(venus2, &wstate2);
  printf(fmt, "EW Wind:", wstate.ewWind, wstate2.ewWind);
  printf(fmt, "NS Wind:", wstate.nsWind, wstate2.nsWind);
  printf(fmt, "Perturbed EW Wind:", wstate.perturbedEWWind, wstate2.perturbedEWWind);
  printf(fmt, "Perturbed NS Wind:", wstate.perturbedNSWind, wstate2.perturbedNSWind);

  // Memory is freed up in the deleteAtmosphere call.
  deleteAtmosphere_V(venus);
  deleteAtmosphere_V(venus2);

  return 0;
}