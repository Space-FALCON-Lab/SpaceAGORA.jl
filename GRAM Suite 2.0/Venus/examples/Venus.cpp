//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "VenusAtmosphere.h" 
#include "VenusNamelistReader.h" 

using namespace std;
using namespace GRAM;

// Helper function for printing
void print(const string& label, greal value1, greal value2)
{
  cout << left << setw(32) << label << right << setw(13) << value1 << "   " << setw(13) << value2 << endl;
}

int main(int argc, char** argv)
{
  cout << endl;
  cout << "=============" << endl;
  cout << "A C++ Example" << endl;
  cout << "=============" << endl;
  cout << endl;

  // Set the Spice data path (this is critical).
  VenusInputParameters inputParameters;
  VenusNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  // Create a Venus model and set default parameters.
  VenusAtmosphere venus;
  venus.setInputParameters(inputParameters);

  // More than one VenusAtmosphere can be created.
  VenusAtmosphere venus2;
  venus2.setInputParameters(inputParameters);

  // Get and print the version string.
  string version = venus.getVersionString();
  cout << version << endl;
  cout << endl;

  // Set the perturbation scale factor, relative step size, and the initial seed.
  venus.setPerturbationScales(1.5);
  venus.setMinRelativeStepSize(0.5);
  venus.setSeed(1001);
  venus2.setPerturbationScales(1.5);
  venus2.setMinRelativeStepSize(0.5);
  venus2.setSeed(1001);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  venus.setStartTime(ttime);
  venus2.setStartTime(ttime);

  // Set the position
  Position position;
  position.height = 50.0_km;
  position.latitude = 22.0_deg;
  position.longitude = 48.0_deg;
  position.elapsedTime = 100.0_sec;
  venus.setPosition(position);
  position.height = 200.0_km;
  venus2.setPosition(position);

  // Update the atmosphere data
  venus.update();
  venus2.update();

  // The position and other output is returned in the Position class.
  const Position& pos = venus.getPosition();
  const Position& pos2 = venus2.getPosition();
  cout << "                                  Venus #1        Venus #2" << scientific << endl;
  print("Total Radius: ", pos.totalRadius, pos2.totalRadius);
  print("Gravity: ", pos.gravity, pos2.gravity);
  cout << endl;

  // The ephemeris values are returned in the EphemerisState class.
  const EphemerisState& estate = venus.getEphemerisState();
  const EphemerisState& estate2 = venus2.getEphemerisState();
  print("Solar Time:", estate.solarTime, estate2.solarTime);
  print("Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  cout << endl;

  // The atmosphere values are returned in the AtmosphereState class.
  const AtmosphereState& atmos = venus.getAtmosphereState();
  const AtmosphereState& atmos2 = venus2.getAtmosphereState();
  print("Temperature:", atmos.temperature, atmos2.temperature);
  print("Pressure:", atmos.pressure, atmos2.pressure);
  print("Density:", atmos.density, atmos2.density);
  print("Pressure Scale Height:", atmos.pressureScaleHeight, atmos2.pressureScaleHeight);
  print("Density Scale Height:", atmos.densityScaleHeight, atmos2.densityScaleHeight);
  cout << endl;

  // Print perturbed density
  print("Mean Density:", atmos.density, atmos2.density);
  print("Perturbed Density:", atmos.perturbedDensity, atmos2.perturbedDensity);
  print("Perturbation Percent:", atmos.densityPerturbation * 100.0, atmos2.densityPerturbation * 100.0);
  cout << endl;

  // Get and print gases.
  print("Average Molecular Weight:", atmos.averageMolecularWeight, atmos2.averageMolecularWeight);
  print("Hydrogen Mole Fraction:", atmos.hydrogen.moleFraction * 100.0, atmos2.hydrogen.moleFraction * 100.0);
  print("Helium Mole Fraction:", atmos.helium.moleFraction * 100.0, atmos2.helium.moleFraction * 100.0);
  print("Oxygen Mole Fraction:", atmos.oxygen.moleFraction * 100.0, atmos2.oxygen.moleFraction * 100.0);
  print("Nitrogen Mole Fraction:", atmos.nitrogen.moleFraction * 100.0, atmos2.nitrogen.moleFraction * 100.0);
  print("Dinitrogen Mole Fraction:", atmos.dinitrogen.moleFraction * 100.0, atmos2.dinitrogen.moleFraction * 100.0);
  print("Carbon Monoxide Mole Fraction:", atmos.carbonMonoxide.moleFraction * 100.0, atmos2.carbonMonoxide.moleFraction * 100.0);
  print("Carbon Dioxide Mole Fraction:", atmos.carbonDioxide.moleFraction * 100.0, atmos2.carbonDioxide.moleFraction * 100.0);
  cout << endl;

  // Print wind data
  print("EW Wind:", atmos.ewWind, atmos2.ewWind);
  print("NS Wind:", atmos.nsWind, atmos2.nsWind);
  print("Perturbed EW Wind:", atmos.perturbedEWWind, atmos2.perturbedEWWind);
  print("Perturbed NS Wind:", atmos.perturbedNSWind, atmos2.perturbedNSWind);

  return 0;
}

