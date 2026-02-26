//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "TitanAtmosphere.h"
#include "TitanNamelistReader.h" 

using namespace std;
using namespace GRAM;

// Helper function for printing
void print(const string& label, greal value1, greal value2)
{
  cout << left << setw(27) << label << right << setw(13) << value1 << "   " << setw(13) << value2 << endl;
}

int main(int argc, char** argv)
{
  cout << endl;
  cout << "=============" << endl;
  cout << "A C++ Example" << endl;
  cout << "=============" << endl;
  cout << endl;

  // Set the Spice data path (this is critical).
  TitanInputParameters inputParameters;
  TitanNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  // Create a Titan model and set default parameters.
  TitanAtmosphere titan;
  titan.setInputParameters(inputParameters);

  // More than one TitanAtmosphere can be created.
  TitanAtmosphere titan2;
  titan2.setInputParameters(inputParameters);

  // Get and print the version string.
  string version = titan.getVersionString();
  cout << version << endl;
  cout << endl;

  // Compare two different models
  titan.setModelType(Yelle97);
  titan.setMinMaxFactor(0.0, true);
  titan2.setModelType(GCM95);

  // Set the perturbation scale factor, relative step size, and the initial seed.
  titan.setPerturbationScales(1.5);
  titan.setMinRelativeStepSize(0.5);
  titan.setSeed(1001);
  titan2.setPerturbationScales(1.5);
  titan2.setMinRelativeStepSize(0.5);
  titan2.setSeed(1001);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  titan.setStartTime(ttime);
  titan2.setStartTime(ttime);

  // Set the position
  Position position;
  position.height = 50.0_km;
  position.latitude = 22.0_deg;
  position.longitude = 48.0_deg;
  position.elapsedTime = 100.0_sec;
  titan.setPosition(position);
  titan2.setPosition(position);

  // Update the atmosphere data
  titan.update();
  titan2.update();

  // The position and other output is returned in the Position class.
  const Position& pos = titan.getPosition();
  const Position& pos2 = titan2.getPosition();
  cout << "                             Titan #1      Titan #2" << scientific << endl;
  print("Total Radius: ", pos.totalRadius, pos2.totalRadius);
  print("Gravity: ", pos.gravity, pos2.gravity);
  cout << endl;

  // The ephemeris values are returned in the EphemerisState class.
  const EphemerisState& estate = titan.getEphemerisState();
  const EphemerisState& estate2 = titan2.getEphemerisState();
  print("Solar Time:", estate.solarTime, estate2.solarTime);
  print("Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  cout << endl;

  // The atmosphere values are returned in the AtmosphereState class.
  const AtmosphereState& atmos = titan.getAtmosphereState();
  const AtmosphereState& atmos2 = titan2.getAtmosphereState();
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
  print("Argon Mole Fraction:", atmos.argon.moleFraction * 100.0, atmos2.argon.moleFraction * 100.0);
  print("Methane Mole Fraction:", atmos.methane.moleFraction * 100.0, atmos2.methane.moleFraction * 100.0);
  print("Dinitrogen Mole Fraction:", atmos.dinitrogen.moleFraction * 100.0, atmos2.dinitrogen.moleFraction * 100.0);
  cout << endl;

  // Print wind data
  print("EW Wind:", atmos.ewWind, atmos2.ewWind);
  print("Perturbed EW Wind:", atmos.perturbedEWWind, atmos2.perturbedEWWind);

  return 0;
}

