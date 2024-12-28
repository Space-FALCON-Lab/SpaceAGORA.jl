//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "JupiterAtmosphere.h" 
#include "JupiterNamelistReader.h" 

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
  JupiterInputParameters inputParameters;
  JupiterNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  // Create a Jupiter model and set default parameters.
  JupiterAtmosphere jupiter;
  jupiter.setInputParameters(inputParameters);

  // More than one JupiterAtmosphere can be created.
  JupiterAtmosphere jupiter2;
  jupiter2.setInputParameters(inputParameters);

  // Get and print the version string.
  string version = jupiter.getVersionString();
  cout << version << endl;
  cout << endl;

  // Set the perturbation scale factor, relative step size, and the initial seed.
  jupiter.setPerturbationScales(1.5);
  jupiter.setMinRelativeStepSize(0.5);
  jupiter.setSeed(1001);
  jupiter2.setPerturbationScales(1.5);
  jupiter2.setMinRelativeStepSize(0.5);
  jupiter2.setSeed(1001);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  jupiter.setStartTime(ttime);
  jupiter2.setStartTime(ttime);

  // Set the position
  Position position;
  position.height = 50.0_km;
  position.latitude = 22.0_deg;
  position.longitude = 48.0_deg;
  position.elapsedTime = 100.0_sec;
  jupiter.setPosition(position);
  position.height = 1000.0_km;
  jupiter2.setPosition(position);

  // Update the atmosphere data
  jupiter.update();
  jupiter2.update();

  // The position and other output is returned in the Position class.
  const Position& pos = jupiter.getPosition();
  const Position& pos2 = jupiter2.getPosition();
  cout << "                             Jupiter #1      Jupiter #2" << scientific << endl;
  print("Total Radius: ", pos.totalRadius, pos2.totalRadius);
  print("Gravity: ", pos.gravity, pos2.gravity);
  cout << endl;

  // The ephemeris values are returned in the EphemerisState class.
  const EphemerisState& estate = jupiter.getEphemerisState();
  const EphemerisState& estate2 = jupiter2.getEphemerisState();
  print("Solar Time:", estate.solarTime, estate2.solarTime);
  print("Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
  cout << endl;

  // The atmosphere values are returned in the AtmosphereState class.
  const AtmosphereState& atmos = jupiter.getAtmosphereState();
  const AtmosphereState& atmos2 = jupiter2.getAtmosphereState();
  print("Temperature:", atmos.temperature, atmos2.temperature);
  print("Pressure:", atmos.pressure, atmos2.pressure);
  print("Density:", atmos.density, atmos2.density);
  print("Pressure Scale Height:", atmos.pressureScaleHeight, atmos2.pressureScaleHeight);
  print("Density Scale Height:", atmos.densityScaleHeight, atmos2.densityScaleHeight);
  cout << endl;

  // Print perturbed density
  //print("Mean Density:", atmos.density, atmos2.density);
  //print("Perturbed Density:", atmos.perturbedDensity, atmos2.perturbedDensity);
  //print("Perturbation Percent:", atmos.densityPerturbation * 100.0, atmos2.densityPerturbation * 100.0);
  //cout << endl;

  // Get and print gases.
  //print("Average Molecular Weight:", atmos.averageMolecularWeight, atmos2.averageMolecularWeight);
  //print("Dihydrogen Mole Fraction:", atmos.dihydrogen.moleFraction * 100.0, atmos2.dihydrogen.moleFraction * 100.0);
  //print("Methane Mole Fraction:", atmos.methane.moleFraction * 100.0, atmos2.methane.moleFraction * 100.0);
  //print("Helium Mole Fraction:", atmos.helium.moleFraction * 100.0, atmos2.helium.moleFraction * 100.0);
  //print("Dinitrogen Mole Fraction:", atmos.dinitrogen.moleFraction * 100.0, atmos2.dinitrogen.moleFraction * 100.0);
  //cout << endl;

  // Print wind data
  //print("EW Wind:", atmos.ewWind, atmos2.ewWind);
  //print("Perturbed EW Wind:", atmos.perturbedEWWind, atmos2.perturbedEWWind);

  return 0;
}

