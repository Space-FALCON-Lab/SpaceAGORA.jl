//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "PerturbedAtmosphere.h"
#include "NamelistReader.h"

using namespace std;
using namespace GRAM;

void atmosphereExample(PerturbedAtmosphere& atmosphere, const string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "An example of the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  atmosphere.setStartTime(ttime);

  // Set the position
  Position position;
  position.height = 50.0_km;
  position.latitude = 22.0_deg;
  position.longitude = 48.0_deg;
  position.longitude = 0.0_deg;
  position.elapsedTime = 100.0_sec;
  position.elapsedTime = 00.0_sec;
  atmosphere.setPosition(position);

  // Update the atmosphere data
  atmosphere.update();

  // Print the results
  AtmosphereState atmos = atmosphere.getAtmosphereState();
  EphemerisState ephem = atmosphere.getEphemerisState();

  cout << "Pressure: " << atmos.pressure << endl;
  cout << "Temperature: " << atmos.temperature << endl;
  cout << "Density: " << atmos.density << endl;
  if (atmos.argon.isPresent) {
    cout << "Argon mole fraction: " << atmos.argon.moleFraction << endl;
  }
  if (atmos.carbonDioxide.isPresent) {
    cout << "Carbon Dioxide mole fraction: " << atmos.carbonDioxide.moleFraction << endl;
  }
  if (atmos.carbonMonoxide.isPresent) {
    cout << "Carbon Monoxide mole fraction: " << atmos.carbonMonoxide.moleFraction << endl;
  }
  if (atmos.helium.isPresent) {
    cout << "Helium mole fraction: " << atmos.helium.moleFraction << endl;
  }
  if (atmos.dihydrogen.isPresent) {
    cout << "Molecular Hydrogen (H2) mole fraction: " << atmos.dihydrogen.moleFraction << endl;
  }
  if (atmos.hydrogen.isPresent) {
    cout << "Atomic Hydrogen (H) mole fraction: " << atmos.hydrogen.moleFraction << endl;
  }
  if (atmos.dinitrogen.isPresent) {
    cout << "Molecular Nitrogen (N2) mole fraction: " << atmos.dinitrogen.moleFraction << endl;
  }
  if (atmos.nitrogen.isPresent) {
    cout << "Atomic Nitrogen (N) mole fraction: " << atmos.nitrogen.moleFraction << endl;
  }
  if (atmos.methane.isPresent) {
    cout << "Methane mole fraction: " << atmos.methane.moleFraction << endl;
  }
  if (atmos.dioxygen.isPresent) {
    cout << "Molecular Oxygen (O2) mole fraction: " << atmos.dioxygen.moleFraction << endl;
  }
  if (atmos.oxygen.isPresent) {
    cout << "Atomic Oxygen (O) mole fraction: " << atmos.oxygen.moleFraction << endl;
  }
  cout << endl;
  cout << "Solar Time: " << ephem.solarTime << endl;
  cout << "Longitude of the Sun: " << setprecision(9) << ephem.longitudeSun << endl;
  cout << endl;
}
