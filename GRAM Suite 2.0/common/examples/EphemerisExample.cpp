//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "EphemerisProfile.h"
#include "ProfilePrinter.h"

using namespace std;
using namespace GRAM;

void ephemerisExample(PerturbedAtmosphere& atmosphere, const std::string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "An Ephemeris Profile example using the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;

  // Create a trajectory profile
  EphemerisProfile ephemerisProfile;

  // Set the atmosphere model
  atmosphere.setSeed(1001);     // Choose a seed for random pertubations
  atmosphere.setMinRelativeStepSize(0.5);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  atmosphere.setStartTime(ttime);

  ephemerisProfile.setAtmosphere(atmosphere);

  // Set the position.
  Position position;
  position.height = 5.0_km;
  position.latitude = 82.0_deg;
  position.setLongitude(48.0_deg);
  position.elapsedTime = 0.0_sec;
  ephemerisProfile.setPosition(position);

  // Set the initial ephemeris values.
  EphemerisState initialEphemeris;
  initialEphemeris.longitudeSun = 0.0_deg;
  initialEphemeris.solarTime = 12.0_hr;
  ephemerisProfile.setInitialEphemeris(initialEphemeris);

  // Set the change in ephemeris values.
  EphemerisState deltaEphemeris;
  deltaEphemeris.longitudeSun = 10.0_deg;
  deltaEphemeris.solarTime = 0.0_hr;
  ephemerisProfile.setDeltaEphemeris(deltaEphemeris);

  // Choose the number of datasteps to generate
  ephemerisProfile.setNumberOfPoints(37);

  // Generate the trajectory
  ephemerisProfile.generate();

  // Print the trajectory
  ProfilePrinter printer;
  printer.setStyle(printer.GRAM_CSV_STYLE | printer.GRAM_MD_STYLE | printer.GRAM_EPHEMERIS_STYLE);
  printer.print(atmosphere, ephemerisProfile.getProfile());

  cout << "Files output: " << printer.getOutputFileNames() << endl;
  cout << endl;

}

