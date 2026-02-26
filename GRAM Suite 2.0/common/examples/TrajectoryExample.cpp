//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "Trajectory.h"
#include "ProfilePrinter.h"

using namespace std;
using namespace GRAM;

void trajectoryExample(PerturbedAtmosphere& atmosphere, const std::string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "A trajectory example using the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;

  // Create a trajectory profile
  Trajectory trajProfile;

  // Set the atmosphere model
  atmosphere.setSeed(1001);     // Choose a seed for random pertubations
  atmosphere.setMinRelativeStepSize(0.5);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  atmosphere.setStartTime(ttime);

  // The trajectory needs an atmosphere model.
  trajProfile.setAtmosphere(atmosphere);

  // Set the initial position of the trajectory
  Position initialPosition;
  initialPosition.height = 0.0_km;
  initialPosition.latitude = 22.0_deg;
  initialPosition.setLongitude(48.0_deg, WEST_POSITIVE);
  initialPosition.elapsedTime = 0.0_sec;
  trajProfile.setInitialPosition(initialPosition);

  // Set the change in position for each update.
  Position deltaPosition;
  deltaPosition.height = 2.0_km;
  deltaPosition.latitude = 0.30_deg;
  deltaPosition.longitude = -0.50_deg;
  deltaPosition.elapsedTime = 500.0_sec;
  trajProfile.setDeltaPosition(deltaPosition);

  // Choose the number of datasteps to generate
  trajProfile.setNumberOfPoints(121);

  // Generate the trajectory
  trajProfile.generate();

  // Print the trajectory
  ProfilePrinter printer;
  printer.setStyle(printer.GRAM_MD_STYLE | printer.GRAM_CSV_STYLE);
  printer.setWestLongitudePositive();
  printer.print(atmosphere, trajProfile.getProfile());

  cout << "Files output: " << printer.getOutputFileNames() << endl;
  cout << endl;

}

