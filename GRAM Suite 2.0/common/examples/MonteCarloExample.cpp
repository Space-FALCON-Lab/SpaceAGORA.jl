//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "MonteCarlo.h"
#include "Trajectory.h"
#include "ProfilePrinter.h"

using namespace std;
using namespace GRAM;

void monteCarloExample(PerturbedAtmosphere& atmosphere, const std::string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "A Monte Carlo example using the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;
  cout << "This may take a while..." << endl;

  // Start with a Monte Carlo object
  MonteCarlo monte;
  
  // The Monte Carlo requires a profile generator
  // Here we use a trajectory profile
  Trajectory trajProfile;

  // The profile requires an atmosphere
  trajProfile.setAtmosphere(atmosphere);

  // A start time is needed
  GramTime ttime;
//  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  ttime.setStartTime(1980, 11, 12, 0, 0, 0.0, UTC, ERT);
  atmosphere.setStartTime(ttime);

  // Also a starting position
  Position initialPosition;
  initialPosition.height = 0.0;
  initialPosition.latitude = 22.0_deg;
  initialPosition.longitude = 360.0_deg - 48.0_deg;
  initialPosition.elapsedTime = 0.0;
  trajProfile.setInitialPosition(initialPosition);

  // And a change in position
  Position deltaPosition;
  deltaPosition.longitude = -0.50_deg;
  deltaPosition.latitude = 0.30_deg;
  deltaPosition.height = 10.0;
  deltaPosition.elapsedTime = 500.00;
  trajProfile.setDeltaPosition(deltaPosition);

  // The number of data steps in each trajectory
  trajProfile.setNumberOfPoints(201);

  // The profile is ready, give it to the Monte Carlo object
  monte.setProfile(trajProfile);

  // The Monte Carlo also needs a profile printer
  ProfilePrinter printer;
  printer.setStyle(printer.GRAM_CSV_STYLE | printer.GRAM_MD_STYLE);
  printer.setWestLongitudePositive();
  monte.setProfilePrinter(printer);

  // You can control the initial seed
  monte.setInitialSeed(1001);

  // Set the number of samples to generate
  monte.setSampleSize(10);

  // Now run the Monte Carlo
  monte.generate();

  cout << "Files output: " << printer.getOutputFileNames() << endl;
  cout << endl;

}
