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

using namespace std;
using namespace GRAM;

void perturbationExample(PerturbedAtmosphere& atmosphere, const string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "An example of Perturbations in the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;

  vector<Position> positions;
  vector<AtmosphereState> firstRun;
  vector<AtmosphereState> secondRun;
  vector<PerturbationState> pertStates;

  // Choose a seed for random pertubations
  atmosphere.setSeed(1111);

  // Set the start time of the trajectory
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  atmosphere.setStartTime(ttime);

  // Set the position
  Position position;
  position.height = 50.0_km;
  position.latitude = 22.0_deg;
  position.longitude = 48.0_deg;
  position.elapsedTime = 100.0_sec;
  atmosphere.setPosition(position);

  // Update the atmosphere data
  atmosphere.update();

  // Save the results
  positions.push_back(position);
  firstRun.push_back(atmosphere.getAtmosphereState());
  pertStates.push_back(atmosphere.getPerturbationState());

  for (int i = 0; i < 4; ++i) {
    // Change Position
    position.height += 1.0_km;
    position.elapsedTime += 1.0;
    atmosphere.setPosition(position);

    // Update the atmosphere data
    atmosphere.update();

    // Save the results
    positions.push_back(position);
    firstRun.push_back(atmosphere.getAtmosphereState());
    pertStates.push_back(atmosphere.getPerturbationState());
  }

  // Now do it all again, but supply the perturbation states
  // Reset the perturbation model by setting the seed.
  // To be certain that the override is working, use a different seed!
  atmosphere.setSeed(9999);
  for (int i = 0; i < 5; ++i) {
    // Set the position
    atmosphere.setPosition(positions[i]);
    // This time, set the perturbation state instead of computing it
    atmosphere.setPerturbationState(pertStates[i]);
    // Update the atmosphere
    atmosphere.update();
    // Save the results
    secondRun.push_back(atmosphere.getAtmosphereState());
  }

  // Print the results to the screen
  cout << "            First Run  Second Run   (Random Number, Rho)" << endl;
  for (int i = 0; i < 5; ++i) {
    cout << "Density:  " << setw(10) << firstRun[i].perturbedDensity << "  " << setw(10) << secondRun[i].perturbedDensity << "    (" << pertStates[i].densityRandomNumber << ", " << pertStates[i].densityRho << ")" << endl;
    cout << "EW Wind:  " << setw(10) << firstRun[i].perturbedEWWind << "  " << setw(10) << secondRun[i].perturbedEWWind << "    (" << pertStates[i].ewWindRandomNumber << ", " << pertStates[i].ewWindRho << ")" << endl;
    cout << "NS Wind:  " << setw(10) << firstRun[i].perturbedNSWind << "  " << setw(10) << secondRun[i].perturbedNSWind << "    (" << pertStates[i].nsWindRandomNumber << ", " << pertStates[i].nsWindRho << ")" << endl;
    cout << endl;
  }
}

