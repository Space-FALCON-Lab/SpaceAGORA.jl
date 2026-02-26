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

void namelistExample(PerturbedAtmosphere& atmosphere, const std::string& name)
{
  cout << endl;
  cout << "======================================================================" << endl;
  cout << "A Namelist example using the " << name << " atmosphere model." << endl;
  cout << "======================================================================" << endl;
  cout << endl;
  cout << "This may take a while..." << endl;

  // Get the input parameters (namelist) from the atmosphere
  const InputParameters& params = atmosphere.getInputParameters();

  // Start with a Monte Carlo object
  MonteCarlo monte;
  
  // The Monte Carlo requires a profile generator
  // Here we use a trajectory profile
  Trajectory trajProfile;
  // The trajectory requires an atmosphere
  trajProfile.setAtmosphere(atmosphere);

  // The profile is ready, give it to the Monte Carlo object
  monte.setProfile(trajProfile);

  // The Monte Carlo also needs a profile printer
  ProfilePrinter printer;
  printer.setStyle(printer.GRAM_MD_STYLE | printer.GRAM_CSV_STYLE);
  monte.setProfilePrinter(printer);

  // Now we pass the input parameters to the Monte Carlo object.
  // They will get passed on to the profile, atmosphere, and printer
  monte.setInputParameters(params);

  // Now run the Monte Carlo
  monte.generate();

  cout << "Files output: " << printer.getOutputFileNames() << endl;
  cout << endl;

}

