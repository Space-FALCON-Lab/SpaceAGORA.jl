//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "gtest/gtest.h"
#include "EarthAtmosphere.h"
#include "EarthNamelistReader.h"
#include "MonteCarlo.h"
#include "BallisticTrajectory.h"
#include "EarthProfilePrinter.h"
#include "SpiceLoader.h"
#include "EarthCorrelator.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  cout << "===========" << endl;
  cout << "The EarthCorrTraj Program" << endl;
  cout << EarthAtmosphere::getVersionString() << endl;
  cout << "===========\n" << endl;

  // Parse the command line.
  string namelistFileName;
  if (argc > 1) {
    string arg1(argv[1]);
    // The namelist filename is supplied on the command line.
    if (arg1 == "-file") {
      if (argc > 2) {
        namelistFileName = argv[2];
      }
    }
    // Command line arg to run the unit tests.
    else if (arg1 == "-test") {
      // Set the SPICE data path
      InputParameters inputParameters;
      NamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      SpiceLoader spiceLoader;
      spiceLoader.setSpiceDataPath(inputParameters.spicePath);

      // Run the unit tests.
      testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
    }
    else {
      cout << "Unknown option." << endl;
      return 0;
    }
  }
  else {
    // No command line options. So prompt the user.
    cout << "Enter the namelist file name: ";
    cin >> namelistFileName;
  }

  // Process timer.
  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  // Feedback.
  cout << "Reading " << namelistFileName << endl;
  std::string outFiles;

  try {
    // Load parameters from the namelist file.
    EarthInputParameters inputParameters;
    EarthNamelistReader reader;
    // See if there is a SPICE path override first.
    reader.tryGetSpicePath(inputParameters);
    reader.getParameters(namelistFileName, inputParameters);

    cout << "Starting simulation." << endl;

    // Start with a Monte Carlo object
    MonteCarlo monte;

    // The Monte Carlo requires a profile generator
    // Here we use a trajectory profile
    BallisticTrajectory trajProfile;

    // The trajectory requires an atmosphere
    EarthAtmosphere earth;
    earth.setInputParameters(inputParameters);
    trajProfile.setAtmosphere(earth);

    EarthCorrelator correlator;
    correlator.setMean(false);
    trajProfile.setCorrelator(correlator);

    // The profile is ready, give it to the Monte Carlo object
    monte.setProfile(trajProfile);

    // The Monte Carlo also needs a profile printer
    EarthProfilePrinter printer;
    if (inputParameters.useLegacyOutputs) {
      // User wants legacy outputs.
      printer.setStyle(printer.EARTH_TRAJ_STYLE);
    }
    else if (inputParameters.listFileName.empty() || inputParameters.listFileName == "null") {
      // List output suppressed.
      printer.setStyle(printer.GRAM_CSV_STYLE);
    }
    else {
      // Default is LIST and CSV output.
      printer.setStyle(printer.GRAM_CSV_STYLE | printer.GRAM_MD_STYLE);
    }
    monte.setProfilePrinter(printer);

    // Now we pass the input parameters to the Monte Carlo object.
    // They will get passed on to the profile, atmosphere, and printer
    monte.setInputParameters(inputParameters);

    // Now run the Monte Carlo
    monte.generate();

    // Save the output file names.
    outFiles = printer.getOutputFileNames();
  }
  catch (const string& msg) {
    cerr << msg << endl;
  }
  catch (const std::runtime_error& err) {
    cerr << "An unanticipated error occurred." << endl;
    cerr << err.what() << endl;
  }

  // Report the processor time.
  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0 << " seconds." << endl;

  // Report the output files
  if (!outFiles.empty()) {
    cout << "Files output: " << outFiles << endl;
  }
  cout << endl;

  return 0;
}

