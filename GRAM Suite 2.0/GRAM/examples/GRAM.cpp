//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "gtest/gtest.h"
#include "MonteCarlo.h"
#include "Trajectory.h"
#include "ProfilePrinter.h"
#include "SpiceLoader.h"

#ifdef VENUS_GRAM
#include "VenusAtmosphere.h"
#include "VenusNamelistReader.h"
#endif
#ifdef EARTH_GRAM
#include "EarthAtmosphere.h"
#include "EarthNamelistReader.h"
#endif
#ifdef MARS_GRAM
#include "MarsAtmosphere.h"
#include "MarsNamelistReader.h"
#endif
#ifdef JUPITER_GRAM
#include "JupiterAtmosphere.h"
#include "JupiterNamelistReader.h"
#endif
#ifdef NEPTUNE_GRAM
#include "NeptuneAtmosphere.h"
#include "NeptuneNamelistReader.h"
#endif
#ifdef SATURN_GRAM
#include "SaturnAtmosphere.h"
#include "SaturnNamelistReader.h"
#endif
#ifdef TITAN_GRAM
#include "TitanAtmosphere.h"
#include "TitanNamelistReader.h"
#endif
#ifdef URANUS_GRAM
#include "UranusAtmosphere.h"
#include "UranusNamelistReader.h"
#endif

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  cout << "===========" << endl;
  cout << "All GRAMs  " << endl;
  cout << "===========\n" << endl;

  // Parse the command line.
  if (argc > 1) {
    string arg1(argv[1]);
    if (arg1 == "-test") {
      // Set the SPICE data path
      InputParameters inputParameters;
      NamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      SpiceLoader spiceLoader;
      spiceLoader.setInputParameters(inputParameters);

      // Run the unit tests.
      testing::InitGoogleTest(&argc, argv);
      return RUN_ALL_TESTS();
    }
    else {
      cout << "Unknown option." << endl;
      return 0;
    }
  }

  int option = 0;
  cout << "0. Quit Program" << endl;
#ifdef VENUS_GRAM
  VenusAtmosphere venus;
  cout << "2. Venus" << endl;
#endif
#ifdef EARTH_GRAM
  EarthAtmosphere earth;
  cout << "3. Earth" << endl;
#endif
#ifdef MARS_GRAM
  MarsAtmosphere mars;
  cout << "4. Mars" << endl;
#endif
#ifdef JUPITER_GRAM
  JupiterAtmosphere jupiter;
  cout << "5. Jupiter" << endl;
#endif
#ifdef SATURN_GRAM
  SaturnAtmosphere saturn;
  cout << "6. Saturn" << endl;
#endif
#ifdef URANUS_GRAM
  UranusAtmosphere uranus;
  cout << "7. Uranus" << endl;
#endif
#ifdef NEPTUNE_GRAM
  NeptuneAtmosphere neptune;
  cout << "8. Neptune" << endl;
#endif
#ifdef TITAN_GRAM
  TitanAtmosphere titan;
  cout << "10. Titan" << endl;
#endif
  cout << "\nSelect one: ";
  cin >> option;

  if (option == 0) {
    exit(0);
  }

  string namelistFileName;
  cout << "Enter the namelist file name: ";
  cin >> namelistFileName;

  PerturbedAtmosphere* body;
  string bodyName;

  // Process timer.
  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();
  std::string outFiles;

  try {

    switch (option) {
#ifdef VENUS_GRAM
    case 2: {
      VenusInputParameters inputParameters;
      VenusNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      venus.setInputParameters(inputParameters);
      body = &venus;
      bodyName = "Venus";
      break;
    }
#endif
#ifdef EARTH_GRAM
    case 3: {
      EarthInputParameters inputParameters;
      EarthNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      earth.setInputParameters(inputParameters);
      body = &earth;
      bodyName = "Earth";
      break;
    }
#endif
#ifdef MARS_GRAM
    case 4: {
      MarsInputParameters inputParameters;
      MarsNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      mars.setInputParameters(inputParameters);
      body = &mars;
      bodyName = "Mars";
      break;
    }
#endif
#ifdef JUPITER_GRAM
    case 5: {
      JupiterInputParameters inputParameters;
      JupiterNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      jupiter.setInputParameters(inputParameters);
      body = &jupiter;
      bodyName = "Jupiter";
      break;
    }
#endif
#ifdef SATURN_GRAM
    case 6: {
      SaturnInputParameters inputParameters;
      SaturnNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      saturn.setInputParameters(inputParameters);
      body = &saturn;
      bodyName = "Saturn";
      break;
    }
#endif
#ifdef URANUS_GRAM
    case 7: {
      UranusInputParameters inputParameters;
      UranusNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      uranus.setInputParameters(inputParameters);
      body = &uranus;
      bodyName = "Uranus";
      break;
    }
#endif
#ifdef NEPTUNE_GRAM
    case 8: {
      NeptuneInputParameters inputParameters;
      NeptuneNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      neptune.setInputParameters(inputParameters);
      body = &neptune;
      bodyName = "Neptune";
      break;
    }
#endif
#ifdef TITAN_GRAM
    case 10: {
      TitanInputParameters inputParameters;
      TitanNamelistReader reader;
      reader.tryGetSpicePath(inputParameters);
      reader.getParameters(namelistFileName, inputParameters);
      titan.setInputParameters(inputParameters);
      body = &titan;
      bodyName = "Titan";
      break;
    }
#endif
    default:
      cout << "Done." << endl;
      exit(0);
    }

    // Process timer.
    t1 = Clock::now();

    // Feedback.
    cout << "Reading " << namelistFileName << endl;

    // Get the input parameters (namelist) from the atmosphere
    const InputParameters& params = body->getInputParameters();

    if (!params.findDates) {
      cout << "Starting simulation for " << bodyName << "." << endl;

      // Start with a Monte Carlo object
      MonteCarlo monte;

      // The Monte Carlo requires a profile generator
      // Here we use a trajectory profile
      Trajectory trajProfile;
      // The trajectory requires an atmosphere
      trajProfile.setAtmosphere(*body);

      // The profile is ready, give it to the Monte Carlo object
      monte.setProfile(trajProfile);

      // The Monte Carlo also needs a profile printer
      ProfilePrinter printer;
      if (params.useLegacyOutputs) {
        printer.setStyle(printer.GRAM_MONTE_CARLO_STYLE);
      }
      else if (params.listFileName.empty() || params.listFileName == "null") {
        printer.setStyle(printer.GRAM_CSV_STYLE);
      }
      else {
        printer.setStyle(printer.GRAM_CSV_STYLE | printer.GRAM_MD_STYLE);
      }
      monte.setProfilePrinter(printer);

      // Now we pass the input parameters to the Monte Carlo object.
      // They will get passed on to the profile, atmosphere, and printer
      monte.setInputParameters(params);

      // Now run the Monte Carlo
      monte.generate();

      // Save the output file names.
      outFiles = printer.getOutputFileNames();
    }
    else {
      cout << "Finding dates for " << bodyName << "." << endl;

      // The start date has been set via the input paramters.
      // Now set the longitude (via position).
      Position position;
      position.setLongitudeDegrees(params.initialLongitude, params.isEastLongitudePositiveOnInput);
      body->setPosition(position);

      // Declare arrays for return data.
      GramTime gramTime[3];
      greal lonSun[3], tlst[3];

      // Search for the target dates.
      switch (option) {
#ifdef VENUS_GRAM
      case 2: {
        venus.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef EARTH_GRAM
      case 3: {
        earth.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef MARS_GRAM
      case 4: {
        mars.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef JUPITER_GRAM
      case 5: {
        jupiter.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef SATURN_GRAM
      case 6: {
        saturn.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef URANUS_GRAM
      case 7: {
        uranus.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef NEPTUNE_GRAM
      case 8: {
        neptune.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
#ifdef TITAN_GRAM
      case 10: {
        titan.findDates(params.targetLongitudeSun, params.targetSolarTime, gramTime, lonSun, tlst);
        break;
      }
#endif
      }
      // Print the longitude
      cout << endl;
      cout << "Longitude: " << params.initialLongitude << " degrees " << (params.isEastLongitudePositiveOnInput ? "E" : "W") << endl;

      // Print the target data.
      ProfilePrinter printer;
      printer.printDates(gramTime, lonSun, tlst);

    }
  }
  catch (const string& msg) {
    cerr << msg << endl;
    exit(1);
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

