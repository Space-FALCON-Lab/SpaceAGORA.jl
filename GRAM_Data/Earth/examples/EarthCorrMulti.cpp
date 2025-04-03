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
#include "EarthProfilePrinter.h"
#include "SpiceLoader.h"
#include "EarthCorrelator.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  cout << "===========" << endl;
  cout << "The EarthCorrMulti Example Program" << endl;
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

    ofstream v1out("Vehicle_1.txt");
    ofstream v2out_un("Vehicle_2_uncorrelated.txt");
    ofstream v2out_corr("Vehicle_2_correlated.txt");

    // The trajectory requires an atmosphere
    EarthAtmosphere earth1, earth2;
    earth1.setInputParameters(inputParameters);
    earth2.setInputParameters(inputParameters);

    Position position1, position2;
    position1.height = 0.0;
    position1.latitude = 28.31;
    position1.longitude = -80.55 + 360.0;
    position1.elapsedTime = 0.0;

    earth1.setSeed(inputParameters.initialRandomSeed);

    earth1.setPosition(position1);
    earth1.update();
    const AtmosphereState& atm = earth1.getAtmosphereState();
    const EarthAtmosphereState& eatm = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();

    v1out << position1.height << " " << position1.latitude << " " << position1.longitude - 360.0 << " "
      << atm.perturbedEWWind << " " << atm.perturbedNSWind << " " << eatm.perturbedTemperature << endl;

    EarthCorrelator correlator;
    correlator.setMean(false);

    for (int i = 0; i < 50; ++i) {
      position1.height += 0.5;
      position1.latitude += 0.1;
      position1.longitude += 0.1;
      position1.elapsedTime += 30.0;

      earth1.setPosition(position1);
      earth1.update();
      const AtmosphereState& atm1 = earth1.getAtmosphereState();
      const EarthAtmosphereState& eatm1 = atm1.getPlanetSpecificMetrics<EarthAtmosphereState>();

      v1out << position1.height << " " << position1.latitude << " " << position1.longitude - 360.0 << " "
        << atm.perturbedEWWind << " " << atm.perturbedNSWind << " " << eatm1.perturbedTemperature << endl;

      if (position1.height == 10.0) {
        position2 = position1;
        // Reseed to get different randomness from vehicle 1.
        earth2.setSeed(inputParameters.initialRandomSeed);
        earth2.setPosition(position2);
        earth2.update();
        AtmosphereState atm2 = earth2.getAtmosphereState();
        const EarthAtmosphereState& eatm2 = atm2.getPlanetSpecificMetrics<EarthAtmosphereState>();

        v2out_un << position2.height << " " << position2.latitude << " " << position2.longitude - 360.0 << " "
          << atm2.perturbedEWWind << " " << atm2.perturbedNSWind << " " << eatm2.perturbedTemperature << endl;

        correlator.updateCoefficients(position1, position2);
        correlator.correlate(atm1, atm2);

        v2out_corr << position2.height << " " << position2.latitude << " " << position2.longitude - 360.0 << " "
          << atm2.perturbedEWWind << " " << atm2.perturbedNSWind << " " << eatm2.perturbedTemperature << endl;
      }
      else if (position1.height > 10.0) {
        if (position1.height > 17.0) {
          position2.height -= 0.5;
        }
        else {
          position2.height += 0.5;
        }
        position2.latitude -= 0.01;
        position2.longitude -= 0.01;
        position2.elapsedTime += 30.0;
        earth2.setPosition(position2);
        earth2.update();
        position2 = earth2.getPosition();
        AtmosphereState atm2 = earth2.getAtmosphereState();
        const EarthAtmosphereState& eatm2 = atm2.getPlanetSpecificMetrics<EarthAtmosphereState>();

        v2out_un << position2.height << " " << position2.latitude << " " << position2.longitude - 360.0 << " "
          << atm2.perturbedEWWind << " " << atm2.perturbedNSWind << " " << eatm2.perturbedTemperature << endl;

        correlator.updateCoefficients(position1, position2);
        correlator.correlate(atm1, atm2);

        v2out_corr << position2.height << " " << position2.latitude << " " << position2.longitude - 360.0 << " "
          << atm2.perturbedEWWind << " " << atm2.perturbedNSWind << " " << eatm2.perturbedTemperature << endl;
      }
    }

    v1out.close();
    v2out_un.close();
    v2out_corr.close();
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

