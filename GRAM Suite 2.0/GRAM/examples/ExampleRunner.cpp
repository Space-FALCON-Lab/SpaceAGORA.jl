//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "GramExamples.h"
#include "NamelistReader.h"

#include "JupiterAtmosphere.h"
#include "EarthAtmosphere.h"
#include "MarsAtmosphere.h"
#include "NeptuneAtmosphere.h"
#include "TitanAtmosphere.h"
#include "UranusAtmosphere.h"
#include "VenusAtmosphere.h"
#include "VenusNamelistReader.h"
#include "EarthNamelistReader.h"
#include "MarsNamelistReader.h"
#include "JupiterNamelistReader.h"
#include "UranusNamelistReader.h"
#include "NeptuneNamelistReader.h"
#include "TitanNamelistReader.h"

const std::string earthDataPath = "../../Earth/data";
const std::string marsDataPath = "../../Mars/data";

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  string spicePath;
  NamelistReader reader;
  reader.tryGetSpicePath(spicePath);

  int bodyOption = 0;
  cout << "====================" << endl;
  cout << "GRAM examples" << endl;
  cout << "====================" << endl;

  while (true) {
    cout << "0. Quit Program" << endl;
    cout << "2. Venus" << endl;
    cout << "3. Earth" << endl;
    cout << "4. Mars" << endl;
    cout << "5. Jupiter" << endl;
    cout << "7. Uranus" << endl;
    cout << "8. Neptune" << endl;
    cout << "10. Titan" << endl;
    cout << "\nSelect one: ";
    cin >> bodyOption;

    PerturbedAtmosphere* body;
    string bodyName;
    switch (bodyOption) {
    case 2: {
      VenusInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      VenusAtmosphere* venus = new VenusAtmosphere();
      venus->setInputParameters(inputParameters);
      body = venus;
      bodyName = "Venus";
      break;
    }
    case 3: {
      EarthInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      inputParameters.dataPath = earthDataPath;
      EarthAtmosphere* earth = new EarthAtmosphere();
      earth->setInputParameters(inputParameters);
      body = earth;
      bodyName = "Earth";
      break;
    }
    case 4: {
      MarsInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      inputParameters.dataPath = marsDataPath;
      MarsAtmosphere* mars = new MarsAtmosphere();
      mars->setInputParameters(inputParameters);
      body = mars;
      bodyName = "Mars";
      break;
    }
    case 5: {
      JupiterInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      JupiterAtmosphere* jupiter = new JupiterAtmosphere();
      jupiter->setInputParameters(inputParameters);
      body = jupiter;
      bodyName = "Jupiter";
      break;
    }
    case 7: {
      UranusInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      UranusAtmosphere* uranus = new UranusAtmosphere();
      uranus->setInputParameters(inputParameters);
      body = uranus;
      bodyName = "Uranus";
      break;
    }
    case 8: {
      NeptuneInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      NeptuneAtmosphere* neptune = new NeptuneAtmosphere();
      neptune->setInputParameters(inputParameters);
      body = neptune;
      bodyName = "Neptune";
      break;
    }
    case 10: {
      TitanInputParameters inputParameters;
      inputParameters.spicePath = spicePath;
      TitanAtmosphere* titan = new TitanAtmosphere();
      titan->setInputParameters(inputParameters);
      body = titan;
      bodyName = "Titan";
      break;
    }
    default:
      cout << "Done." << endl;
      exit(0);
    }

    cout << endl;
    cout << "0. Unit testing" << endl;
    cout << "1. Atmosphere examples" << endl;
    cout << "2. Trajectory example" << endl;
    cout << "3. Monte Carlo example" << endl;
    cout << "4. Namelist example" << endl;
    cout << "5. Ephemeris example" << endl;
    cout << "6. Perturbation example" << endl;
    cout << "\nSelect one: ";
    int option = 0;
    cin >> option;

    typedef std::chrono::high_resolution_clock Clock;
    auto t1 = Clock::now();

    switch (option) {
    case 0: {
      // This kludge is needed to run ALL unit tests.
      // Earth and Mars must have a data path set first.
      EarthInputParameters eparams;
      eparams.spicePath = spicePath;
      eparams.dataPath = earthDataPath;
      EarthAtmosphere* earth = new EarthAtmosphere();
      earth->setInputParameters(eparams);

      MarsInputParameters mparams;
      mparams.spicePath = spicePath;
      mparams.dataPath = marsDataPath;
      MarsAtmosphere* mars = new MarsAtmosphere();
      mars->setInputParameters(mparams);

      unitTesting(argc, argv);
      break;
    }
    case 1:
      atmosphereExample(*body, bodyName);
      break;
    case 2:
      trajectoryExample(*body, bodyName);
      break;
    case 3:
      monteCarloExample(*body, bodyName);
      break;
    case 4:
    {
      string namelistFileName;
      cout << "Enter the namelist file name: ";
      cin >> namelistFileName;
      // Every body has its own custom namelist reader.
      switch (bodyOption) {
      case 2: {
        VenusNamelistReader reader;
        VenusInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 3: {
        EarthNamelistReader reader;
        EarthInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 4: {
        MarsNamelistReader reader;
        MarsInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 5: {
        JupiterNamelistReader reader;
        JupiterInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 7: {
        UranusNamelistReader reader;
        UranusInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 8: {
        NeptuneNamelistReader reader;
        NeptuneInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      case 10: {
        TitanNamelistReader reader;
        TitanInputParameters inputParameters;
        reader.getParameters(namelistFileName, inputParameters);
        body->setInputParameters(inputParameters);
        break;
      }
      default:
        break;
      }
      namelistExample(*body, bodyName);
      break;
    }
    case 5:
      ephemerisExample(*body, bodyName);
      break;
    case 6:
      perturbationExample(*body, bodyName);
      break;
    default:
      cout << "Sorry.  Try again." << endl;
    }

    auto t2 = Clock::now();
    cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0 << " seconds." << endl;
    cout << endl;
  }
  return 0;
}

