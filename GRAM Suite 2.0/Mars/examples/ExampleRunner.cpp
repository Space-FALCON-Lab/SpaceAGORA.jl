//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "GramExamples.h"
#include "MarsAtmosphere.h"
#include "MarsNamelistReader.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  int option = 0;
  cout << "====================" << endl;
  cout << "MarsGRAM examples" << endl;
  cout << "====================" << endl;
  cout << "1. Atmosphere examples" << endl;
  cout << "2. Trajectory example" << endl;
  cout << "3. Monte Carlo example" << endl;
  cout << "4. Namelist example" << endl;
  cout << "Select one: ";
  cin >> option;

  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  try {
    MarsInputParameters inputParameters;
    inputParameters.dataPath = "./";
    inputParameters.mapYear = 0;

    MarsNamelistReader reader;
    reader.tryGetSpicePath(inputParameters);

    MarsAtmosphere mars;
    mars.setInputParameters(inputParameters);

    switch (option) {
    case 1:
      atmosphereExample(mars, "Mars");
      break;
    case 2:
      trajectoryExample(mars, "Mars");
      break;
    case 3:
      monteCarloExample(mars, "Mars");
      break;
    case 4:
    {
      string namelistFileName;
      cout << "Enter the namelist file name: ";
      cin >> namelistFileName;
      MarsNamelistReader reader;
      reader.getParameters(namelistFileName, inputParameters);
      mars.setInputParameters(inputParameters);
      namelistExample(mars, "Mars");
      break;
    }
    default:
      cout << "Sorry.  Try again." << endl;
    }
  }
  catch (const string& msg) {
    cerr << msg << endl;
  }
  catch (const std::runtime_error& err) {
    cerr << "An unanticipated error occurred." << endl;
    cerr << err.what() << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()/1000.0 << " seconds." << endl;
  return 0;
}

