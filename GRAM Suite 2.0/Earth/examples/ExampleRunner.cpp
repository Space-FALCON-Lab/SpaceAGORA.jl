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
#include <string>
#include "GramExamples.h"
#include "EarthAtmosphere.h"
#include "EarthNamelistReader.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  int option = 0;
  cout << "====================" << endl;
  cout << "EarthGRAM examples" << endl;
  cout << "====================" << endl;
  cout << "1. Atmosphere examples" << endl;
  cout << "2. Trajectory example" << endl;
  cout << "3. Monte Carlo example" << endl;
  cout << "4. Namelist example" << endl;
  cout << "Select one: ";
  cin >> option;

  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  EarthInputParameters inputParameters;
  EarthNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  EarthAtmosphere earth;
  earth.setInputParameters(inputParameters);
  auto t3 = Clock::now();
  double init = 0.0;

  switch (option) {
  case 1:
    atmosphereExample(earth, "Earth");
    break;
  case 2:
    trajectoryExample(earth, "Earth");
    break;
  case 3:
    monteCarloExample(earth, "Earth");
    break;
  case 4:
  {
    string namelistFileName;
    cout << "Enter the namelist file name: ";
    cin >> namelistFileName;
    init = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t1).count() / 1000.0;
    t1 = Clock::now();
    reader.getParameters(namelistFileName, inputParameters);
    earth.setInputParameters(inputParameters);
    namelistExample(earth, "Earth");
    break;
  }
  default:
    cout << "Sorry.  Try again." << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << init << " seconds." << endl;
  cout << "Elapsed time: " << init + std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0 << " seconds." << endl;
  return 0;
}

