//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "GramExamples.h"
#include "VenusAtmosphere.h"
#include "VenusNamelistReader.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  int option = 0;
  cout << "====================" << endl;
  cout << "VenusGRAM examples" << endl;
  cout << "====================" << endl;
  cout << "1. Atmosphere examples" << endl;
  cout << "2. Trajectory example" << endl;
  cout << "3. Monte Carlo example" << endl;
  cout << "4. Namelist example" << endl;
  cout << "Select one: ";
  cin >> option;

  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  VenusInputParameters inputParameters;
  VenusNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  VenusAtmosphere venus;
  venus.setInputParameters(inputParameters);

  switch (option) {
  case 1:
    atmosphereExample(venus, "Venus");
    break;
  case 2:
    trajectoryExample(venus, "Venus");
    break;
  case 3:
    monteCarloExample(venus, "Venus");
    break;
  case 4:
  {
    string namelistFileName;
    cout << "Enter the namelist file name: ";
    cin >> namelistFileName;
    VenusNamelistReader reader;
    reader.getParameters(namelistFileName, inputParameters);
    venus.setInputParameters(inputParameters);
    namelistExample(venus, "Venus");
    break;
  }
  default:
    cout << "Sorry.  Try again." << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()/1000.0 << " seconds." << endl;
  return 0;
}

