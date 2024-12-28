//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "GramExamples.h"
#include "UranusAtmosphere.h"
#include "UranusNamelistReader.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  int option = 0;
  cout << "====================" << endl;
  cout << "UranusGRAM examples" << endl;
  cout << "====================" << endl;
  cout << "1. Atmosphere examples" << endl;
  cout << "2. Trajectory example" << endl;
  cout << "3. Monte Carlo example" << endl;
  cout << "4. Namelist example" << endl;
  cout << "Select one: ";
  cin >> option;

  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  UranusInputParameters inputParameters;
  UranusNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  UranusAtmosphere uranus;
  uranus.setInputParameters(inputParameters);

  switch (option) {
  case 1:
    atmosphereExample(uranus, "Uranus");
    break;
  case 2:
    trajectoryExample(uranus, "Uranus");
    break;
  case 3:
    monteCarloExample(uranus, "Uranus");
    break;
  case 4:
  {
    string namelistFileName;
    cout << "Enter the namelist file name: ";
    cin >> namelistFileName;
    UranusNamelistReader reader;
    reader.getParameters(namelistFileName, inputParameters);
    uranus.setInputParameters(inputParameters);
    namelistExample(uranus, "Uranus");
    break;
  }
  default:
    cout << "Sorry.  Try again." << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()/1000.0 << " seconds." << endl;
  return 0;
}

