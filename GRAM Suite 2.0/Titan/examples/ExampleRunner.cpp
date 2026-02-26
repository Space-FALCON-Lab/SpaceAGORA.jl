//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include "GramExamples.h"
#include "TitanAtmosphere.h"
#include "TitanNamelistReader.h"

using namespace std;
using namespace GRAM;

int main(int argc, char** argv)
{
  int option = 0;
  cout << "====================" << endl;
  cout << "TitanGRAM examples" << endl;
  cout << "====================" << endl;
  cout << "1. Atmosphere examples" << endl;
  cout << "2. Trajectory example" << endl;
  cout << "3. Monte Carlo example" << endl;
  cout << "4. Namelist example" << endl;
  cout << "5. Ephemeris example" << endl;
  cout << "6. Perturbation example" << endl;
  cout << "Select one: ";
  cin >> option;

  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  TitanInputParameters inputParameters;
  TitanNamelistReader reader;
  reader.tryGetSpicePath(inputParameters);

  TitanAtmosphere titan;
  titan.setInputParameters(inputParameters);
  titan.setModelType(Yelle97);

  switch (option) {
  case 1:
    atmosphereExample(titan, "Titan Yelle");
    titan.setModelType(GCM95);
    atmosphereExample(titan, "Titan GCM");
    break;
  case 2:
    trajectoryExample(titan, "Titan Yelle");
    break;
  case 3:
    monteCarloExample(titan, "Titan Yelle");
    break;
  case 4:
  {
    string namelistFileName;
    cout << "Enter the namelist file name: ";
    cin >> namelistFileName;
    reader.getParameters(namelistFileName, inputParameters);
    titan.setInputParameters(inputParameters);
    namelistExample(titan, "Titan");
    break;
  }
  case 5:
    ephemerisExample(titan, "Titan");
    break;
  case 6:
    perturbationExample(titan, "Titan");
    break;
  default:
    cout << "Sorry.  Try again." << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()/1000.0 << " seconds." << endl;
  return 0;
}

