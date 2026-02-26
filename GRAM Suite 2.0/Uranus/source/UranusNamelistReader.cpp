//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include "UranusNamelistReader.h"

using namespace std;

namespace GRAM {

UranusNamelistReader::UranusNamelistReader()
{
}

void UranusNamelistReader::getParameters(const std::string& fileName, UranusInputParameters& inputParameters)
{
  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters.
  NamelistReader::getParameters(inputParameters);

  checkForErrors();
}

} // namespace
