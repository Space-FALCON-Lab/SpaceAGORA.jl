//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include "VenusNamelistReader.h"

using namespace std;

namespace GRAM {

//! \copydoc NamelistReader::NamelistReader()
VenusNamelistReader::VenusNamelistReader()
{
}

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param fileName Full or relative path to a namelist file.
//! \param inputParameters An VenusInputParameters object.
void VenusNamelistReader::getParameters(const std::string& fileName, VenusInputParameters& inputParameters)
{
  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters
  NamelistReader::getParameters(inputParameters);

  // Get the Venus input parameters
  // At this point there are no Venus specific parameters

  checkForErrors();
}

} // namespace
