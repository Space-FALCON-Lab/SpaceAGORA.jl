//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#include "NeptuneNamelistReader.h"
#include <algorithm>
#include <fstream>

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
NeptuneNamelistReader::NeptuneNamelistReader() {}

//! \fn  NeptuneNamelistReader::NeptuneNamelistReader(const NeptuneNamelistReader& orig)
//! \brief Copying this object is discouraged.

//! \fn  NeptuneNamelistReader::~NeptuneNamelistReader()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param fileName Full or relative path to a namelist file.
//! \param inputParameters An NeptuneInputParameters object.
void NeptuneNamelistReader::getParameters(const std::string &fileName,
                                          NeptuneInputParameters &inputParameters)
{
  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters
  NamelistReader::getParameters(inputParameters);

  // Get the Neptune input parameters
  getItem("FMINMAX", inputParameters.minMaxFactor);
  getItem("MINMAXFACTOR", inputParameters.minMaxFactor);

  getItem("IFMM", inputParameters.computeMinMaxFactor);
  getItem("COMPUTEMINMAXFACTOR", inputParameters.computeMinMaxFactor);

  getItem("FMOLNITRO", inputParameters.dinitrogenMoleFraction);
  getItem("DINITROGENMOLEFRACTION", inputParameters.dinitrogenMoleFraction);

  checkForErrors();
}

} // namespace GRAM
