//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include "TitanNamelistReader.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanNamelistReader::TitanNamelistReader()
{
}

//! \fn  TitanNamelistReader::TitanNamelistReader(const TitanNamelistReader& orig)
//! \brief Copying this object is discouraged.

//! \fn  TitanNamelistReader::~TitanNamelistReader()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param fileName Full or relative path to a namelist file.
//! \param inputParameters An TitanInputParameters object.
void TitanNamelistReader::getParameters(const std::string& fileName, TitanInputParameters& inputParameters)
{
  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters
  NamelistReader::getParameters(inputParameters);

  // Get the Titan input parameters
  getItem("FMINMAX", inputParameters.minMaxFactor);
  getItem("MINMAXFACTOR", inputParameters.minMaxFactor);

  getItem("COMPUTEMINMAXFACTOR", inputParameters.computeMinMaxFactor);

  getItem("FMOLMETH", inputParameters.userMethaneMoleFraction);
  getItem("METHANEMOLEFRACTION", inputParameters.userMethaneMoleFraction);

  // Legacy form
  int ifmm;
  if (getItem("IFMM", ifmm)) {
    if (ifmm == 2) {
      inputParameters.modelType = GCM95;
    }
    else {
      inputParameters.modelType = Yelle97;
      inputParameters.computeMinMaxFactor = (ifmm == 1);
    }
  }

  // New form
  int type;
  if (getItem("MODELTYPE", type)) {
    switch (type) {
    case 1:
      inputParameters.modelType = Yelle97;
      break;
    case 2:
      inputParameters.modelType = GCM95;
      break;
    default:
      inputParameters.modelType = Yelle97;
      break;
    }
  }

  checkForErrors();
}

} // namespace
