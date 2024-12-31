//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include "JupiterProfilePrinter.h"
#include "JupiterAtmosphere.h"

using namespace std;

namespace GRAM {

JupiterProfilePrinter::JupiterProfilePrinter()
{
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void JupiterProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const JupiterAtmosphere& jupiter = static_cast<const JupiterAtmosphere&>(atmos);

  gramMDFile << "# " << jupiter.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

} // namespace
