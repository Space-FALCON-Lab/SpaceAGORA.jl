//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include "UranusProfilePrinter.h"
#include "UranusAtmosphere.h"

using namespace std;

namespace GRAM {

UranusProfilePrinter::UranusProfilePrinter()
{
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void UranusProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const UranusAtmosphere& uranus = static_cast<const UranusAtmosphere&>(atmos);

  gramMDFile << "# " << uranus.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

} // namespace
