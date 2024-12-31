//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Uranus ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print the
//! UranusGRAM version string in the markdown output.
//! \ingroup UranusGRAM
class UranusProfilePrinter : public ProfilePrinter
{
public:
  UranusProfilePrinter();
  UranusProfilePrinter(const UranusProfilePrinter& orig) = default;
  virtual ~UranusProfilePrinter() = default;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;

};

} // namespace
