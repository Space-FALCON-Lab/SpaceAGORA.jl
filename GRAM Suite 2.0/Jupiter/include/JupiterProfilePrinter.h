//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Jupiter ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print the
//! JupiterGRAM version string in the markdown output.
//! \ingroup JupiterGRAM
class JupiterProfilePrinter : public ProfilePrinter
{
public:
  JupiterProfilePrinter();
  JupiterProfilePrinter(const JupiterProfilePrinter& orig) = default;
  virtual ~JupiterProfilePrinter() = default;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;

};

} // namespace
