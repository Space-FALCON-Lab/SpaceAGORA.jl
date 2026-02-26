//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Titan ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print Column and List styles
//! in a format matching the legacy TitanGRAM output.
//! \ingroup TitanGRAM
class TitanProfilePrinter : public ProfilePrinter
{
public:
  TitanProfilePrinter();
  TitanProfilePrinter(const TitanProfilePrinter& orig) = default;
  virtual ~TitanProfilePrinter() = default;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;

  void printPlanetCSVHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetCSVStyle(const ProfileData& data) override;

  void printGramColumnHeader(const PerturbedAtmosphere& atmos) override;
  void printGramColumnStyle(const ProfileData& data) override;

  void printGramListHeader(const PerturbedAtmosphere& atmos) override;
  void printGramListStyle(const ProfileData& data) override;

  void printGramTPresHgtHeader(const PerturbedAtmosphere& atmos) override;
  void printGramTPresHgtStyle(const ProfileData& data) override;


};

} // namespace
