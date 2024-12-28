//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Neptune ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print Column and List styles
//! in a format matching the legacy NetuneGRAM output.
//! \ingroup NeptuneGRAM
class NeptuneProfilePrinter : public ProfilePrinter
{
public:
  NeptuneProfilePrinter();
  NeptuneProfilePrinter(const NeptuneProfilePrinter& orig) = default;
  virtual ~NeptuneProfilePrinter() = default;

  void printPlanetCSVHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetCSVStyle(const ProfileData& data) override;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetListMDStyle(const ProfileData& data) override;

  void printGramColumnHeader(const PerturbedAtmosphere& atmos) override;
  void printGramColumnStyle(const ProfileData& data) override;

  void printGramListHeader(const PerturbedAtmosphere& atmos) override;
  void printGramListStyle(const ProfileData& data) override;

};

} // namespace
