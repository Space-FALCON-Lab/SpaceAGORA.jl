//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Mars ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print Column and List styles
//! in a format matching the legacy MarsGRAM output.
//! \ingroup MarsGRAM
class MarsProfilePrinter : public ProfilePrinter
{
public:
  MarsProfilePrinter();
  MarsProfilePrinter(const MarsProfilePrinter& orig) = default;
  virtual ~MarsProfilePrinter() = default;

  void printPlanetCSVHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetCSVStyle(const ProfileData& data) override;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetListMDHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetListMDStyle(const ProfileData& data) override;
  void printHeightMDStyle(const ProfileData& data) override;

  void openOutput() override;
  void printFileHeader(const PerturbedAtmosphere& atmos) override;
  void printData(const std::vector<ProfileData>& profile) override;
  void closeOutput() override;

  void printGramColumnHeader(const PerturbedAtmosphere& atmos) override;
  void printGramColumnStyle(const ProfileData& data) override;

  void printGramListHeader(const PerturbedAtmosphere& atmos) override;
  void printGramListStyle(const ProfileData& data) override;

  void printGramTPresHgtHeader(const PerturbedAtmosphere& atmos) override;
  void printGramTPresHgtStyle(const ProfileData& data) override;

  void printGramDensityHeader(const PerturbedAtmosphere& atmos) override;
  void printGramDensityStyle(const ProfileData& data) override;

  void printGramWindsHeader(const PerturbedAtmosphere& atmos) override;
  void printGramWindsStyle(const ProfileData& data) override;

  virtual void printGramPerturbHeader(const PerturbedAtmosphere& atmos);
  virtual void printGramPerturbStyle(const ProfileData& data);

  void printDayDataHeader(const PerturbedAtmosphere& atmos);
  void printDayDataStyle(const ProfileData& data);

  void printMarsRadHeader(const PerturbedAtmosphere& atmos);
  void printMarsRadStyle(const ProfileData& data);

  void printThrmDataHeader(const PerturbedAtmosphere& atmos);
  void printThrmDataStyle(const ProfileData& data);

  int ibougher = 0;  //!< Flag for the selected offset model.
  int mapyear = 0;   //!< Flag for the selected map year.

  static constexpr size_t MARS_DAY_STYLE = 1 << 16;   //!< Custom Mars print style.
  static constexpr size_t MARS_RAD_STYLE = 1 << 17;   //!< Custom Mars print style.
  static constexpr size_t MARS_THRM_STYLE = 1 << 18;  //!< Custom Mars print style.

protected:
  std::ofstream dayDataFile;   //!< Output stream.
  std::ofstream thrmDataFile;  //!< Output stream.
  std::ofstream marsRadFile;   //!< Output stream.

  bool isMolaHeights = true;   //!< Flag for height output.
};

} // namespace
