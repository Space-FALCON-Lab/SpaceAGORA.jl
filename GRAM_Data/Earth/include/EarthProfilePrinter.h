//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"

namespace GRAM {

//! \brief This is a customized Earth ProfilePrinter.
//!
//! This ProfilePrinter has been customized to print styles
//! in a format matching the legacy EarthGRAM output.
//! \ingroup EarthGRAM
class EarthProfilePrinter : public ProfilePrinter
{
public:
  EarthProfilePrinter();
  EarthProfilePrinter(const EarthProfilePrinter& orig) = default;
  virtual ~EarthProfilePrinter() = default;

  void setInputParameters(const InputParameters& params) override;

  void printSectionHeader(const PerturbedAtmosphere& atmos) override;

  void printPlanetCSVHeader(const PerturbedAtmosphere& atmos) override;
  void printPlanetCSVStyle(const ProfileData& data) override;

  void printGramListMDHeader(const PerturbedAtmosphere& atmos) override;
  void printGramListMDStyle(const ProfileData& data) override;

  void openOutput() override;
  void printFileHeader(const PerturbedAtmosphere& atmos) override;
  void printData(const std::vector<ProfileData>& profile) override;
  void closeOutput() override;

  void printGramColumnHeader(const PerturbedAtmosphere& atmos) override;
  void printGramColumnStyle(const ProfileData& data) override;

  void printGramListHeader(const PerturbedAtmosphere& atmos) override;
  void printGramListStyle(const ProfileData& data) override;

  void printSpeciesHeader(const PerturbedAtmosphere& atmos);
  void printSpeciesStyle(const ProfileData& data);

  void printBoundaryLayerHeader(const PerturbedAtmosphere& atmos);
  void printBoundaryLayerStyle(const ProfileData& data);

  void printCorrHeader(const PerturbedAtmosphere& atmos);
  void printCorrStyle(const ProfileData& data);

  void printTrajHeader(const PerturbedAtmosphere& atmos);
  void printTrajStyle(const ProfileData& data);

  static constexpr size_t EARTH_SPECIES_STYLE = 1 << 16;
  static constexpr size_t EARTH_BLTEST_STYLE = 1 << 17;
  static constexpr size_t EARTH_CORR_STYLE = 1 << 18;
  static constexpr size_t EARTH_CORR1_STYLE = 1 << 19;
  static constexpr size_t EARTH_TRAJ_STYLE = 1 << 20;

protected:
  std::ofstream speciesFile;              //!< Output stream.
  std::ofstream bltestFile;               //!< Output stream.
  std::ofstream corrFile;                 //!< Output stream.
  std::ofstream trajFile;                 //!< Output stream.

  std::string speciesFileName = "species";
  std::string bltestFileName = "BLTest";
  std::string corrFileName = "Profiles";
  std::string trajFileName = "CorrTraj";
};

} // namespace
