//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include <fstream>
#include "gram.h"
#include "Profile.h"

namespace GRAM {

//! \brief The ProfilePrinter class is an output utility class for Profile data.
//!
//! This class will write a list of Profile data to a file in various formats.  A new print format
//! can be defined by:
//!    \li Defining a style constant, such as GRAM_COLUMN_STYLE.
//!    \li Defining a header function, such as printGramColumnHeader().
//!    \li Defining a Profile print function, such as printGramColumnStyle().
//!    \li Inserting the style into openOutput(), printFileHeader(), printSectionHeader(), printData(), and closeOutput().
//!
//! This class is easily subclassed for customizing output.  Existing print styles my be overridden, or a new
//! print style can be created.  For custom print styles, it is recommended that the style constant bit mask
//! use a shift value of 16 or greater.
//! \ingroup CommonGRAM
class ProfilePrinter
{
public:
  ProfilePrinter();
  ProfilePrinter(const ProfilePrinter& orig);
  virtual ~ProfilePrinter() = default;

  virtual void setInputParameters(const InputParameters& params);

  virtual void setEastLongitudePositive(bool flag = true) { eastLongitudePositive = flag; }
  virtual void setWestLongitudePositive() { eastLongitudePositive = false; }
  virtual void setStyle(size_t styleFlag) { outputStyle = styleFlag; }
  virtual void setDensityPrintScale(DensityPrintScale scale) { densityPrintScale = scale; }
  virtual void setExtraPrecision(size_t precision) { extra = precision; }

  virtual void print(const PerturbedAtmosphere& atmos, const std::vector<ProfileData>& profile);

  virtual void openOutput();
  virtual void printFileHeader(const PerturbedAtmosphere& atmos);
  virtual void printSectionHeader(const PerturbedAtmosphere& atmos);
  virtual void printData(const std::vector<ProfileData>& profile);
  virtual void closeOutput();

  void setListFileName(const std::string& fileName) { if (!fileName.empty()) listFileName = fileName; }
  void setColumnFileName(const std::string& fileName) { if (!fileName.empty()) columnFileName = fileName; }
  void setFileNamePrefix(const std::string& prefix) { fileNamePrefix = prefix; }

  virtual void printGramCSVHeader(const PerturbedAtmosphere& atmos);
  virtual void printGramCSVStyle(const ProfileData& data);
  virtual void printPlanetCSVHeader(const PerturbedAtmosphere& atmos) {}
  virtual void printPlanetCSVStyle(const ProfileData& data) {}

  virtual void printGramListMDHeader(const PerturbedAtmosphere& atmos);
  virtual void printGramListMDStyle(const ProfileData& data);
  virtual void printHeightMDStyle(const ProfileData& data);
  virtual void printPlanetListMDHeader(const PerturbedAtmosphere& atmos) {}
  virtual void printPlanetListMDStyle(const ProfileData& data) {}

  virtual void printGramColumnHeader(const PerturbedAtmosphere& atmos);
  virtual void printGramColumnStyle(const ProfileData& data);

  virtual void printGramListHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramListStyle(const ProfileData& data);

	virtual void printGramEphemerisHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramEphemerisStyle(const ProfileData& data);

	virtual void printGramDensityHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramDensityStyle(const ProfileData& data);

	virtual void printGramPerturbHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramPerturbStyle(const ProfileData& data);
  
	virtual void printGramWindsHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramWindsStyle(const ProfileData& data);
  
	virtual void printGramTPresHgtHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramTPresHgtStyle(const ProfileData& data);
  
	virtual void printGramSoundHeader(const PerturbedAtmosphere& atmos);
	virtual void printGramSoundStyle(const ProfileData& data);

  const std::string& getOutputFileNames() const { return outputFiles; }

  static constexpr size_t GRAM_COLUMN_STYLE    = 1 << 0;
  static constexpr size_t GRAM_LIST_STYLE      = 1 << 1;
  static constexpr size_t GRAM_EPHEMERIS_STYLE = 1 << 2;
  static constexpr size_t GRAM_DENSITY_STYLE   = 1 << 3;
  static constexpr size_t GRAM_PERTURB_STYLE   = 1 << 4;
  static constexpr size_t GRAM_WINDS_STYLE     = 1 << 5;
  static constexpr size_t GRAM_TPRESHGT_STYLE  = 1 << 6;
  static constexpr size_t GRAM_SOUND_STYLE     = 1 << 7;
  static constexpr size_t GRAM_CSV_STYLE       = 1 << 8;
  static constexpr size_t GRAM_MD_STYLE        = 1 << 9;
  static constexpr size_t GRAM_MONTE_CARLO_STYLE =
    GRAM_COLUMN_STYLE | GRAM_LIST_STYLE | GRAM_DENSITY_STYLE | GRAM_PERTURB_STYLE | GRAM_WINDS_STYLE | GRAM_TPRESHGT_STYLE;
  static constexpr size_t GRAM_ALL_STYLES = ~0;

  void printDates(GramTime gramTime[3], greal lonSun[3], greal tlst[3]);

protected:
  void scaleDensityForOutput(AtmosphereState& atmos);
  const char* eastWestName(bool flip = false) { return (eastLongitudePositive ^ flip) ? "East" : "West"; }
  char eastWestChar(bool flip = false) { return (eastLongitudePositive ^ flip) ? 'E' : 'W'; }
  void printGasHeaders(std::ostream& stream, const PerturbedAtmosphere& atmos, const char units[4]);
  void printGasHeaders(std::ostream& stream, const PerturbedAtmosphere& atmos);
  void openFile(std::ofstream& outputStream, const std::string& fileName);
  void printCSVGasData(const ConstituentGas& gas);
  void printListGasData(const ConstituentGas& gas, const char* name);
  void appendExtension(std::string& fileName, const std::string& ext);

  static const size_t bufferSize = 1200;
  bool hasVerticalWinds = false;

  const char comma = ',';
  const char newline = '\n';
  size_t extra = 0;
  size_t mdCount = 0;
  std::string outputFiles;

  greal massTotal = 0.0;
  greal moleTotal = 0.0;

  std::string densUnits = "kg/m^3";                //!< The units string for density.
  DensityPrintScale densityPrintScale = STANDARD;  //!< The selected density scale.
  bool eastLongitudePositive = true;               //!< The selected convention for longitude output.
  std::string listFileName = "LIST";               //!< The file name for LIST style output
  std::string columnFileName = "OUTPUT";           //!< The file name for COLUMN style output
  std::string fileNamePrefix = "";                 //!< A prefix to go before every output file name.
  size_t outputStyle = 0;                          //!< Bitfield mask for selected print styles.
  std::ofstream gramCSVFile;                       //!< Output stream.
  std::ofstream gramMDFile;                        //!< Output stream.
  std::ofstream gramColumnFile;                    //!< Output stream.
  std::ofstream gramListFile;                      //!< Output stream.
  std::ofstream gramEphemerisFile;                 //!< Output stream.
  std::ofstream gramDensityFile;                   //!< Output stream.
  std::ofstream gramPerturbFile;                   //!< Output stream.
  std::ofstream gramWindsFile;                     //!< Output stream.
  std::ofstream gramTPresHgtFile;                  //!< Output stream.
  std::ofstream gramSoundFile;                     //!< Output stream.
};

} // namespace
