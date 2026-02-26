//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <set>
#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

//! \brief This is a utility class for loading Spice data files into memory.
//!
//! This utility class wraps the NAIF Spice furnsh function for loading data into memory.
//! Consecutive calls of the same file to furnsh will result in memory bloat.
//! So the SpiceLoader will only load unique files into memory.  The methods in 
//! this class are static and may be called directly.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class SpiceLoader
{
public:
  SpiceLoader();
  SpiceLoader(const SpiceLoader& orig) = delete;
  virtual ~SpiceLoader() = default;

  static void load(GRAM_BODY body);
  static void loadFile(const std::string& path);
  static void setSpiceLsk(const std::string& lsk);
  static void setSpicePck(const std::string& pck);
  static void setSpiceKernel(GRAM_BODY body, const std::string& bsp);
  static void setSpiceDataPath(const std::string& path);
  static const std::string& getSpiceDataPath() { return spiceDataPath; }
  static void setInputParameters(const InputParameters& params);

private:
  static std::set<std::string> spiceFiles;  //!< The set of Spice files that have already been loaded into memory.
  static std::string spiceDataPath;         //!< The path to the root Spice data folder.

  static std::string spiceLsk;
  static std::string spicePck;
  static std::string spiceVenus;
  static std::string spiceEarth;
  static std::string spiceMars;
  static std::string spiceJupiter;
  static std::string spiceSaturn;
  static std::string spiceUranus;
  static std::string spiceNeptune;
  static std::string spiceTitan;

};

} // namespace

