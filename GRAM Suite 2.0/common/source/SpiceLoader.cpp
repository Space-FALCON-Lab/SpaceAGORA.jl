//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "SpiceUsr.h"
#include "SpiceLoader.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

std::set<std::string> SpiceLoader::spiceFiles;
std::string SpiceLoader::spiceDataPath = DEFAULT_SPICE_PATH;
std::string SpiceLoader::spiceLsk     = "/lsk/naif0012.tls";
std::string SpiceLoader::spicePck     = "/pck/pck00011.tpc";
std::string SpiceLoader::spiceVenus   = "/spk/planets/de440_GRAM.bsp";
std::string SpiceLoader::spiceEarth   = "/spk/planets/de440_GRAM.bsp";
std::string SpiceLoader::spiceMars    = "/spk/satellites/mar097_GRAM.bsp";
std::string SpiceLoader::spiceJupiter = "/spk/satellites/jup365_GRAM.bsp";
std::string SpiceLoader::spiceSaturn  = "/spk/satellites/sat441_GRAM.bsp";
std::string SpiceLoader::spiceUranus  = "/spk/satellites/ura116_GRAM.bsp";
std::string SpiceLoader::spiceNeptune = "/spk/satellites/nep101_GRAM.bsp";
std::string SpiceLoader::spiceTitan   = "/spk/satellites/sat441_GRAM.bsp";

//! \copydoc Atmosphere::Atmosphere()
SpiceLoader::SpiceLoader()
{
  // SPICE call: set the error action to REPORT
  SpiceInt len = 0;
  SpiceChar action[10] = "REPORT";
  SpiceChar device[10] = "NULL";
//  SpiceChar device[10] = "SCREEN";
  erract_c("SET", len, action);
  errdev_c("SET", len, device);
}

//! \fn SpiceLoader::SpiceLoader(const SpiceLoader& orig)
//! \brief Copying this object is discouraged.

//! \fn SpiceLoader::~SpiceLoader()
//! \copydoc Atmosphere::~Atmosphere()

void SpiceLoader::setSpiceLsk(const std::string& lsk) {
  if (!lsk.empty()) {
    spiceLsk = lsk;
  }
}

void SpiceLoader::setSpicePck(const std::string& pck) {
  if (!pck.empty()) {
    spicePck = pck;
  }
}

void SpiceLoader::setSpiceKernel(GRAM_BODY body, const std::string& bsp) {
  if (!bsp.empty()) {
    switch (body) {
    case VENUS:
      spiceVenus = bsp;
      break;
    case EARTH:
      spiceEarth = bsp;
      break;
    case MARS:
      spiceMars = bsp;
      break;
    case JUPITER:
      spiceJupiter = bsp;
      break;
    case URANUS:
      spiceUranus = bsp;
      break;
    case NEPTUNE:
      spiceNeptune = bsp;
      break;
    case SATURN:
      spiceSaturn = bsp;
      break;
    case TITAN:
      spiceTitan = bsp;
      break;
    default:
      break;
    }
  }
}
 
//! \brief Set the path to the Spice data files.
//!
//! Spice data is distributed in a folder containing lsk, pck, and spk subfolders.
//! The spice data path should be set to point to this root Spice folder.
//! This static method may be called directly with SpiceLoader::setSpiceDataPath("mypath").
//! \param path The path to the root folder of the spice data.
void SpiceLoader::setSpiceDataPath(const std::string& path)
{ 
  spiceDataPath = path; 
  loadFile(spiceLsk);
}

void SpiceLoader::setInputParameters(const InputParameters& params) {

  if (!params.spiceLsk.empty()) {
    spiceLsk = params.spiceLsk;
  }
  if (!params.spicePck.empty()) {
    spicePck = params.spicePck;
  }
  if (!params.spiceVenus.empty()) {
    spiceVenus = params.spiceVenus;
  }
  if (!params.spiceEarth.empty()) {
    spiceEarth = params.spiceEarth;
  }
  if (!params.spiceMars.empty()) {
    spiceMars = params.spiceMars;
  }
  if (!params.spiceJupiter.empty()) {
    spiceJupiter = params.spiceJupiter;
  }
  if (!params.spiceSaturn.empty()) {
    spiceSaturn = params.spiceSaturn;
  }
  if (!params.spiceUranus.empty()) {
    spiceUranus = params.spiceUranus;
  }
  if (!params.spiceNeptune.empty()) {
    spiceNeptune = params.spiceNeptune;
  }
  if (!params.spiceTitan.empty()) {
    spiceTitan = params.spiceTitan;
  }
  setSpiceDataPath(params.spicePath);
}

void SpiceLoader::load(GRAM_BODY body)
{
  // The LSK file is a leap-second database.
  loadFile(spiceLsk);

  // The PCK file defines body size, shape, orientation
  loadFile(spicePck);

  switch (body) {
  case VENUS:
    loadFile(spiceVenus);
    break;
  case EARTH:
    loadFile(spiceEarth);
    break;
  case MARS:
    loadFile(spiceMars);
    break;
  case JUPITER:
    loadFile(spiceJupiter);
    break;
  case URANUS:
    loadFile(spiceUranus);
    break;
  case NEPTUNE:
    loadFile(spiceNeptune);
    break;
  case SATURN:
    loadFile(spiceSaturn);
    break;
  case TITAN:
    loadFile(spiceTitan);
    break;
  default:
    break;
  }
}

//! \brief Loads spice data into memory ensuring uniqueness.
//!
//! After checking for uniqueness of the (relative path) file name, this method
//! loads the Spice data into memory via the Spice furnsh function.
//! \param path A path to the spice file relative to the root Spice folder.
void SpiceLoader::loadFile(const std::string & path)
{
  auto result = spiceFiles.insert(path);
  if (result.second) {
    furnsh_c((spiceDataPath + path).c_str());
    if (failed_c()) {
      throw(SPICE_ERROR_MESSAGE
      + "Unable to load a SPICE file. Check the SPICE data path."
      + "\n       Attempting to load: " + spiceDataPath + path
      + "\n       This is an unrecoverable error.");
    }
  }
}

//! \fn SpiceLoader::getSpiceDataPath()
//! \brief Get the current Spice data path.
//!
//! \returns The current Spice data path.

} // namespace