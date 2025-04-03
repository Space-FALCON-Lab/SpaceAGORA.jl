//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include "MarsNamelistReader.h"
#include "MarsCommon.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MarsNamelistReader::MarsNamelistReader()
{
}

//! \fn  MarsNamelistReader::MarsNamelistReader(const MarsNamelistReader& orig)
//! \brief Copying this object is discouraged.

//! \fn  MarsNamelistReader::~MarsNamelistReader()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param fileName Full or relative path to a namelist file.
//! \param inputParameters A MarsInputParameters object.
void MarsNamelistReader::getParameters(const std::string& fileName, MarsInputParameters& inputParameters)
{
  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters
  NamelistReader::getParameters(inputParameters);

  // Get the Mars input parameters
  getItem("MOLAHGTS", inputParameters.isMolaHeights);
  getItem("ISMOLAHEIGHTS", inputParameters.isMolaHeights);

  if (inputParameters.isMolaHeights == true) {
    inputParameters.isPlanetoCentric = true;
  }

  getItem("MAPYEAR", inputParameters.mapYear);

  getItem("IDAYDATA", inputParameters.computeMinMax);
  getItem("COMPUTEMINMAX", inputParameters.computeMinMax);

  // Bougher model parameters
  int offsetModel = inputParameters.offsetModel;
  getItem("IBOUGHER", offsetModel);
  getItem("OFFSETMODEL", offsetModel);
  switch (offsetModel) {
  case 0:
    inputParameters.offsetModel = MARS_CONSTANT;
    break;
  case 1:
    inputParameters.offsetModel = MARS_SEASONAL;
    break;
  case 2:
    inputParameters.offsetModel = MARS_GLOBAL_MEAN;
    break;
  case 3:
    inputParameters.offsetModel = MARS_DAILY_AVERAGE;
    break;
  case 4:
    inputParameters.offsetModel = MARS_CURRENT;
    break;
  default:
    inputParameters.offsetModel = MARS_GLOBAL_MEAN;
    break;
  }

  getItem("ZOFFSET", inputParameters.constantHeightOffset);
  getItem("CONSTANTHEIGHTOFFSET", inputParameters.constantHeightOffset);

  // Stewart Model Parameters
  getItem("DELTATEX", inputParameters.exosphericTemperatureOffset);
  getItem("EXOSPHERICTEMPERATUREOFFSET", inputParameters.exosphericTemperatureOffset);

  getItem("STDL", inputParameters.exosphericTemperatureFactor);
  getItem("EXOSPHERICTEMPERATUREFACTOR", inputParameters.exosphericTemperatureFactor);

  getItem("F107", inputParameters.F107);

  greal inMeters = 0.0;
  getItem("HGTASFCM", inMeters);
  getItem("HEIGHTABOVESURFACE", inMeters);
  inputParameters.heightAboveSurface = inMeters * 0.001;

  greal rwScale;
  if (getItem("RWSCALE", rwScale)) {
    rwScale = max(0.0, rwScale);
    inputParameters.ewWindPerturbationScale = rwScale;
    inputParameters.nsWindPerturbationScale = rwScale;
  }

  getItem("WLSCALE", inputParameters.perturbationWaveLengthScale);
  getItem("PERTURBATIONWAVELENGTHSCALE", inputParameters.perturbationWaveLengthScale);

  // Slope Winds Parameters
  getItem("WMSCALE", inputParameters.meanWindsScale);
  getItem("MEANWINDSSCALE", inputParameters.meanWindsScale);

  getItem("BLWINFAC", inputParameters.boundaryLayerWindsScale);
  getItem("BOUNDARYLAYERWINDSSCALE", inputParameters.boundaryLayerWindsScale);

  // Dust parameters
  getItem("DUSTTAU", inputParameters.mgcmConstantDustLevel);
  getItem("MGCMCONSTANTDUSTLEVEL", inputParameters.mgcmConstantDustLevel);

  getItem("DUSTMAX", inputParameters.mgcmMaxDustLevel);
  getItem("MGCMMAXDUSTLEVEL", inputParameters.mgcmMaxDustLevel);

  getItem("DUSTMIN", inputParameters.mgcmMinDustLevel);
  getItem("MGCMMINDUSTLEVEL", inputParameters.mgcmMinDustLevel);

  getItem("DUSTNU", inputParameters.dustNu);

  getItem("DUSTDIAM", inputParameters.dustDiameter);
  getItem("DUSTDIAMETER", inputParameters.dustDiameter);

  getItem("DUSTDENS", inputParameters.dustDensity);
  getItem("DUSTDENSITY", inputParameters.dustDensity);

  // Dust Storm Parameters
  getItem("ALS0", inputParameters.stormLongitudeSun);
  getItem("STORMLONGITUDESUN", inputParameters.stormLongitudeSun);

  getItem("ALSDUR", inputParameters.stormDuration);
  getItem("STORMDURATION", inputParameters.stormDuration);

  getItem("INTENS", inputParameters.stormIntensity);
  getItem("STORMINTENSITY", inputParameters.stormIntensity);

  getItem("RADMAX", inputParameters.stormMaxRadius);
  getItem("STORMMAXRADIUS", inputParameters.stormMaxRadius);

  getItem("DUSTLAT", inputParameters.stormLatitude);
  getItem("STORMLATITUDE", inputParameters.stormLatitude);

  getItem("DUSTLON", inputParameters.stormLongitude);
  getItem("STORMLONGITUDE", inputParameters.stormLongitude);

  // Wave Parameters
  getItem("WAVEFILE", inputParameters.waveFile);

  getItem("WSCALE", inputParameters.waveScale);
  getItem("WAVESCALE", inputParameters.waveScale);

  getItem("WAVEDATE", inputParameters.waveDate);

  getItem("WAVEA0", inputParameters.waveMeanOffset);
  getItem("WAVEMEANOFFSET", inputParameters.waveMeanOffset);

  getItem("WAVEA1", inputParameters.waveAmplitude1);
  getItem("WAVEAMPLITUDE1", inputParameters.waveAmplitude1);

  getItem("WAVEA2", inputParameters.waveAmplitude2);
  getItem("WAVEAMPLITUDE2", inputParameters.waveAmplitude2);

  getItem("WAVEA3", inputParameters.waveAmplitude3);
  getItem("WAVEAMPLITUDE3", inputParameters.waveAmplitude3);

  getItem("WAVEPHI1", inputParameters.wavePhase1);
  getItem("WAVEPHASE1", inputParameters.wavePhase1);

  getItem("WAVEPHI2", inputParameters.wavePhase2);
  getItem("WAVEPHASE2", inputParameters.wavePhase2);

  getItem("WAVEPHI3", inputParameters.wavePhase3);
  getItem("WAVEPHASE3", inputParameters.wavePhase3);

  getItem("PHI1DOT", inputParameters.wavePhase1Rate);
  getItem("WAVEPHASE1RATE", inputParameters.wavePhase1Rate);

  getItem("PHI2DOT", inputParameters.wavePhase2Rate);
  getItem("WAVEPHASE2RATE", inputParameters.wavePhase2Rate);

  getItem("PHI3DOT", inputParameters.wavePhase3Rate);
  getItem("WAVEPHASE3RATE", inputParameters.wavePhase3Rate);

  getItem("REQUA", inputParameters.equatorialRadius);
  getItem("EQUATORIALRADIUS", inputParameters.equatorialRadius);

  getItem("RPOLE", inputParameters.polarRadius);
  getItem("POLARRADIUS", inputParameters.polarRadius);

  int IUWAVE;
  getItem("IUWAVE", IUWAVE); // Not used.

  std::string gcmdir;
  getItem("GCMDIR", gcmdir); // Not used.

  checkForErrors();
}

} // namespace
