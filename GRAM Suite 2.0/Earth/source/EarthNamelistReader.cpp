//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include "EarthNamelistReader.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EarthNamelistReader::EarthNamelistReader()
{
}

//! \fn  EarthNamelistReader::EarthNamelistReader(const EarthNamelistReader& orig)
//! \brief Copying this object is discouraged.

//! \fn  EarthNamelistReader::~EarthNamelistReader()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param fileName Full or relative path to a namelist file.
//! \param inputParameters An EarthInputParameters object.
void EarthNamelistReader::getParameters(const std::string& fileName, EarthInputParameters& inputParameters)
{
  // Set the default time frame and scale.
  inputParameters.timeFrame = PET;
  inputParameters.timeScale = UTC;

  // Load up the itemMap from the namelist file.
  readFile(fileName);

  // Get the common input parameters
  NamelistReader::getParameters(inputParameters);

  // Get the duplicate named Earth input parameters
  // New name entries override old names.
  getItem("MN", inputParameters.month);
  getItem("MONTH", inputParameters.month);

  getItem("IDA", inputParameters.day);
  getItem("DAY", inputParameters.day);

  getItem("IYR", inputParameters.year);
  getItem("YEAR", inputParameters.year);

  getItem("IHRO", inputParameters.hour);
  getItem("HOUR", inputParameters.hour);

  getItem("MINO", inputParameters.minute);
  getItem("MINUTE", inputParameters.minute);

  getItem("SECO", inputParameters.seconds);
  getItem("SECONDS", inputParameters.seconds);

  getItem("NMAX", inputParameters.numberOfPositions);
  getItem("NUMBEROFPOSITIONS", inputParameters.numberOfPositions);

  getItem("MC", inputParameters.numberOfMonteCarloRuns);
  getItem("NUMBEROFMONTECARLORUNS", inputParameters.numberOfMonteCarloRuns);
  if (inputParameters.numberOfMonteCarloRuns == 0) {
    inputParameters.numberOfMonteCarloRuns = 1;
  }

  getItem("H1", inputParameters.initialHeight);
  getItem("INITIALHEIGHT", inputParameters.initialHeight);

  getItem("PHI1", inputParameters.initialLatitude);
  getItem("INITIALLATITUDE", inputParameters.initialLatitude);

  getItem("THET1", inputParameters.initialLongitude);
  getItem("INITIALLONGITUDE", inputParameters.initialLongitude);

  getItem("DPHI", inputParameters.deltaLatitude);
  getItem("DELTALATITUDE", inputParameters.deltaLatitude);

  getItem("DTHET", inputParameters.deltaLongitude);
  getItem("DELTALONGITUDE", inputParameters.deltaLongitude);

  getItem("DHGT", inputParameters.deltaHeight);
  getItem("DELTAHEIGHT", inputParameters.deltaHeight);

  getItem("DELT", inputParameters.deltaTime);
  getItem("DELTATIME", inputParameters.deltaTime);

  getItem("RPSCALE", inputParameters.randomPerturbationScale);
  getItem("RANDOMPERTURBATIONSCALE", inputParameters.randomPerturbationScale);

  getItem("RUSCALE", inputParameters.horizontalWindPerturbationScale);
  getItem("HORIZONTALWINDPERTURBATIONSCALE", inputParameters.horizontalWindPerturbationScale);

  getItem("RWSCALE", inputParameters.verticalWindPerturbationScale);
  getItem("VERTICALWINDPERTURBATIONSCALE", inputParameters.verticalWindPerturbationScale);

  if (getItem("ATMPATH", inputParameters.atmPath)) {
    processPath(inputParameters.atmPath);
  }

  if (getItem("NCEPPATH", inputParameters.NCEPPath)) {
    processPath(inputParameters.NCEPPath);
  }

  if (getItem("M2PATH", inputParameters.M2Path)) {
    processPath(inputParameters.M2Path);
  }

  getItem("PROFILE", inputParameters.auxiliaryAtmosphereFileName[0]);
  int iaux = 0;
  if (getItem("IAUX", iaux) && iaux == 0) {
    inputParameters.auxiliaryAtmosphereFileName[0].clear();
  }
  if (getItem("USEAUXILIARYATMOSPHERE", iaux) && iaux == 0) {
    inputParameters.auxiliaryAtmosphereFileName[0].clear();
  }

  getItem("TRAPATH", inputParameters.trajectoryFileName);
  getItem("TRAJECTORYFILENAME", inputParameters.trajectoryFileName);

  int iopt = 0;
  if (getItem("IOPT", iopt) && iopt == 0) {
    inputParameters.trajectoryFileName.clear();
  }
  if (getItem("USETRAJECTORYFILE", iopt) && iopt == 0) {
    inputParameters.trajectoryFileName.clear();
  }

  if (inputParameters.useLegacyOutputs == true) {
    // Legacy list output
    getItem("PRTPATH", inputParameters.listFileName);
    int iprt = 0;
    if (getItem("IPRT", iprt) && iprt == 0) {
      inputParameters.listFileName.clear();
    }

    // Legacy column output
    getItem("NPRPATH", inputParameters.columnFileName);
    int inpr = 0;
    if (getItem("INPR", inpr) && inpr == 0) {
      inputParameters.columnFileName.clear();
    }
  }
  else {
    string dummy;
    int fake;
    getItem("PRTPATH", dummy);
    getItem("IPRT", fake);
    getItem("NPRPATH", dummy);
    getItem("INPR", fake);
  }

  // Legacy species output
  getItem("CONPATH", inputParameters.speciesPath);
  int icon = 0;
  if (getItem("ICON", icon) && icon == 0) {
    inputParameters.speciesPath.clear();
  }

  // Legacy boundary layer output
  getItem("BLPATH", inputParameters.boundaryLayerPath);
  int ibltest = 0;
  if (getItem("IBLTEST", ibltest) && ibltest == 0) {
    inputParameters.boundaryLayerPath.clear();
  }

  int itherm;
  getItem("ITHERM", itherm);
  getItem("THERMOSPHEREMODEL", itherm);
  if (itherm == 2) {
    inputParameters.thermosphereModel = EG_MSIS;
  }
  else if (itherm == 3) {
    inputParameters.thermosphereModel = EG_JB2008;
  }
  else {
    inputParameters.thermosphereModel = EG_MET;
  }

  getItem("USENCEP", inputParameters.useNCEP);
  getItem("NCEPYR", inputParameters.NCEPYear);
  getItem("NCEPYEAR", inputParameters.NCEPYear);
  getItem("NCEPHR", inputParameters.NCEPHour);
  getItem("NCEPHOUR", inputParameters.NCEPHour);

  getItem("M2HOUR", inputParameters.M2Hour);
  getItem("M2MINIMUMLATITUDE", inputParameters.minimumLatitude);
  getItem("M2MAXIMUMLATITUDE", inputParameters.maximumLatitude);
  getItem("M2MINIMUMLONGITUDE", inputParameters.minimumLongitude);
  getItem("M2MAXIMUMLONGITUDE", inputParameters.maximumLongitude);

  // RRA Parameters
  int iurra = 0;
  if (getItem("IURRA", iurra)) {
    inputParameters.useRRA = (iurra != 0);
  }
  getItem("USERRA", inputParameters.useRRA);
  if (inputParameters.useRRA) {
    // Aux profiles with RRA is not allowed.
    inputParameters.auxiliaryAtmosphereFileName[0].clear();
  }

  if (getItem("RRAPATH", inputParameters.rraPath)) {
    processPath(inputParameters.rraPath);
  }

  getItem("RRALIST", inputParameters.rraSiteList);
  getItem("RRASITELIST", inputParameters.rraSiteList);

  int iyrrra = 0;
  if (getItem("IYRRRA", iyrrra)) {
    switch (iyrrra) {
    case 1:
      inputParameters.rraYear = RRA_1983;
      break;
    case 2:
      inputParameters.rraYear = RRA_2006;
      break;
    case 3:
      inputParameters.rraYear = RRA_2013;
      break;
    case 4:
      inputParameters.rraYear = RRA_2019;
      break;
    default:
      inputParameters.rraYear = RRA_2019;
      break;
    }
  }
  int rraYear = 0;
  if (getItem("RRAYEAR", rraYear)) {
    switch (rraYear) {
    case 1983:
      inputParameters.rraYear = RRA_1983;
      break;
    case 2006:
      inputParameters.rraYear = RRA_2006;
      break;
    case 2013:
      inputParameters.rraYear = RRA_2013;
      break;
    case 2019:
      inputParameters.rraYear = RRA_2019;
      break;
    default:
      inputParameters.rraYear = RRA_2019;
      break;
    }
  }

  if (getItem("SITENEAR", inputParameters.rraInnerRadius)) {
    inputParameters.innerRadius[0] = inputParameters.rraInnerRadius;
  }
  getItem("RRAINNERRADIUS", inputParameters.rraInnerRadius);

  if (getItem("SITELIM", inputParameters.rraOuterRadius)) {
    inputParameters.outerRadius[0] = inputParameters.rraOuterRadius;
  }
  getItem("RRAOUTERRADIUS", inputParameters.rraOuterRadius);

  getItem("F10", inputParameters.dailyF10);
  getItem("DAILYF10", inputParameters.dailyF10);

  getItem("F10B", inputParameters.meanF10);
  getItem("MEANF10", inputParameters.meanF10);

  getItem("S10", inputParameters.dailyS10);
  getItem("DAILYS10", inputParameters.dailyS10);

  getItem("S10B", inputParameters.meanS10);
  getItem("MEANS10", inputParameters.meanS10);

  getItem("XM10", inputParameters.dailyXM10);
  getItem("DAILYXM10", inputParameters.dailyXM10);

  getItem("XM10B", inputParameters.meanXM10);
  getItem("MEANXM10", inputParameters.meanXM10);

  getItem("Y10", inputParameters.dailyY10);
  getItem("DAILYY10", inputParameters.dailyY10);

  getItem("Y10B", inputParameters.meanY10);
  getItem("MEANY10", inputParameters.meanY10);

  getItem("AP", inputParameters.ap);

  getItem("DSTTEMPERATURECHANGE", inputParameters.dstdtc);
  getItem("DSTDTC", inputParameters.dstdtc);

  getItem("Z0IN", inputParameters.surfaceRoughness);
  getItem("SURFACEROUGHNESS", inputParameters.surfaceRoughness);

  getItem("INITPERT", inputParameters.initializePerturbations);
  getItem("INITIALIZEPERTURBATIONS", inputParameters.initializePerturbations);

  getItem("RDINIT", inputParameters.initialDensityPerturbation);
  getItem("INITIALDENSITYPERTURBATION", inputParameters.initialDensityPerturbation);

  getItem("RTINIT", inputParameters.initialTemperaturePerturbation);
  getItem("INITIALTEMPERATUREPERTURBATION", inputParameters.initialTemperaturePerturbation);

  getItem("RUINIT", inputParameters.initialEWWindPerturbation);
  getItem("INITIALEWWINDPERTURBATION", inputParameters.initialEWWindPerturbation);

  getItem("RVINIT", inputParameters.initialNSWindPerturbation);
  getItem("INITIALNSWINDPERTURBATION", inputParameters.initialNSWindPerturbation);

  getItem("RWINIT", inputParameters.initialVerticalWindPerturbation);
  getItem("INITIALVERTICALWINDPERTURBATION", inputParameters.initialVerticalWindPerturbation);

  double patchy;
  if (getItem("PATCHY", patchy)) {
    inputParameters.patchy = (patchy != 0.0);
  }

  getItem("CORRMONTE", inputParameters.corrMonte);
  getItem("CORRDELTAHOURS", inputParameters.corrDeltaHours);
  getItem("CORRMEAN", inputParameters.corrMean);

  getItem("MAXIMUMHEIGHT", inputParameters.maximumHeight);

  checkForErrors();
}

} // namespace
