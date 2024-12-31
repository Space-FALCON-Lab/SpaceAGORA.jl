//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "EarthAtmosphere.h"
#include "Interpolator.h"
#include "EarthAtmosphereState.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EarthAtmosphere::EarthAtmosphere()
  : EarthCommon(this)
{
  // Allocate new objects for pointers.
  earthModelPtr = new EarthModel();
  rraPtr = new RRA();

  // Add the Earth specific metrics to the AtmosphereState.
  atmos.setPlanetSpecificMetrics(earthAtmos);

  // Give the RRA object a reference to the atmosphere state.
  rraPtr->setInputState(atmos);

  // This model has vertical winds.
  setVerticalWinds(true);

  // This model requires standard deviations in the auxiliary profiles.
  setSDFlag(true);

  // Identity (for ephemeris).
  setGramBody(EARTH);
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
EarthAtmosphere::EarthAtmosphere(const EarthAtmosphere& orig) 
  : PerturbedAtmosphere(orig), EarthCommon(this)
{
  // Allocate new objects for pointers
  earthModelPtr = new EarthModel(*orig.earthModelPtr);
  rraPtr = new RRA(*orig.rraPtr);

  // Copy over data. 
  // (Some of this may not be necessary, but it is easier to just copy everything.)
  inputParameters = orig.inputParameters;
  earthAtmos = orig.earthAtmos;
  corrRandom = orig.corrRandom;
  initializing = orig.initializing;
  udCorrSmall = orig.udCorrSmall;
  vdCorrSmall = orig.vdCorrSmall;
  pdCorrLarge = orig.pdCorrLarge;
  waveRand = orig.waveRand;
  phi_q = orig.phi_q;
  A_Q = orig.A_Q;
  a_v = orig.a_v;
  b_v = orig.b_v;
  hwn = orig.hwn;
  T = orig.T;
  T_v = orig.T_v;
  ewWaveAngle = orig.ewWaveAngle;
  nsWaveAngle = orig.nsWaveAngle;
  ewCorrPhase = orig.ewCorrPhase;
  nsCorrPhase = orig.nsCorrPhase;
  xMeanSmall = orig.xMeanSmall;
  zMeanSmall = orig.zMeanSmall;
  xLengthSmall = orig.xLengthSmall;
  zLengthSmall = orig.zLengthSmall;
  xSDSmall = orig.xSDSmall;
  zSDSmall = orig.zSDSmall;
  latitudePrev = orig.latitudePrev;
  longitudePrev = orig.longitudePrev;
  heightPrev = orig.heightPrev;
  timePrev = orig.timePrev;
  baseDayOfYear = orig.baseDayOfYear;

  // Constructor settings.
  atmos.setPlanetSpecificMetrics(earthAtmos);
  rraPtr->setInputState(atmos);
  setVerticalWinds(true);
  setSDFlag(true);
  setGramBody(EARTH);
}

//! \copydoc Atmosphere::~Atmosphere()
EarthAtmosphere::~EarthAtmosphere()
{
  if (earthModelPtr != nullptr) {
    delete earthModelPtr;
  }
  if (rraPtr != nullptr) {
    delete rraPtr;
  }
}

//! \fn EarthAtmosphereState& EarthAtmosphere::getEarthAtmosphereState()
//! \brief Get Earth specific metrics.
//!
//! Metrics that are specific to the EarthGRAM model are stored in the EarthAtmosphereState. After a call to update(), 
//! one gets the common metrics with getAtmosphereState() and the Earth specific metrics with this method.
//! \returns An EarthAtmosphereState object.

//! \copydoc Atmosphere::getVersionString()
const std::string& EarthAtmosphere::getVersionString()
{
  static const string version = "EarthGRAM 2023 :: " + Atmosphere::getVersionString();
  return version;
}

//! \copydoc PerturbedAtmosphere::setInputParameters()
void EarthAtmosphere::setInputParameters(const EarthInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;
  dataPath = inputParameters.dataPath;

  if (inputParameters.atmPath.empty()) {
    atmPath = inputParameters.atmPath = inputParameters.dataPath + "/modeldata/";
  }

  if (inputParameters.NCEPPath.empty()) {
    inputParameters.NCEPPath = inputParameters.dataPath + "/NCEPdata/FixedBin/";
  }
  NCEPPath = inputParameters.NCEPPath;

  if (inputParameters.M2Path.empty()) {
    M2Path = inputParameters.M2Path = inputParameters.dataPath + "/MERRA2data/";
  }

  if (inputParameters.rraPath.empty()) {
    rraPath = inputParameters.rraPath = inputParameters.dataPath + "/RRAdata/";
  }

  // Apply input parameters pertinent to the base Atmosphere.
  PerturbedAtmosphere::setInputParameters(inputParameters);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(inputParameters);

  assert(earthModelPtr != nullptr);
  earthModelPtr->setInputParameters(inputParameters);

  assert(rraPtr != nullptr);
  rraPtr->setInputParameters(inputParameters);

  initializeData();

  setSeed(params.initialRandomSeed);

  // Set the flag denoting first pass of a trajectory.
  initializing = true;
}

//! \brief Set the thermosphere model.
//!
//! This parameter selects the model for thermosphere calculations.
//! \param model      A ThermosphereModelType (EG_MET, EG_MSIS, or EG_JB2008).
void EarthAtmosphere::setThermosphereModel(ThermosphereModelType model)
{
  inputParameters.thermosphereModel = model;
  assert(earthModelPtr != nullptr);
  earthModelPtr->setThermosphereModel(model);
}

//! \brief Set the default data path.
//!
//! The default data path must point to the folder containing the Earth atmosphere, NCEP, and RRA data files.
//! This folder should contain folders named "modeldata", "NCEPdata/FixedBin", and "RRAdata".
//! \param path      An absolute or relative path to the data folder.
void EarthAtmosphere::setDataPath(const std::string& path)
{
  inputParameters.dataPath = path;
  setAtmosPath(path + "/modeldata");
  setNCEPPath(path + "/NCEPdata/FixedBin/");
  setMERRA2Path(path + "/MERRA2data/");
  setRRAPath(path + "/RRAdata");
}

//! \brief Set the Earth atmosphere model data path.
//!
//! The Earth atmosphere data path must point to the folder containing the "sdata.txt, "topo.txt", and "zdata.txt" data files.
//! \param path      An absolute or relative path to the atmosphere data.
void EarthAtmosphere::setAtmosPath(const std::string& path)
{
  atmPath = inputParameters.atmPath = path;
}

//! \brief Determine use of NCEP model over default MERRA-2 model.
//!
//! \param useNCEP   Set to \p true to use the NCEP model.
void EarthAtmosphere::setUseNCEP(bool useNCEP) {
  inputParameters.useNCEP = useNCEP;
  assert(earthModelPtr != nullptr);
  earthModelPtr->setUseNCEP(useNCEP);
}

//! \brief Set the NCEP data path.
//!
//! The NCEP data path must point to the folder containing the fixed binary NCEP data files.
//! \param path      An absolute or relative path to the NCEP data.
void EarthAtmosphere::setNCEPPath(const std::string& path)
{
  NCEPPath = inputParameters.NCEPPath = path;
}

//! \brief Identify the desired NCEP data.
//!
//! The NCEP data is identified by specifying two codes, a year and hour code.  NCEP monthly 
//! climatology is determined by the value of month in the initial time input.
//! 
//! \param NCEPYear   The NCEPYear code y1y2 specifies NCEP climatology for period-of-record (POR) from year
//!                   y1 through year y2 (e.g. NCEPYear = 9008 for POR = 1990 through 2008). 
//! \param NCEPHour   The NCEPHour code specifies the UT hour of day for NCEP climatology:
//!                   1 = 00 UT, 2 = 06 UT, 3 = 12 UT, 4 = 18 UT, 5 = all times of day combined,
//!                   or 0 to use NCEP time-of-day based on input UTC hour.
void EarthAtmosphere::setNCEPParameters(int NCEPYear, int NCEPHour)
{
  inputParameters.NCEPYear = NCEPYear;
  inputParameters.NCEPHour = NCEPHour;
  assert(earthModelPtr != nullptr);
  earthModelPtr->setNCEPParameters(NCEPYear, NCEPHour);
}

//! \brief Set the MERRA2 data path.
//!
//! The MERRA2 data path must point to the folder containing the binary MERRA2 data files.
//! \param path      An absolute or relative path to the MERRA2 data.
void EarthAtmosphere::setMERRA2Path(const std::string& path)
{
  M2Path = inputParameters.M2Path = path;
}

//! \brief Identify the desired MERRA-2 data.
//!
//! The MERRA-2 data is identified by specifying an hour code.  MERRA-2 monthly 
//! climatology is determined by the value of month in the initial time input.
//! One can also reduce the amount of data read from the MERRA-2 data file by specifying
//! bounds on latitude and longitude.  Data analysis must, of course, fall in this range.
//! 
//! \param M2Hour   The M2Hour code specifies the UT hour of day for MERRA-2 climatology:
//!                 1 = 00 UT, 2 = 03 UT, 3 = 06 UT, 4 = 09 UT, 5 = 12 UT, 6 = 15 UT,
//!                 7 = 18 UT, 8 = 21 UT, 9 = all times of day combined,
//!                 or 0 to use time-of-day code based on input UTC hour.
//! \param latMin   The lower latitude boundary for reading MERRA-2 data. 
//! \param latMax   The upper latitude boundary for reading MERRA-2 data. 
//! \param lonMin   The lower longitude boundary for reading MERRA-2 data. 
//! \param lonMax   The upper longitude boundary for reading MERRA-2 data. 
void EarthAtmosphere::setMERRA2Parameters(int M2Hour, greal latMin, greal latMax, greal lonMin, greal lonMax) {

    inputParameters.M2Hour = M2Hour;
    inputParameters.minimumLatitude = latMin;
    inputParameters.maximumLatitude = latMax;
    inputParameters.minimumLongitude = lonMin;
    inputParameters.maximumLongitude = lonMax;
    assert(earthModelPtr != nullptr);
    earthModelPtr->setMERRA2Parameters(M2Hour, latMin, latMax, lonMin, lonMax);
}

//! \brief Set the RRA data path.
//!
//! The RRA data path must point to the folder containing the RRA sites, and RRA data files.
//! \param path      An absolute or relative path to the RRA data.
void EarthAtmosphere::setRRAPath(const std::string& path)
{
  rraPath = inputParameters.rraPath = path;
  assert(rraPtr != nullptr);
}

//! \brief Set the RRA site list file name.
//!
//! Use the method to override the default name of the RRA site list file. This file
//! is expected to be located in the folder specified by the RRA data path.
//! \param fileName      The name of the RRA site list file.
void EarthAtmosphere::setRRASiteList(const std::string& fileName)
{
  inputParameters.rraSiteList = fileName;
  assert(rraPtr != nullptr);
  rraPtr->setRRASiteList(fileName);
}

//! \brief Set the RRA parameters.
//!
//! This method selects the year of the RRA data set and defines a region in with RRA data is used.  There
//! is a smooth transition between the \p outerRadius (no RRA) and the \p innerRadius (all RRA) of the region.
//! \param year             RRA_1983, RRA_2006, or RRA_2013.
//! \param innerRadius      Lat-lon radius from RRA site, inside which RRA data is used with full weight of 1 \units{\text{degrees}}.
//! \param outerRadius      Lat-lon radius from RRA site, outside which RRA data are NOT used \units{\text{degrees}}.
void EarthAtmosphere::setRRAParameters(RRAYearType year, greal innerRadius, greal outerRadius)
{
  inputParameters.rraYear = year;
  inputParameters.rraInnerRadius = innerRadius;
  inputParameters.rraOuterRadius = outerRadius;
  assert(rraPtr != nullptr);
  rraPtr->setRRAParameters(year, innerRadius, outerRadius);
}

//! \brief Enable or disable the use of RRA data.
//!
//! \param useFlag     Set to \p true to enable the RRA data regions or \p false to disable.
void EarthAtmosphere::setUseRRA(bool useFlag)
{
  inputParameters.useRRA = useFlag;
  assert(rraPtr != nullptr);
  rraPtr->setUseRRA(useFlag);
}

//! \brief Set solar flux and index parameters.
//!
//! \param dailyF10      Daily 10.7-cm flux \units{sfu}.
//! \param meanF10       Mean 10.7-cm flux \units{sfu}.
//! \param ap            Geomagnetic index.
void EarthAtmosphere::setSolarParameters(greal dailyF10, greal meanF10, greal ap)
{
  inputParameters.dailyF10 = dailyF10;
  inputParameters.meanF10 = meanF10;
  inputParameters.ap = ap;
  assert(earthModelPtr != nullptr);
  earthModelPtr->setSolarParameters(inputParameters);
}

//! \brief Set parameters for the JB2008 thermosphere model.
//!
//! This method is not required when other thermosphere models are selected.
//! \param dailyS10      EUV index (26-34 nm) scaled to F10 units (0.0 -> dailyS10 = dailyF10).
//! \param meanS10       EUV 81-day center-averaged index (0.0 -> meanS10 = meanF10).
//! \param dailyXM10     MG2 index scaled to F10 units (0.0 -> dailyXM10 = dailyF10).
//! \param meanXM10      MG2 81-day center-averaged index (0.0 -> meanXM10 = meanF10).
//! \param dailyY10      Solar X-Ray & Lya index scaled to F10 (0.0 -> dailyY10 = dailyF10).
//! \param meanY10       Solar X-Ray & Lya 81-day avg. centered index (0.0 -> meanY10 = meanF10).
//! \param dstdtc        Temperature change computed from Dst index.
void EarthAtmosphere::setJB2008Parameters(greal dailyS10, greal meanS10, greal dailyXM10, greal meanXM10,
                                          greal dailyY10, greal meanY10, greal dstdtc)
{
  inputParameters.dailyS10 = dailyS10;
  inputParameters.meanS10 = meanS10;
  inputParameters.dailyXM10 = dailyXM10;
  inputParameters.meanXM10 = meanXM10;
  inputParameters.dailyY10 = dailyY10;
  inputParameters.meanY10 = meanY10;
  inputParameters.dstdtc = dstdtc;
  assert(earthModelPtr != nullptr);
  earthModelPtr->setJB2008Parameters(inputParameters);
}

//! \brief Set user defined the initial perturbations.
//!
//! Typically, the initial perturbations are randomly derived in the Earth model.  Calling this method
//! allows the user to set those initial perturbations.
//! \param densPert      Initial density perturbation value \units{\text{\% of mean}}.
//! \param tempPert      Initial temperature perturbation value \units{\text{\% of mean}}.
//! \param ewPert        Initial eastward velocity perturbation value \units{m/s}.
//! \param nsPert        Initial northward velocity perturbation value \units{m/s}.
//! \param verticalPert  Initial upward velocity perturbation value \units{m/s}.
void EarthAtmosphere::setInitialPerturbations(greal densPert, greal tempPert, greal ewPert, greal nsPert, greal verticalPert)
{
  inputParameters.initializePerturbations = true;
  inputParameters.initialDensityPerturbation = densPert;
  inputParameters.initialTemperaturePerturbation = tempPert;
  inputParameters.initialEWWindPerturbation = ewPert;
  inputParameters.initialNSWindPerturbation = nsPert;
  inputParameters.initialVerticalWindPerturbation = verticalPert;
}

//! \brief Set the perturbation scale factor for density, pressure, and temperature.
//!
//! \param scaleFactor Between 0 and 2.
void EarthAtmosphere::setRandomPerturbationScale(greal scaleFactor)
{
  PerturbedAtmosphere::setDensityPerturbationScale(scaleFactor);
  inputParameters.randomPerturbationScale = densityPerturbationScale;
  earthModelPtr->setPerturbationScales(inputParameters);
  rraPtr->setInputParameters(inputParameters);
}

//! \brief Set the perturbation scale factor for E/W and N/S winds.
//!
//! \param scaleFactor Between 0 and 2.
void EarthAtmosphere::setHorizontalWindPerturbationScale(greal scaleFactor)
{
  PerturbedAtmosphere::setEWWindPerturbationScale(scaleFactor);
  inputParameters.horizontalWindPerturbationScale = ewWindPerturbationScale;
  earthModelPtr->setPerturbationScales(inputParameters);
  rraPtr->setInputParameters(inputParameters);
}

//! \copydoc PerturbedAtmosphere::setVerticalWindPerturbationScale
void EarthAtmosphere::setVerticalWindPerturbationScale(greal scaleFactor)
{
  PerturbedAtmosphere::setVerticalWindPerturbationScale(scaleFactor);
  inputParameters.verticalWindPerturbationScale = verticalWindPerturbationScale;
  earthModelPtr->setPerturbationScales(inputParameters);
  rraPtr->setInputParameters(inputParameters);
}

//! \copydoc PerturbedAtmosphere::setStartTime
void EarthAtmosphere::setStartTime(const GramTime& gramTime)
{ 
  PerturbedAtmosphere::setStartTime(gramTime); 
  inputParameters.timeScale = UTC;
  inputParameters.timeFrame = PET;
  time.getStartTime(inputParameters.timeScale, inputParameters.timeFrame, 
    inputParameters.year, inputParameters.month, inputParameters.day,
    inputParameters.hour, inputParameters.minute, inputParameters.seconds);
  earthModelPtr->setInputParameters(inputParameters);
}

//! \fn  EarthAtmosphere::setPatchiness(bool usePatchiness)
//! \brief Enable patchiness in the perturbation model.
//! \param usePatchiness   Set to \p true to enable patchiness.

//! \fn  EarthAtmosphere::setSurfaceRoughness(greal z0)
//! \brief Set the surface roughness for the sigma-w model.
//! \param z0   Surface roughness for the sigma-w model 
//! ( < 0 to use 1-by-1 deg lat-lon surface data; 
//!   = 0 for speed-dependent z0 over water;
//!   or enter a value between 1.0e-5 and 3 for a user-specified z0 value).

//! \brief Interface for the primary atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//!
//! \b Inputs
//! \arg #position      
//!
//! \returns The AtmosphereState and EarthAtmosphereState are populated.
void EarthAtmosphere::update()
{
  if (initializing) {
    rraPtr->setMonth(month);
  }

  // Calculate suface height from topographic map
  updateSurfaceHeight();

  // Update radius and gravity
  updatePosition();

  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Throw error if height is too negative.
  if (height < -0.031_km) {
    // For external sims, only throw this error once.
    if (heightErrorThrown) {
      return;
    }
    heightErrorThrown = true;
    throw string("Error: Height below -31 meters.\n"
                 "       This is an unrecoverable error.");
  }

  // Update ephemeris values.
  updateEphemeris();
  // On the first pass of a trajectory, perturbations get special treatment.
  if (initializing) {
    initializePerturbations();
    time.setElapsedTime(0.0);
    baseDayOfYear = time.getDayOfYear();
    int y, mm, d, h, m;
    greal s;
    time.getStartTime(UTC, PET, y, mm, d, h, m, s);
    leapYear = y % 4;
    time.setElapsedTime(elapsedTime);
  }

  // Set up the earth model
  EarthModel& earthModel = *earthModelPtr;
  earthModel.setPosition(position);
  // Set the current ephemeris state.
  earthModel.setEphemerisState(ephem);
  greal jd = 0.0;
  if (inputParameters.thermosphereModel == EG_JB2008) {
    time.getTime(time.getTimeScale(), time.getTimeFrame(), jd);
  }
  const greal daysPerYear = 365.0;
  greal dayOfYear = baseDayOfYear + elapsedTime / ephem.secondsPerSol;
  int leap = leapYear;
  while (dayOfYear > daysPerYear + (leap == 0 ? 1 : 0)) {
    dayOfYear -= daysPerYear + (leap == 0 ? 1 : 0);
    leap = (leap + 1) % 4;
  }
  earthModel.setDayOfYear(dayOfYear, jd);

  // Compute the atmosphere state.
  earthModel.update();

  // Get the metrics from the earth model.
  const AtmosphereState& em = earthModel.getAtmosphereState();
  getStateVariables(em);

  // Apply perturbations.
  if (!initializing) {
    updatePerturbations();
  }
  else {
    updateRRA();
    updateStatus = NO_UPDATES;
  }
  previousPosition = savedPosition;

  // This section ensures that random number generation aligns correctly
  // when running the CorrMonte program.
  if (inputParameters.corrMonte) {
    if (initializing && elapsedTime != 0.0) {
      if (!corrRandom.initialized()) {
        corrRandom = randomNumberGenerator;
      }
      vector<double> r(6);
      corrRandom.getRandomNumbers(r);
      randomNumberGenerator = corrRandom;
    }
  }

  // First perturbation pass is complete.  Turn off the flag.
  initializing = false;

  // Update the psychrometrics with any RRA data.
  rraPtr->setEpsilon(earthModel.getEpsilon());
  rraPtr->setPpmToND(earthModel.getPpmToND());
  rraPtr->updatePsychrometrics();

  // Update gases after the (possible) RRA update.
  updateGases();

  // Update the atmosphere state with any auxilary atmosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  //Calculate total perturbed values added to mean data.
  perturbedPressure = pressure * (1.0 + pressurePerturbation);
  perturbedDensity = density * (1.0 + densityPerturbation);
  perturbedTemperature = temperature * (1.0 + temperaturePerturbation);
  perturbedEWWind = ewWind + ewWindPerturbation;
  perturbedNSWind = nsWind + nsWindPerturbation;
  perturbedVerticalWind = verticalWind + verticalWindPerturbation;
  lowDensity = density * (1.0 - densityStandardDeviation);
  highDensity = density * (1.0 + densityStandardDeviation);

  // Update perturbed winds and wind speed.
  updateWindSpeed();

  // Pressure scale height (km), density scale height (km)
  earthModel.getScaleHeights(pressure, density, temperature, pressureScaleHeight, densityScaleHeight);

  // Update speed of sound based on constituent gases.
  specificGasConstant = pressure / (density*temperature);
  updateSpecificHeatRatio();
  updateSpeedOfSound();

  // Update metrics that depend on atmosphere state.
  updateMetrics();

  // Find geodetic coordinates.  First copy the geocentric position.
  planetodeticPosition = position;
  // Convert geocentric to geodetic.
  planetodeticPosition.convertToPlanetodetic(polarRadius, equatorialRadius);
  // Save for output.
  earthAtmos.geodeticLatitude = planetodeticPosition.latitude;
}

//! \brief Computes the speed of sound.
//!
//! The routine computes the speed of sound and the perturbed speed of sound
//! given the current atmosphere state.
//!
//! \b Inputs
//! \arg #pressure          
//! \arg #density           
//! \arg #temperature 
//! \arg #perturbedTemperature 
//!
//! \retval #speedOfSound     
void EarthAtmosphere::updateSpeedOfSound()
{
  Atmosphere::updateSpeedOfSound();
  perturbedSpeedOfSound = sqrt(specificHeatRatio * perturbedPressure / perturbedDensity);
}

//! \brief Computes the surface height for the current latitude and longitude.
//!
//! This routine interpolates the surface height array ztopo to geocentric latitude and longitude.
//!
//! \b Inputs
//! \arg #latitude          
//! \arg #longitude           
//!
//! \retval #surfaceHeight     
void EarthAtmosphere::updateSurfaceHeight()
{
  // Lower longitude index, (-179.5 -> 0, -178.5 -> 1, ..., 179.5 -> 359)
  // OR (0.5 -> 180, 1.5 -> 181, ..., 179.5 -> 359, 180.5 -> 0, 181.5 -> 1, ..., 359.5 -> 179)
  int iLon = int(longitude + 179.5) % 360;

  // Upper longitude index (use modulo 360 near prime meridian)
  int iLonUp = (iLon + 1) % 360;

  // relative longitude deviation from corner reference point
  greal dLon = longitude - (iLon + 180.5);
  if (dLon < 0.0) {
    dLon += 360.0;
  }

  // Lower latitude index (-89.5 -> 0, -88.5 -> 1, ..., 89.5 -> 179)
  int iLat = int(latitude + 89.5);

  // Upper latitude index
  int iLatUp = iLat + 1;

  // relative latitude deviation from corner reference location
  greal dLat = latitude - (iLat - 89.5);

  // Handle the cases near the poles (abs(latitude) > 89.5).
  // Use values at +-89.5 (constant extrapolation)
  if (latitude < -89.5) {
    iLat = 0;
    iLatUp = 1;
    dLat = 0.0;
  }
  else if (latitude > 89.5) {
    iLat = 178;
    iLatUp = 179;
    dLat = 1.0;
  }

  //lat-lon interpolation of topographic height
  Interpolator interp(dLon, dLat);
  surfaceHeight = interp.linear(ztopo[iLon][iLat], ztopo[iLon][iLatUp], ztopo[iLonUp][iLat], ztopo[iLonUp][iLatUp]);
}

//! \brief Copies select state variables.
//!
//! This routine copies select state variables from the input AtmosphereState to the current
//! AtmosphereState (atmos).
//!
//! \param state The AtmosphereState to copy.           
void EarthAtmosphere::getStateVariables(const AtmosphereState& state)
{
  // Standard metrics.
  temperature = state.temperature;
  pressure = state.pressure;
  density = state.density;
  ewWind = state.ewWind;
  nsWind = state.nsWind;
  verticalWind = state.verticalWind;

  // Gases.
  nitrogen = state.nitrogen;
  dinitrogen = state.dinitrogen;
  dioxygen = state.dioxygen;
  oxygen = state.oxygen;
  argon = state.argon;
  helium = state.helium;
  hydrogen = state.hydrogen;
  ozone = state.ozone;
  nitrousOxide = state.nitrousOxide;
  carbonMonoxide = state.carbonMonoxide;
  methane = state.methane;
  carbonDioxide = state.carbonDioxide;
  water = state.water;

  // Earth specific metrics.
  const EarthAtmosphereState& earth = state.getPlanetSpecificMetrics<EarthAtmosphereState>();
  vaporPressure = earth.vaporPressure;
  vaporPressureSD = earth.vaporPressureSD;
  vaporDensity = earth.vaporDensity;
  vaporDensitySD = earth.vaporDensitySD;
  dewPoint = earth.dewPoint;
  dewPointSD = earth.dewPointSD;
  relativeHumidity = earth.relativeHumidity;
  relativeHumiditySD = earth.relativeHumiditySD;
  windSpeed = earth.windSpeed;
  windSpeedStandardDeviation = earth.windSpeedStandardDeviation;

  if (!initializing) {
    temperatureStandardDeviation = state.temperatureStandardDeviation;
    pressureStandardDeviation = state.pressureStandardDeviation;
    densityStandardDeviation = state.densityStandardDeviation;
    ewStandardDeviation = state.ewStandardDeviation;
    nsStandardDeviation = state.nsStandardDeviation;
  }

  pressureAtSurface = state.pressureAtSurface;
  earthAtmos.densityAtSurface = earth.densityAtSurface;
  ewWindAtSurface = earth.ewWindAtSurface;
  nsWindAtSurface = earth.nsWindAtSurface;
  ewWindSDAtSurface = earth.ewWindSDAtSurface;
  nsWindSDAtSurface = earth.nsWindSDAtSurface;
  temperatureAtSurface = earth.temperatureAtSurface;
  windCorrelationAtSurface = earth.windCorrelationAtSurface;
  windSpeedAtSurface = earth.windSpeedAtSurface;
  windSpeedSDAtSurface = earth.windSpeedSDAtSurface;
  if (!initializing) {
    temperatureSDAtSurface = earth.temperatureSDAtSurface;
    windCorrelation = earth.windCorrelation;
  }
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void EarthAtmosphere::updatePressureAtSurface()
{
  // Nothing to do here.  This value was computed in the NCEP updates.
}

//! \copydoc Atmosphere::updateReferenceValues()
void EarthAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the US-76 model.
  getUS76StandardAtmosphere(referenceTemperature, referencePressure, referenceDensity);
}

//! \brief Update perturbed winds and wind speed.
//!
//! Perturbed winds are limited by the speed of sound.  Wind speed and standard deviation are
//! weighted with RRA or profile values.
//!
//! \b Inputs
//! \arg #perturbedPressure          
//! \arg #perturbedDensity           
//! \arg #perturbedEWWind 
//! \arg #perturbedNSWind 
//! \arg #ewWind 
//! \arg #nsWind 
//! \arg #ewStandardDeviation 
//! \arg #nsStandardDeviation 
//! \arg #windSpeed     
//! \arg #windSpeedStandardDeviation     
//!
//! \retval #perturbedEWWind     
//! \retval #perturbedNSWind     
//! \retval #windSpeed     
//! \retval #windSpeedStandardDeviation     
void EarthAtmosphere::updateWindSpeed()
{
  // Ensure wind component magnitudes < 0.7 times sound speed
  greal speedLimit = 0.7 * sqrt(1.4 * perturbedPressure / perturbedDensity);
  if (abs(perturbedEWWind) > speedLimit) {
    perturbedEWWind = copysign(speedLimit, perturbedEWWind);
  }
  if (abs(perturbedNSWind) > speedLimit) {
    perturbedNSWind = copysign(speedLimit, perturbedNSWind);
  }

  //Weight RRA wind speed data when it is available.
  rraPtr->updateWindSpeed();

  // Update when auxiliary atmospheres are in effect.
  greal profileWeight = atmos.profileWeight[0];
  if (hasAuxiliaryAtmopheres() && profileWeight > 0.0) {
    greal FS = square(ewStandardDeviation) + square(nsStandardDeviation);
    windSpeed = (sqrt(square(ewWind) + square(nsWind) + 0.605 * FS) * profileWeight)
      + windSpeed * (1 - profileWeight);
    windSpeedStandardDeviation = (sqrt(0.3950 * FS) * profileWeight) + windSpeedStandardDeviation * (1 - profileWeight);
  }
}

//! \brief Rescales gases to correct for water vapor.
//!
//! Rescale non-water-vapor concentrations to correct for the effects of 
//! interpolation, fairing, and water vapor.  Updates total number density and
//! average molecular weight.
//!
//! \b Inputs
//! \arg #pressure          
//! \arg #temperature          
//! \arg #gases          
//!
//! \retval #gases     
//! \retval #totalNumberDensity     
//! \retval #averageMolecularWeight     
void EarthAtmosphere::updateGases() {
  // Iterate over the set of all activated gases to find total moles.
  // This macro defines "gas" and "gasType" for each activated gas.
  greal totalMole = 0.0;
  FOR_ALL_GASES_PRESENT(
    if (gasType != WATER) {
      totalMole += gas.moleFraction;
    }
  );

  // Rescale non-water-vapor concentrations to correct for the effects of 
  // interpolation, fairing, and water vapor.
  greal rescale = (1.0 - water.moleFraction) / totalMole;
  
  greal moleFractionToNumberDensity = AVOGADRO * pressure / (UNIVERSAL_GAS * temperature);

  // Iterate over the set of all activated gases to rescale and compute total number density.
  // This macro defines "gas" and "gasType" for each activated gas.
  totalNumberDensity = 0.0;
  FOR_ALL_GASES_PRESENT(
    if (gasType != WATER) {
      gas.moleFraction *= rescale;
    }
    gas.numberDensity = gas.moleFraction * moleFractionToNumberDensity;
    totalNumberDensity += gas.numberDensity;
  );

  // Update the average molecular weight and mass fractions.
  averageMolecularWeight = getTotalMass() / totalNumberDensity;
  updateMassFractions();
}

//********************************* Perturbation Factor Methods ********************************//

////! \brief Get high and low perturbation factors.
////!
////! Not used in this perturbation model
////!
////! \retval pertLow The low perturbation factor.
////! \retval pertHigh The high perturbation factor.
//void EarthAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
//{
//  pertHigh = 1.0;
//  pertLow = 1.0;
//}
//
////! \brief Get vertical and horizontal scale parameters.
////!
////! Not used in this perturbation model
////!
////! \retval verticalScale The vertical scale parameter.
////! \retval horizontalScale The horizontal scale parameter.
//void EarthAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
//{
//  verticalScale = 1.0;
//  horizontalScale = 1.0;
//}
//
////! \brief Get east/west and north/south wind deviations.
////!
////! Not used in this perturbation model
////!
////! At this time, there is no winds model for Earth.
////! \retval ewStdDev The east/west wind deviation.
////! \retval nsStdDev The north/south wind deviation.
//void EarthAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
//{
//  ewStdDev = 0.0;
//  nsStdDev = 0.0;
//  vertStdDev = 0.0;
//}
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


//! \brief The initialization method for perturbation calculations.
//!
//! This method is called at the first time step of a trajectory.  It calculates
//! the initial perturbations for pressure, density, temperature, east-west wind,
//! north-south wind, and vertical wind.
void EarthAtmosphere::initializePerturbations()
{
  // Do this so we can use dots instead of arrows.
  EarthModel& earthModel = *earthModelPtr;

  // Save the lat/lon for next round computation of dx.
  latitudePrev = latitude;
  longitudePrev = longitude;
  heightPrev = height;
  timePrev = elapsedTime;

  // Get u-density and u-v cross correlations from global data
  interpolatePCData(udl, uvt, height, latitude, udCorrLarge, windCorrelation);
  interpolatePCData(uds, vds, height, latitude, udCorrSmall, vdCorrSmall);

  // Initialize standard deviations using NCEP data.
  earthModel.setPosition(position);
  earthModel.getStandardDeviations(atmos, true);
  verticalStandardDeviation = interpolateRSData(wr, height);

  // Fair current state with RRA state.
  updateRRA(true);

  pressureStandardDeviation = sqrt(abs(pressureStandardDeviation));
  densityStandardDeviation = sqrt(abs(densityStandardDeviation));
  temperatureStandardDeviation = sqrt(abs(temperatureStandardDeviation));

  //If present, use auxiliary profile data
  if (hasAuxiliaryAtmopheres()) {
    updateAuxiliaryStandardDeviations(position, atmos);
  }

  // Update large and small scale standard deviations
  updateStandardDeviations();

  // Update the pressure-density correlation
  updatePDCorrelation();

  // Get 7 random numbers
  vector<double> randomNumber(10);    // (was rvec)
  randomNumberGenerator.getRandomNumbers(randomNumber);

  // Evaluate initial small-scale perturbations values from the initial standard deviations
  densPertSmall = densSDSmall * normalPercentagePoint(randomNumber[0]);

  presPertSmall = presSDSmall * (pdCorrLarge * densPertSmall / densSDSmall
           + sqrt(1 - pow(pdCorrLarge, 2)) * normalPercentagePoint(randomNumber[1]));

  // Compute new small-scale perturbations for temperature
  // Adjust temperature perturbation assuming small-scale = large-scale p-d correlation
  greal tempScaleSmall = sqrt(abs(square(presSDSmall) + square(densSDSmall)
                                  - 2.0 * pdCorrLarge * presSDSmall * densSDSmall)); // (was stsc)
  tempPertSmall = (presPertSmall - densPertSmall) * tempSDSmall / tempScaleSmall;

  // Get initial small-scale velocity perturbations
  ewWindPertSmall = ewWindSDSmall * normalPercentagePoint(randomNumber[2]);
  nsWindPertSmall = nsWindSDSmall * normalPercentagePoint(randomNumber[3]);

  // Compute large-scale wave amplitude factor (equation 49) and wave phase (in equation 48)
  A_Q = 0.4808 + 0.96 * randomNumber[6];
  phi_q = TWO_PI * randomNumber[4];

  // Compute randomized wave factors
  waveRand = randomNumber[5];
  greal Q_nm = normalPercentagePoint(waveRand);  // (was dphid)
  if (abs(Q_nm) > 3.0) {
    Q_nm = copysign(3.0, Q_nm);
  }
  // Randomized factor in equation 51.
  a_v = 25.0 + 3.0 * Q_nm;
  // Horizontal wave numbers for wave perturbations (m and n in equation 50).
  // Note that here it is assumed that m = n.
  hwn = round(4.0 + 0.833 * Q_nm);

  // Period and wavelengt parameters for large-scale perturbations (T in equation 48)
  T = 2.0 + 60.0 * waveRand;             // (was wper)
  T_v = 2.0 + 40.0 * pow(waveRand, 2);   // (was wperv)

  // Update the large scale perturbations and the boundary layer model.
  updateLargeScalePerturbations(1.0);

  // Calculate sigma-w from boundary layer (BL) model
  updateBoundaryLayerModel(true);

  // Randomize the vertical wind perturbation.
  verticalWindPerturbation = verticalStandardDeviation * normalPercentagePoint(randomNumber[7]);

  // Substitute user-selected initial perturbations if provided
  if (initPerturbations) {
    presPertSmall = (initialDensityPerturbation + initialTemperaturePerturbation) / 100.0 - presPertLarge;
    densPertSmall = initialDensityPerturbation / 100.0 - densPertLarge;
    tempPertSmall = initialTemperaturePerturbation / 100.0 - tempPertLarge;
    ewWindPertSmall = initialEWWindPerturbation - ewWindPertLarge;
    nsWindPertSmall = initialNSWindPerturbation - nsWindPertLarge;
    verticalWindPerturbation = initialVerticalWindPerturbation;
  }

  //Compute variable small-scale lengths x = horizontal, z = vertical
  xSDSmall = interpolateRSData(xsigl, height);
  zSDSmall = interpolateRSData(zsigl, height);
  xMeanSmall = interpolateRSData(xlbar, height);
  zMeanSmall = interpolateRSData(zlbar, height);
  xLengthSmall = xMeanSmall + xSDSmall * normalPercentagePoint(randomNumber[8]);
  zLengthSmall = zMeanSmall + zSDSmall * normalPercentagePoint(randomNumber[9]);

  //Total perturbation
  densityPerturbation = densPertSmall + densPertLarge;
  temperaturePerturbation = tempPertSmall + tempPertLarge;
  pressurePerturbation = presPertSmall + presPertLarge;
  ewWindPerturbation = ewWindPertSmall + ewWindPertLarge;
  nsWindPerturbation = nsWindPertSmall + nsWindPertLarge;
}

//! \brief The primary interface for perturbation calculations.
//!
//! This method should be called after the atmospheric state has been updated
//! with EarthModel::update().  It computes the large-scale and small-scale perturbation 
//! values in pressure, density, temperature, horizontal and vertical wind components.
void EarthAtmosphere::updatePerturbations(greal meanDensity, greal meanEWWind, greal meanNSWind)
{
  greal xLengthSmallPrev = xLengthSmall;
  greal xSDSmallPrev = xSDSmall;
  greal xMeanSmallPrev = xMeanSmall;
  greal zLengthSmallPrev = zLengthSmall;
  greal zSDSmallPrev = zSDSmall;
  greal zMeanSmallPrev = zMeanSmall;

  //Compute parameters for variable small scale lengths
  greal xMin, zMin;
  greal xScaleSmall, zScaleSmall;  // (was xlen zlen)
  xMeanSmall = interpolateRSData(xlbar, height);
  zMeanSmall = interpolateRSData(zlbar, height);
  xSDSmall = interpolateRSData(xsigl, height);
  zSDSmall = interpolateRSData(zsigl, height);
  xMin = interpolateRSData(xlmin, height);
  zMin = interpolateRSData(zlmin, height);
  xScaleSmall = interpolateRSData(xscale, height);
  zScaleSmall = interpolateRSData(zscale, height);

  greal vertSDPrev = verticalStandardDeviation;

  // horizontal great circle distance between successive positions
  greal dx = totalRadius * toRadians(getArcAngle(latitude, longitude, latitudePrev, longitudePrev));
  greal deltaHeight = height - heightPrev;
  greal deltaTime = elapsedTime - timePrev;

  //If certain conditions are met new perturbations are not calculated
  if ((dx < nearlyZero && deltaHeight == 0.0 && deltaTime == 0.0) || perturbationAction == DO_NOT_UPDATE_PERTS) {
    xLengthSmall = xMeanSmall + (xLengthSmallPrev - xMeanSmallPrev) / xSDSmallPrev;
    zLengthSmall = zMeanSmall + (zLengthSmallPrev - zMeanSmallPrev) / zSDSmallPrev;
    updateRRA();
    updateStatus = NO_UPDATES;
    return;
  }

  // Save the lat/lon for next round computation of dx.
  latitudePrev = latitude;
  longitudePrev = longitude;
  heightPrev = height;
  timePrev = elapsedTime;

  greal presPertSmallPrev = presPertSmall;  // previous cycle prhs
  greal densPertSmallPrev = densPertSmall;
  greal ewPertSmallPrev = ewWindPertSmall;
  greal nsPertSmallPrev = nsWindPertSmall;
  greal presSDSmallPrev = presSDSmall;  // previous cycle sphs
  greal densSDSmallPrev = densSDSmall;
  greal tempSDSmallPrev = tempSDSmall;
  greal ewSDSmallPrev = ewWindSDSmall;
  greal nsSDSmallPrev = nsWindSDSmall;
  greal vertPertPrev = verticalWindPerturbation;

  // Get 7 random numbers
  vector<double> randomNumber(7); 
  randomNumberGenerator.getRandomNumbers(randomNumber);

  // Compute coefficients for small scale length model
  // Time scale (sec) for small-scale perturbations
  greal timeScaleSmall = 86400.0 * 0.735 * pow(abs(height), 0.116);   // (was tpers)
  timeScaleSmall = max(10800.0, timeScaleSmall);
  greal xPrevCoef = getCorrelation(abs(dx / xScaleSmall)  // (was as)
                     + abs(deltaHeight / zScaleSmall)
                     + abs(deltaTime / timeScaleSmall));
  greal xRandCoef = sqrt(1.0 - xPrevCoef * xPrevCoef);  // (was bs)

  // Compute small horizontal scale
  xLengthSmall = xMeanSmall + xSDSmall * (xPrevCoef * (xLengthSmallPrev - xMeanSmallPrev)
    / xSDSmallPrev + xRandCoef * normalPercentagePoint(randomNumber[5]));

  // Evaluate xzCorr = correlation between horizontal and vertical small scales
  greal xzCorr;   // (was rcor)
  if (height >= 200.0) {
    xzCorr = 0.9;
  }
  else {
    xzCorr = 0.5 + 0.002 * height;
  }
  greal denom = 1.0 - square((xzCorr * xPrevCoef));
  greal zPrevCoef = xPrevCoef * (1.0 - square(xzCorr)) / denom;  // (was alph)
  greal zCorrCoef = xzCorr * (1.0 - square(xPrevCoef)) / denom;  // (was beta)
  greal zRandCoef = sqrt(abs(1.0 - square(zPrevCoef) - square(zCorrCoef)  // (was gamm)
                              - 2.0 * zPrevCoef * zCorrCoef * xzCorr * xPrevCoef));

  //Compute small vertical scale (was zlh)
  zLengthSmall = zMeanSmall + zSDSmall * (zPrevCoef * (zLengthSmallPrev - zMeanSmallPrev) / zSDSmallPrev
                        + zCorrCoef * (xLengthSmall - xMeanSmall) / xSDSmall
                        + zRandCoef * normalPercentagePoint(randomNumber[6]));

  // Compute severeFactor for severe perturbed standard deviations
  greal severeFactor;  // (was fbig)
  if (height <= 10.0) {
    severeFactor = 6.0 + 0.6 * height;
  }
  else if (height >= 16.0) {
    severeFactor = 6.0;
  }
  else {
    severeFactor = 22.0 - height;
  }

  // Compute probability for severe perturbed conditions
  greal xTail = (xMeanSmall - xMin) / xSDSmall;
  greal zTail = (zMeanSmall - zMin) / zSDSmall;
  xTail = min(zTail, xTail);
  greal xProb = getProbabilityTail(xTail);  // (was ptailx)

  // Compute factor for perturbed length scales
  greal scaleFactor = xProb * xTail / (1.0 - xProb);

  // Compute normalFactor for normal perturbation conditions
  greal normalFactor = (1.0 - severeFactor * xProb) / (1.0 - xProb);  // (was fsmall)
  if (normalFactor < 0.0005) {
    normalFactor = 0.0005;
    severeFactor = (0.9995 + 0.0005 * xProb) / xProb;
  }

  // Compute scale sizes for normal perturbed conditions
  greal xMax = xMeanSmall + scaleFactor * xSDSmall;
  greal zMax = zMeanSmall + scaleFactor * zSDSmall;

  // Small-scale multiplier for better comparison with KSC shears
  greal zfactor = 1.0;
  if (height <= 30.0) {
    zfactor = exp(1.4 * sin(TWO_PI * height / 30.0));
  }

  // If severe perturbed conditions, set scales to min and standard deviation to big
  greal zScale, xScale, densFactor;  // (was vdsm=vts=vus, hls, densfact)
  if (patchy && ((zLengthSmall < zMin) || (xLengthSmall < xMin))) {
    zScale = zMin * zfactor;
    xScale = xMin * zfactor * 2.0;
    densFactor = sqrt(severeFactor);
  }
  // If normal perturbed conditions set scales to max and standard deviation to small
  else {
    zScale = zMax * zfactor / 2.0;
    xScale = xMax * zfactor * 2.0;
    densFactor = sqrt(normalFactor);
  }

  // Compute relative separations (x/Lx and z/Lz) for small and large-scale models, and compute the
  // correlations across the separation distance with the correl function. Include small-scale time
  // effect with explicit time scale tpers
  greal horizAndTime = abs(dx / xScale) + abs(deltaTime / timeScaleSmall);
  greal densDistCorr = getCorrelation(horizAndTime + abs(deltaHeight / zScale));  // (was rds)
  greal presDistCorr = getCorrelation((horizAndTime + abs(deltaHeight / zScale)) / 3.0);  // (was rps)
  greal windDistCorr = getCorrelation(horizAndTime + abs(deltaHeight / zScale));  // (was rvs)

  // Save density-velocity correlations
  greal udCorrSmallPrev = udCorrSmall;
  greal vdCorrSmallPrev = vdCorrSmall;

  // Get u-density and u-v cross correlations from global data
  double wc;
  interpolatePCData(udl, uvt, height, latitude, udCorrLarge, wc);
  interpolatePCData(uds, vds, height, latitude, udCorrSmall, vdCorrSmall);

  if (windCorrelation == 0.0) {
    windCorrelation = wc;
  }

  // Calculate standard deviations for calculating perturbations 
  verticalStandardDeviation = interpolateRSData(wr, height);

  // Fair current state with RRA state.
  updateRRA();

  // Use auxiliary profile data if present
  updateAuxiliaryStandardDeviations(position, atmos);

  // Check for Buell constraint violations
  if (densityStandardDeviation >= BuellLimit * (pressureStandardDeviation + temperatureStandardDeviation)) {
    densityStandardDeviation = BuellLimit * (pressureStandardDeviation + temperatureStandardDeviation);
  }
  else if (temperatureStandardDeviation >= BuellLimit * (pressureStandardDeviation + densityStandardDeviation)) {
    temperatureStandardDeviation = BuellLimit * (pressureStandardDeviation + densityStandardDeviation);
  }
  else if (pressureStandardDeviation >= BuellLimit * (densityStandardDeviation + temperatureStandardDeviation)) {
    pressureStandardDeviation = BuellLimit * (densityStandardDeviation + temperatureStandardDeviation);
  }

  // Update large and small scale standard deviations
  updateStandardDeviations();

  // Update the pressure-density correlation
  updatePDCorrelation();

  // Evaluate the perturbation model coefficients (*Coef) for the small scale model
  greal densPrevCoef, densRandCoef, presPrevCoef, presDensCoef, presRandCoef;
  greal ewPrevCoef, ewDensCoef, ewRandCoef, nsPrevCoef, nsDensCoef, nsRandCoef;
  getPerturbationCoefficients(
      densSDSmallPrev, densSDSmall, presSDSmallPrev, presSDSmall, tempSDSmallPrev, tempSDSmall,
      ewSDSmallPrev, ewWindSDSmall, nsSDSmallPrev, nsWindSDSmall, udCorrSmallPrev, udCorrSmall,
      vdCorrSmallPrev, vdCorrSmall, pdCorrLarge, presDistCorr, densDistCorr, windDistCorr,
      densPrevCoef, densRandCoef, presPrevCoef, presDensCoef, presRandCoef, ewPrevCoef, ewDensCoef,
      ewRandCoef, nsPrevCoef, nsDensCoef, nsRandCoef);

  // Get normally distrubuted (Gaussian) random numbers.
  greal densGaussSmall = normalPercentagePoint(randomNumber[0]);  // (was zd)
  greal presGaussSmall = normalPercentagePoint(randomNumber[1]);  // (was zt)

  // Compute new small-scale perturbations for density.
  // Set densFactor to 1 if patchiness is off
  if (!patchy) {
    densFactor = 1.0;
  }
  greal densNewSmall = (densPrevCoef - 1.0) * densPertSmallPrev
                       + densRandCoef * densGaussSmall * densFactor; // (was z2)
  // Assure change in density during severe perturbation is not large compared to normal
  // perturbation standard deviation
  if (patchy && (abs(densNewSmall) > 3.0 * densSDSmall)) {
    densNewSmall = copysign(3.0 * densSDSmall, densNewSmall);
  }
  densPertSmall = densPertSmallPrev + densNewSmall;

  // Compute new small-scale perturbations for pressure
  greal presNewSmall = (presPrevCoef - 1.0) * presPertSmallPrev + presDensCoef * densPertSmall
                       + presRandCoef * presGaussSmall * densFactor; // (was z2)
  // Assure change in pressure during severe perturbations is not large compared to normal
  // perturbation standard deviation
  if (patchy && (abs(presNewSmall) > 3.0 * presSDSmall)) {
    presNewSmall = copysign(3.0 * presSDSmall, presNewSmall);
  }
  presPertSmall = presPertSmallPrev + presNewSmall;
  if (abs(presPertSmall / presSDSmall) > 4.0) {
    presPertSmall = copysign(abs(presPertSmall) - 0.1 * presSDSmall, presPertSmall);
  }

  // Compute new small-scale perturbations for temperature
  // Adjust temperature perturbation assuming small-scale = large-scale p-d correlation
  greal tempScaleSmall = sqrt(abs(square(presSDSmall) + square(densSDSmall)
                                  - 2.0 * pdCorrLarge * presSDSmall * densSDSmall));
  tempPertSmall = (presPertSmall - densPertSmall) * tempSDSmall / tempScaleSmall;

  // Get normally distrubuted (Gaussian) random numbers.
  greal ewGaussSmall = normalPercentagePoint(randomNumber[2]);  // (was zd)
  greal nsGaussSmall = normalPercentagePoint(randomNumber[3]);  // (was zt)

  // Apply scalfact for normal or severe perturbed conditions
  scaleFactor = densFactor;
  // Assure wind perturbation standard deviations are not large compared with severe turbulence
  // sigmas from Figure 2, NASA TM 4168 (1990)
  if (patchy) {
    scaleFactor = min((4.0 + 0.7 * height) / ewWindSDSmall, scaleFactor);
    scaleFactor = min((4.0 + 0.7 * height) / nsWindSDSmall, scaleFactor);
    scaleFactor = max(densFactor, scaleFactor);
  }

  // Compute new small-scale ew wind perturbations
  greal ewNewSmall = (ewPrevCoef - 1.0) * ewPertSmallPrev + ewDensCoef * densPertSmall + ewRandCoef * ewGaussSmall * scaleFactor;
  // Assure change in E-W wind during severe perturbations is not large compared to normal
  // perturbation standard deviation
  if (patchy && (abs(ewNewSmall) > 3.0 * ewWindSDSmall)) {
    ewNewSmall = copysign(3.0 * ewWindSDSmall, ewNewSmall);
  }
  ewWindPertSmall = ewPertSmallPrev + ewNewSmall;

  // Compute new small-scale ns wind perturbations 
  greal nsNewSmall = (nsPrevCoef - 1.0) * nsPertSmallPrev + nsDensCoef * densPertSmall + nsRandCoef * nsGaussSmall * scaleFactor;
  // Assure change in N-S wind during severe perturbations is not large compared to normal
  // perturbation standard deviation
  if (patchy && (abs(nsNewSmall) > 3.0 * nsWindSDSmall)) {
    nsNewSmall = copysign(3.0 * nsWindSDSmall, nsNewSmall);
  }
  nsWindPertSmall = nsPertSmallPrev + nsNewSmall;

  greal minSD = min(400.0, 200.0 + height);
  greal ewMaxSmall = min(minSD, 4.0 * ewWindSDSmall * scaleFactor);
  ewMaxSmall = min(ewMaxSmall, 6.0 * ewWindSDSmall);

  greal nsMaxSmall = min(minSD, 4.0 * nsWindSDSmall * scaleFactor);
  nsMaxSmall = min(nsMaxSmall, 6.0 * nsWindSDSmall);

  if (abs(ewWindPertSmall) > ewMaxSmall)
    ewWindPertSmall = copysign(abs(ewWindPertSmall) - 0.1 * ewMaxSmall, ewWindPertSmall);
  if (abs(nsWindPertSmall) > nsMaxSmall)
    nsWindPertSmall = copysign(abs(nsWindPertSmall) - 0.1 * nsMaxSmall, nsWindPertSmall);

  // Update the large scale perturbations and the boundary layer model.
  updateLargeScalePerturbations(scaleFactor);

  // Calculate sigma-w from boundary layer (BL) model
  updateBoundaryLayerModel();

  if (vertSDPrev <= 0.0) {
    vertSDPrev = nearlyZero;
  }
  if (verticalStandardDeviation <= 0.0) {
    verticalStandardDeviation = nearlyZero;
  }

  // Compute coefficients needed for vertical wind perturbations
  greal vertPrevCoef = windDistCorr * verticalStandardDeviation / vertSDPrev;                // (was aw)
  greal vertRandCoef = verticalStandardDeviation * sqrt(1.0 - windDistCorr * windDistCorr);  // (was bw)
  greal vertGaussSmall = normalPercentagePoint(randomNumber[4]);  // (was zd)

  // Compute new vertical wind (small-scale) perturbation (verticalWindPerturbation)
  greal vertNewSmall = (vertPrevCoef - 1.0) * vertPertPrev
                       + vertRandCoef * vertGaussSmall * densFactor;          // (was z2)
  if (abs(vertNewSmall) > 3.0 * verticalStandardDeviation) {
    vertNewSmall = copysign(3.0 * verticalStandardDeviation, vertNewSmall);
  }
  verticalWindPerturbation = vertPertPrev + vertNewSmall;

  // Total (small + large-scale) perturbations
  densityPerturbation = densPertSmall + densPertLarge;
  temperaturePerturbation = tempPertSmall + tempPertLarge;

  // Use 2nd order gas law for total and large-scale pressure perturbations
  pressurePerturbation = densityPerturbation + temperaturePerturbation
                         + densityPerturbation * temperaturePerturbation;
  presPertLarge = pressurePerturbation - presPertSmall;

  // Avoid perturbation less than -90% of mean value
  densityPerturbation = max(-0.9, densityPerturbation);
  temperaturePerturbation = max(-0.9, temperaturePerturbation);
  pressurePerturbation = max(-0.9, pressurePerturbation);

  // Total (small + large-scale) wind perturbations
  ewWindPerturbation = ewWindPertSmall + ewWindPertLarge;
  nsWindPerturbation = nsWindPertSmall + nsWindPertLarge;

  greal ewMax = min(minSD, 4.0 * ewStandardDeviation * scaleFactor);  // (was umax)
  ewMax = min(ewMax, 6.0 * ewStandardDeviation);

  greal nsMax = min(minSD, 4.0 * nsStandardDeviation * scaleFactor);  // (was vmax)
  nsMax = min(nsMax, 6.0 * nsStandardDeviation);

  if (abs(ewWindPerturbation) > ewMax) {
    ewWindPerturbation = copysign(abs(ewWindPerturbation) - 0.1 * ewMax, ewWindPerturbation);
    ewWindPertSmall = ewWindPerturbation - ewWindPertLarge;
  }
  if (abs(nsWindPerturbation) > nsMax) {
    nsWindPerturbation = copysign(abs(nsWindPerturbation) - 0.1 * nsMax, nsWindPerturbation);
    nsWindPertSmall = nsWindPerturbation - nsWindPertLarge;
  }

  if (densFactor > 1.0) {
    severityLevel = 1;
  }
  else {
    severityLevel = 0;
  }
  updateStatus = PERTS_UPDATED;
}

//! \brief Updates the large and small scale standard deviations.
//!
//! Fractional variances are interpolated by height and latitude.  These variances are
//! applied to the current standard deviations to derive the large and small scale
//! standard deviations.
//!
//! \b Inputs
//! \arg #height      
//! \arg #latitude      
//! \arg #pressureStandardDeviation      
//! \arg #densityStandardDeviation      
//! \arg #temperatureStandardDeviation      
//! \arg #ewStandardDeviation      
//! \arg #nsStandardDeviation      
//!
//! \returns Large and small scale standard deviations are populated.
void EarthAtmosphere::updateStandardDeviations()
{
  //Fractional variances
  greal presVar, densVar, ewVar, nsVar;    // (was plph, dlph, ulph, vlph)
  interpolatePCData(plp, dlp, height, latitude, presVar, densVar);
  interpolatePCData(ulp, vlp, height, latitude, ewVar, nsVar);

  presVar = min(maximumPressureVariance, presVar);
  densVar = min(maximumVariance, densVar);
  greal tempVar = densVar;  // (was tlph)
  nsVar = min(maximumVariance, nsVar);
  ewVar = nsVar;

  // Compute large scale and small scale standard deviations
  // from the total standard deviations
  presSDLarge = sqrt(presVar) * pressureStandardDeviation;
  presSDSmall = sqrt(1.0 - presVar) * pressureStandardDeviation;
  densSDLarge = sqrt(densVar) * densityStandardDeviation;
  densSDSmall = sqrt(1.0 - densVar) * densityStandardDeviation;
  tempSDLarge = sqrt(tempVar) * temperatureStandardDeviation;
  tempSDSmall = sqrt(1.0 - tempVar) * temperatureStandardDeviation;
  ewWindSDLarge = sqrt(ewVar) * ewStandardDeviation;
  ewWindSDSmall = sqrt(1.0 - ewVar) * ewStandardDeviation;
  nsWindSDLarge = sqrt(nsVar) * nsStandardDeviation;
  nsWindSDSmall = sqrt(1.0 - nsVar) * nsStandardDeviation;
}

//! \brief Updates the large-scale pressure-density correlation.
//!
//! This routine computes the pressure-density correlation from the Buell relation.  It then applies an
//! upper limit to the correlation value.  Finally the correlation value is scaled using the current
//! standard deviations and large and small scale pressure and density standard deviations.
//!
//! \b Inputs
//! \arg #height      
//! \arg #pressureStandardDeviation      
//! \arg #densityStandardDeviation      
//! \arg #temperatureStandardDeviation      
//! \arg #presSDLarge      
//! \arg #densSDLarge      
//! \arg #presSDSmall      
//! \arg #densSDSmall      
//!
//! \returns #pdCorrLarge
void EarthAtmosphere::updatePDCorrelation()
{
  // Will use shorter names here for readability of formulas.
  greal& presSD = pressureStandardDeviation;
  greal& densSD = densityStandardDeviation;
  greal& tempSD = temperatureStandardDeviation;

  // Compute pressure-density correlation from the Buell relation
  greal pdCorr = (square(presSD) + square(densSD) - square(tempSD)) / (2.0 * presSD * densSD);  // (was rpd)

  // Get upper limit on dt Correlation (which should be negative)
  greal dtCorrLimit;  // (was rdtlim)
  if (height < 10.0) {
    dtCorrLimit = -0.9 + 0.05 * height;
  }
  else if (height < 15.0) {
    dtCorrLimit = -0.4 - 0.08 * (height - 10.0);
  }
  else if (height < 30.0) {
    dtCorrLimit = -0.8 + 0.04 * (height - 15.0);
  }
  else {
    dtCorrLimit = -0.2;
  }

  // Get upper limit on pd Correlation from dt Correlation limit
  greal pdCorrLimit = (densSD + dtCorrLimit * tempSD) / presSD;  // (was rpdlim)
  if (abs(pdCorrLimit) > 0.99) {
    pdCorrLimit = copysign(0.99, pdCorrLimit);
  }
  if (pdCorr > pdCorrLimit) {
    pdCorr = pdCorrLimit;
  }

  // Assume large=small for perturbation p-d correlation
  pdCorrLarge = pdCorr * presSD * densSD  / (presSDLarge * densSDLarge + presSDSmall * densSDSmall);
  if (abs(pdCorrLarge) > 0.99) {
    pdCorrLarge = copysign(0.99, pdCorrLarge);
  }
}

//! \brief Updates the large-scale perturbations.
//!
//! This routine computes large-scale pertubations using the wave model described in Section 2.7.3 of the 
//! GRAM2010 Users Guide.  Some wave values are computed in this routine for later use in the boundary layer model.
//!
//! \b Inputs
//! \arg #position      
//! \arg #surfaceHeight      
//! \arg #presSDLarge      
//! \arg #densSDLarge      
//! \arg #tempSDLarge      
//! \arg #densSDSmall      
//! \arg #T      
//! \arg #T_v      
//! \arg #a_v      
//! \arg #b_v      
//! \arg #phi_q      
//! \arg #A_Q      
//! \arg #pdCorrLarge      
//! \arg #udCorrLarge      
//! \arg #windCorrelation      
//! \arg #ewWindSDLarge      
//! \arg #nsWindSDLarge      
//! \arg #ewStandardDeviation      
//! \arg #waveRand      
//!
//! \param scaleFactor Perturbation scaling factor. Used in setting wind perturbation limits.
//!
//! \returns The large scale perturbations are computed.
void EarthAtmosphere::updateLargeScalePerturbations(greal scaleFactor)
{
  // Compute vertical wavelenghth (from equation 51)
  greal lambda_z = a_v + b_v * pow(abs(height), 3.0 / 2.0); // (was vll)
  lambda_z = min(200.0, lambda_z);

  // Compute new large-scale perturbations for density, temperature, and pressure
  // Get the wave angle (equation 48, argument to cos())
  greal nTheta = hwn * toRadians(wrapDegrees180(longitude));
  greal mPhi = hwn * toRadians(latitude);
  greal zOverLamda_z = TWO_PI * height / lambda_z;                      // (was phidenz)
  greal tOverT = TWO_PI * elapsedTime / (T * 86400.0);
  greal waveAngle = nTheta + mPhi + zOverLamda_z + tOverT + phi_q;

  greal lambda_zSurface = a_v + b_v * pow(abs(surfaceHeight), 3.0 / 2.0); // (was vlls)
  greal zOverLamda_zSurface = TWO_PI * surfaceHeight / lambda_zSurface;    // (was phidenz)
  ewWaveAngle = nTheta + mPhi + zOverLamda_zSurface + tOverT + phi_q;

  greal tOverT_v = TWO_PI * elapsedTime / (T_v * 86400.0);
  greal waveAngle_v = nTheta + mPhi + zOverLamda_z + tOverT_v + phi_q;
  nsWaveAngle = nTheta + mPhi + zOverLamda_zSurface + tOverT_v + phi_q;

  // This is 1/sigma_cos (in equation 48)
  const static greal overSigma_cos = sqrt(2);

  // Large-scale perturbations (equation 48)
  densPertLarge = densSDLarge * A_Q * overSigma_cos * cos(waveAngle);
  // As explained below equation 51, a phase shift is used in the pressure perturbations
  // so that it properly correlates with density.
  greal presCorrPhase = acos(pdCorrLarge);   // (was dphiu)
  presPertLarge = presSDLarge * A_Q * overSigma_cos * cos(waveAngle + presCorrPhase);

  // Adjust temperature perturbation 
  tempPertLarge = (presPertLarge - densPertLarge) * tempSDLarge
    / sqrt(abs(square(presSDLarge) + square(densSDLarge) - 2 * pdCorrLarge * presSDLarge * densSDLarge));

  // Adjust large-scale pressure using 2nd-order gas law
  if (initializing) {
    presPertLarge = densPertLarge + densPertSmall + tempPertLarge + tempPertSmall
      + (densPertLarge + densPertSmall) * (tempPertLarge + tempPertSmall)
      - presPertSmall;
  }

  // Get phase angles for large-scale u and v components
  ewCorrPhase = acos(udCorrLarge);    // (was dphiu)
  greal windCorrLarge = windCorrelation / (0.02 + 0.98 * square(ewWindSDLarge / ewStandardDeviation));
  if (abs(windCorrLarge) > 1.0) {
    windCorrLarge = copysign(1.0, windCorrLarge);
  }
  nsCorrPhase = ewCorrPhase + copysign(acos(windCorrLarge), waveRand - 0.5);  // (was dphiv)

  // Compute new large-scale wind perturbations (equation 48, discussion under equation 51)
  ewWindPertLarge = ewWindSDLarge * A_Q * overSigma_cos * cos(waveAngle + ewCorrPhase);
  nsWindPertLarge = nsWindSDLarge * A_Q * overSigma_cos * cos(waveAngle_v + nsCorrPhase);

  greal minSD = min(400.0, 200.0 + height);
  greal ewMaxLarge = min(minSD, 4.0 * ewWindSDLarge * scaleFactor);  // (was umax)
  ewMaxLarge = min(ewMaxLarge, 6.0 * ewWindSDLarge);

  greal nsMaxLarge = min(minSD, 4.0 * nsWindSDLarge * scaleFactor);  // (was vmax)
  nsMaxLarge = min(nsMaxLarge, 6.0 * nsWindSDLarge);

  if (abs(ewWindPertLarge) > ewMaxLarge) {
    ewWindPertLarge = copysign(ewMaxLarge, ewWindPertLarge);
  }
  if (abs(nsWindPertLarge) > nsMaxLarge) {
    nsWindPertLarge = copysign(nsMaxLarge, nsWindPertLarge);
  }

  // Ensure large-scale perturbations approach zero at the poles
  if (abs(latitude) >= 85.0) {
    greal poleFactor = cos(PI * (abs(latitude) - 85.0) / 10.0);
    densPertLarge *= poleFactor;
    presPertLarge *= poleFactor;
    tempPertLarge *= poleFactor;
    ewWindPertLarge *= poleFactor;
    nsWindPertLarge *= poleFactor;
  }
}

//! \brief The boundary layer model provides a vertical winds standard deviation.
//!
//! Uses the boundary layer (BL) model to get sigma-w as a function of height, lat, lon,
//! wind speed at 10-m, surface roughness (z0), and atmospheric stability (which depends on time of day).
//!
//! \b Inputs
//! \arg time          
//! \arg #position
//! \arg #A_Q
//! \arg #windCorrelationAtSurface      
//! \arg #ewWindSDAtSurface
//! \arg #nsWindSDAtSurface
//! \arg #verticalWindPerturbationScale
//! \arg #ewWindAtSurface
//! \arg #nsWindAtSurface
//! \arg #temperatureAtSurface
//! \arg #ewWaveAngle 
//! \arg #nsWaveAngle 
//!
//! \param init Set to true if initializing.
//!
//! \retval #verticalStandardDeviation
//! \retval BLmetrics
void EarthAtmosphere::updateBoundaryLayerModel(bool init)
{
  // See page 16 of GRAM2010 manual.
  static const greal srfRuf[14] = { 3.2e-04, 0.6,   0.48,  0.42,  0.0056, 0.45,    0.12,
                                     0.046,   0.015, 0.042, 0.065, 0.45,   1.0e-04, 3.2e-04 };
  static const greal overSigma_cos = sqrt(2.0);
  constexpr greal third = 1.0 / 3.0;

  // Calculate the change in the large-scale wave phase.
  greal ewPhaseDelta, nsPhaseDelta;  // (was dphiu, dphiv)
  if (height <= 10.0 || init) {
    // Get phase angles for large-scale U and V components at surface
    greal udCorrSurface, nsVarSurface;  // (was udlsrf, vlphsrf)
    interpolatePCData(udl, vlp, surfaceHeight, latitude, udCorrSurface, nsVarSurface);
    ewPhaseDelta = acos(udCorrSurface);
    greal windCorrSurface = windCorrelationAtSurface / (0.02 + 0.98 * nsVarSurface);
    if (abs(windCorrSurface) > 1.0) {
      windCorrSurface = copysign(1.0, windCorrSurface);
    }
    nsPhaseDelta = ewPhaseDelta + copysign(acos(windCorrSurface), waveRand - 0.5);
  }
  else {
    ewPhaseDelta = ewCorrPhase;
    nsPhaseDelta = nsCorrPhase;
  }

  // Use sigm-w from atmosdat table if height > 10 km
  if (height > 10.0) {
    verticalStandardDeviation *= verticalWindPerturbationScale;
    return;
  }

  int imin, ihr;
  greal sec;
  if (elapsedTime == 0.0) {
    imin = minute;
    sec = seconds;
    ihr = hour;
  }
  else {
    imin = minute + int(elapsedTime) / 60;
    sec = seconds + elapsedTime;
    sec = int(sec) % 60;
    ihr = hour + imin / 60;
    imin = imin % 60;
  }

  // Default landcode and z0
  landCode = 99;
  surfaceRoughness = 9.99999;

  // Get large-scale wind perturbations at surface for boundary-layer model
  greal ewVarSurface, nsVarSurface;   // (was ulphsrf, vlphsrf)
  interpolatePCData(ulp, vlp, surfaceHeight, latitude, ewVarSurface, nsVarSurface);
  greal ewSD = ewWindSDAtSurface * sqrt(nsVarSurface);  // (was sul2b)
  greal nsSD = nsWindSDAtSurface * sqrt(nsVarSurface);  // (was svl2b)
  greal ewPertSurface = ewSD * A_Q * overSigma_cos * cos(ewWaveAngle + ewPhaseDelta);
  greal nsPertSurface = nsSD * A_Q * overSigma_cos * cos(nsWaveAngle + nsPhaseDelta);

  // Compute surface (10 m) wind speed
  perturbedWindSpeedAtSurface = sqrt(square(ewWindAtSurface + ewPertSurface)
    + square(nsWindAtSurface + nsPertSurface));
  if (perturbedWindSpeedAtSurface < 0.1) {
    perturbedWindSpeedAtSurface = 0.1;
  }

  // Index for first sigma-w above boundary layer
  int indexAboveBoundaryLayer = 1;  // (was iwb)
  if (surfaceHeight > 2.0) {
    indexAboveBoundaryLayer = 2;
  }

  // Get solar declination
  greal solarDec = toRadians(solarDeclination); // (was sda)
  greal xmin = imin + sec / 60.0;

  // Get some sines and cosines
  greal slat = sin(toRadians(latitude));
  greal clat = cos(toRadians(latitude));
  greal sdec = sin(solarDec);
  greal cdec = cos(solarDec);

  // Coriolis parameter (absolute value) with minimum near 20 degrees latitude
  greal coriolis = max(1.4584e-04 * abs(slat), 5.0e-05);

  // Get solar elevation angle (deg.)
  solarElevation = toDegrees(asin(sdec * slat + cdec * clat * cos(toRadians(solarHourAngle))));

  // Solar elevation at midnight (deg.)
  elevationAtMidnight = toDegrees(asin(sdec * slat - cdec * clat));
  if (elevationAtMidnight == 0.0) {
    elevationAtMidnight = 1.0e-03;
  }

  // Solar elevation at midday (deg.)
  elevationAtNoon = toDegrees(asin(sdec * slat + cdec * clat));
  if (elevationAtNoon == 0.0) {
    elevationAtNoon = 1.0e-03;
  }

  // Height factor for unstable B.L. during early daytime
  unstableBLFactor = 1.0;
  if ((elevationAtMidnight < 0.0) && (solarElevation > 0.0) && (solarElevation <= elevationAtNoon)
    && (solarHourAngle < 0.0)) {
    unstableBLFactor = 0.3 + 0.7 * solarElevation / elevationAtNoon;
  }

  // Get (modified) Net Radiation Index from simplified versions
  // of Table 4-7 and 4-8 of Justus (1978)
  // Case when dark all 24 hours
  if (sdec * slat + cdec * clat < 0.0) {
    netRadiationIndex = -3.5;
  }
  else {
    if ((solarHourAngle < 0.0) & (solarElevation > elevationAtMidnight) & (solarElevation < 0.0)) {
      // Special interpolated nri from midnight to dawn
      greal nrimn = 0.5 + 6.1154e-02 * elevationAtMidnight - 2.6390e-06 * pow(elevationAtMidnight, 3);
      netRadiationIndex = nrimn + (3.5 + nrimn) * (solarElevation / elevationAtMidnight - 1.0);
    }
    else {
      // Usual nri versus elevation (-90 to +90)
      netRadiationIndex = 0.5 + 6.1154e-02 * solarElevation - 2.6390e-06 * pow(solarElevation, 3);
    }
  }
  netRadiationIndex = clamp(netRadiationIndex, -3.5, 4.5);

  // Use input surface roughness value or compute from model
  if (userSurfaceRoughness > 0.0) {
    surfaceRoughness = userSurfaceRoughness;
    landCode = 99;
  }
  else {
    if (userSurfaceRoughness == 0.0) {
      landCode = 0;
    }
    else {
      // Get surface type from landcode array
      int i = int(wrapDegrees180(longitude) + 180.0);
      if (i > 359) {
        i = 0;
      }
      int j = clampSize(int(latitude + 90.0), 180);
      landCode = landcd[i][j];
    }

    // Get surface type from landcode value
    if (landCode == 0) {
      // use iterative solution if water type
      surfaceRoughness = srfRuf[0];

      for (int i = 0; i < 4; i++) {
        surfaceRoughness = 2.612e-04
          * square(perturbedWindSpeedAtSurface / log(10.0 / surfaceRoughness));
      }
    }
    else {
      surfaceRoughness = srfRuf[landCode];
    }

    // Model to increase surface roughness z0 when in high-altitude (mountainous) terrain
    if (surfaceHeight >= 1.5) {
      surfaceRoughness += surfaceHeight - 1.5;
    }

    // Bound surface roughness values.
    surfaceRoughness = clamp(surfaceRoughness, 1.0e-5, 3.0);
  }

  // Get stability category S and inverse Monin-Obukhov length 1/L = inverseLength, 
  // versus wind speed u, and Net Radiation Index (nri).
  // Wind speed factor for stability category [for simplified version
  // of Table 4-7 of Justus (1978)]
  greal windSpeedFactor;  // (was f)
  if (perturbedWindSpeedAtSurface <= 6.0) {
    windSpeedFactor = (1.0 - perturbedWindSpeedAtSurface / 7.5);
  }
  else {
    windSpeedFactor = 0.2 * exp(12.0 - 2.0 * perturbedWindSpeedAtSurface);
  }

  // Stability category for wind speed and net radiation index
  //[simplified version of Justus (1978) Pages 57-58, neglecting
  // cloud cover and ceiling information, unavailable in GRAM]
  stability = 4.229 - netRadiationIndex * windSpeedFactor;
  stability = clamp(stability, 0.5, 7.5);

  // Inverse of Monin-Obukhov length (inverseLength=1/L, in m**-1) vs stability
  // and surface roughness from Figure 4-9 of Justus (1978)
  inverseLength = 0.25 * (-0.2161 + 0.0511 * stability) * log10(10.0 / surfaceRoughness);

  greal psi;
  if (inverseLength >= 0.0) {
    // Stable case psi(10/L) [psi(z/L) at z = 10 m]
    psi = -50.0 * inverseLength;
  }
  else {
    // Simplified version of Paulson unstable relation for psi(10/L),
    // from equation (5) of Hsu et al. (1999)
    psi = 1.0496 * pow((-10.0 * inverseLength), 0.4591);
  }

  // Friction velocity ustar (m/s) versus surface roughness and psi(10/L)
  greal denom = max(log(10.0 / surfaceRoughness) - psi, 0.1);
  frictionVelocity = 0.4 * perturbedWindSpeedAtSurface / denom;
  frictionVelocity = clamp(frictionVelocity, 0.04 / denom, 4.0);

  // Calculate hbl = height of boundary layer
  // Method for stable-to-neutral, Section 2.1, Sugiyama and Nasstrom,
  // 1999, with neutral BL height changed from hN = 0.2*ustar/coriol to
  // hN = ustar*(4.0*B/(coriol*BVfsq))**third.  As expressed in Seibert
  // (2000), with time dependent term d(hbl)/dt replaced by hbl*coriol/2.
  constexpr greal a = 0.4;
  constexpr greal b = 20.0;
  constexpr greal gamma = 3.3e-03;

  BVFrequencySquare = (9.8 / temperatureAtSurface) * gamma;
  neutralBoundaryLayerDepth = frictionVelocity * pow(((4.0 * b) / (coriolis * BVFrequencySquare)), third);
  if (inverseLength >= 0.0) {
    boundaryLayerDepth = 2.0 * neutralBoundaryLayerDepth / (1.0 + sqrt(1.0 + 4.0 * neutralBoundaryLayerDepth * inverseLength));
  }
  else {
    greal ha = neutralBoundaryLayerDepth;
    greal hb = 0.0;
    greal r = -(1.0 + 2.0 * a) * inverseLength / (2.0 * b * 0.4);

    for (int m = 0; m < 7; m++) {
      hb = neutralBoundaryLayerDepth * pow((1.0 + r * ha), third);
      if (hb < hb) {
        hb = neutralBoundaryLayerDepth;
      }
      if (abs(hb - ha) <= 5.0) {
        break;
      }
      ha = hb;
    }

    // BL height hbl in meters
    boundaryLayerDepth = hb;
  }

  // Factor unstableBLFactor for unstable conditions, dawn to mid-morning,
  // to lessen sudden increase in B.L. height at dawn
  if (boundaryLayerDepth > 3000.0)
    boundaryLayerDepth = 3000.0;
  if (inverseLength < -0.0001)
    boundaryLayerDepth = unstableBLFactor * boundaryLayerDepth;
  if (boundaryLayerDepth < 200.0)
    boundaryLayerDepth = 200.0;

  // Height (above MSL) for top of BL (in km)
  greal topOfBoundaryLayer = surfaceHeight + (boundaryLayerDepth / 1000.0);

  if (height < 5.0 * indexAboveBoundaryLayer) {
    // Use height (AGL) at top of BL (in meters) if current height above BL
    metersAboveSurface = boundaryLayerDepth;
    // Use actual height (AGL in meters) if current height within BL
    if (height < topOfBoundaryLayer) {
      metersAboveSurface = (height - surfaceHeight) * 1000.0;
    }
    if (metersAboveSurface < 0.0) {
      metersAboveSurface = 0.0;
    }
    // Revised BL sigma-w from z0 and ustar
    // Ratio sigma-w/ustar
    greal sigmaRatioLimit;  // (was siglim)
    if (inverseLength >= 0.0) {
      // Eq. 1.33, Kaimel and Finnigan (1994); larger range from
      // Pahlow et al. (2001)
      sigmaRatio = 1.25 * (1.0 + 0.2 * metersAboveSurface * inverseLength);
      sigmaRatioLimit = 3.0;
    }
    else {
      // Limit sigma-w/ustar with convective velocity wstar: Eq. 1.48,
      // 1.51 at z/zi = 0.5, and Fig 1.10 of Kaimal and Finnigan (1994)
      sigmaRatioLimit = 0.62 * pow((-boundaryLayerDepth * inverseLength / 0.4), third);
      if (sigmaRatioLimit < 1.25) {
        sigmaRatioLimit = 1.25;
      }
      // Equation (2), Page 161, Panofsky and Dutton (1984)
      sigmaRatio = 1.25 * pow((1.0 - 3.0 * metersAboveSurface * inverseLength), third);
    }
    if (sigmaRatio > sigmaRatioLimit) {
      sigmaRatio = sigmaRatioLimit;
    }

    // sigma-w values from ustar and ratios
    sigmaW = frictionVelocity * sigmaRatio;
    sigmaW = clamp(sigmaW, 0.1, 5.0);
    // If height above top of BL, interpolate between sigma-w at top
    // of BL and next highest level
    if (height > topOfBoundaryLayer) {
      greal fact = (height - topOfBoundaryLayer) / (5.0 * indexAboveBoundaryLayer - topOfBoundaryLayer);
      verticalStandardDeviation = sigmaW + (wr[indexAboveBoundaryLayer] - sigmaW) * fact;
    }
    else {
      // fact = 0.0;
      verticalStandardDeviation = sigmaW;
    }

    // Get solar days (local solar time from start)
    greal days = day + (60.0 * (ihr - hour) + xmin - minute) / 1440.0 - solarTime / 24.0;
    solarDays = 1 + int(days) + solarTime / 24.0;
  }

  verticalStandardDeviation *= verticalWindPerturbationScale;
}

//! \brief Computes small-scale perturbation model coefficients.
//!
//! Computes perturbation model coefficients (*Coef) at previous 
//! positions (*Prev) and current position from the standard deviations and correlations.
//! \param densSDPrev    Density standard deviation at the previous position.
//! \param densSD        Density standard deviation at the current position.
//! \param presSDPrev    Pressure standard deviation at the previous position.
//! \param presSD        Pressure standard deviation at the current position.
//! \param tempSDPrev    Temperature standard deviation at the previous position.
//! \param tempSD        Temperature standard deviation at the current position.
//! \param ewSDPrev      East/west wind standard deviation at the previous position.
//! \param ewSD          East/west wind standard deviation at the current position.
//! \param nsSDPrev      North/south wind standard deviation at the previous position.
//! \param nsSD          North/south wind standard deviation at the current position.
//! \param udCorrPrev    Correlation between EW wind and density at the previous position.
//! \param udCorr        Correlation between EW wind and density at the current position.
//! \param vdCorrPrev    Correlation between NS wind and density at the previous position.
//! \param vdCorr        Correlation between NS wind and density at the current position.
//! \param pdCorr        Correlation between pressure and density at the current position.
//! \param presDistCorr  Correlation between pressure and distance at the current position.
//! \param densDistCorr  Correlation between density and distance at the current position.
//! \param windDistCorr  Correlation between winds and distance at the current position.
//! \param[out] densPrevCoef  Density coefficient for the previous density perturbations.
//! \param[out] densRandCoef  Density coefficient for the current density randomization.
//! \param[out] presPrevCoef  Pressure coefficient for the previous pressure perturbations.
//! \param[out] presDensCoef  Pressure coefficient for the current density perturbation.
//! \param[out] presRandCoef  Pressure coefficient for the current pressure randomization.
//! \param[out] ewPrevCoef    East/west wind coefficient for the previous east/west wind perturbations.
//! \param[out] ewDensCoef    East/west wind coefficient for the current density perturbation.
//! \param[out] ewRandCoef    East/west wind coefficient for the current east/west wind randomization.
//! \param[out] nsPrevCoef    North/south wind coefficient for the previous north/south wind perturbations.
//! \param[out] nsDensCoef    North/south wind coefficient for the current density perturbation.
//! \param[out] nsRandCoef    North/south wind coefficient for the current north/south wind randomization.
void EarthAtmosphere::getPerturbationCoefficients(greal densSDPrev, greal densSD, greal presSDPrev, greal presSD,
  greal tempSDPrev, greal tempSD, greal ewSDPrev, greal ewSD, greal nsSDPrev, greal nsSD,
  greal udCorrPrev, greal udCorr, greal vdCorrPrev, greal vdCorr, 
  greal pdCorr,  greal presDistCorr, greal densDistCorr, greal windDistCorr, 
  greal &densPrevCoef, greal &densRandCoef, greal &presPrevCoef, greal &presDensCoef, greal &presRandCoef,
  greal &ewPrevCoef, greal &ewDensCoef, greal &ewRandCoef, greal &nsPrevCoef, greal &nsDensCoef, greal &nsRandCoef)
{
  constexpr greal nearlyOne = 0.99999;

  // Default values avoid division by zero
  if (densSDPrev <= 0.0) densSDPrev = nearlyZero;
  if (tempSDPrev <= 0.0) tempSDPrev = nearlyZero;
  if (densSD <= 0.0) densSD = nearlyZero;
  if (tempSD <= 0.0) tempSD = nearlyZero;
  if (presDistCorr <= 0.0) presDistCorr = nearlyZero;
  if (densDistCorr <= 0.0) densDistCorr = nearlyZero;
  if (windDistCorr <= 0.0) windDistCorr = nearlyZero;
  if (abs(udCorrPrev) <= 0.0) udCorrPrev = nearlyZero;
  if (abs(vdCorrPrev) <= 0.0) vdCorrPrev = nearlyZero;
  if (abs(ewSDPrev) <= 0.0) ewSDPrev = nearlyZero;
  if (abs(nsSDPrev) <= 0.0) nsSDPrev = nearlyZero;
  if (abs(udCorrPrev) >= nearlyOne) udCorrPrev = nearlyOne * udCorrPrev / abs(udCorrPrev);
  if (abs(vdCorrPrev) >= nearlyOne) vdCorrPrev = nearlyOne * vdCorrPrev / abs(vdCorrPrev);
  if (abs(udCorr) <= 0.0) udCorr = nearlyZero;
  if (abs(vdCorr) <= 0.0) vdCorr = nearlyZero;
  if (abs(ewSD) <= 0.0) ewSD = nearlyZero;
  if (abs(nsSD) <= 0.0) nsSD = nearlyZero;
  if (abs(udCorr) >= nearlyOne) udCorr = nearlyOne * udCorr / abs(udCorr);
  if (abs(vdCorr) >= nearlyOne) vdCorr = nearlyOne * vdCorr / abs(vdCorr);
  if (abs(presDistCorr) >= nearlyOne) presDistCorr = nearlyOne * presDistCorr / abs(presDistCorr);
  if (abs(densDistCorr) >= nearlyOne) densDistCorr = nearlyOne * densDistCorr / abs(densDistCorr);
  if (abs(windDistCorr) >= nearlyOne) windDistCorr = nearlyOne * windDistCorr / abs(windDistCorr);
  if (abs(pdCorr) <= 0.0) pdCorr = nearlyZero;
  if (abs(pdCorr) >= nearlyOne) pdCorr = nearlyOne * pdCorr / abs(pdCorr);

  // Compute perturbation model coeffcients
  densPrevCoef = densDistCorr * densSD / densSDPrev;
  densRandCoef = densSD * sqrt(1 - square(densDistCorr));

  presPrevCoef = (presSD / presSDPrev) * (presDistCorr - (densDistCorr * square(pdCorr)))
                 / (1.0 - (square(densDistCorr) * square(pdCorr)));
  presDensCoef = pdCorr * (presSD - (presPrevCoef * densDistCorr * presSDPrev)) / densSD;
  presRandCoef = square(presSD)
                 - square(presPrevCoef) * square(presSDPrev)
                 - square(presDensCoef) * square(densSD)
                 - 2 * presPrevCoef * presDensCoef * densDistCorr * pdCorr * presSDPrev * densSD;
  presRandCoef = sqrt(max(0.0, presRandCoef));

  ewPrevCoef = (ewSD / ewSDPrev) * (windDistCorr - densDistCorr * udCorr * udCorrPrev)
               / (1.0 - square(densDistCorr) * square(udCorrPrev));
  ewDensCoef = (udCorr * ewSD - ewPrevCoef * densDistCorr * udCorrPrev * ewSDPrev) / densSD;
  ewRandCoef = square(ewSD)
               - square(ewPrevCoef) * square(ewSDPrev)
               - square(ewDensCoef) * square(densSD)
               - 2.0 * ewPrevCoef * ewDensCoef * densDistCorr * udCorrPrev * densSD * ewSDPrev;
  ewRandCoef = sqrt(max(0.0, ewRandCoef));

  nsPrevCoef = (nsSD / nsSDPrev) * (windDistCorr - densDistCorr * vdCorr * vdCorrPrev)
               / (1.0 - square(densDistCorr) * square(vdCorrPrev));
  nsDensCoef = (vdCorr * nsSD - nsPrevCoef * densDistCorr * vdCorrPrev * nsSDPrev) / densSD;
  nsRandCoef = square(nsSD)
               - square(nsPrevCoef) * square(nsSDPrev)
               - square(nsDensCoef) * square(densSD)
               - 2.0 * nsPrevCoef * nsDensCoef * densDistCorr * vdCorrPrev * densSD * nsSDPrev;
  nsRandCoef = sqrt(max(0.0, nsRandCoef));
}

//! \brief Interpolates over two arrays by height and latitude.
//!
//! This method performs linear interpolation on two arrays at the specified height (km) and geocentric latitude (degrees). 
//! 
//! \param xarray  A two-dimensional array.
//! \param yarray  A two-dimensional array.
//! \param hgt     The interpolation height in km.
//! \param phi     The interpolation latitude in degrees.
//! \param[out] xint  A greal.
//! \param[out] yint  A greal.
//!
//! \retval xint The interpolated x value.
//! \retval yint The interpolated y value.
void EarthAtmosphere::interpolatePCData(greal xarray[pcHeightSize][pcLatSize], greal yarray[pcHeightSize][pcLatSize],
                             greal hgt, greal phi, greal &xint, greal &yint)
{
//  assert(hgt >= 0.0);

  // i - lower height index
  int i;
  if (hgt < 95.0) {
    i = int(hgt) / 5;
  }
  else {
    i = 18 + (int(hgt) - 80) / 20;
  }
  i = min(24, i);

  // Upper height index
  int ip = i + 1;
  ip = min(24, ip);

  assert(phi <= 90.0);
  assert(phi >= -90.0);

  // Lower latitude index
  int j = int(phi + 110.0) / 20 - 1;
  j = min(8, j);

  // Upper latitude index
  int jp = j + 1;
  jp = min(9, jp);

  // phi1 - lower latitude index for ur and vr arrays
  greal phi1 = (-110.0) + 20.0 * (j + 1);
  // phi2 - upper latitude index for ur and vr arrays
  greal phi2 = (-110.0) + 20.0 * (jp + 1);

  // lower height for ur and vr array values
  greal z1;
  if (i > 18) {
    z1 = 20.0 * (i - 14);
  }
  else {
    z1 = 5.0 * i;
  }

  // Upper height for ur and vr array values
  greal z2;
  if (ip > 18) {
    z2 = 20.0 * (ip - 14);
  }
  else {
    z2 = 5.0 * ip;
  }

  // The interpolation fraction for height has a special case.
  greal heightFraction;
  if (i == ip) {
    // We are above the height range.  Use the upper bound.
    heightFraction = 1.0;
  }
  else {
    heightFraction = (hgt - z1) / (z2 - z1);
  }

  // Interpolate on geocentric latitude and height
  Interpolator interp2d(heightFraction, (phi - phi1) / (phi2 - phi1));
  xint = interp2d.linear(xarray[i][j], xarray[i][jp], xarray[ip][j], xarray[ip][jp]);
  yint = interp2d.linear(yarray[i][j], yarray[i][jp], yarray[ip][j], yarray[ip][jp]);
}

//! \brief Interpolates variable-scale perturbation data arrays based on the current height.
//!
//! This function returns the linear interpolation of the input array at the current height.
//! Arrays are variable-scale perturbation data of size rsHeightSize.
//!
//! \param data    An RS data array indexed by height.
//! \param hgt     The interpolation height in km.
//!
//! \returns The interpolated value.
greal EarthAtmosphere::interpolateRSData(greal data[rsHeightSize], greal hgt)
{
  // lower height index
  int lowerHeightIndex = int(hgt) / 5;
  if (hgt >= 125.0) {
    lowerHeightIndex = 24 + (int(hgt) - 120) / 20;
  }
  lowerHeightIndex = min(28, lowerHeightIndex);

  // Upper height index
  int upperHeightIndex = lowerHeightIndex + 1;
  upperHeightIndex = min(28, upperHeightIndex);

  greal z1;
  if (lowerHeightIndex > 24) {
    // Lower height for wr array values
    z1 = 20 * (lowerHeightIndex - 18);
  }
  else {
    z1 = 5.0 * lowerHeightIndex;
  }

  greal z2;
  if (upperHeightIndex > 24) {
    // Upper height for wr array values
    z2 = 20.0 * (upperHeightIndex - 18);
  }
  else {
    z2 = 5.0 * upperHeightIndex;
  }

  // Set the interpolation fraction based on height
  Interpolator interpHeight;
  if (lowerHeightIndex == upperHeightIndex) {
    // use the upper bound
    interpHeight.setFraction(1.0);
  }
  else {
    interpHeight.makeFraction(z1, z2, hgt);
  }

  //Interpolate on height
  return interpHeight.linear(data[lowerHeightIndex], data[upperHeightIndex]);
}


//! \brief Compute the correlation factor.
//!
//! Correlation function for relative separation distance x (r/L).
//!
//! \param x  A relative separation distance.
//!
//! returns The correlation factor.
greal EarthAtmosphere::getCorrelation(greal x)
{
  greal correl = 0.0;

  // Exponential function, with defaults for small and large x
  if (abs(x) < 1.0e-7) {
    correl = 0.9999999;
  }
  else if (abs(x) > 23.0) {
    correl = exp(-23.0);
  }
  else {
    correl = exp(-abs(x));
  }
  return correl;
}

//! \brief Tail probability from normal distribution.
//!
//! Tail probability from normal distribution (Integral x to infinity of Gaussian),
//! where x is normalized deviate (deviation from mean divided by standard deviation
//! about mean).  Note \f$\text{Ptail} = (1/2)*(1 - \text{erf}(x/\sqrt{2})\f$, where \p erf is the error function,
//! or \f$\text{Ptail} = (1/2)\text{erfcomp}(x/\sqrt{2})\f$, where \p erfcomp is the complementary error function.
//! This implementation of erfcomp is based on Chebyshev fitting.  Coefficient values
//! are from "Numerical Recipes in Fortran", 2nd Edition, page 214.
//!
//! \param x A normalized deviate.
//!
//! returns The probability tail.
greal EarthAtmosphere::getProbabilityTail(greal x)
{
  greal y = x / sqrt(2.0);
  greal q = abs(y);
  greal t = 1.0 / (1.0 + q / 2.0);

  greal erfcomp = t * exp((-q * q) - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
    t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
      t * (1.48851587 + t * (-0.82215223 + t * (0.17087277))))))))));

  if (x < 0.0) {
    erfcomp = 2.0 - erfcomp;
  }

  greal ptail = erfcomp / 2.0;
  return ptail;
}

//! \brief Fair atmosphere state with RRA data.
//!
//! This method uses the RRA class to determine if the current position lies within an RRA site limits.
//! If it does, then the current atmosphere state is faired with the RRA values.
//!
//! \b Inputs:
//! \arg #atmos
//!
//! \param init True if values are being initialized.  False by default.
//! 
//! \returns atmos The updated atmosphere state.
void EarthAtmosphere::updateRRA(bool init)
{
  RRA& rra = *rraPtr;
  rra.setPosition(position);

  // If conditions are met use RRA data
  if (height <= 70.0 && rra.isActive()) {
    // Save the incoming values.
    greal presMean = pressure;
    greal densMean = density;
    greal tempMean = temperature;
    greal ewMean = ewWind;
    greal nsMean = nsWind;
    greal presSD = pressureStandardDeviation;
    greal densSD = densityStandardDeviation;
    greal tempSD = temperatureStandardDeviation;
    greal ewSD = ewStandardDeviation;
    greal nsSD = nsStandardDeviation;

    greal densitySDRRA = 0.0;     // (was sdhr)
    greal pressureSDRRA = 0.0;    // (was sphr)
    greal temperatureSDRRA = 0.0; // (was sthr)
    greal ewWindSDRRA = 0.0;      // (was suhr)
    greal nsWindSDRRA = 0.0;      // (was svhr)
    if (init) {
      pressureStandardDeviation = sqrt(pressureStandardDeviation);
      densityStandardDeviation = sqrt(densityStandardDeviation);
      temperatureStandardDeviation = sqrt(temperatureStandardDeviation);

      // Call RRA class to get weighted standard deviations
      // This call updates SDs as a weighted blend of current SDs and RRA SDs.
      rra.update();
      surfaceHeight = rra.getPosition().surfaceHeight;

      if (rra.isWeighted()) {
        pressureSDRRA = square(pressureStandardDeviation);
        densitySDRRA = square(densityStandardDeviation);
        temperatureSDRRA = square(temperatureStandardDeviation);
        ewWindSDRRA = ewStandardDeviation;
        nsWindSDRRA = nsStandardDeviation;
      }
      else {
        // No RRA weighting, use the incoming SDs
        pressureSDRRA = presSD;
        densitySDRRA = densSD;
        temperatureSDRRA = tempSD;
        ewWindSDRRA = ewSD;
        nsWindSDRRA = nsSD;
      }
    }
    else {
      // Call RRA class to get weighted standard deviations
      // This call updates SDs as a weighted blend of current SDs and RRA SDs.
      rra.update();
      surfaceHeight = rra.getPosition().surfaceHeight;

      pressureSDRRA = pressureStandardDeviation;
      densitySDRRA = densityStandardDeviation;
      temperatureSDRRA = temperatureStandardDeviation;
      ewWindSDRRA = ewStandardDeviation;
      nsWindSDRRA = nsStandardDeviation;
    }

    greal lowerHeight, upperHeight;
    rra.getFairingHeights(lowerHeight, upperHeight);

    // Taller height fairing region for 2013 data.
    if (height >= lowerHeight && height <= upperHeight 
      && (rra.getRRAYear() == RRA_2013 || rra.getRRAYear() == RRA_2019)) {

      //RRA Fair
      greal wdum = 0.0;
      fair(2, 
        lowerHeight, pressureSDRRA, densitySDRRA, temperatureSDRRA, ewWindSDRRA, nsWindSDRRA, wdum,
        upperHeight, presSD, densSD, tempSD, ewSD, nsSD, wdum, 
        height, pressureStandardDeviation, densityStandardDeviation, temperatureStandardDeviation, 
        ewStandardDeviation, nsStandardDeviation, wdum);

      greal presRRA = pressure;
      greal densRRA = density;
      greal tempRRA = temperature;
      greal ewRRA = ewWind;
      greal nsRRA = nsWind;

      //RRA Fair
      fair(1, 
        lowerHeight, presRRA, densRRA, tempRRA, ewRRA, nsRRA, wdum,
        upperHeight, presMean, densMean, tempMean, ewMean, nsMean, wdum,
        height, pressure, density, temperature, ewWind, nsWind, wdum);
    }
    else {
      pressureStandardDeviation = pressureSDRRA;
      densityStandardDeviation = densitySDRRA;
      temperatureStandardDeviation = temperatureSDRRA;
      ewStandardDeviation = ewWindSDRRA;
      nsStandardDeviation = nsWindSDRRA;
    }
  }
}

//! \brief Compute the US-76 standard atmosphere components.
//!
//! US-76 Standard Atmosphere values of the temperature (K), pressure (N/m^2),
//! and density (kg/m^3) at height (km).  Uses vertical interpolation between
//! values in data table s (zs = height, tms = temperature, wms = molecular 
//! weight, ps = pressure).
//!
//! \b Inputs:
//! \arg #height
//!
//! \param[out] temp A greal.
//! \param[out] pres A greal.
//! \param[out] dens A greal.
//!
//! \retval temp The standared temperature (K).
//! \retval pres The standared pressure (Pa).
//! \retval dens The standared density (kg/m^3).
void EarthAtmosphere::getUS76StandardAtmosphere(greal &temp, greal &pres, greal &dens)
{
  // DATA
  static const greal zs[49] =
  {
    0.0,   11.019, 20.063, 32.162, 47.35, 51.413, 71.802, 86.0,  91.0,  94.0,
    97.0,  100.0,  103.0,  106.0,  108.0, 110.0,  112.0,  115.0, 120.0, 125.0,
    130.0, 135.0,  140.0,  145.0,  150.0, 155.0,  160.0,  165.0, 170.0, 180.0,
    190.0, 210.0,  230.0,  265.0,  300.0, 350.0,  400.0,  450.0, 500.0, 550.0,
    600.0, 650.0,  700.0,  750.0,  800.0, 850.0,  900.0,  950.0, 1000.0
  };

  static const greal tms[49] =
  {
    288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95, 186.87, 187.74,
    190.40, 195.08, 202.23, 212.89, 223.29, 240.0,  264.0,  300.00, 360.00, 417.23,
    469.27, 516.59, 559.63, 598.78, 634.39, 666.80, 696.29, 723.13, 747.57, 790.07,
    825.31, 878.84, 915.78, 955.20, 976.01, 990.06, 995.83, 998.22, 999.24, 999.67,
    999.85, 999.93, 999.97, 999.99, 999.99, 1000.0, 1000.0, 1000.0, 1000.0
  };

  static const greal wms[49] =
  {
    28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9522, 28.889, 28.783,
    28.620,  28.395,  28.104,  27.765,  27.521,  27.268,  27.020,  26.680,  26.205, 25.803,
    25.436,  25.087,  24.749,  24.422,  24.103,  23.792,  23.488,  23.192,  22.902, 22.342,
    21.809,  20.825,  19.952,  18.688,  17.726,  16.735,  15.984,  15.247,  14.330, 13.092,
    11.505,  9.718,   7.998,   6.579,   5.543,   4.849,   4.404,   4.122,   3.940
  };

  static const greal ps[49] =
  {
    1013.25,    226.32,     54.7487,    8.68014,    1.10905,    0.66938,    0.039564,
    3.7338e-3,  1.5381e-3,  9.0560e-4,  5.3571e-4,  3.2011e-4,  1.9742e-4,  1.2454e-4,
    9.3188e-5,  7.1042e-5,  5.5547e-5,  4.0096e-5,  2.5382e-5,  1.7354e-5,  1.25054e-5,
    9.3568e-6,  7.2028e-6,  5.6691e-6,  4.5422e-6,  3.6930e-6,  3.0395e-6,  2.5278e-6,
    2.1210e-6,  1.5271e-6,  1.1266e-6,  6.4756e-7,  3.9276e-7,  1.7874e-7,  8.7704e-8,
    3.4498e-8,  1.4518e-8,  6.4468e-9,  3.0236e-9,  1.5137e-9,  8.2130e-10, 4.8865e-10,
    3.1908e-10, 2.2599e-10, 1.7036e-10, 1.3415e-10, 1.0873e-10, 8.9816e-11, 7.5138e-11
  };

  // Set values to 0 if z < 0 (or z > 1000 km)
  if (height < 0.0 || height > 1000.0) {
    temp = 0.0;
    pres = 0.0;
    dens = 0.0;
    return;
  }

  // Constants used in the 76 standard.
  constexpr greal ro = 6356.766;  // Average earth radius in km
  constexpr greal go = 9.80665;   // Average earth gravity in m/s^2
  constexpr greal wmo = 28.9644;  // Molar mass of air g/mol
  constexpr greal rs = 8314.472;  // universal gas constant in J/(kmol.K)

  // Do height interpolation for temperature t (K)
  if (height < 86.0_km) {
    int i;
    for (i = 0; i < 7; ++i) {
      if ((zs[i] <= height) && (height < zs[i + 1])) {
        break;
      }
    }

    // Geopotential heights
    greal zl = ro * zs[i] / (ro + zs[i]);
    greal zu = ro * zs[i + 1] / (ro + zs[i + 1]);
    greal ht = (ro * height) / (ro + height);
    // in meters
    greal zlm = zl * 1000.0;
    greal hm = ht * 1000.0;

    // Compute temperature gradient (deg/km)
    greal g = (tms[i + 1] - tms[i]) / (zu - zl);
    // in deg/m
    greal gm = g * 0.001;

    if (abs(g) <= 0.0001) {
      pres = ps[i] * exp(-(go * wmo * (hm - zlm)) / (rs * tms[i])) * 100.0;
    }
    else {
      pres = ps[i] * (pow(tms[i] / (tms[i] + g * (ht - zl)), (go * wmo) / (rs * gm))) * 100.0;
    }

    temp = tms[i] + g * (ht - zl);

    // Compute density from perfect gas law and molecular weight
    dens = (wmo * pres) / (rs * temp);

    return;
  }

  // For heights greater than 86 km.
  int i;
  for (i = 7; i < 48; ++i) {
    if (zs[i] <= height && height < zs[i + 1]) {
      break;
    }
  }
  if (i > 47) {
    i = 47;
  }

  if (i == 7) {
    temp = tms[8];
  }
  else if (i >= 15 && i < 18) {
    temp = 240.0 + 12.0 * (height - 110.0);
  }
  else if (i < 15) {
    temp = 263.1905 - 76.3232 * sqrt(1.0 - pow((height - 91.0) / 19.9429, 2.0));
  }
  else {
    greal xi = (height - 120.0) * (ro + 120.0) / (ro + height);
    temp = (1000.0 - 640.0 * exp(-0.01875 * xi));
  }

  int j = min(i, 46);

  // Do height interpolation for molecular weight, wm
  greal z0 = zs[j];
  greal z1 = zs[j + 1];
  greal z2 = zs[j + 2];
  greal fraction1 = (height - z1) * (height - z2) / ((z0 - z1) * (z0 - z2));
  greal fraction2 = (height - z0) * (height - z2) / ((z1 - z0) * (z1 - z2));
  greal fraction3 = (height - z0) * (height - z1) / ((z2 - z0) * (z2 - z1));

  greal wma = wms[j] * fraction1 + wms[j + 1] * fraction2 + wms[j + 2] * fraction3;

  greal alp0 = log(ps[j]);
  greal alp1 = log(ps[j + 1]);
  greal alp2 = log(ps[j + 2]);
  greal alpa = alp0 * fraction1 + alp1 * fraction2 + alp2 * fraction3;

  greal alpb = alpa;
  greal wmb = wma;

  if (i != 7 && i != 47) {
    j = j - 1;

    z0 = zs[j];
    z1 = zs[j + 1];
    z2 = zs[j + 2];
    greal fraction1 = (height - z1) * (height - z2) / ((z0 - z1) * (z0 - z2));
    greal fraction2 = (height - z0) * (height - z2) / ((z1 - z0) * (z1 - z2));
    greal fraction3 = (height - z0) * (height - z1) / ((z2 - z0) * (z2 - z1));

    alp0 = log(ps[j]);
    alp1 = log(ps[j + 1]);
    alp2 = log(ps[j + 2]);

    alpb = alp0 * fraction1 + alp1 * fraction2 + alp2 * fraction3;

    wmb = wms[j] * fraction1 + wms[j + 1] * fraction2 + wms[j + 2] * fraction3;
  }

  // Convert pressure in mb to N/m^2
  pres = 100.0 * exp((alpa + alpb) / 2.0);
  greal wm = (wma + wmb) / 2.0;

  // Compute density from perfect gas law and molecular weight
  dens = (wm * pres) / (rs * temp);
}

//! \brief Generate correlation coefficients between two positions.
//!
//! Data corresponding to two positions, such as a previous and a current position, must
//! be correlated based on time, horizontal position, and vertical position.  This routine
//! generates the coefficients that would be used as in this equation:
//! \f$ \text{correlatedData} = \text{baseCoefficient} * \text{baseData} + \text{targetCoefficient} * \text{targetData} \f$
//!
//! \param base    The base position.
//! \param target  The target position.
//! \param[out] baseCoefficient   A greal.
//! \param[out] targetCoefficient A greal.
//!
//! \retval baseCoefficient    The coefficient to apply to the base data.
//! \retval targetCoefficient  The coefficient to apply to the target data.
void EarthAtmosphere::getCorrelationCoefficients(const Position& base, const Position& target, 
                                                 greal& baseCoefficient, greal& targetCoefficient)
{
  greal timeScaleSmall = 86400.0 * 0.735 * pow(abs(target.height), 0.116);
  timeScaleSmall = max(10800.0, timeScaleSmall);
  greal xScaleSmall = interpolateRSData(EarthAtmosphere::xscale, target.height);
  greal zScaleSmall = interpolateRSData(EarthAtmosphere::zscale, target.height);

  greal dx = toRadians(target.totalRadius * getArcAngle(base.latitude, base.longitude, target.latitude, target.longitude));
  greal dz = abs(base.height - target.height);
  greal dtime = abs(base.elapsedTime - target.elapsedTime);

  baseCoefficient = getCorrelation(dx / xScaleSmall + dz / zScaleSmall + dtime / timeScaleSmall);
  targetCoefficient = sqrt(1.0 - square(baseCoefficient));
}


} // namespace
