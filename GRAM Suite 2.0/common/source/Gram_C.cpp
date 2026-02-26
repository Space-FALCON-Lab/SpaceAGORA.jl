//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include "Gram_C.h"
#include "SpiceLoader.h"
#include "NamelistReader.h"

namespace GRAM {

static std::map<GRAM::GRAM_BODY, GRAM::CreateFunction> createMap;
static std::map<GRAM::GRAM_BODY, GRAM::CopyFunction> copyMap;
static std::map<GRAM::GRAM_BODY, GRAM::UpdateFunction> updateMap;
static std::map<GRAM::GRAM_BODY, GRAM::DeleteFunction> deleteMap;
std::string errorMessage;

void registerBody_C(GRAM_BODY gramBody, CreateFunction createFunction, CopyFunction copyFunction, UpdateFunction updateFunction, DeleteFunction deleteFunction)
{
  createMap[gramBody] = createFunction;
  copyMap[gramBody] = copyFunction;
  updateMap[gramBody] = updateFunction;
  deleteMap[gramBody] = deleteFunction;
}

//! \brief Sets the start time (epoch) of the simulation.
//!
//! \param atmos An atmosphere handle.
//! \param time A GramTime_C structure containing date/time components.
void setStartTime_C(PerturbedAtmosphere* atmos, const GramTime_C* time)
{
  if (atmos != nullptr && time != nullptr) {
    GramTime gtime;
    GRAM_TIME_SCALE timeScale;
    switch (time->scale) {
    case 0:
      timeScale = TDT;
      break;
    case 1:
      timeScale = UTC;
      break;
    case 2:
      timeScale = TDB;
      break;
    default:
      timeScale = UTC;
      break;
    }
    GRAM_TIME_FRAME timeFrame;
    switch (time->frame) {
    case 0:
      timeFrame = PET;
      break;
    case 1:
      timeFrame = ERT;
      break;
    default:
      timeFrame = ERT;
      break;
    }
    gtime.setStartTime(time->year, time->month, time->day, time->hour, time->minute, time->seconds, timeScale, timeFrame);
    atmos->setStartTime(gtime);
  }
}

//! \copydoc GRAM::Atmosphere::setPosition()
//! \param atmos An atmosphere handle.
void setPosition_C(PerturbedAtmosphere* atmos, const Position_C* pos)
{
  if (atmos != nullptr && pos != nullptr) {
    Position position;
    position.height = pos->height;
    position.latitude = pos->latitude;
    position.longitude = pos->longitude;
    position.elapsedTime = pos->elapsedTime;
    position.isPlanetoCentric = (pos->isPlanetoCentric != 0);
    position.latitudeRadius = pos->latitudeRadius;
    position.totalRadius = pos->totalRadius;
    position.surfaceHeight = pos->surfaceHeight;
    position.gravity = pos->gravity;
    atmos->setPosition(position);
  }
}

//! \copydoc GRAM::PerturbedAtmosphere::setDelta()
//! \param atmos An atmosphere handle.
void setDelta_C(PerturbedAtmosphere* atmos, const Position_C* delta)
{
  if (atmos != nullptr && delta != nullptr) {
    Position pos;
    pos.height = delta->height;
    pos.latitude = delta->latitude;
    pos.longitude = delta->longitude;
    pos.elapsedTime = delta->elapsedTime;
    pos.isPlanetoCentric = true;
    atmos->setDelta(pos);
  }
}

//! \brief Set the pertubation random number seed.
//!
//! \param atmos An atmosphere handle.
//! \param seed An integer between 0 and 30,000.
void setSeed_C(PerturbedAtmosphere* atmos, int seed)
{
  if (atmos != nullptr && seed != 0) {
    atmos->setSeed(seed);
  }
}

//! \brief Set the minimum relative step size for perturbations.
//!
//! \param atmos An atmosphere handle.
//! \param minRelativeStepSize Between 0 and 1.
void setMinRelativeStepSize_C(PerturbedAtmosphere* atmos, greal minRelativeStepSize)
{
  if (atmos != nullptr) {
    atmos->setMinRelativeStepSize(minRelativeStepSize);
  }
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param densityScale Between 0 and 2.
//! \param ewWindScale Between 0 and 2.
//! \param nsWindScale Between 0 and 2.
//! \param verticalWindScale Between 0 and 2.
void setPerturbationScales_C(PerturbedAtmosphere* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale)
{
  if (atmos != nullptr) {
    atmos->setDensityPerturbationScale(densityScale);
    atmos->setEWWindPerturbationScale(ewWindScale);
    atmos->setNSWindPerturbationScale(nsWindScale);
    atmos->setVerticalWindPerturbationScale(verticalWindScale);
  }
}

//! \brief Adds an auxiliary atmosphere to the list.
//!
//! Use this function to append to the list of auxiliary atmospheres.
//! Fairing between the current atmospheric state and auxiliary atmospheres
//! is performed in the order in which the atmospheres are added to the list.
//! \param atmos An atmosphere handle.
//! \param fileName The name of the file containing the auxiliary atmosphere data.
//! \param innerRadius The inner radius of the fairing region.
//! \param outerRadius The outer radius of the fairing region.
//! \param isEastLongitudePositive If input data file uses east longitude positive convention set to 1. Set to 0 for west longitude positive.
void addAuxiliaryAtmosphere_C(PerturbedAtmosphere* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  if (atmos != nullptr) {
    atmos->addAuxiliaryAtmosphere(fileName, innerRadius, outerRadius, (isEastLongitudePositive != 0));
  }
}

//! \brief Set auxiliary atmosphere values for the next update.
//!
//! \param atmos An atmosphere handle.
//! \param dens Density \units{kg/m^3}.
//! \param pres Pressure \units{Pa}.
//! \param temp Temperuture \units{K}.
//! \param ew East/West windswith east postive \units{m/s}.
//! \param ns North/South winds with north positive \units{m/s}.
void setAuxiliaryValues_C(PerturbedAtmosphere* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  if (atmos != nullptr) {
    atmos->getAuxiliaryAtmosphere(0).setValues(dens, pres, temp, ew, ns);
  }
}

//! \brief Controls the computation of perturbations.
//!
//! \param atmos An atmosphere handle.
//! \param action Set to 0 for no perturbations, 1 for perturbations.
void setPerturbationAction_C(PerturbedAtmosphere* atmos, int action)
{
  if (atmos != nullptr) {
    PerturbationAction pertAction = (action == 0) ? DO_NOT_UPDATE_PERTS : UPDATE_PERTS;
    atmos->setPerturbationAction(pertAction);
  }
}

//! \copydoc GRAM::PerturbedAtmosphere::setEphemerisState()
//! \param atmos An atmosphere handle.
void setEphemerisState_C(PerturbedAtmosphere* atmos, const EphemerisState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    EphemerisState ephem;
    ephem.solarTime = state->solarTime;
    ephem.longitudeSun = state->longitudeSun;
    ephem.subsolarLatitude = state->subsolarLatitude;
    ephem.subsolarLongitude = state->subsolarLongitude;
    ephem.orbitalRadius = state->orbitalRadius;
    ephem.oneWayLightTime = state->oneWayLightTime;
    ephem.solarZenithAngle = state->solarZenithAngle;
    ephem.secondsPerSol = state->secondsPerSol;
    atmos->setEphemerisState(ephem);
  }
}

//! \copydoc GRAM::Ephemeris::setFastModeOn()
//! \param atmos An atmosphere handle.
void setEphemerisFastModeOn_C(PerturbedAtmosphere* atmos, int flag)
{
  if (atmos != NULL) {
    atmos->setEphemerisFastModeOn(flag != 0);
  }
}

//! \copydoc GRAM::Ephemeris::setSubsolarUpdateTime()
//! \param atmos An atmosphere handle.
void setSubsolarUpdateTime_C(PerturbedAtmosphere* atmos, greal utime)
{
  if (atmos != NULL) {
    atmos->setSubsolarUpdateTime(utime);
  }
}


//! \copybrief GRAM::PerturbedAtmosphere::getPosition()
//! \param atmos An atmosphere handle.
//! \param[out] position A non-null pointer.
//! \retval position A structure containing the current position. It also contains
//! a few computed parameters which are valid after update.
void getPosition_C(PerturbedAtmosphere* atmos, Position_C* position)
{
  if (atmos != nullptr && position != nullptr) {
    const Position& pos = atmos->getPosition();
    position->height = pos.height;
    position->latitude = pos.latitude;
    position->longitude = pos.longitude;
    position->elapsedTime = pos.elapsedTime;
    position->totalRadius = pos.totalRadius;
    position->latitudeRadius = pos.latitudeRadius;
    position->surfaceHeight = pos.surfaceHeight;
    position->gravity = pos.gravity;
    position->isPlanetoCentric = pos.isPlanetoCentric ? 1 : 0;
  }
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The dynamics state computed in the last update.
void getDynamicsState_C(PerturbedAtmosphere* atmos, DynamicsState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    state->temperature = astate.temperature;
    state->pressure = astate.pressure;
    state->density = astate.density;
    state->pressureScaleHeight = astate.pressureScaleHeight;
    state->densityScaleHeight = astate.densityScaleHeight;
    state->speedOfSound = astate.speedOfSound;
    state->referenceDensity = astate.referenceDensity;
    state->referenceTemperature = astate.referenceTemperature;
    state->referencePressure = astate.referencePressure;
    state->sigmaLevel = astate.sigmaLevel;
    state->pressureAtSurface = astate.pressureAtSurface;
    state->pressureAltitude = astate.pressureAltitude;
  }
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The density state computed in the last update.
void getDensityState_C(PerturbedAtmosphere* atmos, DensityState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    state->density = astate.density;
    state->lowDensity = astate.lowDensity;
    state->highDensity = astate.highDensity;
    state->densityPerturbation = astate.densityPerturbation;
    state->perturbedDensity = astate.perturbedDensity;
    state->densityStandardDeviation = astate.densityStandardDeviation;
    state->perturbedSpeedOfSound = astate.perturbedSpeedOfSound;
    state->relativeStepSize = astate.relativeStepSize;
    state->referenceDensity = astate.referenceDensity;
    state->densityDeviation = astate.densityDeviation;
    state->lowDensityDeviation = astate.lowDensityDeviation;
    state->highDensityDeviation = astate.highDensityDeviation;
    state->perturbedDensityDeviation = astate.perturbedDensityDeviation;
    switch (astate.updateStatus) {
    case STEP_UPDATED:
      state->updateStatus = 1;
      break;
    case PERTS_UPDATED:
      state->updateStatus = 2;
      break;
    default:
      state->updateStatus = 0;
      break;
    }
  }
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The winds state computed in the last update.
void getWindsState_C(PerturbedAtmosphere* atmos, WindsState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    state->ewWind = astate.ewWind;
    state->nsWind = astate.nsWind;
    state->verticalWind = astate.verticalWind;
    state->ewWindPerturbation = astate.ewWindPerturbation;
    state->nsWindPerturbation = astate.nsWindPerturbation;
    state->verticalWindPerturbation = astate.verticalWindPerturbation;
    state->perturbedEWWind = astate.perturbedEWWind;
    state->perturbedNSWind = astate.perturbedNSWind;
    state->perturbedVerticalWind = astate.perturbedVerticalWind;
    state->ewStandardDeviation = astate.ewStandardDeviation;
    state->nsStandardDeviation = astate.nsStandardDeviation;
    state->verticalStandardDeviation = astate.verticalStandardDeviation;
  }
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The gases state computed in the last update.
void getGasesState_C(PerturbedAtmosphere* atmos, GasesState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    state->totalNumberDensity = astate.totalNumberDensity;
    state->averageMolecularWeight = astate.averageMolecularWeight;
    state->compressibilityFactor = astate.compressibilityFactor;
    state->specificGasConstant = astate.specificGasConstant;
    state->specificHeatRatio = astate.specificHeatRatio;
  }
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param type The type of gas being requested (see GasType).
//! \param[out] gas A non-NULL pointer.
//! \retval gas The gas state computed in the last update.
void getConstituentGas_C(PerturbedAtmosphere* atmos, ConstituentGas_C* gas, GasType type)
{
  if (atmos != nullptr && gas != nullptr) {
    ConstituentGas g = atmos->getConstituentGas(type);
    gas->moleFraction = g.moleFraction;
    gas->massFraction = g.massFraction;
    gas->numberDensity = g.numberDensity;
    gas->averageMolecularWeight = g.averageMolecularWeight;
    gas->specificHeatCapacity = g.specificHeatCapacity;
  }
}

//! \brief Get the current ephemeris values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The ephemeris state computed in the last update.
void getEphemerisState_C(PerturbedAtmosphere* atmos, EphemerisState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const EphemerisState& ephem = atmos->getEphemerisState();
    state->solarTime = ephem.solarTime;
    state->longitudeSun = ephem.longitudeSun;
    state->subsolarLatitude = ephem.subsolarLatitude;
    state->subsolarLongitude = ephem.subsolarLongitude;
    state->orbitalRadius = ephem.orbitalRadius;
    state->oneWayLightTime = ephem.oneWayLightTime;
    state->solarZenithAngle = ephem.solarZenithAngle;
    state->secondsPerSol = ephem.secondsPerSol;
  }
}

//! \brief Get the current perturbation values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The perturbation state computed in the last update.
void getPerturbationState_C(PerturbedAtmosphere* atmos, PerturbationState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const PerturbationState& pert = atmos->getPerturbationState();
    state->densityRandomNumber = pert.densityRandomNumber;
    state->ewWindRandomNumber = pert.ewWindRandomNumber;
    state->nsWindRandomNumber = pert.nsWindRandomNumber;
    state->verticalWindRandomNumber = pert.verticalWindRandomNumber;
  }
}

//! \brief Gets the start time (epoch) of the simulation.
//!
//! \param atmos An atmosphere handle.
//! \param[out] time A non-NULL pointer.
//! \retval time A GramTime_C structure containing date/time components.
void getStartTime_C(PerturbedAtmosphere* atmos, GramTime_C* time)
{
  if (atmos != nullptr && time != nullptr) {
    const GramTime& stime = atmos->getStartTime();
    GRAM_TIME_SCALE timeScale;
    switch (time->scale) {
    case 0:
      timeScale = TDT;
      break;
    case 1:
      timeScale = UTC;
      break;
    case 2:
      timeScale = TDB;
      break;
    default:
      timeScale = UTC;
      break;
    }
    GRAM_TIME_FRAME timeFrame;
    switch (time->frame) {
    case 0:
      timeFrame = PET;
      break;
    case 1:
      timeFrame = ERT;
      break;
    default:
      timeFrame = ERT;
      break;
    }
    stime.getStartTime(timeScale, timeFrame, time->year, time->month,
      time->day, time->hour, time->minute, time->seconds);
  }
}


}

///////////////////////////////////////////////////////////////////////////////////////////
using namespace GRAM;

//! \defgroup C_GRAM The C (Generic) Interface for the GRAM Suite 
//! @{
//! \brief A comlete listing of the C interface functions for all GRAMs.
//!
//! The C interface for Mars is declared in the following:
//! \code #include "Gram_C.h" \endcode 
//! An example using the generic C interface can be found  "here" (TODO).
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. 

//! \brief Control the output to the console.
//!
//! By default, console output is muted when using the GRAM Suite API.
//! Use this function to enable/disable console output.
//! \param flag Set to true to enable console output.
void enableConsoleOutput_C(int flag)
{
  enableConsoleOutput(flag != 0);
}

//! \brief Check if there is a namelist Spice data path.
//!
//! This routine will return the SPICE data path located in a namelist
//! file called "spice.txt" in the current folder.  If such a file is
//! not found, then the default SPICE path will be returned.  The 
//! character buffer must be large enough to contain the path.
//! \param bufferSize The size of the spicePath character buffer.
//! \param[out] spicePath A non-NULL pointer.
//! \retval spicePath A string containing the path to the spice data
void tryGetSpicePath_C(char* spicePath, int bufferSize)
{
  NamelistReader reader;
  std::string path;
  if (bufferSize > 1) {
    reader.tryGetSpicePath(path);
    path.copy(spicePath, (size_t)bufferSize - 1);
    spicePath[std::min((int)path.length(), bufferSize - 1)] = '\0';
  }
}

//! \brief Check if there is a namelist Spice data path and planetary data path.
//!
//! This routine will return the SPICE data path located in a namelist
//! file called "spice.txt" in the current folder.  If such a file is
//! not found, then the default SPICE path will be returned.  The 
//! character buffer must be large enough to contain the path.
//! The namelist file will also be scanned for a dataPath pointing to the
//! folder containing planetary input data.
//! \param bufferSize The size of the spicePath and dataPath character buffers.
//! \param[out] spicePath A non-NULL pointer.
//! \param[out] dataPath A non-NULL pointer.
//! \retval spicePath A string containing the path to the spice data
void tryGetDataPaths_C(char* spicePath, char* dataPath, int bufferSize)
{
  NamelistReader reader;
  std::string sPath, dPath;
  if (bufferSize > 1) {
    reader.tryGetSpicePath(sPath, dPath);
    sPath.copy(spicePath, (size_t)bufferSize - 1);
    spicePath[std::min((int)sPath.length(), bufferSize - 1)] = '\0';
    dPath.copy(dataPath, (size_t)bufferSize - 1);
    dataPath[std::min((int)dPath.length(), bufferSize - 1)] = '\0';
  }
}

//! \brief Override the default SPICE leapsecond file name.
//!
//! Critical to timekeeping is the SPICE leapsecond file.  The location of this
//! file within the SPICE data folder can be specified with this function.
//! The relative path/name is typically of the form "/lsk/naif0012.tls".
//! \param lsk  A string containing the LSK file path within the SPICE data folder
void setSpiceLsk_C(const char* lsk) {
  try {
    SpiceLoader::setSpiceLsk(lsk);
  }
  catch (const std::string& msg) {
    errorMessage = msg; 
  }
}

//! \brief Override the default SPICE planetary constants file name.
//!
//! The NAIF SPICE library use a file containing a number of planetary constants.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The relative path/name is typically of the form "/pck/pck00011.tpc".
//! \param pck  A string containing the PCK file path within the SPICE data folder
void setSpicePck_C(const char* pck) {
  try {
    SpiceLoader::setSpicePck(pck);
  }
  catch (const std::string& msg) {
    errorMessage = msg; 
  }
}

//! \brief Override the default SPICE planetary kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The Earth and Venus kernel path/name is typically of the form "/spk/planets/de440.bsp".
//! The other planets kernel path/name is typically of the form "/spk/satellites/mar097.bsp" where
//! the first three letters name the planet and the digits provide for versioning.
//! \param body  An enumeration designating desired body
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_C(enum GRAM_BODY_C body, const char* bsp) {
  GRAM_BODY gramBody = NO_BODY;
  switch (body) {
  case VENUS_C:
    gramBody = VENUS;
    break;
  case EARTH_C:
    gramBody = EARTH;
    break;
  case MARS_C:
    gramBody = MARS;
    break;
  case JUPITER_C:
    gramBody = JUPITER;
    break;
  case SATURN_C:
    gramBody = SATURN;
    break;
  case URANUS_C:
    gramBody = URANUS;
    break;
  case NEPTUNE_C:
    gramBody = NEPTUNE;
    break;
  case TITAN_C:
    gramBody = TITAN;
    break;
  default:
    gramBody = NO_BODY;
    break;
  }
  try {
    SpiceLoader::setSpiceKernel(gramBody, bsp);
  }
  catch (const std::string& msg) {
    errorMessage = msg; 
  }
}

//! \brief Initialize the Spice data path.
//!
//! Ensure that the GRAM software can locate the Spice data.
//! The function should be called before creating an atmosphere.
//! \param spicePath A string containing the path to the spice data
void initialize_C(const char* spicePath)
{
  try {
    SpiceLoader::setSpiceDataPath(spicePath);
  }
  catch (const std::string& msg) {
    errorMessage = msg; 
  }
}

//! \brief Loads a Spice library.
//!
//! This is a convenience function for accessing the FURNSH in Spice.
//! The function will store the file name and not re-FURNSH should the same
//! file name appear as an argument.  The file name should be specified as
//! a path relative to the initialized Spice data path.
//! \param fileName A string containing the relative path to the spice data
void loadSpiceFile_C(const char* fileName)
{
  try {
    SpiceLoader::loadFile(std::string(fileName));
  }
  catch (const std::string& msg) {
    errorMessage = msg;
  }
}

//! \brief Sets the start time (epoch) of the simulation.
//!
//! \param atmos An atmosphere handle.
//! \param time A GramTime_C structure containing date/time components.
void setStartTime_C(void* atmos, const GramTime_C* time)
{
  GRAM::setStartTime_C(static_cast<PerturbedAtmosphere*>(atmos), time);
}

//! \copydoc GRAM::Atmosphere::setPosition()
//! \param atmos An atmosphere handle.
void setPosition_C(void* atmos, const Position_C* pos)
{
  GRAM::setPosition_C(static_cast<PerturbedAtmosphere*>(atmos), pos);
}

//! \copydoc GRAM::PerturbedAtmosphere::setDelta()
//! \param atmos An atmosphere handle.
void setDelta_C(void* atmos, const Position_C* delta)
{
  GRAM::setDelta_C(static_cast<PerturbedAtmosphere*>(atmos), delta);
}

//! \brief Set the pertubation random number seed.
//!
//! \param atmos An atmosphere handle.
//! \param seed An integer between 0 and 30,000.
void setSeed_C(void* atmos, int seed)
{
  GRAM::setSeed_C(static_cast<PerturbedAtmosphere*>(atmos), seed);
}

//! \brief Set the minimum relative step size for perturbations.
//!
//! \param atmos An atmosphere handle.
//! \param minRelativeStepSize Between 0 and 1.
void setMinRelativeStepSize_C(void* atmos, greal minRelativeStepSize)
{
  GRAM::setMinRelativeStepSize_C(static_cast<PerturbedAtmosphere*>(atmos), minRelativeStepSize);
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param densityScale Between 0 and 2.
//! \param ewWindScale Between 0 and 2.
//! \param nsWindScale Between 0 and 2.
//! \param verticalWindScale Between 0 and 2.
void setPerturbationScales_C(void* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale)
{
  GRAM::setPerturbationScales_C(static_cast<PerturbedAtmosphere*>(atmos), densityScale, ewWindScale, nsWindScale, verticalWindScale);
}

//! \brief Adds an auxiliary atmosphere to the list.
//!
//! Use this function to append to the list of auxiliary atmospheres.
//! Fairing between the current atmospheric state and auxiliary atmospheres
//! is performed in the order in which the atmospheres are added to the list.
//! \param atmos An atmosphere handle.
//! \param fileName The name of the file containing the auxiliary atmosphere data.
//! \param innerRadius The inner radius of the fairing region.
//! \param outerRadius The outer radius of the fairing region.
//! \param isEastLongitudePositive If input data file uses east longitude positive convention set to 1. Set to 0 for west longitude positive.
void addAuxiliaryAtmosphere_C(void* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  GRAM::addAuxiliaryAtmosphere_C(static_cast<PerturbedAtmosphere*>(atmos), fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \brief Set auxiliary atmosphere values for the next update.
//!
//! \param atmos An atmosphere handle.
//! \param dens Density \units{kg/m^3}.
//! \param pres Pressure \units{Pa}.
//! \param temp Temperuture \units{K}.
//! \param ew East/West windswith east postive \units{m/s}.
//! \param ns North/South winds with north positive \units{m/s}.
void setAuxiliaryValues_C(void* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  GRAM::setAuxiliaryValues_C(static_cast<PerturbedAtmosphere*>(atmos), dens, pres, temp, ew, ns);
}

//! \brief Controls the computation of perturbations.
//!
//! \param atmos An atmosphere handle.
//! \param action Set to 0 for no perturbations, 1 for perturbations.
void setPerturbationAction_C(void* atmos, int action)
{
  GRAM::setPerturbationAction_C(static_cast<PerturbedAtmosphere*>(atmos), action);
}

//! \copydoc GRAM::PerturbedAtmosphere::setEphemerisState()
//! \param atmos An atmosphere handle.
void setEphemerisState_C(void* atmos, const EphemerisState_C* state)
{
  GRAM::setEphemerisState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \copydoc GRAM::Ephemeris::setFastModeOn()
//! \param atmos An atmosphere handle.
void setEphemerisFastModeOn_C(void* atmos, int flag)
{
  GRAM::setEphemerisFastModeOn_C(static_cast<PerturbedAtmosphere*>(atmos), flag);
}

//! \copydoc GRAM::Ephemeris::setSubsolarUpdateTime()
//! \param atmos An atmosphere handle.
void setSubsolarUpdateTime_C(void* atmos, greal utime)
{
  GRAM::setSubsolarUpdateTime_C(static_cast<PerturbedAtmosphere*>(atmos), utime);
}


//! \brief Create a handle to an atmosphere model.
//!
//! This function creates and initializes an atmosphere model.
//! The desired model is specified via the \ref GRAM_BODY_C enumeration.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \param body The desired model.
//! \returns A handle to the new model.
void* createAtmosphere_C(enum GRAM_BODY_C body)
{
  GRAM_BODY gramBody = NO_BODY;
  switch (body) {
  case VENUS_C:
    gramBody = VENUS;
    break;
  case EARTH_C:
    gramBody = EARTH;
    break;
  case MARS_C:
    gramBody = MARS;
    break;
  case JUPITER_C:
    gramBody = JUPITER;
    break;
  case SATURN_C:
    gramBody = SATURN;
    break;
  case URANUS_C:
    gramBody = URANUS;
    break;
  case NEPTUNE_C:
    gramBody = NEPTUNE;
    break;
  case TITAN_C:
    gramBody = TITAN;
    break;
  default:
    gramBody = NO_BODY;
    break;
  }
  auto iter = createMap.find(gramBody);
  if (iter != createMap.end()) {
    CreateFunction createBody = iter->second;
    return createBody();
  }
  else {
    return nullptr;
  }
}

//! \brief Creates a copy of an atmosphere model.
//!
//! During a simulation it may be necessary to duplicate an existing atmosphere 
//! model, say with a vehicle separation.  This function will duplicate the model
//! along with its current state.
//! \param atmos An atmosphere handle.
//! \returns A handle to the copy.
void* copyAtmosphere_C(void* atmos)
{
  GRAM_BODY gramBody = static_cast<PerturbedAtmosphere*>(atmos)->getBody();
  auto iter = copyMap.find(gramBody);
  if (iter != copyMap.end()) {
    CopyFunction copyBody = iter->second;
    return copyBody(atmos);
  }
  else {
    return nullptr;
  }

}

//! \brief Releases memory associated with an atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos An atmosphere handle.
void deleteAtmosphere_C(void* atmos)
{
  GRAM_BODY gramBody = static_cast<PerturbedAtmosphere*>(atmos)->getBody();
  auto iter = deleteMap.find(gramBody);
  if (iter != deleteMap.end()) {
    DeleteFunction deleteBody = iter->second;
    deleteBody(atmos);
  }
}

//! \brief Performs the atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated.
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos An atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_C(void* atmos)
{
  GRAM_BODY gramBody = static_cast<PerturbedAtmosphere*>(atmos)->getBody();
  auto iter = updateMap.find(gramBody);
  if (iter != updateMap.end()) {
    UpdateFunction updateBody = iter->second;
    return updateBody(atmos);
  }
  errorMessage = "GRAM body not found. Check atmosphere pointer.";
  return 1;
}

//! \brief Retrieve the current error message.
//!
//! If an update call results in a non-zero error code, then a call to 
//! this function will capture the error message.  The message is copied
//! into the supplied character buffer.  This buffer must be large enough
//! to contain the error message.  Attempts are made to keep messages
//! brief. However, a buffer size of 1000 is suggested to avoid overruns.
//! \param[out] message A non-NULL pointer.
//! \param bufferSize The size of the pre-allocated message buffer.
//! \retval message The error message will be copied into this buffer.
size_t getErrorMessage_C(char *message, size_t bufferSize)
{
  if (message != nullptr) {
    size_t len = errorMessage.copy(message, bufferSize - 2);
    message[len] = '\0';
    return len;
  }
  return 0;
}

//! \copybrief GRAM::PerturbedAtmosphere::getPosition()
//! \param atmos An atmosphere handle.
//! \param[out] position A non-NULL pointer.
//! \retval position A structure containing the current position. It also contains
//! a few computed parameters which are valid after update.
void getPosition_C(void* atmos, Position_C* position)
{
  GRAM::getPosition_C(static_cast<PerturbedAtmosphere*>(atmos), position);
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The dynamics state computed in the last update.
void getDynamicsState_C(void* atmos, DynamicsState_C* state)
{
  GRAM::getDynamicsState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The density state computed in the last update.
void getDensityState_C(void* atmos, DensityState_C* state)
{
  GRAM::getDensityState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The winds state computed in the last update.
void getWindsState_C(void* atmos, WindsState_C* state)
{
  GRAM::getWindsState_C(static_cast<PerturbedAtmosphere*>(atmos), state);

}

//! \brief Get the current atmosphere values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The gases state computed in the last update.
void getGasesState_C(void* atmos, GasesState_C* state)
{
  GRAM::getGasesState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \brief Get the current constituent gas values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param gasTypeC The type of gas being requested (see GasType_C).
//! \param[out] gas A non-NULL pointer.
//! \retval gas The gas state computed in the last update.
void getConstituentGas_C(void* atmos, enum GasType_C gasTypeC, ConstituentGas_C* gas)
{
  GasType gasType;
  switch (gasTypeC) {
  case ARGON_C:
    gasType = ARGON;
    break;
  case CARBON_DIOXIDE_C:
    gasType = CARBON_DIOXIDE;
    break;
  case CARBON_MONOXIDE_C:
    gasType = CARBON_MONOXIDE;
    break;
  case HELIUM_C:
    gasType = HELIUM;
    break;
  case HYDROGEN_C:
    gasType = HYDROGEN;
    break;
  case DIHYDROGEN_C:
    gasType = DIHYDROGEN;
    break;
  case NITROGEN_C:
    gasType = NITROGEN;
    break;
  case DINITROGEN_C:
    gasType = DINITROGEN;
    break;
  case OXYGEN_C:
    gasType = OXYGEN;
    break;
  case DIOXYGEN_C:
    gasType = DIOXYGEN;
    break;
  case METHANE_C:
    gasType = METHANE;
    break;
  case OZONE_C:
    gasType = OZONE;
    break;
  case NITROUS_OXIDE_C:
    gasType = NITROUS_OXIDE;
    break;
  case WATER_C:
    gasType = WATER;
    break;
  default:
    gasType = ARGON;
    break;
  }
  GRAM::getConstituentGas_C(static_cast<PerturbedAtmosphere*>(atmos), gas, gasType);
}

//! \brief Get the current ephemeris values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The ephemeris state computed in the last update.
void getEphemerisState_C(void* atmos, EphemerisState_C* state)
{
  GRAM::getEphemerisState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \brief Get the current perturbation values after on update.
//!
//! \param atmos An atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The perturbation state computed in the last update.
void getPerturbationState_C(void* atmos, PerturbationState_C* state)
{
  GRAM::getPerturbationState_C(static_cast<PerturbedAtmosphere*>(atmos), state);
}

//! \brief Gets the start time (epoch) of the simulation.
//!
//! \param atmos An atmosphere handle.
//! \param[out] time A non-NULL pointer.
//! \retval time A GramTime_C structure containing date/time components.
void getStartTime_C(void* atmos, GramTime_C* time)
{
  GRAM::getStartTime_C(static_cast<PerturbedAtmosphere*>(atmos), time);
}

//! @}
