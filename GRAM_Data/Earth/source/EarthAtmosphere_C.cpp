//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cstring>
#include "gram.h"
#include "EarthAtmosphere_C.h"
#include "SpiceLoader.h"
#include "GramTime.h"

using namespace GRAM;

//! \defgroup C_Earth The C Interface for Earth 
//! @{
//! The C interface for Earth is declared in the following:
//! \code #include "EarthAtmosphere_C.h" \endcode 
//! An example using the C interface can be found \ref Earth/examples/Earth_C.c "here".
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
//! C and C++ code below.
//! \code
//! // In C
//! EarthAtmosphere_C* earth = createAtmosphere_E();
//! setPosition_E(earth, &pos);
//! update_E(earth); 
//! 
//! // In C++
//! EarthAtmosphere earth;
//! earth.setPosition(pos);
//! earth.update(); 
//! \endcode
//! \example Earth/examples/Earth_C.c

//! \copydoc GRAM::tryGetSpicePath_C
void tryGetSpicePath_E(char* spicePath, int bufferSize)
{
  tryGetSpicePath_C(spicePath, bufferSize);
}

//! \copydoc GRAM::tryGetDataPaths_C
void tryGetDataPaths_E(char* spicePath, char* dataPath, int bufferSize)
{
  tryGetDataPaths_C(spicePath, dataPath, bufferSize);
}

//! \copydoc GRAM::setSpiceLsk_C
void setSpiceLsk_E(const char* lsk)
{
  setSpiceLsk_C(lsk);
}

//! \copydoc GRAM::setSpicePck_C
void setSpicePck_E(const char* pck)
{
  setSpicePck_C(pck);
}

//! \brief Override the default SPICE Earth kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The Earth and Venus kernel path/name is typically of the form "/spk/planets/de440.bsp".
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_E(const char* bsp)
{
  setSpiceKernel_C(EARTH_C, bsp);
}

//! \copydoc GRAM::initialize_C
void initialize_E(const char* spicePath)
{
  initialize_C(spicePath);
  registerBody_C(EARTH,
    reinterpret_cast<CreateFunction>(createAtmosphere_E),
    reinterpret_cast<CopyFunction>(copyAtmosphere_E),
    reinterpret_cast<UpdateFunction>(update_E),
    reinterpret_cast<DeleteFunction>(deleteAtmosphere_E));
}

//! \copydoc GRAM::loadSpiceFile_C
void loadSpiceFile_E(const char* fileName)
{
  loadSpiceFile_C(fileName);
}

//! \brief Create a handle to a Earth atmosphere model.
//!
//! This function creates and initializes a Earth atmosphere model.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \returns A handle to the new model.
EarthAtmosphere_C* createAtmosphere_E(const char* dataPath)
{
  // Load default settings (defined in EarthInputParameters)
  EarthInputParameters params;
  // Except let's not override the SpiceDataPath
  params.spicePath = SpiceLoader::getSpiceDataPath();
  params.dataPath = dataPath;
  EarthAtmosphere* atmos = new EarthAtmosphere();
  atmos->setInputParameters(params);
  return atmos;
}

//! \brief Copies an Earth atmosphere model.
//!
//! This function creates an Earth atmosphere model and copies the state of
//! the supplied Earth atmosphere model. It returns a handle (or id) to the copied model.
//! \param atmos A handle to the source model.
//! \returns A handle to the copied model.
EarthAtmosphere_C* copyAtmosphere_E(EarthAtmosphere_C* atmos)
{
  if (atmos != nullptr) {
    EarthAtmosphere* atmosCopy = new EarthAtmosphere(*atmos);
    return atmosCopy;
  }
  else {
    return nullptr;
  }
}

//! \brief Releases memory associated with a Earth atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos A Earth atmosphere handle.
void deleteAtmosphere_E(EarthAtmosphere_C* atmos)
{
  if (atmos != nullptr) {
    delete atmos;
  }
}

//! \brief Set the thermosphere model.
//!
//! This parameter selects the model for thermosphere calculations.
//! \param model      A ThermosphereModelType (1 = MET, 2 = MSIS, or 3 = JB2008).
//! \param atmos      An Earth atmosphere handle.
void setThermosphereModel_E(EarthAtmosphere_C* atmos, int model)
{
  if (atmos != nullptr) {
    ThermosphereModelType modelType = EG_MET;
    if (model == 2) {
      modelType = EG_MSIS;
    }
    else if (model == 3) {
      modelType = EG_JB2008;
    }
    atmos->setThermosphereModel(modelType);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setSurfaceRoughness
//! \param atmos      An Earth atmosphere handle.
void setSurfaceRoughness_E(EarthAtmosphere_C* atmos, greal z0)
{
  if (atmos != nullptr) {
    atmos->setSurfaceRoughness(z0);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setAtmosPath
//! \param atmos      An Earth atmosphere handle.
void setAtmosPath_E(EarthAtmosphere_C* atmos, const char* path)
{
  if (atmos != nullptr && path != nullptr) {
    atmos->setAtmosPath(std::string(path));
  }
}

//! \brief Determine use of NCEP model over default MERRA-2 model.
//!
//! \param atmos      An Earth atmosphere handle.
//! \param useNCEP    Set to 1 to use the NCEP model, 0 for the MERRA-2 model.
void setUseNCEP_E(EarthAtmosphere_C* atmos, int useNCEP)
{
  if (atmos != nullptr) {
    atmos->setUseNCEP(useNCEP!=0);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setNCEPPath
//! \param atmos      An Earth atmosphere handle.
void setNCEPPath_E(EarthAtmosphere_C* atmos, const char* path)
{
  if (atmos != nullptr && path != nullptr) {
    atmos->setNCEPPath(std::string(path));
  }
}

//! \copydoc GRAM::EarthAtmosphere::setNCEPParameters
//! \param atmos      An Earth atmosphere handle.
void setNCEPParameters_E(EarthAtmosphere_C* atmos, int NCEPYear, int NCEPHour)
{
  if (atmos != nullptr) {
    atmos->setNCEPParameters(NCEPYear, NCEPHour);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setMERRA2Path
//! \param atmos      An Earth atmosphere handle.
void setMERRA2Path_E(EarthAtmosphere_C* atmos, const char* path)
{
  if (atmos != nullptr && path != nullptr) {
    atmos->setMERRA2Path(std::string(path));
  }
}

//! \copydoc GRAM::EarthAtmosphere::setMERRA2Parameters
//! \param atmos      An Earth atmosphere handle.
void setMERRA2Parameters_E(EarthAtmosphere_C* atmos, int M2Hour, greal latMin, greal latMax, greal lonMin, greal lonMax)
{
  if (atmos != nullptr) {
    atmos->setMERRA2Parameters(M2Hour, latMin, latMax, lonMin, lonMax);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setRRAPath
//! \param atmos      An Earth atmosphere handle.
void setRRAPath_E(EarthAtmosphere_C* atmos, const char* path)
{
  if (atmos != nullptr && path != nullptr) {
    atmos->setRRAPath(std::string(path));
  }
}

//! \copydoc GRAM::EarthAtmosphere::setRRASiteList
//! \param atmos      An Earth atmosphere handle.
void setRRASiteList_E(EarthAtmosphere_C* atmos, const char* fileName)
{
  if (atmos != nullptr && fileName != nullptr) {
    atmos->setRRASiteList(std::string(fileName));
  }
}

//! \brief Set the RRA parameters.
//!
//! This method selects the year of the RRA data set and defines a region in which RRA data is used.  There
//! is a smooth transition between the \p outerRadius (no RRA) and the \p innerRadius (all RRA) of the region.
//! \param year             Year code: 1 for 1983 RRAs, 2 for 2006 RRAs, 3 for 2013 RRAs, 4 for 2019 RRAs.
//! \param innerRadius      Lat-lon radius from RRA site, inside which RRA data is used with full weight of 1 \units{degrees}.
//! \param outerRadius      Lat-lon radius from RRA site, outside which RRA data are NOT used \units{degrees}.
//! \param atmos      An Earth atmosphere handle.
void setRRAParameters_E(EarthAtmosphere_C* atmos, int year, greal innerRadius, greal outerRadius)
{
  if (atmos != nullptr) {
    RRAYearType rraYear = RRA_2013;
    switch (year) {
    case 1:
      rraYear = RRA_1983;
      break;
    case 2:
      rraYear = RRA_2006;
      break;
    case 3:
      rraYear = RRA_2013;
      break;
    default:
      rraYear = RRA_2019;
      break;
    }
    atmos->setRRAParameters(rraYear, innerRadius, outerRadius);
  }
}

//! \brief Enable or disable the use of RRA data.
//!
//! \param useFlag    Set to 1 to enable the RRA data regions or 0 to disable.
//! \param atmos      An Earth atmosphere handle.
void setUseRRA_E(EarthAtmosphere_C* atmos, int useFlag)
{
  if (atmos != nullptr) {
    atmos->setUseRRA(useFlag != 0);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setSolarParameters
//! \param atmos      An Earth atmosphere handle.
void setSolarParameters_E(EarthAtmosphere_C* atmos, greal dailyF10, greal meanF10, greal ap)
{
  if (atmos != nullptr) {
    atmos->setSolarParameters(dailyF10, meanF10, ap);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setJB2008Parameters
//! \param atmos      An Earth atmosphere handle.
void setJB2008Parameters_E(EarthAtmosphere_C* atmos, greal dailyS10, greal meanS10, greal dailyXM10, greal meanXM10, greal dailyY10, greal meanY10, greal dstdtc)
{
  if (atmos != nullptr) {
    atmos->setJB2008Parameters(dailyS10, meanS10, dailyXM10, meanXM10, dailyY10, meanY10, dstdtc);
  }
}

//! \copydoc GRAM::EarthAtmosphere::setInitialPerturbations
//! \param atmos      An Earth atmosphere handle.
void setInitialPerturbations_E(EarthAtmosphere_C* atmos, greal densPert, greal tempPert, greal ewPert, greal nsPert, greal verticalPert)
{
  if (atmos != nullptr) {
    atmos->setInitialPerturbations(densPert, tempPert, ewPert, nsPert, verticalPert);
  }
}

//! \brief Enable patchiness in the perturbation model.
//! \param usePatchiness   Set to 1 to enable patchiness or 0 to disable.
//! \param atmos           An Earth atmosphere handle.
void setPatchiness_E(EarthAtmosphere_C* atmos, int usePatchiness)
{
  if (atmos != nullptr) {
    bool patchy = (usePatchiness != 0);
    atmos->setPatchiness(patchy);
  }
}

//! \copydoc GRAM::setStartTime_C
void setStartTime_E(EarthAtmosphere_C* atmos, const GramTime_C* time)
{
  setStartTime_C(atmos, time);
}

//! \copydoc GRAM::setSeed_C()
void setSeed_E(EarthAtmosphere_C* atmos, int seed)
{
  setSeed_C(atmos, seed);
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param randomScale For density, pressure, and temperature. Between 0 and 2.
//! \param horizontalWindScale For horizontal (E/W, N/W) winds.  Between 0 and 2.
//! \param verticalWindScale For vertical winds. Between 0 and 2.
void setPerturbationScales_E(EarthAtmosphere_C* atmos, greal randomScale, greal horizontalWindScale, greal verticalWindScale)
{
  if (atmos != nullptr) {
    atmos->setRandomPerturbationScale(randomScale);
    atmos->setHorizontalWindPerturbationScale(horizontalWindScale);
    atmos->setVerticalWindPerturbationScale(verticalWindScale);
  }
}

//! \copydoc GRAM::addAuxiliaryAtmosphere_C()
void addAuxiliaryAtmosphere_E(EarthAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  addAuxiliaryAtmosphere_C(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \copydoc GRAM::setAuxiliaryValues_C()
void setAuxiliaryValues_E(EarthAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  setAuxiliaryValues_C(atmos, dens, pres, temp, ew, ns);
}

//! \copydoc GRAM::setPosition_C()
void setPosition_E(EarthAtmosphere_C* atmos, const Position_C* pos)
{
  setPosition_C(atmos, pos);
}

//! \copydoc GRAM::setDelta_C()
void setDelta_E(EarthAtmosphere_C* atmos, const Position_C* delta)
{
  setDelta_C(atmos, delta);
}

//! \copydoc GRAM::setPerturbationAction_C()
void setPerturbationAction_E(EarthAtmosphere_C* atmos, int action)
{
  setPerturbationAction_C(atmos, action);
}


//! \copydoc GRAM::setEphemerisState_C()
void setEphemerisState_E(EarthAtmosphere_C* atmos, const EphemerisState_C* state)
{
  setEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::setEphemerisFastModeOn_C()
void setEphemerisFastModeOn_E(EarthAtmosphere_C* atmos, int flag)
{
  setEphemerisFastModeOn_C(atmos, flag);
}

//! \copydoc GRAM::setSubsolarUpdateTime_C()
void setSubsolarUpdateTime_E(EarthAtmosphere_C* atmos, greal utime)
{
  setSubsolarUpdateTime_C(atmos, utime);
}

//! \brief Performs the Earth atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos A Earth atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_E(EarthAtmosphere_C* atmos)
{
  int error = 0;
  errorMessage.clear();
  try {
    if (atmos != nullptr) {
      atmos->update();
    }
  }
  catch (const std::string& msg) {
    errorMessage = msg;
    error = 1;
  }
  return error;
}

//! \copydoc GRAM::getErrorMessage_C()
size_t getErrorMessage_E(char *message, size_t bufferSize)
{
  return getErrorMessage_C(message, bufferSize);
}

//! \brief Get the current version as a string.
//!
//! \param atmos A Earth atmosphere handle.
//! \param bufferSize The size of the versionString buffer.
//! \param[out] versionString A non-NULL pointer.
//! \retval versionString A character array large enough to hold the version string.
//! \returns The length of the version string.
size_t getVersionString_E(EarthAtmosphere_C* atmos, char* versionString, size_t bufferSize)
{
  if (atmos != nullptr && versionString != nullptr) {
    std::string tempString = atmos->getVersionString();
    size_t len = tempString.copy(versionString, bufferSize - 2);
    versionString[len] = '\0';
    return len;
  }
  return 0;
}

//! \copydoc GRAM::getPosition_C()
void getPosition_E(EarthAtmosphere_C* atmos, Position_C* position)
{
  getPosition_C(atmos, position);
}

//! \copydoc GRAM::getDynamicsState_C()
void getDynamicsState_E(EarthAtmosphere_C* atmos, DynamicsState_C* state)
{
  getDynamicsState_C(atmos, state);
}

//! \copydoc GRAM::getDensityState_C()
void getDensityState_E(EarthAtmosphere_C* atmos, DensityState_C* state)
{
  getDensityState_C(atmos, state);
}

//! \copydoc GRAM::getWindsState_C()
void getWindsState_E(EarthAtmosphere_C* atmos, WindsState_C* state)
{
  getWindsState_C(atmos, state);
}

//! \copydoc GRAM::getGasesState_C()
//! \param[out] argon A non-NULL pointer.
//! \param[out] carbonDioxide A non-NULL pointer.
//! \param[out] carbonMonoxide A non-NULL pointer.
//! \param[out] dinitrogen A non-NULL pointer.
//! \param[out] dioxygen A non-NULL pointer.
//! \param[out] helium A non-NULL pointer.
//! \param[out] hydrogen A non-NULL pointer.
//! \param[out] methane A non-NULL pointer.
//! \param[out] nitrogen A non-NULL pointer.
//! \param[out] nitrousOxide A non-NULL pointer.
//! \param[out] oxygen A non-NULL pointer.
//! \param[out] ozone A non-NULL pointer.
//! \param[out] water A non-NULL pointer.
//! \retval argon The argon (Ar) components.
//! \retval carbonDioxide The molecular carbon dioxide (CO2) components.
//! \retval carbonMonoxide The molecular carbon monoxide (CO) components.
//! \retval dinitrogen The molecular nitrogen (N2) components.
//! \retval dioxygen The molecular oxygen (O2) components.
//! \retval helium The atomic helium (He) components.
//! \retval hydrogen The atomic hydrogen (H) components.
//! \retval methane The methane (CH4) components.
//! \retval nitrogen The atomic nitrogen (N) components.
//! \retval nitrousOxide The nitrousOxide (N2O) components.
//! \retval oxygen The atomic oxygen (O) components.
//! \retval ozone The ozone (O3) components.
//! \retval water The water vapor (H2O) components.
void getGasesState_E(EarthAtmosphere_C* atmos, GasesState_C* state,
  ConstituentGas_C* argon, ConstituentGas_C* carbonDioxide, ConstituentGas_C* carbonMonoxide, ConstituentGas_C* dinitrogen,
  ConstituentGas_C* dioxygen, ConstituentGas_C* helium, ConstituentGas_C* hydrogen, ConstituentGas_C* methane, 
  ConstituentGas_C* nitrogen, ConstituentGas_C* nitrousOxide, ConstituentGas_C* oxygen, ConstituentGas_C* ozone, ConstituentGas_C* water)
{
  getGasesState_C(atmos, state);
  getConstituentGas_C(atmos, argon, ARGON);
  getConstituentGas_C(atmos, carbonDioxide, CARBON_DIOXIDE);
  getConstituentGas_C(atmos, carbonMonoxide, CARBON_MONOXIDE);
  getConstituentGas_C(atmos, dinitrogen, DINITROGEN);
  getConstituentGas_C(atmos, dioxygen, DIOXYGEN);
  getConstituentGas_C(atmos, helium, HELIUM);
  getConstituentGas_C(atmos, hydrogen, HYDROGEN);
  getConstituentGas_C(atmos, methane, METHANE);
  getConstituentGas_C(atmos, nitrogen, NITROGEN);
  getConstituentGas_C(atmos, nitrousOxide, NITROUS_OXIDE);
  getConstituentGas_C(atmos, oxygen, OXYGEN);
  getConstituentGas_C(atmos, ozone, OZONE);
  getConstituentGas_C(atmos, water, WATER);
}

//! \copydoc GRAM::getEphemerisState_C()
void getEphemerisState_E(EarthAtmosphere_C* atmos, EphemerisState_C* state)
{
  getEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::getPerturbationState_C()
void getPerturbationState_E(EarthAtmosphere_C* atmos, PerturbationState_C* state)
{
  getPerturbationState_C(atmos, state);
}

//! \brief Get the current Earth specific values after on update.
//!
//! \param atmos An Earth atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The Earth state computed in the last update.
void getEarthState_E(EarthAtmosphere_C* atmos, EarthState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const EarthAtmosphereState& estate = astate.getPlanetSpecificMetrics<EarthAtmosphereState>();

    state->perturbedTemperature = estate.perturbedTemperature;
    state->temperaturePerturbation = estate.temperaturePerturbation;
    state->temperatureStandardDeviation = astate.temperatureStandardDeviation;
    state->perturbedPressure = estate.perturbedPressure;
    state->pressurePerturbation = estate.pressurePerturbation;
    state->pressureStandardDeviation = astate.pressureStandardDeviation;

    state->vaporPressure = estate.vaporPressure;
    state->vaporPressureSD = estate.vaporPressureSD;
    state->vaporDensity = estate.vaporDensity;
    state->vaporDensitySD = estate.vaporDensitySD;
    state->dewPoint = estate.dewPoint;
    state->dewPointSD = estate.dewPointSD;
    state->relativeHumidity = estate.relativeHumidity;
    state->relativeHumiditySD = estate.relativeHumiditySD;

    state->geodeticLatitude = estate.geodeticLatitude;
    state->rraWeight = estate.rraWeight;
    memcpy(state->rraSiteName, estate.rraSiteName, 6);

    state->windSpeed = estate.windSpeed;
    state->windSpeedStandardDeviation = estate.windSpeedStandardDeviation;
    state->windCorrelation = estate.windCorrelation;
    state->severityLevel = estate.severityLevel;
  }
}

//! \brief Get the current Earth specific perturbations after on update.
//!
//! \param atmos An Earth atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The Earth specific perturbations computed in the last update.
void getEarthPerts_E(EarthAtmosphere_C* atmos, EarthPerts_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const EarthAtmosphereState& estate = astate.getPlanetSpecificMetrics<EarthAtmosphereState>();

    state->presPertSmall = estate.presPertSmall;
    state->densPertSmall = estate.densPertSmall;
    state->tempPertSmall = estate.tempPertSmall;
    state->ewWindPertSmall = estate.ewWindPertSmall;
    state->nsWindPertSmall = estate.nsWindPertSmall;

    state->presSDSmall = estate.presSDSmall;
    state->densSDSmall = estate.densSDSmall;
    state->tempSDSmall = estate.tempSDSmall;
    state->ewWindSDSmall = estate.ewWindSDSmall;
    state->nsWindSDSmall = estate.nsWindSDSmall;

    state->presPertLarge = estate.presPertLarge;
    state->densPertLarge = estate.densPertLarge;
    state->tempPertLarge = estate.tempPertLarge;
    state->ewWindPertLarge = estate.ewWindPertLarge;
    state->nsWindPertLarge = estate.nsWindPertLarge;

    state->presSDLarge = estate.presSDLarge;
    state->densSDLarge = estate.densSDLarge;
    state->tempSDLarge = estate.tempSDLarge;
    state->ewWindSDLarge = estate.ewWindSDLarge;
    state->nsWindSDLarge = estate.nsWindSDLarge;
  }
}

//! \brief Get the current Earth surface values after on update.
//!
//! \param atmos An Earth atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The Earth surface values computed in the last update.
void getEarthSurface_E(EarthAtmosphere_C* atmos, EarthSurface_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const EarthAtmosphereState& estate = astate.getPlanetSpecificMetrics<EarthAtmosphereState>();

    state->windSpeedAtSurface = estate.windSpeedAtSurface;
    state->windSpeedSDAtSurface = estate.windSpeedSDAtSurface;
    state->temperatureAtSurface = estate.temperatureAtSurface;
    state->temperatureSDAtSurface = estate.temperatureSDAtSurface;
    state->pressureSDAtSurface = estate.pressureSDAtSurface;
    state->densitySDAtSurface = estate.densitySDAtSurface;
    state->densityAtSurface = estate.densityAtSurface;
    state->ewWindAtSurface = estate.ewWindAtSurface;
    state->nsWindAtSurface = estate.nsWindAtSurface;
    state->ewWindSDAtSurface = estate.ewWindSDAtSurface;
    state->nsWindSDAtSurface = estate.nsWindSDAtSurface;
    state->windCorrelationAtSurface = estate.windCorrelationAtSurface;
  }
}

//! \brief Get the current Earth boundary layer values after on update.
//!
//! \param atmos An Earth atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The Earth boundary layer values computed in the last update.
void getEarthBoundaryLayer_E(EarthAtmosphere_C* atmos, EarthBoundaryLayer_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const EarthAtmosphereState& estate = astate.getPlanetSpecificMetrics<EarthAtmosphereState>();

    state->landCode = estate.landCode;
    state->surfaceRoughness = estate.surfaceRoughness;
    state->netRadiationIndex = estate.netRadiationIndex;
    state->stability = estate.stability;
    state->inverseLength = estate.inverseLength;
    state->frictionVelocity = estate.frictionVelocity;
    state->BVFrequencySquare = estate.BVFrequencySquare;
    state->boundaryLayerDepth = estate.boundaryLayerDepth;
    state->neutralBoundaryLayerDepth = estate.neutralBoundaryLayerDepth;
    state->unstableBLFactor = estate.unstableBLFactor;
    state->sigmaRatio = estate.sigmaRatio;
    state->sigmaW = estate.sigmaW;
    state->metersAboveSurface = estate.metersAboveSurface;
    state->perturbedWindSpeedAtSurface = estate.perturbedWindSpeedAtSurface;
  }
}


//! \copydoc GRAM::getStartTime_C()
void getStartTime_E(EarthAtmosphere_C* atmos, GramTime_C* time)
{
  getStartTime_C(atmos, time);
}

//! @}


