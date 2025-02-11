//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include "gram.h"
#include "MarsAtmosphere_C.h"
#include "SpiceLoader.h"
#include "GramTime.h"

using namespace GRAM;

//! \defgroup C_Mars The C Interface for Mars 
//! @{
//! \brief A listing of the C interface functions for MarsGRAM.
//!
//! The C interface for Mars is declared in the following:
//! \code #include "MarsAtmosphere_C.h" \endcode 
//! An example using the C interface can be found \ref Mars/examples/Mars_C.c "here".
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
//! C and C++ code below.
//! \code
//! // In C
//! MarsAtmosphere_C* mars = createAtmosphere_M(dataPath);
//! setPosition_M(mars, &pos);
//! update_M(mars); 
//! 
//! // In C++
//! MarsAtmosphere mars(dataPath);
//! mars.setPosition(pos);
//! mars.update(); 
//! \endcode
//! \example Mars/examples/Mars_C.c

//! \copydoc GRAM::tryGetSpicePath_C
void tryGetSpicePath_M(char* spicePath, int bufferSize)
{
  tryGetSpicePath_C(spicePath, bufferSize);
}

//! \copydoc GRAM::tryGetDataPaths_C
void tryGetDataPaths_M(char* spicePath, char* dataPath, int bufferSize)
{
  tryGetDataPaths_C(spicePath, dataPath, bufferSize);
}

//! \copydoc GRAM::setSpiceLsk_C
void setSpiceLsk_M(const char* lsk)
{
  setSpiceLsk_C(lsk);
}

//! \copydoc GRAM::setSpicePck_C
void setSpicePck_M(const char* pck)
{
  setSpicePck_C(pck);
}

//! \brief Override the default SPICE Mars kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The Mars kernel path/name is typically of the form "/spk/satellites/mar097.bsp" where
//! the first three letters name the planet and the digits provide for versioning.
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_M(const char* bsp)
{
  setSpiceKernel_C(MARS_C, bsp);
}

//! \copydoc GRAM::initialize_C
void initialize_M(const char* spicePath)
{
  initialize_C(spicePath);
  registerBody_C(MARS,
    reinterpret_cast<CreateFunction>(createAtmosphere_M),
    reinterpret_cast<CopyFunction>(copyAtmosphere_M),
    reinterpret_cast<UpdateFunction>(update_M),
    reinterpret_cast<DeleteFunction>(deleteAtmosphere_M));
}

//! \copydoc GRAM::loadSpiceFile_C
void loadSpiceFile_M(const char* fileName)
{
  loadSpiceFile_C(fileName);
}

//! \brief Create a handle to a Mars atmosphere model.
//!
//! This function creates and initializes a Mars atmosphere model.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \returns A handle to the new model.
MarsAtmosphere_C* createAtmosphere_M(const char* dataPath)
{
  MarsAtmosphere* atmos = nullptr;
  try {
    // Load default settings (defined in MarsInputParameters)
    MarsInputParameters params;
    // Except let's not override the SpiceDataPath
    params.spicePath = SpiceLoader::getSpiceDataPath();
    if (dataPath != NULL) {
      params.dataPath = dataPath;
    }
    else {
      params.dataPath.clear();
    }
    atmos = new MarsAtmosphere();
    atmos->setInputParameters(params);
  }
  catch (const std::string& msg) {
    errorMessage = msg;
    if (atmos != nullptr) {
      delete atmos;
      atmos = nullptr;
    }
  }
  return atmos;
}

//! \brief Copies a Mars atmosphere model.
//!
//! This function creates a Mars atmosphere model and copies the state of
//! the supplied Mars atmosphere model. It returns a handle (or id) to the copied model.
//! \param atmos A handle to the source model.
//! \returns A handle to the copied model.
MarsAtmosphere_C* copyAtmosphere_M(MarsAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    MarsAtmosphere* atmosCopy = new MarsAtmosphere(*atmos);
    return atmosCopy;
  }
  else {
    return NULL;
  }
}

//! \brief Releases memory associated with a Mars atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos A Mars atmosphere handle.
void deleteAtmosphere_M(MarsAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    delete atmos;
  }
}

//! \brief Set the radii for the reference ellipsoid.
//!
//! Set the equatorial and polar radii for the reference ellipsiod.  The values will override
//! the internal constants, \p equatorialRadiusDefault and ]p polarRadiusDefault, set in MarsCommon.h.
//! \param atmos      A Mars atmosphere handle.
//! \param equatorial The equatorial radius \units{km}.
//! \param polar      The polar radius      \units{km}.
void setPlanetaryRadii_M(MarsAtmosphere_C* atmos, greal equatorial, greal polar)
{
  if (atmos != NULL) {
    atmos->setPlanetaryRadii(equatorial, polar);
  }
}

//! \copydoc GRAM::MarsGCM::setMapYear(int year)
//! \param atmos      A Mars atmosphere handle.
void setMapYear_M(MarsAtmosphere_C* atmos, int year)
{
  if (atmos != NULL) {
    atmos->setMapYear(year);
  }
}

//! \brief Set the height offset model.
//!
//! Height offsets can be used for two purposes: (1) to control the smoothness of 
//! the transition at 80 km altitude between MGCM data and MTGCM data, or (2) as
//! another means (besides wave multipliers) to adjust MTGCM data for better
//! agreement with observations, such as obtained during aerobraking operations.
//! \param atmos     A Mars atmosphere handle.
//! \param model     The Mars offset model option (0=MARS_CONSTANT, 1=MARS_SEASONAL, 2=MARS_GLOBAL_MEAN, 3=MARS_DAILY_AVERAGE, 4=MARS_CURRENT).
//! \param hgtOffset The height offset required by models 0 and 1.
void setHeightOffsetModel_M(MarsAtmosphere_C* atmos, int model, greal hgtOffset)
{
  if (atmos != NULL) {
    MarsOffsetModel offsetModel;
    switch (model) {
    case 0: 
      offsetModel = MARS_CONSTANT; 
      break;
    case 1: 
      offsetModel = MARS_SEASONAL; 
      break;
    case 2: 
      offsetModel = MARS_GLOBAL_MEAN; 
      break;
    case 3: 
      offsetModel = MARS_DAILY_AVERAGE; 
      break;
    case 4: 
      offsetModel = MARS_CURRENT;
      break;
    default: 
      offsetModel = MARS_GLOBAL_MEAN;
      break;
    }
    atmos->setHeightOffsetModel(offsetModel, hgtOffset);
  }
}

//! \copydoc  GRAM::MarsAtmosphere::setHeightAboveSurface(greal heightAboveSurface)
//! \param atmos      A Mars atmosphere handle.
void setHeightAboveSurface_M(MarsAtmosphere_C* atmos, greal heightAboveSurface)
{
  if (atmos != NULL) {
    atmos->setHeightAboveSurface(heightAboveSurface);
  }
}

//! \copydoc  GRAM::MarsAtmosphere::setMGCMDustLevels(greal constantDustLevel, greal minDustLevel, greal maxDustLevel)
//! \param atmos      A Mars atmosphere handle.
void setMGCMDustLevels_M(MarsAtmosphere_C* atmos, greal constantDustLevel, greal minDustLevel, greal maxDustLevel)
{
  if (atmos != NULL) {
    atmos->setMGCMDustLevels(constantDustLevel, minDustLevel, maxDustLevel);
  }
}

//! \copydoc  GRAM::MarsAtmosphere::setF107(greal f107)
//! \param atmos      A Mars atmosphere handle.
void setF107_M(MarsAtmosphere_C* atmos, greal f107)
{
  if (atmos != NULL) {
    atmos->setF107(f107);
  }
}

//! \copydoc  GRAM::MarsAtmosphere::setPerturbationWaveLengthScale(greal scale)
//! \param atmos      A Mars atmosphere handle.
void setPerturbationWaveLengthScale_M(MarsAtmosphere_C* atmos, greal scale)
{
  if (atmos != NULL) {
    atmos->setPerturbationWaveLengthScale(scale);
  }
}

//! \brief Set the use of MOLA topographic heights.
//!
//! Set to true if input heights are relative to the MOLA areoid. 
//! Otherwise, input heights are relative to reference ellipsoid.
//! \param isMola     Set to 1 (true) if input heights are relative to the MOLA areoid.
//! \param atmos      A Mars atmosphere handle.
void setMOLAHeights_M(MarsAtmosphere_C* atmos, int isMola)
{
  if (atmos != NULL) {
    atmos->setMOLAHeights(isMola != 0);
  }
}

//! \brief Controls the computation of daily min/max values.
//!
//! Set to true if daily min/max values are to be computed.
//! \param minMax     Set to 1 (true) if daily min/max values are to be computed.
//! \param atmos      A Mars atmosphere handle.
void setMinMax_M(MarsAtmosphere_C* atmos, int minMax)
{
  if (atmos != NULL) {
    atmos->setMinMax(minMax != 0);
  }
}

//! \copydoc  GRAM::MarsAtmosphere::setDustStorm(greal lonSun, greal duration, greal intensity, greal maxRadius, greal lat, greal lon)
//! \param atmos      A Mars atmosphere handle.
void setDustStorm_M(MarsAtmosphere_C* atmos, greal lonSun, greal duration, greal intensity, greal maxRadius, greal lat, greal lon)
{
  if (atmos != NULL) {
    atmos->setDustStorm(lonSun, duration, intensity, maxRadius, lat, lon);
  }
}

//! \copydoc GRAM::MarsAtmosphere::setDustDensity(greal nu, greal diameter, greal dens)
//! \param atmos      A Mars atmosphere handle.
void setDustDensity_M(MarsAtmosphere_C* atmos, greal nu, greal diameter, greal dens)
{
  if (atmos != NULL) {
    atmos->setDustDensity(nu, diameter, dens);
  }
}

//! \copydoc GRAM::MarsAtmosphere::setExosphericTemperature(greal offset, greal factor)
//! \param atmos      A Mars atmosphere handle.
void setExosphericTemperature_M(MarsAtmosphere_C* atmos, greal offset, greal factor)
{
  if (atmos != NULL) {
    atmos->setExosphericTemperature(offset, factor);
  }
}

//! \copydoc GRAM::MarsAtmosphere::setWaveDefaults(greal date, greal scale, greal mean, greal a1, greal p1, greal r1, greal a2, greal p2, greal r2, greal a3, greal p3, greal r3)
//! \param atmos      A Mars atmosphere handle.
void setWaveDefaults_M(MarsAtmosphere_C* atmos, greal date, greal scale, greal mean, greal a1, greal p1, greal r1, greal a2, greal p2, greal r2, greal a3, greal p3, greal r3)
{
  if (atmos != NULL) {
    atmos->setWaveDefaults(date, scale, mean, a1, p1, r1, a2, p2, r2, a3, p3, r3);
  }
}

//! \brief Set the wave perturbation parameter file.
//!
//! Specifies the path to the optional file containing time-dependent wave coefficient data.
//! See setWaveDefaults() for a list of parameters preceded by an elasped time \units{\text{seconds}}.
//! \param atmos      A Mars atmosphere handle.
//! \param waveFile   The full or relative path to the wave file.
void setWaveFile_M(MarsAtmosphere_C* atmos, const char* waveFile)
{
  if (atmos != NULL) {
    atmos->setWaveFile(std::string(waveFile));
  }
}

//! \copydoc GRAM::MarsAtmosphere::setWindScales(greal meanWinds, greal boundaryLayerWinds)
//! \param atmos      A Mars atmosphere handle.
void setWindScales_M(MarsAtmosphere_C* atmos, greal meanWinds, greal boundaryLayerWinds)
{
  if (atmos != NULL) {
    atmos->setWindScales(meanWinds, boundaryLayerWinds);
  }
}

//! \copydoc GRAM::setStartTime_C
void setStartTime_M(MarsAtmosphere_C* atmos, const GramTime_C* time)
{
  setStartTime_C(atmos, time);
}

//! \copydoc GRAM::setSeed_C()
void setSeed_M(MarsAtmosphere_C* atmos, int seed)
{
  setSeed_C(atmos, seed);
}

//! \copydoc GRAM::setMinRelativeStepSize_C()
void setMinRelativeStepSize_M(MarsAtmosphere_C* atmos, greal minRelativeStepSize)
{
  setMinRelativeStepSize_C(atmos, minRelativeStepSize);
}

//! \copydoc GRAM::setPerturbationScales_C()
void setPerturbationScales_M(MarsAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale)
{
  setPerturbationScales_C(atmos, densityScale, ewWindScale, nsWindScale, verticalWindScale);
}

//! \copydoc GRAM::addAuxiliaryAtmosphere_C()
void addAuxiliaryAtmosphere_M(MarsAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  addAuxiliaryAtmosphere_C(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \copydoc GRAM::setAuxiliaryValues_C()
void setAuxiliaryValues_M(MarsAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  setAuxiliaryValues_C(atmos, dens, pres, temp, ew, ns);
}

//! \copydoc GRAM::setPosition_C()
void setPosition_M(MarsAtmosphere_C* atmos, const Position_C* pos)
{
  setPosition_C(atmos, pos);
}

//! \copydoc GRAM::setDelta_C()
void setDelta_M(MarsAtmosphere_C* atmos, const Position_C* delta)
{
  setDelta_C(atmos, delta);
}

//! \copydoc GRAM::setPerturbationAction_C()
void setPerturbationAction_M(MarsAtmosphere_C* atmos, int action)
{
  setPerturbationAction_C(atmos, action);
}


//! \copydoc GRAM::setEphemerisState_C()
void setEphemerisState_M(MarsAtmosphere_C* atmos, const EphemerisState_C* state)
{
  setEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::setEphemerisFastModeOn_C()
void setEphemerisFastModeOn_M(MarsAtmosphere_C* atmos, int flag)
{
  setEphemerisFastModeOn_C(atmos, flag);
}

//! \copydoc GRAM::setSubsolarUpdateTime_C()
void setSubsolarUpdateTime_M(MarsAtmosphere_C* atmos, greal utime)
{
  setSubsolarUpdateTime_C(atmos, utime);
}

//! \copydoc GRAM::MOLATopography::setAreoidRadiusCallback()
//! \param atmos      A Mars atmosphere handle.
void setAreoidRadiusCallback_M(MarsAtmosphere_C* atmos, TopoCallback callback)
{
  if (atmos != nullptr) {
    atmos->setAreoidRadiusCallback(callback);
  }
}

//! \copydoc GRAM::MOLATopography::setTopographicHeightCallback()
//! \param atmos      A Mars atmosphere handle.
void setTopographicHeightCallback_M(MarsAtmosphere_C* atmos, TopoCallback callback)
{
  if (atmos != nullptr) {
    atmos->setTopographicHeightCallback(callback);
  }
}

//! \copydoc GRAM::MOLATopography::setCallbackData()
//! \param atmos      A Mars atmosphere handle.
void setCallbackData_M(MarsAtmosphere_C* atmos, void* dataPointer)
{
  if (atmos != nullptr) {
    atmos->setCallbackData(dataPointer);
  }
}

//! \brief Performs the Mars atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated.
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos A Mars atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_M(MarsAtmosphere_C* atmos)
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
size_t getErrorMessage_M(char *message, size_t bufferSize)
{
  return getErrorMessage_C(message, bufferSize);
}


//! \brief Get the current version as a string.
//!
//! \param atmos A Mars atmosphere handle.
//! \param bufferSize The size of the versionString buffer.
//! \param[out] versionString A non-NULL pointer.
//! \retval versionString A character array large enough to hold the version string.
//! \returns The length of the version string.
size_t getVersionString_M(MarsAtmosphere_C* atmos, char* versionString, size_t bufferSize)
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
void getPosition_M(MarsAtmosphere_C* atmos, Position_C* position)
{
  getPosition_C(atmos, position);
}

//! \copydoc GRAM::getDynamicsState_C()
void getDynamicsState_M(MarsAtmosphere_C* atmos, DynamicsState_C* state)
{
  getDynamicsState_C(atmos, state);
}

//! \copydoc GRAM::getDensityState_C()
void getDensityState_M(MarsAtmosphere_C* atmos, DensityState_C* state)
{
  getDensityState_C(atmos, state);
}

//! \copydoc GRAM::getWindsState_C()
void getWindsState_M(MarsAtmosphere_C* atmos, WindsState_C* state)
{
  getWindsState_C(atmos, state);
}

//! \copydoc GRAM::getGasesState_C()
//! \param[out] hydrogen A non-NULL pointer.
//! \param[out] dihydrogen A non-NULL pointer.
//! \param[out] argon A non-NULL pointer.
//! \param[out] helium A non-NULL pointer.
//! \param[out] oxygen A non-NULL pointer.
//! \param[out] dioxygen A non-NULL pointer.
//! \param[out] dinitrogen A non-NULL pointer.
//! \param[out] carbonMonoxide A non-NULL pointer.
//! \param[out] carbonDioxide A non-NULL pointer.
//! \param[out] water A non-NULL pointer.
//! \retval hydrogen The atomic hydrogen (H) components.
//! \retval dihydrogen The molecular hydrogen (H2) components.
//! \retval argon The atomic argon (Ar) components.
//! \retval helium The atomic helium (He) components.
//! \retval oxygen The atomic oxygen (O) components.
//! \retval dioxygen The molecular oxygen (O2) components.
//! \retval dinitrogen The molecular nitrogen (N2) components.
//! \retval carbonMonoxide The molecular carbon monoxide (CO) components.
//! \retval carbonDioxide The molecular carbon dioxide (CO2) components.
//! \retval water The water vapor (H2O) components.
void getGasesState_M(MarsAtmosphere_C* atmos, GasesState_C* state,
  ConstituentGas_C* argon, ConstituentGas_C* carbonDioxide, ConstituentGas_C* carbonMonoxide,
  ConstituentGas_C* dihydrogen, ConstituentGas_C* dinitrogen, ConstituentGas_C* dioxygen,
  ConstituentGas_C* helium, ConstituentGas_C* hydrogen, ConstituentGas_C* oxygen, ConstituentGas_C* water)
{
  getGasesState_C(atmos, state);
  getConstituentGas_C(atmos, hydrogen, HYDROGEN);
  getConstituentGas_C(atmos, dihydrogen, DIHYDROGEN);
  getConstituentGas_C(atmos, argon, ARGON);
  getConstituentGas_C(atmos, helium, HELIUM);
  getConstituentGas_C(atmos, oxygen, OXYGEN);
  getConstituentGas_C(atmos, dioxygen, DIOXYGEN);
  getConstituentGas_C(atmos, dinitrogen, DINITROGEN);
  getConstituentGas_C(atmos, carbonMonoxide, CARBON_MONOXIDE);
  getConstituentGas_C(atmos, carbonDioxide, CARBON_DIOXIDE);
  getConstituentGas_C(atmos, water, WATER);
}

//! \copydoc GRAM::getEphemerisState_C()
void getEphemerisState_M(MarsAtmosphere_C* atmos, EphemerisState_C* state)
{
  getEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::getPerturbationState_C()
void getPerturbationState_M(MarsAtmosphere_C* atmos, PerturbationState_C* state)
{
  getPerturbationState_C(atmos, state);
}

//! \brief Get the current daily mean values after on update.
//!
//! \param atmos A Mars atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The daily mean state computed in the last update.
void getDailyDynamicsState_M(MarsAtmosphere_C* atmos, DailyDynamicsState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const MarsAtmosphereState& mstate = astate.getPlanetSpecificMetrics<MarsAtmosphereState>();
    state->temperatureDaily = mstate.temperatureDaily;
    state->pressureDaily = mstate.pressureDaily;
    state->densityDaily = mstate.densityDaily;
    state->ewWindDaily = mstate.ewWindDaily;
    state->nsWindDaily = mstate.nsWindDaily;
    state->densityMin = mstate.densityMin;
    state->densityMax = mstate.densityMax;
    state->temperatureMin = mstate.temperatureMin;
    state->temperatureMax = mstate.temperatureMax;
  }
}

//! \brief Get the current Mars specific values after on update.
//!
//! \param atmos A Mars atmosphere handle.
//! \param[out] state A non-NULL pointer.
//! \retval state The ephemeris state computed in the last update.
void getMarsState_M(MarsAtmosphere_C* atmos, MarsState_C* state)
{
  if (atmos != nullptr && state != nullptr) {
    const AtmosphereState& astate = atmos->getAtmosphereState();
    const MarsAtmosphereState& mstate = astate.getPlanetSpecificMetrics<MarsAtmosphereState>();
    state->planetoGraphicHeight = mstate.planetoGraphicHeight;
    state->planetoGraphicLatitude = mstate.planetoGraphicLatitude;
    state->referenceHeight = mstate.referenceHeight;
    state->referenceRadius = mstate.referenceRadius;
    state->groundTemperature = mstate.groundTemperature;
    state->thermosphereBaseHeight = mstate.thermosphereBaseHeight;
    state->thermosphereBaseTemperature = mstate.thermosphereBaseTemperature;
    state->exosphericTemperature = mstate.exosphericTemperature;
    state->f1PeakHeight = mstate.f1PeakHeight;
    state->albedo = mstate.albedo;
    state->heightOffset = mstate.heightOffset;
    state->localHeightOffset = mstate.localHeightOffset;
    state->dustOpticalDepth = mstate.dustOpticalDepth;
    state->dustColumnArealDensity = mstate.dustColumnArealDensity;
    state->dustMixingRatio = mstate.dustMixingRatio;
    state->dustMassDensity = mstate.dustMassDensity;
    state->dustNumberDensity = mstate.dustNumberDensity;
    state->iceIsPresent = mstate.iceIsPresent;
  }
}

//! \copydoc GRAM::getStartTime_C()
void getStartTime_M(MarsAtmosphere_C* atmos, GramTime_C* time)
{
  getStartTime_C(atmos, time);
}


//! @}


