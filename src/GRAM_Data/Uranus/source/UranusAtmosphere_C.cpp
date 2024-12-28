//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "gram.h"
#include "UranusAtmosphere_C.h"
#include "SpiceLoader.h"
#include "GramTime.h"

using namespace GRAM;

//! \defgroup C_Uranus The C Interface for Uranus 
//! @{
//! The C interface for Uranus is declared in the following:
//! \code #include "UranusAtmosphere_C.h" \endcode 
//! An example using the C interface can be found \ref Uranus/examples/Uranus_C.c "here".
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
//! C and C++ code below.
//! \code
//! // In C
//! UranusAtmosphere_C* uranus = createAtmosphere_U();
//! setPosition_U(uranus, &pos);
//! update_U(uranus); 
//! 
//! // In C++
//! UranusAtmosphere uranus;
//! uranus.setPosition(pos);
//! uranus.update(); 
//! \endcode
//! \example Uranus/examples/Uranus_C.c

//! \copydoc GRAM::tryGetSpicePath_C
void tryGetSpicePath_U(char* spicePath, int bufferSize)
{
  tryGetSpicePath_C(spicePath, bufferSize);
}

//! \copydoc GRAM::initialize_C
void initialize_U(const char* spicePath)
{
  initialize_C(spicePath);
  registerBody_C(URANUS,
    reinterpret_cast<CreateFunction>(createAtmosphere_U),
    reinterpret_cast<CopyFunction>(copyAtmosphere_U),
    reinterpret_cast<UpdateFunction>(update_U),
    reinterpret_cast<DeleteFunction>(deleteAtmosphere_U));
}

//! \copydoc GRAM::loadSpiceFile_C
void loadSpiceFile_U(const char* fileName)
{
  loadSpiceFile_C(fileName);
}

//! \copydoc GRAM::setSpiceLsk_C
void setSpiceLsk_U(const char* lsk)
{
  setSpiceLsk_C(lsk);
}

//! \copydoc GRAM::setSpicePck_C
void setSpicePck_U(const char* pck)
{
  setSpicePck_C(pck);
}

//! \brief Override the default SPICE Uranus kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The Uranus kernel path/name is typically of the form "/spk/satellites/ura116.bsp" where
//! the first three letters name the planet and the digits provide for versioning.
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_U(const char* bsp)
{
  setSpiceKernel_C(EARTH_C, bsp);
}

//! \brief Create a handle to a Uranus atmosphere model.
//!
//! This function creates and initializes a Uranus atmosphere model.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \returns A handle to the new model.
UranusAtmosphere_C* createAtmosphere_U()
{
  // Load default settings (defined in UranusInputParameters)
  UranusInputParameters params;
  // Except let's not override the SpiceDataPath
  params.spicePath = SpiceLoader::getSpiceDataPath();
  UranusAtmosphere* atmos = new UranusAtmosphere();
  atmos->setInputParameters(params);
  return atmos;
}

//! \brief Copies a Uranus atmosphere model.
//!
//! This function creates a Uranus atmosphere model and copies the state of
//! the supplied Uranus atmosphere model. It returns a handle (or id) to the copied model.
//! \param atmos A handle to the source model.
//! \returns A handle to the copied model.
UranusAtmosphere_C* copyAtmosphere_U(UranusAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    UranusAtmosphere* atmosCopy = new UranusAtmosphere(*atmos);
    return atmosCopy;
  }
  else {
    return NULL;
  }
}

//! \brief Releases memory associated with a Uranus atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos A Uranus atmosphere handle.
void deleteAtmosphere_U(UranusAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    delete atmos;
  }
}

//! \copydoc GRAM::setStartTime_C
void setStartTime_U(UranusAtmosphere_C* atmos, const GramTime_C* time)
{
  setStartTime_C(atmos, time);
}

//! \copydoc GRAM::setSeed_C()
void setSeed_U(UranusAtmosphere_C* atmos, int seed)
{
  setSeed_C(atmos, seed);
}

//! \copydoc GRAM::setMinRelativeStepSize_C()
void setMinRelativeStepSize_U(UranusAtmosphere_C* atmos, greal minRelativeStepSize)
{
  setMinRelativeStepSize_C(atmos, minRelativeStepSize);
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param densityScale Between 0 and 2.
//! \param ewWindScale Between 0 and 2.
//! \param nsWindScale Between 0 and 2.
void setPerturbationScales_U(UranusAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale)
{
  setPerturbationScales_C(atmos, densityScale, ewWindScale, nsWindScale, 0.0);
}

//! \copydoc GRAM::addAuxiliaryAtmosphere_C()
void addAuxiliaryAtmosphere_U(UranusAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  addAuxiliaryAtmosphere_C(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \copydoc GRAM::setAuxiliaryValues_C()
void setAuxiliaryValues_U(UranusAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  setAuxiliaryValues_C(atmos, dens, pres, temp, ew, ns);
}

//! \copydoc GRAM::setPosition_C()
void setPosition_U(UranusAtmosphere_C* atmos, const Position_C* pos)
{
  setPosition_C(atmos, pos);
}

//! \copydoc GRAM::setDelta_C()
void setDelta_U(UranusAtmosphere_C* atmos, const Position_C* delta)
{
  setDelta_C(atmos, delta);
}

//! \copydoc GRAM::setPerturbationAction_C()
void setPerturbationAction_U(UranusAtmosphere_C* atmos, int action)
{
  setPerturbationAction_C(atmos, action);
}


//! \copydoc GRAM::setEphemerisState_C()
void setEphemerisState_U(UranusAtmosphere_C* atmos, const EphemerisState_C* state)
{
  setEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::setEphemerisFastModeOn_C()
void setEphemerisFastModeOn_U(UranusAtmosphere_C* atmos, int flag)
{
  setEphemerisFastModeOn_C(atmos, flag);
}

//! \copydoc GRAM::setSubsolarUpdateTime_C()
void setSubsolarUpdateTime_U(UranusAtmosphere_C* atmos, greal utime)
{
  setSubsolarUpdateTime_C(atmos, utime);
}

//! \brief Performs the Uranus atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos A Uranus atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_U(UranusAtmosphere_C* atmos)
{
  int error = 0;
  errorMessage.clear();
  try {
    if (atmos != NULL) {
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
size_t getErrorMessage_U(char *message, size_t bufferSize)
{
  return getErrorMessage_C(message, bufferSize);
}


//! \brief Get the current version as a string.
//!
//! \param atmos A Uranus atmosphere handle.
//! \param bufferSize The size of the versionString buffer.
//! \param[out] versionString A non-NULL pointer.
//! \retval versionString A character array large enough to hold the version string.
//! \returns The length of the version string.
size_t getVersionString_U(UranusAtmosphere_C* atmos, char* versionString, size_t bufferSize)
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
void getPosition_U(UranusAtmosphere_C* atmos, Position_C* position)
{
  getPosition_C(atmos, position);
}

//! \copydoc GRAM::getDynamicsState_C()
void getDynamicsState_U(UranusAtmosphere_C* atmos, DynamicsState_C* state)
{
  getDynamicsState_C(atmos, state);
}

//! \copydoc GRAM::getDensityState_C()
void getDensityState_U(UranusAtmosphere_C* atmos, DensityState_C* state)
{
  getDensityState_C(atmos, state);
}

//! \copydoc GRAM::getWindsState_C()
void getWindsState_U(UranusAtmosphere_C* atmos, WindsState_C* state)
{
  getWindsState_C(atmos, state);
}

//! \copydoc GRAM::getGasesState_C()
//! \param[out] dihydrogen A non-NULL pointer.
//! \param[out] methane A non-NULL pointer.
//! \param[out] helium A non-NULL pointer.
//! \retval dihydrogen The molecular hydrogen (H2) components.
//! \retval methane The atomic methane (CH4) components.
//! \retval helium The atomic helium (He) components.
//! \returns state The gases state computed in the last update.
void getGasesState_U(UranusAtmosphere_C* atmos, GasesState_C* state,
  ConstituentGas_C* dihydrogen, ConstituentGas_C* methane, ConstituentGas_C* helium)
{
  getGasesState_C(atmos, state);
  getConstituentGas_C(atmos, dihydrogen, DIHYDROGEN);
  getConstituentGas_C(atmos, methane, METHANE);
  getConstituentGas_C(atmos, helium, HELIUM);
}

//! \copydoc GRAM::getEphemerisState_C()
void getEphemerisState_U(UranusAtmosphere_C* atmos, EphemerisState_C* state)
{
  getEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::getPerturbationState_C()
void getPerturbationState_U(UranusAtmosphere_C* atmos, PerturbationState_C* state)
{
  getPerturbationState_C(atmos, state);
}

//! \copydoc GRAM::getStartTime_C()
void getStartTime_U(UranusAtmosphere_C* atmos, GramTime_C* time)
{
  getStartTime_C(atmos, time);
}

//! @}


