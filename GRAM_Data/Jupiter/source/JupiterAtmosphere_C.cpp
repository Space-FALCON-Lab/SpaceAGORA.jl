//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "gram.h"
#include "JupiterAtmosphere_C.h"
#include "SpiceLoader.h"
#include "GramTime.h"

using namespace GRAM;

//! \defgroup C_Jupiter The C Interface for Jupiter 
//! @{
//! The C interface for Jupiter is declared in the following:
//! \code #include "JupiterAtmosphere_C.h" \endcode 
//! An example using the C interface can be found \ref Jupiter/examples/Jupiter_C.c "here".
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
//! C and C++ code below.
//! \code
//! // In C
//! JupiterAtmosphere_C* jupiter = createAtmosphere_J();
//! setPosition_J(jupiter, &pos);
//! update_J(jupiter); 
//! 
//! // In C++
//! JupiterAtmosphere jupiter;
//! jupiter.setPosition(pos);
//! jupiter.update(); 
//! \endcode
//! \example Jupiter/examples/Jupiter_C.c

//! \copydoc GRAM::tryGetSpicePath_C
void tryGetSpicePath_J(char* spicePath, int bufferSize)
{
  tryGetSpicePath_C(spicePath, bufferSize);
}

//! \copydoc GRAM::setSpiceLsk_C
void setSpiceLsk_J(const char* lsk)
{
  setSpiceLsk_C(lsk);
}

//! \copydoc GRAM::setSpicePck_C
void setSpicePck_J(const char* pck)
{
  setSpicePck_C(pck);
}

//! \brief Override the default SPICE Jupiter kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The jupiter kernel path/name is typically of the form "/spk/satellites/jup365.bsp" where
//! the first three letters name the planet and the digits provide for versioning.
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_J(const char* bsp)
{
  setSpiceKernel_C(JUPITER_C, bsp);
}

//! \copydoc GRAM::initialize_C
void initialize_J(const char* spicePath)
{
  initialize_C(spicePath);
  registerBody_C(JUPITER,
    reinterpret_cast<CreateFunction>(createAtmosphere_J),
    reinterpret_cast<CopyFunction>(copyAtmosphere_J),
    reinterpret_cast<UpdateFunction>(update_J),
    reinterpret_cast<DeleteFunction>(deleteAtmosphere_J));
}

//! \copydoc GRAM::loadSpiceFile_C
void loadSpiceFile_J(const char* fileName)
{
  loadSpiceFile_C(fileName);
}

//! \brief Create a handle to a Jupiter atmosphere model.
//!
//! This function creates and initializes a Jupiter atmosphere model.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \returns A handle to the new model.
JupiterAtmosphere_C* createAtmosphere_J()
{
  // Load default settings (defined in JupiterInputParameters)
  JupiterInputParameters params;
  // Except let's not override the SpiceDataPath
  params.spicePath = SpiceLoader::getSpiceDataPath();
  JupiterAtmosphere* atmos = new JupiterAtmosphere();
  atmos->setInputParameters(params);
  return atmos;
}

//! \brief Copies a Jupiter atmosphere model.
//!
//! This function creates a Jupiter atmosphere model and copies the state of
//! the supplied Jupiter atmosphere model. It returns a handle (or id) to the copied model.
//! \param atmos A handle to the source model.
//! \returns A handle to the copied model.
JupiterAtmosphere_C* copyAtmosphere_J(JupiterAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    JupiterAtmosphere* atmosCopy = new JupiterAtmosphere(*atmos);
    return atmosCopy;
  }
  else {
    return NULL;
  }
}

//! \brief Releases memory associated with a Jupiter atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos A Jupiter atmosphere handle.
void deleteAtmosphere_J(JupiterAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    delete atmos;
  }
}

//! \copydoc GRAM::setStartTime_C
void setStartTime_J(JupiterAtmosphere_C* atmos, const GramTime_C* time)
{
  setStartTime_C(atmos, time);
}

//! \copydoc GRAM::setSeed_C()
void setSeed_J(JupiterAtmosphere_C* atmos, int seed)
{
  setSeed_C(atmos, seed);
}

//! \copydoc GRAM::setMinRelativeStepSize_C()
void setMinRelativeStepSize_J(JupiterAtmosphere_C* atmos, greal minRelativeStepSize)
{
  setMinRelativeStepSize_C(atmos, minRelativeStepSize);
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param densityScale Between 0 and 2.
//! \param ewWindScale Between 0 and 2.
//! \param nsWindScale Between 0 and 2.
void setPerturbationScales_J(JupiterAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale)
{
  setPerturbationScales_C(atmos, densityScale, ewWindScale, nsWindScale, 0.0);
}

//! \copydoc GRAM::addAuxiliaryAtmosphere_C()
void addAuxiliaryAtmosphere_J(JupiterAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  addAuxiliaryAtmosphere_C(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \copydoc GRAM::setAuxiliaryValues_C()
void setAuxiliaryValues_J(JupiterAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  setAuxiliaryValues_C(atmos, dens, pres, temp, ew, ns);
}

//! \copydoc GRAM::setPosition_C()
void setPosition_J(JupiterAtmosphere_C* atmos, const Position_C* pos)
{
  setPosition_C(atmos, pos);
}

//! \copydoc GRAM::setDelta_C()
void setDelta_J(JupiterAtmosphere_C* atmos, const Position_C* delta)
{
  setDelta_C(atmos, delta);
}

//! \copydoc GRAM::setPerturbationAction_C()
void setPerturbationAction_J(JupiterAtmosphere_C* atmos, int action)
{
  setPerturbationAction_C(atmos, action);
}


//! \copydoc GRAM::setEphemerisState_C()
void setEphemerisState_J(JupiterAtmosphere_C* atmos, const EphemerisState_C* state)
{
  setEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::setEphemerisFastModeOn_C()
void setEphemerisFastModeOn_J(JupiterAtmosphere_C* atmos, int flag)
{
  setEphemerisFastModeOn_C(atmos, flag);
}

//! \copydoc GRAM::setSubsolarUpdateTime_C()
void setSubsolarUpdateTime_J(JupiterAtmosphere_C* atmos, greal utime)
{
  setSubsolarUpdateTime_C(atmos, utime);
}


//! \brief Performs the Jupiter atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos A Jupiter atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_J(JupiterAtmosphere_C* atmos)
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
size_t getErrorMessage_J(char *message, size_t bufferSize)
{
  return getErrorMessage_C(message, bufferSize);
}

//! \brief Get the current version as a string.
//!
//! \param atmos A Jupiter atmosphere handle.
//! \param bufferSize The size of the versionString buffer.
//! \param[out] versionString A non-NULL pointer.
//! \retval versionString A character array large enough to hold the version string.
//! \returns The length of the version string.
size_t getVersionString_J(JupiterAtmosphere_C* atmos, char* versionString, size_t bufferSize)
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
void getPosition_J(JupiterAtmosphere_C* atmos, Position_C* position)
{
  getPosition_C(atmos, position);
}

//! \copydoc GRAM::getDynamicsState_C()
void getDynamicsState_J(JupiterAtmosphere_C* atmos, DynamicsState_C* state)
{
  getDynamicsState_C(atmos, state);
}

//! \copydoc GRAM::getDensityState_C()
void getDensityState_J(JupiterAtmosphere_C* atmos, DensityState_C* state)
{
  getDensityState_C(atmos, state);
}

//! \copydoc GRAM::getWindsState_C()
void getWindsState_J(JupiterAtmosphere_C* atmos, WindsState_C* state)
{
  getWindsState_C(atmos, state);
}

//! \copydoc GRAM::getGasesState_C()
// \retval dihydrogen The molecular hydrogen (H2) components.
// \retval helium The atomic helium (He) components.
void getGasesState_J(JupiterAtmosphere_C* atmos, GasesState_C* state)
//  ConstituentGas_C* dihydrogen, ConstituentGas_C* helium)
{
  getGasesState_C(atmos, state);
//  getConstituentGas_C(atmos, dihydrogen, DIHYDROGEN);
//  getConstituentGas_C(atmos, helium, HELIUM);
}

//! \copydoc GRAM::getEphemerisState_C()
void getEphemerisState_J(JupiterAtmosphere_C* atmos, EphemerisState_C* state)
{
  getEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::getPerturbationState_C()
void getPerturbationState_J(JupiterAtmosphere_C* atmos, PerturbationState_C* state)
{
  getPerturbationState_C(atmos, state);
}

//! \copydoc GRAM::getStartTime_C()
void getStartTime_J(JupiterAtmosphere_C* atmos, GramTime_C* time)
{
  getStartTime_C(atmos, time);
}

//! @}


