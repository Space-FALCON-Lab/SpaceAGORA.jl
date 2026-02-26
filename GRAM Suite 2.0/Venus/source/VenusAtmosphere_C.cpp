//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "gram.h"
#include "VenusAtmosphere_C.h"
#include "SpiceLoader.h"
#include "GramTime.h"

using namespace GRAM;

//! \defgroup C_Venus The C Interface for Venus 
//! @{
//! The C interface for Venus is declared in the following:
//! \code #include "VenusAtmosphere_C.h" \endcode 
//! An example using the C interface can be found \ref Venus/examples/Venus_C.c "here".
//! <br><br>The C interface is a wrapper around the C++ GRAM library. The design of the interface
//! intentionally mimics the C++ interface as closely as is possible. As an example, compare the
//! C and C++ code below.
//! \code
//! // In C
//! VenusAtmosphere_C* venus = createAtmosphere_V();
//! setPosition_V(venus, &pos);
//! update_V(venus); 
//! 
//! // In C++
//! VenusAtmosphere venus;
//! venus.setPosition(pos);
//! venus.update(); 
//! \endcode
//! \example Venus/examples/Venus_C.c

//! \copydoc GRAM::tryGetSpicePath_C
void tryGetSpicePath_V(char* spicePath, int bufferSize)
{
  tryGetSpicePath_C(spicePath, bufferSize);
}

//! \copydoc GRAM::initialize_C
void initialize_V(const char* spicePath)
{
  initialize_C(spicePath);
  registerBody_C(VENUS,
    reinterpret_cast<CreateFunction>(createAtmosphere_V),
    reinterpret_cast<CopyFunction>(copyAtmosphere_V),
    reinterpret_cast<UpdateFunction>(update_V),
    reinterpret_cast<DeleteFunction>(deleteAtmosphere_V));
}

//! \copydoc GRAM::loadSpiceFile_C
void loadSpiceFile_V(const char* fileName)
{
  loadSpiceFile_C(fileName);
}

//! \copydoc GRAM::setSpiceLsk_C
void setSpiceLsk_V(const char* lsk)
{
  setSpiceLsk_C(lsk);
}

//! \copydoc GRAM::setSpicePck_C
void setSpicePck_V(const char* pck)
{
  setSpicePck_C(pck);
}

//! \brief Override the default SPICE Venus kernel file name.
//!
//! Planetary data used by the SPICE library is contained within a binary kernel file.
//! The location of this file within the SPICE data folder can be specified with this function.
//! The Earth and Venus kernel path/name is typically of the form "/spk/planets/de440.bsp".
//! \param bsp  A string containing the PCK file path within the SPICE data folder
void setSpiceKernel_V(const char* bsp)
{
  setSpiceKernel_C(EARTH_C, bsp);
}

//! \brief Create a handle to a Venus atmosphere model.
//!
//! This function creates and initializes a Venus atmosphere model.
//! It returns a handle (or id) to the new model.  This handle must be 
//! passed to subsequent calls that utilize this model.  Memory associated
//! with the model is only accessible via the GRAM interface calls.
//! \returns A handle to the new model.
VenusAtmosphere_C* createAtmosphere_V()
{
  // Load default settings (defined in VenusInputParameters)
  VenusInputParameters params;
  // Except let's not override the SpiceDataPath
  params.spicePath = SpiceLoader::getSpiceDataPath();
  VenusAtmosphere* atmos = new VenusAtmosphere();
  atmos->setInputParameters(params);
  return atmos;
}

//! \brief Copies a Venus atmosphere model.
//!
//! This function creates a Venus atmosphere model and copies the state of
//! the supplied Venus atmosphere model. It returns a handle (or id) to the copied model.
//! \param atmos A handle to the source model.
//! \returns A handle to the copied model.
VenusAtmosphere_C* copyAtmosphere_V(VenusAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    VenusAtmosphere* atmosCopy = new VenusAtmosphere(*atmos);
    return atmosCopy;
  }
  else {
    return NULL;
  }
}

//! \brief Releases memory associated with a Venus atmosphere model.
//!
//! When an atmosphere model is no longer needed, the associated memory 
//! should be released by a call to this function.
//! \param atmos A Venus atmosphere handle.
void deleteAtmosphere_V(VenusAtmosphere_C* atmos)
{
  if (atmos != NULL) {
    delete atmos;
  }
}

//! \copydoc GRAM::setStartTime_C
void setStartTime_V(VenusAtmosphere_C* atmos, const GramTime_C* time)
{
  setStartTime_C(atmos, time);
}

//! \copydoc GRAM::setSeed_C()
void setSeed_V(VenusAtmosphere_C* atmos, int seed)
{
  setSeed_C(atmos, seed);
}

//! \copydoc GRAM::setMinRelativeStepSize_C()
void setMinRelativeStepSize_V(VenusAtmosphere_C* atmos, greal minRelativeStepSize)
{
  setMinRelativeStepSize_C(atmos, minRelativeStepSize);
}

//! \brief Set the perturbation scale factors.
//!
//! \param atmos An atmosphere handle.
//! \param densityScale Between 0 and 2.
//! \param ewWindScale Between 0 and 2.
//! \param nsWindScale Between 0 and 2.
void setPerturbationScales_V(VenusAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale)
{
  setPerturbationScales_C(atmos, densityScale, ewWindScale, nsWindScale, 0.0);
}

//! \copydoc GRAM::addAuxiliaryAtmosphere_C()
void addAuxiliaryAtmosphere_V(VenusAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive)
{
  addAuxiliaryAtmosphere_C(atmos, fileName, innerRadius, outerRadius, isEastLongitudePositive);
}

//! \copydoc GRAM::setAuxiliaryValues_C()
void setAuxiliaryValues_V(VenusAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns)
{
  setAuxiliaryValues_C(atmos, dens, pres, temp, ew, ns);
}

//! \copydoc GRAM::setPosition_C()
void setPosition_V(VenusAtmosphere_C* atmos, const Position_C* pos)
{
  setPosition_C(atmos, pos);
}

//! \copydoc GRAM::setDelta_C()
void setDelta_V(VenusAtmosphere_C* atmos, const Position_C* delta)
{
  setDelta_C(atmos, delta);
}

//! \copydoc GRAM::setPerturbationAction_C()
void setPerturbationAction_V(VenusAtmosphere_C* atmos, int action)
{
  setPerturbationAction_C(atmos, action);
}


//! \copydoc GRAM::setEphemerisState_C()
void setEphemerisState_V(VenusAtmosphere_C* atmos, const EphemerisState_C* state)
{
  setEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::setEphemerisFastModeOn_C()
void setEphemerisFastModeOn_V(VenusAtmosphere_C* atmos, int flag)
{
  setEphemerisFastModeOn_C(atmos, flag);
}

//! \copydoc GRAM::setSubsolarUpdateTime_C()
void setSubsolarUpdateTime_V(VenusAtmosphere_C* atmos, greal utime)
{
  setSubsolarUpdateTime_C(atmos, utime);
}


//! \brief Performs the Venus atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Then perturbations are computed prior to computing a few final metrics.
//! The current position must be set prior to calling this function.
//! \param atmos A Venus atmosphere handle.
//! \returns An error code. Zero is no error.  Non-zero signals an error.
int update_V(VenusAtmosphere_C* atmos)
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
size_t getErrorMessage_V(char *message, size_t bufferSize)
{
  return getErrorMessage_C(message, bufferSize);
}


//! \brief Get the current version as a string.
//!
//! \param atmos A Venus atmosphere handle.
//! \param bufferSize The size of the versionString buffer.
//! \param[out] versionString A non-NULL pointer.
//! \retval versionString A character array large enough to hold the version string.
//! \returns The length of the version string.
size_t getVersionString_V(VenusAtmosphere_C* atmos, char* versionString, size_t bufferSize)
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
void getPosition_V(VenusAtmosphere_C* atmos, Position_C* position)
{
  getPosition_C(atmos, position);
}

//! \copydoc GRAM::getDynamicsState_C()
void getDynamicsState_V(VenusAtmosphere_C* atmos, DynamicsState_C* state)
{
  getDynamicsState_C(atmos, state);
}

//! \copydoc GRAM::getDensityState_C()
void getDensityState_V(VenusAtmosphere_C* atmos, DensityState_C* state)
{
  getDensityState_C(atmos, state);
}

//! \copydoc GRAM::getWindsState_C()
void getWindsState_V(VenusAtmosphere_C* atmos, WindsState_C* state)
{
  getWindsState_C(atmos, state);
}

//! \copydoc GRAM::getGasesState_C()
//! \param[out] hydrogen A non-NULL pointer.
//! \param[out] oxygen A non-NULL pointer.
//! \param[out] helium A non-NULL pointer.
//! \param[out] nitrogen A non-NULL pointer.
//! \param[out] dinitrogen A non-NULL pointer.
//! \param[out] carbonMonoxide A non-NULL pointer.
//! \param[out] carbonDioxide A non-NULL pointer.
//! \retval hydrogen The atomic hydrogen (H) components.
//! \retval oxygen The atomic oxygen (O) components.
//! \retval helium The atomic helium (He) components.
//! \retval nitrogen The atomic nitrogen (N) components.
//! \retval dinitrogen The molecular dinitrogen (N2) components.
//! \retval carbonMonoxide The molecular carbon monoxide (CO) components.
//! \retval carbonDioxide The molecular carbon dioxide (CO2) components.
void getGasesState_V(VenusAtmosphere_C* atmos, GasesState_C* state,
  ConstituentGas_C* hydrogen, ConstituentGas_C* helium, ConstituentGas_C* oxygen,
  ConstituentGas_C* nitrogen, ConstituentGas_C* dinitrogen,
  ConstituentGas_C* carbonMonoxide, ConstituentGas_C* carbonDioxide)
{
  getGasesState_C(atmos, state);
  getConstituentGas_C(atmos, hydrogen, HYDROGEN);
  getConstituentGas_C(atmos, oxygen, OXYGEN);
  getConstituentGas_C(atmos, helium, HELIUM);
  getConstituentGas_C(atmos, nitrogen, NITROGEN);
  getConstituentGas_C(atmos, dinitrogen, DINITROGEN);
  getConstituentGas_C(atmos, carbonMonoxide, CARBON_MONOXIDE);
  getConstituentGas_C(atmos, carbonDioxide, CARBON_DIOXIDE);
}

//! \copydoc GRAM::getEphemerisState_C()
void getEphemerisState_V(VenusAtmosphere_C* atmos, EphemerisState_C* state)
{
  getEphemerisState_C(atmos, state);
}

//! \copydoc GRAM::getPerturbationState_C()
void getPerturbationState_V(VenusAtmosphere_C* atmos, PerturbationState_C* state)
{
  getPerturbationState_C(atmos, state);
}

//! \copydoc GRAM::getStartTime_C()
void getStartTime_V(VenusAtmosphere_C* atmos, GramTime_C* time)
{
  getStartTime_C(atmos, time);
}

//! @}


