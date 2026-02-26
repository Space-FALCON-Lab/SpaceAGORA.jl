//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath>
#include "JupiterAtmosphere.h"
#include "JupiterModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
JupiterAtmosphere::JupiterAtmosphere()
: JupiterCommon(this)
{
  gramBody = JUPITER;
  jupiterModelPtr = new JupiterModel();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
JupiterAtmosphere::JupiterAtmosphere(const JupiterAtmosphere& orig)
: PerturbedAtmosphere(orig), JupiterCommon(this)
{
  gramBody = JUPITER;
  inputParameters = orig.inputParameters;
  jupiterModelPtr = new JupiterModel(*static_cast<JupiterModel*>(orig.jupiterModelPtr));
}

//! \copydoc Atmosphere::~Atmosphere()
JupiterAtmosphere::~JupiterAtmosphere()
{
  if (jupiterModelPtr != NULL) {
    delete jupiterModelPtr;
  }
}

//! \copydoc Atmosphere::getVersionString()
const std::string& JupiterAtmosphere::getVersionString()
{
  static const string version = "JupiterGRAM 2021a :: " + Atmosphere::getVersionString();
  return version;
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void JupiterAtmosphere::setInputParameters(const JupiterInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;

  // Apply input parameters pertinent to the base Atmosphere.
  PerturbedAtmosphere::setInputParameters(inputParameters);
  setSeed(params.initialRandomSeed);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(params);
}

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
//! \returns The AtmosphereState is populated.
void JupiterAtmosphere::update()
{
  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Update ephemeris values.
  updateEphemeris();

  // Set the position in the Jupiter model
  jupiterModelPtr->setPosition(position);
  // Currently, no ephemeris data is needed in the Jupiter model.
  //jupiterModelPtr->setEphemerisState(ephem);

  // Compute the atmosphere state.
  jupiterModelPtr->update();
  atmos = jupiterModelPtr->getAtmosphereState();

  // Update the atmosphere state with any auxilary atosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Apply perturbations.
  updatePerturbations();

  // Update metrics that depend on atmosphere state.
  updateMetrics();
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void JupiterAtmosphere::updatePressureAtSurface()
{
  // Get the pressure from the Jupiter height model.
  pressureAtSurface = jupiterModelPtr->getPressureAtSurface();
}

//! \copydoc Atmosphere::updateReferenceValues()
void JupiterAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the Jupiter height model.
  jupiterModelPtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}


//********************************* Perturbation Factor Methods ********************************//

//! \brief Get high and low perturbation factors.
//!
//! \retval pertLow The low perturbation factor.
//! \retval pertHigh The high perturbation factor.
void JupiterAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
{
  pertHigh = 1.0;
  pertLow = 1.0;
}

//! \brief Get vertical and horizontal scale parameters.
//!
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void JupiterAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
{
  verticalScale = 1.0;
  horizontalScale = 1.0;
}

//! \brief Get east/west and north/south wind deviations.
//!
//! At this time, there is no winds model for Jupiter.
//! \retval ewStdDev The east/west wind deviation.
//! \retval nsStdDev The north/south wind deviation.
void JupiterAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  ewStdDev = 0.0;
  nsStdDev = 0.0;
  vertStdDev = 0.0;
}

} // namespace
