//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath>
#include "UranusAtmosphere.h"
#include "UranusModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
UranusAtmosphere::UranusAtmosphere()
: UranusCommon(this)
{
  gramBody = URANUS;
  uranusModelPtr = new UranusModel();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
UranusAtmosphere::UranusAtmosphere(const UranusAtmosphere& orig)
  : PerturbedAtmosphere(orig), UranusCommon(this)
{
  gramBody = URANUS;
  inputParameters = orig.inputParameters;
  uranusModelPtr = new UranusModel(*static_cast<UranusModel*>(orig.uranusModelPtr));
}

//! \copydoc Atmosphere::~Atmosphere()
UranusAtmosphere::~UranusAtmosphere()
{
  if (uranusModelPtr != NULL) {
    delete uranusModelPtr;
  }
}

//! \copydoc Atmosphere::getVersionString()
const std::string& UranusAtmosphere::getVersionString()
{
  static const string version = "UranusGRAM 2021b :: " + Atmosphere::getVersionString();
  return version;
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void UranusAtmosphere::setInputParameters(const UranusInputParameters& params)
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
void UranusAtmosphere::update()
{
  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Update ephemeris values.
  updateEphemeris();

  // Set the position in the Uranus model
  uranusModelPtr->setPosition(position);
  // Currently, no ephemeris data is needed in the Uranus model.
  //uranusModel.setEphemerisState(ephem);

  // Compute the atmosphere state.
  uranusModelPtr->update();
  atmos = uranusModelPtr->getAtmosphereState();

  // Update the atmosphere state with any auxilary atosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Apply perturbations.
  updatePerturbations();

  // Update metrics that depend on atmosphere state.
  updateMetrics();
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void UranusAtmosphere::updatePressureAtSurface()
{
  // Get the pressure from the Uranus height model.
  pressureAtSurface = uranusModelPtr->getPressureAtSurface();
}

//! \copydoc Atmosphere::updateReferenceValues()
void UranusAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the Uranus height model.
  uranusModelPtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}


//********************************* Perturbation Factor Methods ********************************//

//! \brief Get high and low perturbation factors.
//!
//! \retval pertLow The low perturbation factor.
//! \retval pertHigh The high perturbation factor.
void UranusAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
{
  //Compute high and low factors for perturbations
  if (height < 0.0) {
    pertHigh = 1.02;
  }
  else if (height < 80.0) {
    pertHigh = 1.02 + 0.13 * height / 80.0;
  }
  else if (height < 300.0) {
    pertHigh = 1.15 - 0.05 * (height - 80.0) / 220.0;
  }
  else if (height < 500.0) {
    pertHigh = 1.10;
  }
  else if (height < 550.0) {
    pertHigh = 1.10 + 0.035 * (height - 500.0) / 50.0;
  }
  else {
    pertHigh = 1.135;
  }
  pertLow = 1.0 / pertHigh;
}

//! \brief Get vertical and horizontal scale parameters.
//!
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void UranusAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
{
  // New scale model approximates VLS/Hrho = 2
  if (height < 0.0_km) {
    verticalScale = 128.0 - 38.0 * (height + 85.0) / 85.0;
  }
  else if (height < 35.0_km) {
    verticalScale = 90.0 - 56.0 * height / 35.0;
  }
  else if (height < 100.0_km) {
    verticalScale = 34.0;
  }
  else if (height < 150.0_km) {
    verticalScale = 34.0 + 46.0 * (height - 100.0) / 50.0;
  }
  else if (height < 330.0_km) {
    verticalScale = 80.0 + 30.0 * (height - 150.0) / 180.0;
  }
  else if (height < 550.0_km) {
    verticalScale = 110.0;
  }
  else if (height < 980.0_km) {
    verticalScale = 110.0 + 90.0 * (height - 550.0) / 430.0;
  }
  else {
    verticalScale = 200.0;
  }
  horizontalScale = 5.0 * verticalScale;
}

//! \brief Get east/west and north/south wind deviations.
//!
//! At this time, there is no winds model for Uranus.
//! \retval ewStdDev The east/west wind deviation.
//! \retval nsStdDev The north/south wind deviation.
void UranusAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  // We have no winds model for Uranus.  So zero out the deviations.
  ewStdDev = 0.0;
  nsStdDev = 0.0;
  vertStdDev = 0.0;
}

} // namespace
