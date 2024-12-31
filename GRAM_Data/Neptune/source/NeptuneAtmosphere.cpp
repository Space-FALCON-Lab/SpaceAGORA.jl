//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//
// Adapted from Neptune-GRAM 2004 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath>
#include "NeptuneAtmosphere.h"
#include "InputParameters.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
NeptuneAtmosphere::NeptuneAtmosphere()
: NeptuneCommon(this)
{
  gramBody = NEPTUNE;
  neptuneMinMaxPtr = new NeptuneMinMax();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
NeptuneAtmosphere::NeptuneAtmosphere(const NeptuneAtmosphere& orig)
: PerturbedAtmosphere(orig), NeptuneCommon(this)
{
  gramBody = NEPTUNE;
  inputParameters = orig.inputParameters;
  neptuneMinMaxPtr = new NeptuneMinMax(*orig.neptuneMinMaxPtr);
}

//! \copydoc Atmosphere::~Atmosphere()
NeptuneAtmosphere::~NeptuneAtmosphere()
{
  if (neptuneMinMaxPtr != NULL) {
    delete neptuneMinMaxPtr;
  }
}

//! \copydoc Atmosphere::getVersionString()
const std::string& NeptuneAtmosphere::getVersionString()
{
  static const string version = "NeptuneGRAM 2019d :: " + Atmosphere::getVersionString();
  return version;
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void NeptuneAtmosphere::setInputParameters(const NeptuneInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;

  // Set the common perturbation parameters
  PerturbedAtmosphere::setInputParameters(params);
  setSeed(params.initialRandomSeed);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(params);

  ephemeris.setFastModeOn(params.fastModeOn);

  // Set the Neptune specific parameters
  neptuneMinMaxPtr->setMinMaxFactor(params.minMaxFactor, params.computeMinMaxFactor);
  neptuneMinMaxPtr->setDinitrogenMoleFraction(params.dinitrogenMoleFraction/100.0);
}

//! \copydoc MinMaxModel::setMinMaxFactor()
void NeptuneAtmosphere::setMinMaxFactor(greal factor, bool computeFlag)
{
  neptuneMinMaxPtr->setMinMaxFactor(factor, computeFlag);
}

//! \brief Set the molecular nitrogen (N2) mole fraction.
//!
//! Diatomic nitrogen (N2) quantities will not be computed unless the 
//! mole fraction is set with this function prior to an update.
//! \param n2mf The mole fraction.
void NeptuneAtmosphere::setDinitrogenMoleFraction(greal n2mf)
{
  neptuneMinMaxPtr->setDinitrogenMoleFraction(n2mf);
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
void NeptuneAtmosphere::update()
{
  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Update ephemeris values.
  updateEphemeris();

  // Set the current position and ephemeris in the Neptune minmax model.
  neptuneMinMaxPtr->setPosition(position);
  neptuneMinMaxPtr->setEphemerisState(ephem);

  // Compute the atmosphere state.
  neptuneMinMaxPtr->update();
  atmos = neptuneMinMaxPtr->getAtmosphereState();
 
  // Update the atmosphere state with any auxilary atosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Apply perturbations.
  updatePerturbations();

  // Update metrics that depend on atmosphere state.
  updateMetrics();
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void NeptuneAtmosphere::updatePressureAtSurface()
{
  // Get the pressure from the Neptune minmax model.
  pressureAtSurface = neptuneMinMaxPtr->getPressureAtSurface();
}

//! \copydoc Atmosphere::updateReferenceValues()
void NeptuneAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the Neptune minmax model.
  neptuneMinMaxPtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}

//! \copydoc Atmosphere::updateCompressibilityFactor()
void NeptuneAtmosphere::updateCompressibilityFactor()
{
  // Compute as usual using the function in the Atmosphere base class.
  Atmosphere::updateCompressibilityFactor();
  // Force constant value over 500 km.
  if (height > 500.0) {
    compressibilityFactor = 1.0;
  }
}


//********************************* Perturbation Factor Methods ********************************//

//! \brief Get high and low perturbation factors.
//!
//! \param[out] pertLow A greal.
//! \param[out] pertHigh A greal.
//! \retval pertLow The low perturbation factor.
//! \retval pertHigh The high perturbation factor.
void NeptuneAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
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
//! \param[out] verticalScale A greal.
//! \param[out] horizontalScale A greal.
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void NeptuneAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
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

//! \brief Get wind standard deviations.
//!
//! \param[out] ewStdDev A greal.
//! \param[out] nsStdDev A greal.
//! \param[out] vertStdDev A greal.
//! \retval ewStdDev The east/west wind standard deviation.
//! \retval nsStdDev The north/south wind standard deviation.
//! \retval vertStdDev The vertical wind standard deviation.
void NeptuneAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  // Sigma U and Sigma V estimates from "Neptune and Triton", page 625
  greal sigmaU = 16.0 + 0.32 * height;
  if (sigmaU > 112.0) {
    sigmaU = 112.0;
  }
  ewStdDev = sigmaU * ewWindPerturbationScale;
  nsStdDev = 0.5 * sigmaU * nsWindPerturbationScale;
  if (ewWindPerturbationScale > 0.0 && ewStdDev < 25.0) {
    ewStdDev = 25.0;
  }
  if (nsWindPerturbationScale > 0.0 && nsStdDev < 25.0) {
    nsStdDev = 25.0;
  }
  vertStdDev = 0.0;
}

} // namespace
