//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//
// Adapted from Titan-GRAM 2004 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath>
#include "TitanAtmosphere.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanAtmosphere::TitanAtmosphere()
  : TitanCommon(this)
{
  gramBody = TITAN;
  yellePtr = new Yelle();
  titanGCMPtr = new TitanGCM();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
TitanAtmosphere::TitanAtmosphere(const TitanAtmosphere& orig)
  : PerturbedAtmosphere(orig), TitanCommon(this)
{
  gramBody = TITAN;
  inputParameters = orig.inputParameters;
  yellePtr = new Yelle(*orig.yellePtr);
  titanGCMPtr = new TitanGCM(*orig.titanGCMPtr);
  setModelType(orig.modelType);
}

//! \copydoc Atmosphere::~Atmosphere()
TitanAtmosphere::~TitanAtmosphere()
{
  if (yellePtr != NULL) {
    delete yellePtr;
  }
  if (titanGCMPtr != NULL) {
    delete titanGCMPtr;
  }
}

//! \brief Set the active model.
//!
//! The Titan models include Yelle97 and GCM95.  The model type should be set as part of the model
//! initialization process.
//! \param type The Titan model type.
void TitanAtmosphere::setModelType(TitanModelType type)
{
  modelType = type;

  switch (modelType) {
  case Yelle97:
    atmosModelPtr = yellePtr;
    break;
  case GCM95:
    atmosModelPtr = titanGCMPtr;
    break;
  }
}


//! \copydoc Atmosphere::getVersionString()
const std::string& TitanAtmosphere::getVersionString()
{
  static const string version = "TitanGRAM 2020c :: " + Atmosphere::getVersionString();
  return version;
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void TitanAtmosphere::setInputParameters(const TitanInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;

  // Apply input parameters pertinent to the base Atmosphere.
  PerturbedAtmosphere::setInputParameters(params);
  setSeed(params.initialRandomSeed);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(params);

  // Set the atmosModel.
  setModelType(inputParameters.modelType);
  if (inputParameters.modelType == Yelle97) {
    setMinMaxFactor(params.minMaxFactor, params.computeMinMaxFactor);
  }

  setMethaneMoleFraction(params.userMethaneMoleFraction);
}

//! \brief Overrides the models' methane mole fraction.
//!
//! The methane mole fraction computed by the models can be override using this method.
//! \param mmf The override methane mole fraction between 1.0 and 5.0 (percent).  No override if mmf is 0.0.
void TitanAtmosphere::setMethaneMoleFraction(greal mmf)
{
  if (mmf == 0) {
    userSuppliedMethaneMoleFraction = false;
    userMethaneMoleFraction = 0.0;
  }
  else {
    userMethaneMoleFraction = 0.01 * clamp(mmf, 1.0, 5.0);
    userSuppliedMethaneMoleFraction = true;
  }
  titanGCMPtr->setMethaneMoleFraction(userMethaneMoleFraction);
}

//! \copydoc MinMaxModel::setMinMaxFactor()
void TitanAtmosphere::setMinMaxFactor(greal factor, bool computeFlag)
{
  yellePtr->setMinMaxFactor(factor, computeFlag);
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
void TitanAtmosphere::update()
{
  if (atmosModelPtr == NULL) {
    atmosModelPtr = yellePtr;
  }
  Atmosphere& atmosModel = *atmosModelPtr;

  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Update ephemeris values.
  updateEphemeris();

  // Set the current position and ephemeris in the atmosModel.
  atmosModel.setPosition(position);
  atmosModel.setEphemerisState(ephem);

  // Compute and get the atmosphere state from Yelle or TitanGCM atmosModel
  atmosModel.update();
  atmos = atmosModel.getAtmosphereState();

  // User override of methane mole fraction
  if (userMethaneMoleFraction) {
    dinitrogen.moleFraction = (argon.averageMolecularWeight - averageMolecularWeight
           + userMethaneMoleFraction * (methane.averageMolecularWeight - argon.averageMolecularWeight))
          / (argon.averageMolecularWeight - dinitrogen.averageMolecularWeight);
    argon.moleFraction = (averageMolecularWeight - dinitrogen.averageMolecularWeight
           + userMethaneMoleFraction * (dinitrogen.averageMolecularWeight - methane.averageMolecularWeight))
          / (argon.averageMolecularWeight - dinitrogen.averageMolecularWeight);
    methane.moleFraction = userMethaneMoleFraction;
    if (argon.moleFraction <= 0.0) {
      argon.moleFraction = 0.0;
      dinitrogen.moleFraction
          = (methane.averageMolecularWeight - averageMolecularWeight)
            / (methane.averageMolecularWeight - dinitrogen.averageMolecularWeight);
      methane.moleFraction
          = (averageMolecularWeight - dinitrogen.averageMolecularWeight)
            / (methane.averageMolecularWeight - dinitrogen.averageMolecularWeight);
    }
    // Preserves mean molecular weight   
    // AMz = fmolN2*amn2 + fmolCH4 *amch4 + fmolAr *amar and total number density ytot
    dinitrogen.numberDensity = dinitrogen.moleFraction * totalNumberDensity;
    methane.numberDensity = methane.moleFraction * totalNumberDensity;
    argon.numberDensity = argon.moleFraction * totalNumberDensity;
    updateMassFractions();
  }

  // Update the atmosphere state with any auxilary atosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Apply perturbations.
  updatePerturbations();

  // Update metrics that depend on atmosphere state.
  updateMetrics();
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void TitanAtmosphere::updatePressureAtSurface()
{
  // Get the pressure from the Yelle or the TitanGCM atmosModel.
  if (modelType == Yelle97) {
    pressureAtSurface = yellePtr->getPressureAtSurface();
  }
  else {
    pressureAtSurface = titanGCMPtr->getPressureAtSurface();
  }
}

//! \copydoc Atmosphere::updateReferenceValues()
void TitanAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the Yelle atmosModel.
  yellePtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}


//********************************* Perturbation Factor Methods ********************************//

//! \brief Get high and low perturbation factors.
//!
//! \retval pertLow The low perturbation factor.
//! \retval pertHigh The high perturbation factor.
void TitanAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
{
  if (height < 60.0) {
    pertHigh = height * 0.002 + 1.03;
  }
  else if (height < 150.0) {
    pertHigh = 1.15 - (height - 60.0) * 0.05 / 90.0;
  }
  else {
    pertHigh = 1.1;
  }
  pertLow = 1.0 / pertHigh;
}

//! \brief Get vertical and horizontal scale parameters.
//!
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void TitanAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
{
  // New scale atmosModel approximates VLS/Hrho = 2
  if (height <= 350.0_km) {
    verticalScale = 40.0 + 0.2 * height;
  }
  else if (height < 700.0_km) {
    verticalScale = 110.0;
  }
  else {
    verticalScale = 110.0 + 0.2 * (height - 700.0);
  }
  horizontalScale = 5.0 * verticalScale;
}

//! \brief Get east/west and north/south wind deviations.
//!
//! \retval ewStdDev The east/west wind deviation.
//! \retval nsStdDev The north/south wind deviation.
void TitanAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  // Sigma U estimate from ESA SP-1177, page 292
  greal sigmaU = 7.0 + 0.14 * height;
  if (sigmaU > 45.0) {
    sigmaU = 45.0;
  }
  ewStdDev = sigmaU * ewWindPerturbationScale;
  nsStdDev = sigmaU * nsWindPerturbationScale;
  vertStdDev = 0.0;
}

} // namespace
