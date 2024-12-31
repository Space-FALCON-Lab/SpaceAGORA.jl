//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//
// Adapted from Venus-GRAM 2005 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <cmath>
#include <algorithm>
#include "VenusAtmosphere.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
VenusAtmosphere::VenusAtmosphere()
: VenusCommon(this)
{
  gramBody = VENUS;
  venusIRAPtr = new VenusIRA();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
VenusAtmosphere::VenusAtmosphere(const VenusAtmosphere& orig)
  : PerturbedAtmosphere(orig), VenusCommon(this)
{
  gramBody = VENUS;
  inputParameters = orig.inputParameters;
  venusIRAPtr = new VenusIRA(*orig.venusIRAPtr);
}

//! \copydoc Atmosphere::~Atmosphere()
VenusAtmosphere::~VenusAtmosphere()
{
  if (venusIRAPtr != NULL) {
    delete venusIRAPtr;
  }
}

//! \copydoc Atmosphere::getVersionString()
const std::string& VenusAtmosphere::getVersionString()
{
  static const string version = "VenusGRAM 2023a :: " + Atmosphere::getVersionString();
  return version;
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void VenusAtmosphere::setInputParameters(const VenusInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;

  // Apply input parameters pertinent to the base Atmosphere.
  PerturbedAtmosphere::setInputParameters(inputParameters);
  setSeed(params.initialRandomSeed);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(params);

  // Notify ephemeris who we are.
  ephemeris.setBody(VENUS);
  ephemeris.setFastModeOn(params.fastModeOn);
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
void VenusAtmosphere::update()
{
  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Compute the surface height
  surfaceHeight = topography.getTopographicHeight(latitude, longitude);

  // Update ephemeris values.
  updateEphemeris();

  // Set the current position and ephemeris in the VIRA model.
  venusIRAPtr->setPosition(position);
  venusIRAPtr->setEphemerisState(ephem);

  // Compute the atmosphere state.
  venusIRAPtr->update();
  atmos = venusIRAPtr->getAtmosphereState();

  // Update the atmosphere state with any auxilary atosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Apply perturbations.
  updatePerturbations();

  // Update metrics that depend on atmosphere state.
  updateMetrics();
}

//! \copydoc Atmosphere::updatePressureAtSurface()
void VenusAtmosphere::updatePressureAtSurface()
{
  // Get the pressure from the VIRA model.
  pressureAtSurface = venusIRAPtr->getPressureAtSurface();
}

//! \copydoc Atmosphere::updateReferenceValues()
void VenusAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the VIRA model.
  venusIRAPtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}

//! \copydoc Atmosphere::updateCompressibilityFactor()
void VenusAtmosphere::updateCompressibilityFactor()
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
void VenusAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
{
  // Compute high and low factors for perturbations: parameterized    
  // from data in Figs. 1-8(d) and 1-12(a) of Seiff et al. VIRA data, 
  // pp. 247, 259 & 278 of "Venus", pp. 200 & 201 of "The Planet      
  // Venus", p. 283. of "Venus II", and Fig. 6 of Hinson and Jenkins  
  // Icarus 114, 310-327 (1995).                                      
  if (height < 50.0) {
    pertHigh = 1.02;
  }
  else {
    pertHigh = 1.02 + 0.0013 * (height - 50.0);
  }
  pertHigh = min(pertHigh, (greal)1.15);
  pertLow = 1.0 / pertHigh;
}

//! \brief Get vertical and horizontal scale parameters.
//!
//! \param[out] verticalScale A greal.
//! \param[out] horizontalScale A greal.
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void VenusAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
{
  // Set vertical and horizontal scale parameters                     
  // Vertical scale approximates VLS/Hrho = 2                         
  verticalScale = max(10.0, 32.0 - 22.0 * height / 65.0);
  if (height > 150.0) {
    verticalScale = 10.0 + 0.6 * (height - 150.0);
  }
  // Horizontal scale selected to be fairly consistent with wavelength
  // estimates from Fig. 6 of Bougher and Borucki, J. Geophys. Res.   
  // 99(E2), 3759-3776 (1994), Fig. 2 of Mayr et al., J. Geophys. Res.
  // 93(A10), 11247-262 (1988), Kasprzak et al. J. Geophys. Res.      
  // 93(A10), 11237-246 (1988), and Kasprzak et al. Geophys. Res.     
  // Lett. 20 2755-2758 (1993)                                        
  horizontalScale = 60.0 * verticalScale;
}

//! \brief Get wind standard deviations.
//!
//! \param[out] ewStdDev A greal.
//! \param[out] nsStdDev A greal.
//! \param[out] vertStdDev A greal.
//! \retval ewStdDev The east/west wind standard deviation.
//! \retval nsStdDev The north/south wind standard deviation.
//! \retval vertStdDev The vertical wind standard deviation.
void VenusAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  // Standard deviation for wind perturbations, from approximation to
  // VIRA model                                                      
  greal sigma = 0.0;
  if (height <= 30.0) {
    sigma = 1.0 + 7.0 * height / 30.0;
  }
  else if (height <= 45.0) {
    sigma = 8.0;
  }
  else if (height <= 60.0) {
    sigma = 8.0 + 7.0 * (height - 45.0) / 15.0;
  }
  else if (height <= 160.0) {
    sigma = 9.0 + height / 10.0;
  }
  else {
    sigma = 25.0;
  }
  ewStdDev = sigma * ewWindPerturbationScale;
  nsStdDev = 0.5 * sigma * nsWindPerturbationScale;
  vertStdDev = 0.0;
}

} // namespace
