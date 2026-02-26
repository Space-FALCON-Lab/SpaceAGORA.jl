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

#include <cmath>
#include <algorithm>
#include "VenusIRAThermos.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
VenusIRAThermos::VenusIRAThermos()
: Atmosphere(), VenusCommon(this)
{
}

//! \fn  VenusIRAThermos::VenusIRAThermos(const VenusIRAThermos& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  VenusIRAThermos::~VenusIRAThermos()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Set the base height and atmosphere state.
//!
//! This functions sets a base height and atmosphere state needed by the update() function.
//!
//! \param ht The base height \units{km}.
//! \param atm The base AtmosphereState.
void VenusIRAThermos::setBase(greal ht, const AtmosphereState& atm)
{
  baseHeight = ht;
  baseAtmos = atm; 
}

//! \brief Interface for the primary atmosphere computations.
//!
//! This routine computes the atmospheric state for the current position.
//! The height and AtmosphereState at the base of the thermosphere must be set
//! using setBase() prior to calling this function.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//! \arg baseHeight The base height must be set proior to calling this function.
//! \arg baseAtmos The base AtmosphereState must be set proior to calling this function.
//!
//! \returns The AtmosphereState is populated.
void VenusIRAThermos::update()
{
  // Get the average molecular weight at 1 km above the current height.
  greal saveHeight = height;
  height += 1.0;
  updateAtmosphereState();
  greal amwHigh = averageMolecularWeight;

  // Restore the current height and compute the atmospheric state
  height = saveHeight;
  updateAtmosphereState();

  // Get density scale height from pressure scale height and dM/dz  
  densityScaleHeight = pressureScaleHeight 
    / (1.0 - pressureScaleHeight * (amwHigh - averageMolecularWeight) / averageMolecularWeight);

}

//! \brief An internal function for the primary atmosphere computations.
//!
//! This routine computes the atmospheric state for the current position.
//! This is an internal routine that is called by the update() function.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//! \arg baseHeight The base height must be set proior to calling this function.
//! \arg baseAtmos The base AtmosphereState must be set proior to calling this function.
//!
//! \returns The AtmosphereState is populated.
void VenusIRAThermos::updateAtmosphereState()
{
  // Special Venus thermosphere model asssuming constant exospheric   
  // temperature above height zf (250 km)                             
  // Evaluates pressure (presz), density (densz), molecular weight    
  // (AMz), scale height (hscale), and number densities yxxx for 7    
  // species, given conditions at height zf (= 250 km), for altitude  
  // (z) and latitude (clat)  

  // Height above 250 km
  greal deltaHeight = height - baseHeight;

  // Surface radius (Rsfc) and gravity (gf) at height zf
  greal Rsfc = getRadius(latitude);
  gf = getGravity(latitude, Rsfc, baseHeight);

  // Radius at height zf
  greal Rf = Rsfc + baseHeight;
  // Height parameter
  y = -deltaHeight * Rf / (Rf + deltaHeight);

  // Assume constant temperature versus height
  temperature = baseAtmos.temperature;

  // Get scale height, partial pressure, and number density of each constituent
  // Total pressure is also computed as the sum of the partial pressures.
  pressure = 0.0;
  carbonDioxide.numberDensity = getNumberDensity(baseAtmos.carbonDioxide.numberDensity, carbonDioxide.averageMolecularWeight);
  carbonMonoxide.numberDensity = getNumberDensity(baseAtmos.carbonMonoxide.numberDensity, carbonMonoxide.averageMolecularWeight);
  oxygen.numberDensity = getNumberDensity(baseAtmos.oxygen.numberDensity, oxygen.averageMolecularWeight);
  dinitrogen.numberDensity = getNumberDensity(baseAtmos.dinitrogen.numberDensity, dinitrogen.averageMolecularWeight);
  nitrogen.numberDensity = getNumberDensity(baseAtmos.nitrogen.numberDensity, nitrogen.averageMolecularWeight);
  helium.numberDensity = getNumberDensity(baseAtmos.helium.numberDensity, helium.averageMolecularWeight);
  hydrogen.numberDensity = getNumberDensity(baseAtmos.hydrogen.numberDensity, hydrogen.averageMolecularWeight);
  // Total pressure from partial pressures is done incrementally in calls above
  // pressure = ppco2 + ppn2 + ppo + ppco + pphe + ppn + pph;

  // Limit number densities to > 1 per km**3
  greal ymin = 1.0e-9;
  carbonDioxide.numberDensity = max(carbonDioxide.numberDensity, ymin);
  carbonMonoxide.numberDensity = max(carbonMonoxide.numberDensity, ymin);
  dinitrogen.numberDensity = max(dinitrogen.numberDensity, ymin);

  // Get total number density from constituent number densities
  totalNumberDensity = carbonDioxide.numberDensity + dinitrogen.numberDensity + oxygen.numberDensity + carbonMonoxide.numberDensity + helium.numberDensity + nitrogen.numberDensity + hydrogen.numberDensity;

  // Get mean molecular weight
  averageMolecularWeight = (carbonDioxide.numberDensity * carbonDioxide.averageMolecularWeight
    + dinitrogen.numberDensity * dinitrogen.averageMolecularWeight
    + oxygen.numberDensity * oxygen.averageMolecularWeight
    + carbonMonoxide.numberDensity * carbonMonoxide.averageMolecularWeight
    + helium.numberDensity * helium.averageMolecularWeight
    + nitrogen.numberDensity * nitrogen.averageMolecularWeight
    + hydrogen.numberDensity * hydrogen.averageMolecularWeight) / totalNumberDensity;

  // Get density from perfect gas law
  specificGasConstant = UNIVERSAL_GAS * 1000.0 / (averageMolecularWeight);
  density = pressure / (specificGasConstant * temperature);

  // Get gravity at height z
  greal Rz = getRadius(latitude);
  greal gz = getGravity(latitude, Rz, height);

  // Get pressure scale height
  pressureScaleHeight = UNIVERSAL_GAS * temperature / (averageMolecularWeight * gz);
}

//! \brief Computes partial pressure and number density.
//!
//! \param baseNumberDensity The base atmosphere number density for the gas.
//! \param molecularWeight  The molecular weight for th gas.
//! \returns The number density for the gas and adds the partial pressure to the total pressure.
greal VenusIRAThermos::getNumberDensity(greal baseNumberDensity, greal molecularWeight)
{
  // Partial pressures of constituents at height zf
  greal ppatzf = BOLTZMANN * baseNumberDensity * temperature;

  // Scale height
  greal h = UNIVERSAL_GAS * temperature / (molecularWeight * gf);

  // partial pressure
  greal pp = ppatzf * exp(y / h);

  // Add partial pressure total pressure.
  pressure += pp;

  // Return number density
  return pp / (BOLTZMANN * temperature);
}

} // namespace
