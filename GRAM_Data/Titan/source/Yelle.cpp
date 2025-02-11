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

#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "Yelle.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
Yelle::Yelle()
: MinMaxModel(), TitanCommon(this)
{
  initializeGases();
  initializeData();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
Yelle::Yelle(const Yelle& orig)
  : MinMaxModel(orig), TitanCommon(this)
{
  initializeGases();
}

//! \fn  Yelle::~Yelle()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void Yelle::update()
{
  // This model only works within the confines of the data.
  height = clamp(height, mdHeight()[0], mdHeight()[mdHeight().size() - 1]);

  updateMinMaxFactor();
  updateAtmosphereState();
  updateWinds();
  atmos.minMaxFactor = minMaxFactor;
}

//!  \brief Calculate the effects of season, latitude, and local solar time.
//!
//! If the #minMaxFactor has not been set by the user, then the factor is
//! updated based on its current value and input parameters.
//!
//! \b Inputs
//! \arg #latitude     
//! \arg #longitudeSun 
//! \arg #solarTime    
//! \arg #computeMinMaxFactor 
//! \arg #userMinMaxFactor   
//!
//! \retval #minMaxFactor 
void Yelle::updateMinMaxFactor()
{
  // Inputs are computeMinMaxFactor = true to compute longitudeSun, latitude, and solarTime effect,
  // or computeMinMaxFactor = false to use input value (userMinMaxFactor)
  // If evaluated, the input value massDensityFactor
  // is re-scaled by a factor of 1 minus magnitude of Ls-Lat-LST effect,
  // so that final calculated value cannot exceed range +/- 1
  if (computeMinMaxFactor) {
    greal alat = 0.46;  // *** SOURCE UNKNOWN ***
    greal alst = 0.2;   // *** SOURCE UNKNOWN ***
    minMaxFactor = alat * sin(toRadians(longitudeSun)) * sin(toRadians(latitude))
      - alst * cos(toRadians(15.0_deg * solarTime))
      + userMinMaxFactor * (1.0 - alat - alst);
  }
  else {
    minMaxFactor = userMinMaxFactor;
  }
}

//! \brief Compute zonal and meridional winds at height and latitude.
//!
//! Zonal wind model from ESA SP-1177, page 140.
//! Latitude variation of zonal wind from equation (2) page 143 of ESA SP-1177.
//!
//! \b Inputs
//! \arg #height  
//! \arg #latitude  
//!
//! \retval #ewWind 
//! \retval #nsWind 
void Yelle::updateWinds()
{
  // Compute zonal and meridional winds (u and v) at height and latitude

  // Zonal wind model from ESA SP-1177, page 140
  // zonal velocity at the equator (U sub e)
  greal U_e;
  if (height < 45.0) {
    U_e = 40.0 * pow(height / 45.0, 0.8);
  }
  else if (height < 100.0) {
    U_e = 40.0 + 20.0 * (height - 45.0) / 55.0;
  }
  else {
    U_e = 60.0 + 50.0 * (height - 100.0) / 120.0;
  }
  if (U_e > 130.0) {
    U_e = 130.0;
  }

  // Latitude variation of zonal wind from equation (2) page 143 of ESA SP-1177
  greal coslat = cos(toRadians(latitude));
  // Latitude variation exponent (2 / R_i - 1)   *** SOURCE UNKNOWN ***
  const greal xpon = 0.5;
  // Omega_a in m/s (page 147)
  const greal Omega_a = 11.7;
  ewWind = (U_e + Omega_a) * pow(coslat, xpon) - Omega_a * coslat;
  ewWind = max(ewWind, (greal)0.0);
  nsWind = 0.0;
}

} // namespace