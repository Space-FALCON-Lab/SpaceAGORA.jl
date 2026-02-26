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
#include "TitanGCM.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanGCM::TitanGCM()
: TitanCommon(this), gcmTerp(), gcmMid(gcmTerp), gcmThermos(gcmTerp)
{
}

//! \fn  TitanGCM::TitanGCM(const TitanGCM& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TitanGCM::~TitanGCM()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void TitanGCM::update()
{
  // Compute altitude profiles from GCM data (surface to 0.3 mbar)
  gcmTerp.setPosition(position);
  gcmTerp.setEphemerisState(ephem);
  gcmTerp.updateGCMTables();
  greal z03mb = gcmTerp.getHeightAt03mb();

  if (height < z03mb) {
    // Compute atmosphere at height z from GCM data if 0<z<z03mb
    gcmTerp.update();
    atmos = gcmTerp.getAtmosphereState();
  }
  else if (height > 925.0) {
    // Get exospheric temperature Texos for heights >= 925 km
    updateExosphericTemperature();
    // Compute atmosphere in simplified thermosphere (constant
    // exospheric temperature) above 925 km
    gcmThermos.setPosition(position);
    gcmThermos.setEphemerisState(ephem);
    gcmThermos.setExosphericTemperature(exosphericTemperature);
    gcmThermos.update();
    atmos = gcmThermos.getAtmosphereState();
  }
  else {
    // Get exospheric temperature Texos for heights >= 925 km
    updateExosphericTemperature();
    // Compute atmosphere for middle height range z03mb to 925 km
    gcmMid.setPosition(position);
    gcmMid.setEphemerisState(ephem);
    gcmMid.setExosphericTemperature(exosphericTemperature);
    gcmMid.update();
    atmos = gcmMid.getAtmosphereState();
  }

  updateMoleFractions();
  updateMassFractions();
}

//! \brief Titan exospheric temperature
//!
//! Computed by Jacchia function, with
//! parameters fit to exospheric temperature graphs in Mueller-
//! Wodarg "The Application of General Circulation Models to the
//! Atmospheres of Terrestrial-Type Moons of the Giant Planets",
//! chapter in AGU monograph "Comparative Atmospheres in the Solar
//! System", M. Mendillo et al. editors, American Geophysical Union,
//! 2002.
//!
//! \b Inputs
//! \arg #latitude      
//! \arg #solarTime      
//! \arg #longitudeSun      
//!
//! \retval #exosphericTemperature
void TitanGCM::updateExosphericTemperature()
{
  // Jacchia function parameters fit to Mueller-Wodarg exospheric
  // temperature graph for Ls = 0
  greal Tex000 = Jacchia(latitude, solarTime, 0.0_deg, 2.2, 2.4, 0.14, 0.24,
          -5.0_deg, 0.0_deg, -25.0_deg, 157.0);

  // Assume same parameters with lat=-lat for Ls = 180
  greal Tex180 = Jacchia(-latitude, solarTime, 0.0_deg, 2.2, 2.4, 0.14, 0.24,
          -5.0_deg, 0.0_deg, -25.0_deg, 157.0);

  // Jacchia function parameters fit to Mueller-Wodarg exospheric
  // temperature graph for Ls = 270
  greal Tex270 = Jacchia(latitude, solarTime, -26.43_deg, 1.25, 3.7, 0.077,
          0.064, -5.0_deg, 0.0_deg, -25.0_deg, 147.5);

  // Assume same parameters with dec=-dec for Ls = 90
  greal Tex090 = Jacchia(latitude, solarTime, 26.43_deg, 1.25, 3.7, 0.077,
          0.064, -5.0_deg, 0.0_deg, -25.0_deg, 147.5);

  // Amplitudes for sines and cosines of assumed Ls variation
  greal Acoef = 0.25 * (Tex000 + Tex090 + Tex180 + Tex270);
  greal Bcoef = 0.5 * (Tex090 - Tex270);
  greal Ccoef = 0.5 * (Tex000 - Tex180);
  greal Dcoef = 0.25 * (Tex000 + Tex180 - Tex090 - Tex270);

  // Exospheric temperature at given Ls
  greal lsr = toRadians(longitudeSun);
  exosphericTemperature = Acoef + Bcoef * sin(lsr)
          + Ccoef * cos(lsr) + Dcoef * cos(2.0 * lsr);
}

//! \brief Empirical function for exospheric temperature
//!
//! Originally developed for Earth thermosphere by Jacchia (Smithsonian
//! Astrophysical Observatory, Special Report 313, 1970). Input
//! parameters [derived to fit Titan exospheric temperature versus xlat and xlst].
//! \param xlat latitude
//! \param xlst time of day
//! \param dec  Titan solar declination angle (deg)
//! \param exm  latitude term exponent
//! \param exn  time-of-day exponent
//! \param R
//! \param R2
//! \param beta (deg)
//! \param p (deg)
//! \param gamma (deg) time offset and time assymetry parameters for solar hour angle (tau, deg)
//! \param Tc  minimum exospheric temperature (K) over all latitudes and times of day.
//! \returns  exopsheric temperature (K).
greal TitanGCM::Jacchia(greal xlat, greal xlst, greal dec, greal exm, greal exn,
                greal R, greal R2, greal beta, greal gamma, greal p, greal Tc)
{
  // Solar hour angle, with time offset from solar noon and shape
  // parameter with amplitude p.
  greal tau = 15.0_deg * (xlst - 12.0_hr) + beta
            + p * sin(toRadians(15.0_deg * (xlst - 12.0_hr) + gamma));
  tau = wrapDegrees180(tau);

  // Latitude sine term
  greal slat = pow(sin(toRadians(0.5 * fabs(xlat + dec))), exm);
  // Factor for time-of-day dependence
  greal Afact = pow(cos(toRadians(0.5 * fabs(xlat - dec))), exm) - slat;
  // Factor for latitude dependence
  greal Cfact = 1.0 + R2 * slat;
  // Computed exospheric temperature for latitude, local time, and
  // solar declination (time of year)
  greal Tinf = Tc * (Cfact + R * Afact * pow(cos(toRadians(0.5 * tau)), exn));
  return Tinf;
}

} // namespace
