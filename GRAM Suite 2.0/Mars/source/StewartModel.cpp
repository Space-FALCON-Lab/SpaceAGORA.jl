//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include <algorithm>
#include <cmath>
#include "StewartModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
StewartModel::StewartModel()
  : Atmosphere(), MarsCommon(this)
{
}

//! \fn  StewartModel::StewartModel(const StewartModel& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  StewartModel::~StewartModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void StewartModel::setInputParameters(const MarsInputParameters& params)
{
  F107 = params.F107;
  exosphericTemperatureFactor = params.exosphericTemperatureFactor;
}

//! \copydoc Atmosphere::update()
void StewartModel::update()
{
  updateExosphericTemperature();
  updateThermos();
  updateMoleFractions();
}


//! \brief Computes the local exospheric temperature.                                                           //
//!
//! Computes the exosphere temperature for the current latitude, and time.  This value depends
//! on the mean 10.7 cm solar flux at 1 AU which is converted to the Mars orbit.  The computed 
//! value can be adjusted using exosphericTemperatureFactor and exosphericTemperatureOffset.
//!
//! \b Inputs
//! \arg #latitude      
//! \arg #orbitalRadius      
//! \arg #subsolarLatitude      
//! \arg #solarTime
//! \arg #F107
//! \arg #exosphericTemperatureOffset
//! \arg #exosphericTemperatureFactor
//!
//! \returns #exosphericTemperature
void StewartModel::updateExosphericTemperature() {
  // Many undocumented/unsourced constants

  // maximum sub-solar latitude (axial tilt angle)
  constexpr greal LAT_MAX = 25.4_deg;

  // mean F10.7 flux at Mars orbit radius
  greal fbarr = F107 / (orbitalRadius * orbitalRadius);

  // global mean exosphere temperature
  greal tbar = 156.3 + 0.9427 * fbarr;

  // zonally-averaged, latitude, Ls-dependent exospheric temperature
  greal lat_ls = subsolarLatitude * latitude;
  greal tave = tbar * (1.0 + 1.369e-4 * lat_ls);

  // phase angles for local solar time variations of exosphere T
  greal t1 = 13.2 - 0.00119 * lat_ls;
  greal t2 = 9.4 - 0.00231 * lat_ls;

  // adjustment factor for pole proximity
  // factor = 1, not near pole (-85 to 85 degrees latitude)
  // factor = 0 at pole, falls to 0 as latitude approaches poles
  greal poleAmplitudeFactor = 1.0;
  if (abs(latitude) > 85.0_deg) {
    // if within 5 degrees of a pole, adjust poleAmplitudeFactor
    poleAmplitudeFactor -= (abs(latitude) - 85.0_deg) / 5.0_deg;
  }

  // amplitude factor for local solar time variations
  greal cphi = poleAmplitudeFactor * cos(toRadians((subsolarLatitude + latitude) / (1.0 + LAT_MAX / 90.0_deg)));

  // Exospheric temperature un-corrected for Stewart (ES array) variations,
  // pressure and dust effects.
  greal initialExoTemperature = tave * (
    1.0 + 0.22 * cphi * cos(toRadians(15.0_deg * (solarTime - t1)))
    + 0.04 * cphi * cos(toRadians(30.0_deg * (solarTime - t2)))
    );
  exosphericTemperature = initialExoTemperature * exp(0.16 * exosphericTemperatureFactor) + exosphericTemperatureOffset;
}

//! \brief Updates the atmospheric state.
//!
//!  Updates the atmospheric state of the thermosphere.
//!
//! \b Inputs
//! \arg #position
//! \arg #solarTime
//! \arg #orbitalRadius
//! \arg #F107
//! \arg #exosphericTemperature
//! \arg #thermosphereBaseHeight
//! \arg #thermosphereBaseTemperature
//! \arg gas.averageMolecularWeight
//!
//! \returns The AtmosphereState is populated.
void StewartModel::updateThermos() {
  // height above Zf
  greal heightAboveThermoBase = height - thermosphereBaseHeight; 
  // total planetary radius to Zf height
  greal thermosphereBaseRadius = latitudeRadius + thermosphereBaseHeight;  

  // mean F10.7 flux at Mars orbit radius
  greal fbarr = F107 / (orbitalRadius * orbitalRadius);
  // zonal average (latitude-dependent) T scale height (km) in thermosphere
  greal temperatureScaleHeight = (8.38 + 0.09725 * fbarr) * (1.14 - 0.18 * cos(toRadians(latitude)));   

  // local values to avoid duplicate computations
  greal tinf_tf = exosphericTemperature - thermosphereBaseTemperature;
  // scaling parameter for current position
  greal ysc = heightAboveThermoBase * thermosphereBaseRadius / (thermosphereBaseRadius + heightAboveThermoBase);
  greal exps = exp(-ysc / temperatureScaleHeight);
  // scaling parameter for 1 km higher than current position
  greal yscAbove = (heightAboveThermoBase + 1.0) * thermosphereBaseRadius / (thermosphereBaseRadius + heightAboveThermoBase + 1.0);

  // accelleration of gravity at Zf
  greal thermosphereBaseGravity = 100.0 * getGravity(latitude, latitudeRadius, thermosphereBaseHeight); 
  // temperature gradient in thermosphere (K/km)
  greal temperatureGradient = tinf_tf * (thermosphereBaseRadius + ysc) / (temperatureScaleHeight * (thermosphereBaseRadius + heightAboveThermoBase)) * exps;

  // T in thermosphere at current height
  temperature = exosphericTemperature - tinf_tf * exps;  
  // T at 1 km  above current height
  greal temperatureAbove = exosphericTemperature - tinf_tf * exp(-yscAbove / temperatureScaleHeight);

  /////////////////////////////////////////////////////////////////////////////////////////////
  // NOTE: Pressure values are in bar in the following code block.
  //       Mass and density values are in grams and grams/cm^3.
  /////////////////////////////////////////////////////////////////////////////////////////////
  // Pressure at thermosphere base height (ZF) in bar.
  greal thermoBasePressure = 1.26e-9;

  // Partial pressure of constituent gases
  greal pf_he = 3.3e-16 *  exosphericTemperature;     // helium partial pressure at Zf
  greal pf_h2 = 2.4e-15;                              // molecular hydrogen partial pressure at Zf
  greal pf_h;                                         // atomic hydrogen partial pressure at Zf
  if (exosphericTemperature <= 330.0) {
    // Compute pf_h by method 1
    pf_h = 5.2e-16 * exosphericTemperature * exp(-exosphericTemperature / 70.0);
  }
  else {
    // Compute pf_h by method 2
    greal ratio = 1440.0 / exosphericTemperature;
    pf_h = 5.8e-18 * sqrt(exosphericTemperature) * exp(ratio) / (1.0 + ratio);
  }

  // constituent concentrations
  carbonDioxide.moleFraction = 0.932;
  dinitrogen.moleFraction = 0.027;
  argon.moleFraction = 0.016;
  dioxygen.moleFraction = 0.002;
  carbonMonoxide.moleFraction = 0.013;
  greal FO = 0.01 * exp(0.40 * exosphericTemperatureFactor);
  oxygen.moleFraction = FO * (1.0 - 0.18 * sin(toRadians(15.0_deg * solarTime)) * cos(toRadians(latitude)));
  helium.moleFraction = pf_he / thermoBasePressure;
  dihydrogen.moleFraction = pf_h2 / thermoBasePressure;
  hydrogen.moleFraction = pf_h / thermoBasePressure;

 //=================================================================================================//
 //  loop over constituents and compute individual/summed quantities                                //
 //-------------------------------------------------------------------------------------------------//

  // Initialize sums to 0
  // molecular weights at current position weighted by partial pressures
  greal molecularWeightSum = 0.0;
  // total pressure at current position in bars
  greal pressureSum = 0.0;
  // mass density at current position in g/m^3
  greal densitySum = 0.0;
  // molecular weights 1 km higher than current position weighted by partial pressures
  greal molecularWeightSumAbove = 0.0;
  // total pressure 1 km higher than current position in bars
  greal pressureSumAbove = 0.0;

  // Iterate over the set of all activated gases.
  for (auto& gasIter : gases) {
    ConstituentGas& gas = gasIter.second;
    if (gas.isPresent) {
      // Molecular mass of the gas in grams
      greal dm = GRAMS_PER_AMU * gas.averageMolecularWeight;
      // scale height for current constituent
      greal hh = 100.0 *  BOLTZMANN * exosphericTemperature / (thermosphereBaseGravity * dm);

      // partial pressure of the constituent at the current height
      greal partialPressure = thermoBasePressure * gas.moleFraction
        * exp((-ysc / hh) - (temperatureScaleHeight / hh) * log(temperature / thermosphereBaseTemperature));
      // number density for current constituent
      gas.numberDensity = partialPressure / (10.0 * BOLTZMANN * temperature);
      // mass density of current constituent in g/m^3
      greal constituentDensity = gas.numberDensity * dm;

      // add current mass density to running total
      densitySum += constituentDensity;
      // add current partial pressure to running total
      pressureSum += partialPressure;
      // add molecular weight to running total
      molecularWeightSum += partialPressure * gas.averageMolecularWeight;

      // partial pressure 1 km higher for current constituent
      double partialPressureAbove = thermoBasePressure * gas.moleFraction
        * exp((-yscAbove / hh) - (temperatureScaleHeight / hh) * log(temperatureAbove / thermosphereBaseTemperature));
      // add partial pressure 1 km higher to running total
      pressureSumAbove += partialPressureAbove;
      // add molecular weight 1 km higher to running total
      molecularWeightSumAbove += partialPressureAbove * gas.averageMolecularWeight;
    }
  }

  //=================================================================================================//
  //  compute final quantities and return values to calling routine                                  //
  //-------------------------------------------------------------------------------------------------//
  // limit T >= sublimation T for CO2
  temperature = max(temperature, getCO2SublimationTemperature(pressureSum));
  averageMolecularWeight = molecularWeightSum / pressureSum;
  totalNumberDensity = pressureSum / (10.0 * BOLTZMANN * temperature);

  // Convert from bar to Pascal.
  pressure = 1.0e5 * pressureSum;
  // Convert from g/cm^3 to kg/m^3
  density = 1.0e3 * densitySum;

  gravity = getGravity(latitude, latitudeRadius, height);
  // Specific gas constant over gravity  (km/K)
  greal rOverG = UNIVERSAL_GAS / (averageMolecularWeight * gravity);
  // change in molecular weight over 1 km height
  greal dmdz = molecularWeightSumAbove / pressureSumAbove - averageMolecularWeight;

  pressureScaleHeight = rOverG * temperature;
  densityScaleHeight = pressureScaleHeight 
    / (1.0 + temperatureGradient * rOverG - dmdz * pressureScaleHeight / averageMolecularWeight); 

}

//! \brief Computes CO2 sublimation temperature.
//!
//! Computes CO2 sublimation temperature (K) as a funciton of P (mb) 
//! Taken from Kieffer and Jakosky, "Mars", 1992, U. of Ariz. Press, p.459 
//! \param pres Pressure.
//! \retval The CO2 sublimation temperature.
greal StewartModel::getCO2SublimationTemperature(greal pres)
{
  return 3182.48 / (23.3494 - log(pres / 100.0));
}

//! \fn  StewartModel::setThermosphereBase(greal height, greal temp)
//! \brief Set the thermosphere base height and temperature.
//!
//! Set the height and temperature for the thermosphere base occurring at the 1.26 nbar pressure level.
//! \param height The thermosphere base height.
//! \param temp The thermosphere base temperature.

//! \fn  StewartModel::setExosphericTemperatureOffset(greal offset)
//! \brief Set exospheric temperature offset.
//!
//! The exospheric temperature offset will be added to the initial (uncorrected) exospheric temperature.
//! The formula used is \f$  \mathrm{exosphericTemperature} = \mathrm{initialExoTemperature}\: e^{(0.16\ \mathrm{exosphericTemperatureFactor})} + \mathrm{exosphericTemperatureOffset}\f$.
//! \param offset The exospheric temperature offset \units{K}.

//! \fn  StewartModel::setExosphericTemperatureFactor(greal factor)
//! \brief Set exospheric temperature offset.
//!
//! The exospheric temperature factor will be applied to the initial (uncorrected) exospheric temperature.
//! The formula used is \f$  \mathrm{exosphericTemperature} = \mathrm{initialExoTemperature}\: e^{(0.16\ \mathrm{exosphericTemperatureFactor})} + \mathrm{exosphericTemperatureOffset}\f$.
//! \param factor The exospheric temperature factor (-3.0 to 3.0).

//! \fn  StewartModel::getExosphericTemperature()
//! \brief Get the exospheric temperature for the current position.
//!
//! The exospheric temperature is returned for the current latitude and solar time.
//! \returns The exospheric temperature \units{K}.

} // namespace


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//

  // This was computed but never used in the legacy code.  (untested)
  ////=================================================================================================//
  ////  compute local height of 1.26 nbar level, z_f0                                                  //
  ////-------------------------------------------------------------------------------------------------//
  //  greal zbar = 197.94 - 49.058 * radius;                                                           // global average height of Zf pressure level
  //  greal latitudeVariationFactor = (subsolarLatitude / LAT_MAX) * pow(latitude / 77.5, 3);          // intermediate latitude factor
  //  greal z_ave = zbar + 4.3 * latitudeVariationFactor;                                              // zonal average (latitude-dependent) height of Zf pressure level
  //  greal a1 = (1.5 - cos(TO_RADIANS * 4.0 * latitude)) * poleAmplitudeFactor;                       // amplitude A1 for local solar time variations of Zf
  //  greal a2 = (2.3 * pow(cos(TO_RADIANS * (latitude + 0.5 * subsolarLatitude)), 3)) * poleAmplitudeFactor;  // amplitude A2 for local solar time variations of Zf
  //	t1 = 16.2 - (subsolarLatitude / LAT_MAX) * atan(TO_RADIANS * 10.0 * latitude);                   // phase angle T1 (hours) for local solar time variations of Zf
  //	t2 = 11.5;                                                                                       // phase angle T2 (hours) for local solar time variations of Zf
  //	greal ZF0 = z_ave + a1 * cos(TO_RADIANS * 15.0 * (solarTime - t1)) +                             // compute current Zf height (km) for base of thermosphere
  //		a2 * cos(TO_RADIANS * 30.0 * (solarTime - t2));                                                // cont.
  ////=================================================================================================//
  ////  compute local T at base of thermosphere, t_f0                                                  //
  ////-------------------------------------------------------------------------------------------------//
  //	tbar = 113.7 + 0.5791 * fbarr;                                                                   // global average T at base of thermosphere (Zf)
  //	tave = tbar * (1.0 + 0.186 * latitudeVariationFactor);                                           // zonal-average (latitude-dependent) T at base of thermosphere
  //	a1 = (0.06 - 0.05 * cos(TO_RADIANS * 4.0 * latitude)) * poleAmplitudeFactor;                     // amplitude A1 for local solar time variations of T_f0
  //	a2 = (0.1 * pow(cos(TO_RADIANS*(latitude + 0.5* subsolarLatitude)), 3)) * poleAmplitudeFactor;   // amplitude A2 for local solar time variations of T_f0
  //	t1 = 17.5 - 2.5 * (subsolarLatitude / LAT_MAX)*atan(TO_RADIANS*10.0* latitude);                  // phase angle T1 for local solar time variations of T_f0
  //	t2 = 10.0 + 2.0 * pow(latitude / 77.5, 2);                                                       // phase angle T2 for local solar time variations of T_f0
  //	greal TF0 = tave * (1.0 + a1 * cos(TO_RADIANS*15.0 * (solarTime - t1)) +                         // thermospheric base T (K) at current position at LTST
  //		a2 * cos(TO_RADIANS * 30.0 * (solarTime - t2)));                                               // cont.

