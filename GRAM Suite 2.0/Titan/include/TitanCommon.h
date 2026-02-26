//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Titan atmosphere models.
//!
//! Common values for the Titan atmosphere models. All Titan
//! atmosphere models should inherit this class.
//!
//! \anchor TC1 (1) Data from NSSDCA Planetary Fact Sheet for Titan (https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturniansatfact.html).
//! Page last updated on July 18, 2018.
//!
//! \anchor TC2 (2) Data from NAIF Spice files (jup310.cmt).
//!
//! \anchor TC3 (3) "Titan\'s rotation: A 3-dimensional theory", A\&A 478, 959 - 970 (2008).
//! \ingroup TitanGRAM
class TitanCommon
{
public:
  TitanCommon(Atmosphere* atmos);
  virtual ~TitanCommon() = default;

private:
  Atmosphere& atmosphere;                   //!< \brief A reference to the parent Atmosphere object.

  const greal muDefault = 8.97803e12;              //!< \brief Gravitational constant \f$m^3/s^2\f$.  \ref TC2 "(2)"
  const greal equatorialRadiusDefault = 2575.5;    //!< \brief \f$km\f$.  \ref TC1 "(1)"
  const greal polarRadiusDefault = 2575.5;         //!< \brief \f$km\f$.  \ref TC1 "(1)"
  const greal J2Default = 0.0000315;               //!< \brief Gravitational coefficient.  \ref TC3 "(3)"
  const greal periodDefault = 1377684.0;           //!< \brief Sidereal period in seconds.  \ref TC1 "(1)"
  
  // Assume specific heat ratio = 1.4 (nitrogen)
  const greal specificHeatRatioDefault = 1.4;      //!< \brief \f$C_P/C_V\f$ (unitless).

  const greal amwDinitrogen = 28.0;         //!< \brief Average molecular weight of Nitrogen diatom (N2) on Titan.
  const greal amwMethane = 16.0;            //!< \brief Average molecular weight of methane (CH4) on Titan.

  // Values consistent with Yelle model assuming primordial Ar(36) and Ar(38).
  // See page 206 of ESA SP-1177.
  const greal amwArgon = 36.4;              //!< \brief Average molecular weight of argon on Titan.
};

} // namespace

