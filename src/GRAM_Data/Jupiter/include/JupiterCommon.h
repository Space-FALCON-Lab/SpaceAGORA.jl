//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Jupiter atmosphere models.
//!
//! Common values for the Jupiter atmosphere models. All Jupiter
//! atmosphere models should inherit this class.
//!
//! \anchor JC1 (1) Data from NSSDCA Planetary Fact Sheet for Jupiter (https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html).
//! Page last updated on July 18, 2018.
//!
//! \anchor JC2 (2) Data from NAIF Spice files (jup310.cmt).
//! \ingroup JupiterGRAM
class JupiterCommon
{
public:
  JupiterCommon(Atmosphere* atmos);
  virtual ~JupiterCommon() = default;

private:
  Atmosphere& atmosphere;                   //!< \brief A reference to the parent Atmosphere object.

  const greal muDefault = 1.266865341960128e17;    //!< \brief Gravitational constant \f$m^3/s^2\f$.  \ref JC2 "(2)"
  const greal equatorialRadiusDefault = 71492.0;   //!< \brief \f$km\f$.  \ref JC1 "(1)"
  const greal polarRadiusDefault = 66854.0;        //!< \brief \f$km\f$.  \ref JC1 "(1)"
  const greal J2Default = 0.014736;                //!< \brief Gravitational coefficient.  \ref JC1 "(1)"
  const greal periodDefault = 35730.0;             //!< \brief Sidereal period in seconds.  \ref JC1 "(1)"

  // Using hydrogen
  const greal specificHeatRatioDefault = 1.41;     //!< \brief \f$C_P/C_V\f$ (unitless).
};

} // namespace

