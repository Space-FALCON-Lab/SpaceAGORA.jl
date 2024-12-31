//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Venus atmosphere models.
//!
//! Common values for the Venus atmosphere models. All Venus
//! atmosphere models should inherit this class.
//!
//! \anchor VC1 (1) Data from NSSDCA Planetary Fact Sheet for Venus (https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html).
//! Page last updated on Sept. 27, 2018.
//!
//! \anchor VC2 (2) Data from NAIF Spice files (de430_tech-comments.txt).
//! \ingroup VenusGRAM
class VenusCommon
{
public:
  VenusCommon(Atmosphere* atmos);
  virtual ~VenusCommon() = default;

private:
  Atmosphere& atmosphere;                   //!< \brief A reference to the parent Atmosphere object.

  const greal muDefault = 3.24858592e14;           //!< \brief Gravitational constant \f$m^3/s^2\f$.  \ref VC2 "(2)"
  const greal equatorialRadiusDefault = 6051.8;    //!< \brief \f$km\f$.  \ref VC1 "(1)"
  const greal polarRadiusDefault = 6051.8;         //!< \brief \f$km\f$.  \ref VC1 "(1)"
  const greal J2Default = 4.458e-6;                //!< \brief Gravitational coefficient.  \ref VC1 "(1)"
  const greal periodDefault = -20997360.0;         //!< \brief Sidereal period in seconds.  \ref VC1 "(1)"

  // Assume specific heat ratio = 1.286 (CO2 - N2 mixture)
  const greal specificHeatRatioDefault = 1.286;

  const greal amwHydrogen = 1.01;        //!< Average molecular weight of hydrogen atom (H) on Venus.
  const greal amwCarbonDioxide = 44.0;   //!< Average molecular weight of carbon dioxide (CO2) on Venus.
  const greal amwHelium = 4.0;           //!< Average molecular weight of helium (HE) on Venus.
  const greal amwDinitrogen = 28.0;      //!< Average molecular weight of Nitrogen diatom (N2) on Venus.
  const greal amwOxygen = 16.0;          //!< Average molecular weight of Oxygen atom (O) on Venus.
  const greal amwCarbonMonoxide = 28.0;  //!< Average molecular weight of carbon monoxide (CO) on Venus.
  const greal amwNitrogen = 14.0;        //!< Average molecular weight of Nitrogen atom (N) on Venus.

};

} // namespace

