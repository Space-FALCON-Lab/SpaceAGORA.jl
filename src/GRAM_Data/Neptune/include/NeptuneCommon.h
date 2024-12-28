//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Neptune atmosphere models.
//!
//! Common values for the Neptune atmosphere models. All Neptune
//! atmosphere models should inherit this class.
//!
//! \anchor NC1 (1) Data from NSSDCA Planetary Fact Sheet for Neptune (https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html).
//! Page last updated on Sept 27, 2018.
//!
//! \anchor NC2 (2) Data from NAIF Spice files (nep081.cmt).
//! \ingroup NeptuneGRAM
class NeptuneCommon
{
public:
  NeptuneCommon(Atmosphere* atmos);
  virtual ~NeptuneCommon() = default;

  void useDinitrogen();

private:
  Atmosphere& atmosphere;                  //!< \brief A reference to the parent Atmosphere object.

  const greal muDefault = 6.835099502439672e15;   //!< \brief Gravitational constant \f$m^3/s^2\f$.  \ref NC2 "(2)"
  const greal equatorialRadiusDefault = 24764.0;  //!< \brief \f$km\f$.  \ref NC1 "(1)"
  const greal polarRadiusDefault = 24341.0;       //!< \brief \f$km\f$.  \ref NC1 "(1)"
  const greal J2Default = 0.003411;               //!< \brief Gravitational coefficient.  \ref NC1 "(1)"
  const greal periodDefault = 57996.0;            //!< \brief Sidereal period in seconds.  \ref NC1 "(1)"
  
  // Assume specific heat ratio = 1.45 (hydrogen-helium mixture)
  const greal specificHeatRatioDefault = 1.45;    //!< \brief \f$C_P/C_V\f$ (unitless).
  const greal amwDihydrogen = 2.02;        //!< \brief Average molecular weight of hydrogen diatom (H2) on Neptune.
  const greal amwMethane = 16.04;          //!< \brief Average molecular weight of methane (CH4) on Neptune.
  const greal amwHelium = 4.0;             //!< \brief Average molecular weight of helium (He) on Neptune.
  const greal amwDinitrogen = 28.0;        //!< \brief Average molecular weight of Nitrogen diatom (N2) on Neptune.
};

} // namespace

