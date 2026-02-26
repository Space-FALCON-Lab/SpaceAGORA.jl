//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Uranus atmosphere models.
//!
//! Common values for the Uranus atmosphere models. All Uranus
//! atmosphere models should inherit this class.
//!
//! \anchor UC1 (1) Data from 
//!     Lindal, G.F.; Lyons, J.R.; Sweetnam, D.N.; Eshleman, V.R.; Hinson, D.P.; and Tyler, G.L.:
//!     "The Atmosphere of Uranus: Results of Radio Occultation Measurements with Voyager 2",
//!     Journal of Geophysical Research, Vol. 92, No.A13, pp. 14, 987 - 15, 001, December 30, 1987.
//!
//! \ingroup UranusGRAM
class UranusCommon
{
public:
  UranusCommon(Atmosphere* atmos);
  virtual ~UranusCommon() = default;

private:
  Atmosphere& atmosphere;                  //!< \brief A reference to the parent Atmosphere object.

  const greal muDefault = 5.793964e15;            //!< \brief Gravitational constant \f$m^3/s^2\f$.  \ref UC1 "(1)"
  const greal equatorialRadiusDefault = 25559.0;  //!< \brief \f$km\f$.  \ref UC1 "(1)"
  const greal polarRadiusDefault = 24973.0;       //!< \brief \f$km\f$.  \ref UC1 "(1)"
  const greal J2Default = 0.00334129;             //!< \brief Gravitational coefficient.  \ref UC1 "(1)"
  const greal periodDefault = -62063.71199;       //!< \brief Sidereal period in seconds.  \ref UC1 "(1)"

  // Assume specific heat ratio = 1.45 (hydrogen-helium mixture)
  const greal specificHeatRatioDefault = 1.45;    //!< \brief \f$C_P/C_V\f$ (unitless).
  const greal amwDihydrogen = 2.02;        //!< Average molecular weight of hydrogen diatom (H2) on Uranus.
  const greal amwMethane = 16.04;          //!< Average molecular weight of methane (CH4) on Uranus.
  const greal amwHelium = 4.0;             //!< Average molecular weight of helium (HE) on Uranus.

};

} // namespace

