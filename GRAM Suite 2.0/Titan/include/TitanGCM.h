//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "unittest_friend.h"
#include "Atmosphere.h"
#include "TitanCommon.h"
#include "TitanGCMThermos.h"
#include "TitanGCMMid.h"
#include "TitanGCMTerp.h"

namespace GRAM {

//! \brief The General Circulation Model of a Titan Atmosphere.
//!
//! The Titan GCM data used is from the model of Hourdin et al., Icarus, vol. 
//! 117, pp. 358 - 374 (1995). Compressibility effects of nitrogen are included
//! at lower altitude levels of the GCM data.
//! 
//! Upper altitudes for the GCM option are computed using a parameterized fit
//! to Titan exospheric temperatures, taken from graphs in Mueller - Wodarg "The 
//! Application of General Circulation Models to the Atmospheres of Terrestrial -
//! Type Moons of the Giant Planets", chapter in AGU monograph "Comparative
//! Atmospheres in the Solar System", M. Mendillo et al. editors, American 
//! Geophysical Union, 2002.
//! \ingroup TitanGRAM
class TitanGCM : public Atmosphere, public TitanCommon
{
public:
  TitanGCM();
  TitanGCM(const TitanGCM& orig) = default;
  virtual ~TitanGCM() override = default;

  void setMethaneMoleFraction(greal mmf) { 
    gcmTerp.setMethaneMoleFraction(mmf); 
    gcmMid.setMethaneMoleFraction(mmf);
    gcmThermos.setMethaneMoleFraction(mmf);
  }

  void update() override;

  greal getExosphericTemperature() const { return exosphericTemperature; }
  greal getPressureAtSurface() const { return 1.465e5; }

private:
  void updateExosphericTemperature();
  greal Jacchia(greal xlat, greal xlst, greal dec, greal exm, greal exn,
                greal R, greal R2, greal beta, greal gamma, greal p, greal Tc);

  TitanGCMTerp gcmTerp;         //!< \brief  For heights less than 0.3 mbar.
  TitanGCMMid gcmMid;           //!< \brief  For heights above 0.3 mbar and less than 925 km.
  TitanGCMThermos gcmThermos;   //!< \brief  For heights above 925 km.

  greal exosphericTemperature = 0; //!< \brief Depends on latitude, solar time, and longitude of the sun.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TitanGCM, updateExosphericTemperature);
  FRIEND_TEST(TitanGCM, Jacchia);
#endif // GRAM_UNIT_TEST
};

} // namespace

