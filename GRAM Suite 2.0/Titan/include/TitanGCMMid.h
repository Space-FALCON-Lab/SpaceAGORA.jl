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
#include "TitanCommon.h"
#include "TitanGCMTerp.h"

namespace GRAM {

//! \brief The GCM mid model of a Titan Atmosphere.
//!
//! This model is for heights over the 0.3 mbar level and below 925 km.  
//! \ingroup TitanGRAM
class TitanGCMMid : public Atmosphere, public TitanCommon
{
public:
  TitanGCMMid(TitanGCMTerp& gcmTerp);
  TitanGCMMid(const TitanGCMMid& orig) = default;
  virtual ~TitanGCMMid() override = default;

  void setMethaneMoleFraction(greal mmf) { gcmTerp.setMethaneMoleFraction(mmf); }

  void update() override;

  void setExosphericTemperature(greal tex) { exosphericTemperature = tex; }

private:
  void getMoleFractions(greal z, greal z03mb, greal mfMethane03mb, greal mfArgon03mb, greal& mfMethane, greal& mfArgon);
  greal getMiddleAtmosphereTemperature(greal z, greal z03mb, greal T03mb, greal Texos);

  TitanGCMTerp& gcmTerp;  //!< \brief  For heights less than 0.3 mbar.

  // Inputs
  greal exosphericTemperature = 0; //!< \brief Depends on latitude, solar time, and longitude of the sun.

   // Constants
  const greal maxHeight = 925.0;  // Maximum allowed height (km), base of the exosphere.


};

} // namespace
