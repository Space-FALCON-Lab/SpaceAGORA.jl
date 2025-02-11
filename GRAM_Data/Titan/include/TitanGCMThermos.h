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
#include "TitanGCMMid.h"

namespace GRAM {

//! \brief The GCM Thermosphere model of a Titan Atmosphere.
//!
//! Simplified themospheric model.  Following Yelle et al. (ESA
//! SP-1177) temperature is assumed to be constant with height (at
//! exospheric temperature) above height zf=925 km.  Constant
//! temperature is actually reached only at about 1250-1300 km (see
//! Mueller-Wodarg et al., J. Geophys. Res., 105 (A9), 20,833-
//! 20,856, 2000.  Constant exospheric temperatures assumed here
//! and their dependence on latitude, time-of-day, and Ls angle (as
//! described in subroutine TexosGCM_T04) are consistent with
//! Figures 8 and 9 of the Mueller-Wodard JGR 2000 article.
//! \ingroup TitanGRAM
class TitanGCMThermos : public Atmosphere, public TitanCommon
{
public:
  TitanGCMThermos(TitanGCMTerp& gcmTerp);
  TitanGCMThermos(const TitanGCMThermos& orig) = default;
  virtual ~TitanGCMThermos() override = default;

  void setExosphericTemperature(greal tex) { exosphericTemperature = tex; }
  void setMethaneMoleFraction(greal mmf) { gcmMid.setMethaneMoleFraction(mmf); }

  void update() override;

private:
  TitanGCMMid gcmMid;

  // Inputs
  greal exosphericTemperature = 0;

  // Constants
  const greal baseHeight = 925.0;

};

} // namespace
