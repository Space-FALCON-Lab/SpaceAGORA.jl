//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "gram.h"
#include "VenusCommon.h"
#include "Atmosphere.h"
#include "VenusIRALow.h"
#include "VenusIRAMid.h"
#include "VenusIRAHigh.h"
#include "VenusIRAThermos.h"

namespace GRAM {

//! \brief The VIRA model of a Venus Atmosphere.
//!
//! The VIRA model of a Venus Atmosphere.
//! \ingroup VenusGRAM
class VenusIRA : public Atmosphere, public VenusCommon
{
public:
  VenusIRA();
  VenusIRA(const VenusIRA& orig) = default;
  virtual ~VenusIRA() override = default;

  void update() override;

  void getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity);
  greal getPressureAtSurface() const;

private:
  void initializeData();
  void updateAtmosphereState();
  void updateWinds() override;
  void interpolateHeight(greal z, greal height1, greal height2,
    const AtmosphereState& atmos1, const AtmosphereState& atmos2);
  const AtmosphereState averageAtmospheres(const AtmosphereState& atmos1, const AtmosphereState& atmos2);

  VenusIRALow venusLow;             //!< The lower atmosphere model.
  VenusIRAMid venusMid;             //!< The middle atmosphere model.
  VenusIRAHigh venusHigh;           //!< The upper atmosphere model.
  VenusIRAThermos venusThermos;     //!< The thermosphere model.

  static bool isInitialized;        //!< Flags initialization of the reference profile.
  static std::vector<greal> zref;   //!< Heights of the reference profile.
  static std::vector<greal> pref;   //!< Profile of reference pressures.
  static std::vector<greal> tref;   //!< Profile of reference temperatures.
  static std::vector<greal> dref;   //!< Profile of reference densities.
};

} // namespace

