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

namespace GRAM {

//! \brief The VIRA thermosphere model of a Venus Atmosphere.
//!
//! The VIRA thermosphere model of a Venus Atmosphere.
//! \ingroup VenusGRAM
class VenusIRAThermos : public Atmosphere, public VenusCommon
{
public:
  VenusIRAThermos();
  VenusIRAThermos(const VenusIRAThermos& orig) = default;
  virtual ~VenusIRAThermos() override = default;

  void update() override;

  void setBase(greal ht, const AtmosphereState& atm); 

private:
  greal getNumberDensity(greal baseNumberDensity, greal molecularWeight);
  void updateAtmosphereState();

  greal baseHeight = 0.0;      //!< Height \units{km} of the base atmosphere state.
  AtmosphereState baseAtmos;   //!< The base atmosphere state.

  greal y = 0.0;   //!< Height parameter used in partial pressure computations.
  greal gf = 0.0;  //!< Surface gravity.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(VenusIRAThermos, getNumberDensity);
  FRIEND_TEST(VenusIRAThermos, updateAtmosphereState);
#endif // GRAM_UNIT_TEST
};

} // namespace

