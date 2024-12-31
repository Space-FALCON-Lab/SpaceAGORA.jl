//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "unittest_friend.h"
#include "Atmosphere.h"
#include "MarsInputParameters.h"
#include "MarsCommon.h"

namespace GRAM {

//! \brief The StewartModel of the Mars atmosphere.
//!
//! The Stewart Model computes a time-dependent Mars atmosphere model, adapted from Ian Stewart, University of 
//! Colorado. Final Report JPL PO # NQ - 802429 
//!
//! \ingroup MarsGRAM
class StewartModel : public Atmosphere, public MarsCommon
{
public:
  StewartModel();
  StewartModel(const StewartModel& orig) = default;
  virtual ~StewartModel() = default;

  void setInputParameters(const MarsInputParameters& params);
  void setThermosphereBase(greal height, greal temp) { thermosphereBaseHeight = height; thermosphereBaseTemperature = temp; }
  void setExosphericTemperatureOffset(greal offset) { exosphericTemperatureOffset = offset; }
  void setExosphericTemperatureFactor(greal factor) { exosphericTemperatureFactor = factor; }

  void update() override;

  greal getExosphericTemperature() const { return exosphericTemperature; }

protected:
  void updateExosphericTemperature();
  void updateThermos();
  greal getCO2SublimationTemperature(greal pres);

  // inputs
  greal F107 = 68.0;                        //!< Solar flux at 10.7 cm in \units{\text{sfu}}.
  greal exosphericTemperatureFactor = 0.0;  //!< Ranges from -3.0 to 3.0.
  greal thermosphereBaseTemperature = 0.0;  //!< \units{K}.
  greal thermosphereBaseHeight = 0.0;       //!< \units{km}.
  greal exosphericTemperatureOffset = 0.0;  //!< \units{K}.

  // output
  greal exosphericTemperature = 0.0;        //!< \units{K}.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(StewartModel, updateExosphericTemperature);
  FRIEND_TEST(StewartModel, updateThermos);
#endif // GRAM_UNIT_TEST

};

} // namespace

