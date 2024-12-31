//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#include "gram.h"
#include "Position.h"
#include "MarsInputParameters.h"
#include "MarsCommon.h"

namespace GRAM {

//! \brief The base class for Mars dust models.
//!
//! This is the base class for Mars dust models. This model will compute the
//! dust optical depth (tau) based on position, longitude of the sun, and dust
//! storm conditions.  Any subclass of this model must implement updateDustOpticalDepth().
//! After the optical depth is calculated, this base class will add in a dust storm
//! intensity offset.  Dust storms are defined by user inputs.
//!
//! \ingroup MarsGRAM
class MarsDustModelBase : public MarsCommon
{
public:
  MarsDustModelBase();
  MarsDustModelBase(const MarsDustModelBase& orig) = default;
  virtual ~MarsDustModelBase() = default;

  virtual void setInputParameters(const MarsInputParameters& params);

  void setLongitudeSun(greal lonSun) { longitudeSun = lonSun; }
  void setPosition(const Position& pos) { position = pos; }

  virtual void update();

  greal getDustOpticalDepth() const { return dustOpticalDepth; }
  greal getDustOffset() const;

protected:
  virtual void updateDustOpticalDepth() = 0;
  virtual void updateDustIntensity();

  // Inputs
  Position position;             //!< Current position.
  greal longitudeSun = 0;        //!< Current Ls value \units{\text{degrees}}.
  greal stormLongitudeSun = 0.0; //!< Starting Ls value \units{\text{degrees}} for dust storm.
  greal stormDuration = 48.0;    //!< Duration (in Ls degrees) for dust storm.
  greal stormIntensity = 0.0;    //!< Dust storm intensity (0.0 - 3.0).
  greal stormMaxRadius = 0.0;    //!< Maximum radius \units{km} of the dust storm (0 or >10000 = global).
  greal stormLatitude = 0.0;     //!< Latitude \units{\text{degrees}} for center of dust storm.
  greal stormLongitude = 0.0;    //!< East positive longitude \units{\text{degrees}} for center of dust storm.

  // Outputs
  greal dustOpticalDepth = 0.0;  //!< Optical depth (tau)

  // Internal
  greal intensity = 0.0;         //!< Intensity factor at current position.
};

} // namespace
