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
#include "gram.h"
#include "Position.h"
#include "MarsCommon.h"

namespace GRAM {

//! \brief The Mars Slope Winds model.
//!
//! Analytical slope winds solution from Ye, Segal, and Pielke, J.Atmos.Sci., 47, 612, 1990.  (hereafter YSP)
//! To use this class, set all parameters and current winds. Call update() to perform the computations.
//! Then get the modified winds.
//!
//! \ingroup MarsGRAM
class SlopeWindsModel
{
public:
  SlopeWindsModel();
  SlopeWindsModel(const SlopeWindsModel& orig) = delete;
  virtual ~SlopeWindsModel() = default;

  void setPosition(const Position& pos) { position = pos; }
  void setSolarTime(greal time) { solarTime = time; }
  void setTemperature(greal temp) { temperature = temp; }
  void setSpeedOfSound(greal sos) { speedOfSound = sos; }
  void setMeanWindsScale(greal wms) { meanWindsScale = wms; }
  void setBoundaryLayerWindsScale(greal blw) { boundaryLayerWindsScale = blw; }
  void setWinds(greal ew, greal ns) { ewWind = ew; nsWind = ns; }
  void setDailyAverageWinds(greal ew, greal ns) { ewWindDailyAverage = ew; nsWindDailyAverage = ns; }
  void setAreoidRadiusCallback(TopoCallback callback) { getAreoidRadiusCallback = callback; }
  void setTopographicHeightCallback(TopoCallback callback) { getTopographicHeightCallback = callback; }
  void setCallbackData(void* dataPointer) { callbackDataPointer = dataPointer; }

  void update();

  void getWinds(greal& ew, greal& ns, greal& vert) const { ew = ewWind; ns = nsWind; vert = verticalWind; }
  void getDailyAverageWinds(greal& ew, greal& ns) const { ew = ewWindDailyAverage; ns = nsWindDailyAverage; }

private:
  void updateWinds();
  void slopeWinds(greal time, greal& ewSlopeWinds, greal& nsSlopeWinds, greal& verticalSlopeWinds);
  void updateSlopes();
  greal cp(greal t);

  // Outputs
  greal ewWind = 0.0;             //!< The east-west winds.
  greal nsWind = 0.0;             //!< The north-south winds.
  greal verticalWind = 0.0;       //!< The vertical winds.
  greal ewWindDailyAverage = 0.0; //!< The east-west average winds.
  greal nsWindDailyAverage = 0.0; //!< The north-south average winds.

  // Internal
  greal nsSlope = 0.0;  //!< The north-south slope (inclination).
  greal ewSlope = 0.0;  //!< The east-west slope (inclination).

  // Inputs
  greal meanWindsScale = 0.0;          //!< User supplied scaling factor.
  greal boundaryLayerWindsScale = 0.0; //!< User supplied scaling factor.
  greal solarTime = 0.0;               //!< Local solar time in hours.
  greal temperature = 0.0;             //!< Current temperature \units{K}.
  greal speedOfSound = 0.0;            //!< \units{m/s}.
  Position position;                   //!< Current position.

  TopoCallback getAreoidRadiusCallback = nullptr;      //!< Override of MOLATopography::getAreoidRadius.
  TopoCallback getTopographicHeightCallback = nullptr; //!< Override of MOLATopography::getTopographicHeight.
  void* callbackDataPointer = nullptr;                 //!< Data pointer for overrides.

  // Convenience variables (Position)
  greal& height         = position.height;          //!< \copydoc Position::height
  greal& latitude       = position.latitude;        //!< \copydoc Position::latitude
  greal& longitude      = position.longitude;       //!< \copydoc Position::longitude
  greal& latitudeRadius = position.latitudeRadius;  //!< \copydoc Position::latitudeRadius
  greal& surfaceHeight  = position.surfaceHeight;   //!< \copydoc Position::surfaceHeight

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(SlopeWindsModel, updateSlopes);
  FRIEND_TEST(SlopeWindsModel, slopeWinds);
#endif // GRAM_UNIT_TEST
};

} // namespace

