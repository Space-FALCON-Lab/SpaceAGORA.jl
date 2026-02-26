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
#include "MarsCommon.h"
#include "Atmosphere.h"
#include "MarsInputParameters.h"
#include "MarsInterpolatorBase.h"
#include "MarsDustModelBase.h"

namespace GRAM {

//! \brief The base class for the MGCM and TesGCM models.
//!
//! MarsGCMBase is a base class for the MGCM and TesGCM models.  In order to define the models,
//! all pure virtual functions must be implemented in the sub-class.  This class performs the 
//! majority of the calculations for performing interpolation over the height dimension.
//!
//! \ingroup MarsGRAM
class MarsGCMBase : public Atmosphere, public MarsCommon
{
public:
  MarsGCMBase(MarsInterpolatorBase& sfc, MarsInterpolatorBase& low, MarsInterpolatorBase& high, MarsDustModelBase& dust);
  MarsGCMBase(const MarsGCMBase& orig) = delete;
  virtual ~MarsGCMBase() override = default;

  void update();

  greal getSurfaceFairingHeight() const { return marsSurface.getLowerFairingHeight(); }

  // Pure virtual functions (override in sub-class).
  virtual void setInputParameters(const MarsInputParameters& params) = 0;
  virtual void updateDustModel() = 0;
  virtual greal getMaxHeight() const = 0;

protected:
  // Pure virtual functions (override in sub-class).
  virtual greal getGlobalMeanOffset() = 0;
  virtual greal getSeasonalOffset() = 0;

  // Update methods
  void updateHeightOffsets();
  void updateFairingHeights();
  void updateUpperLower();
  void updateScaleHeights();
  void updateGasConstant();
  void updateAtmos();

  // Utility methods
  greal bltp(greal tg, greal z5, greal t5, greal u5, greal v5, greal zeval, greal factor, greal grav);
  greal getCO2SublimationTemperature(greal p);
  greal cp(greal t);

  // Set surface roughness parameter (km). Value 1 cm (1.0e-5 km) is
  // consistent with NASA Ames MGCM boundary layer model
  const greal minSurfaceRoughness = 1.0e-5;  //!< Surface roughness above CO2 sublimation temperature.
  const greal icySurfaceRoughness = 1.0e-7;  //!< Surface roughness below CO2 sublimation temperature.

  // Inputs
  MarsOffsetModel offsetModel = MARS_CONSTANT; //!< Height offset model (was IBOUGHER).
  bool computeMinMax = true;         //!< User input. When true, daily extrema are calculated for temperature and density.
  greal constantHeightOffset = 0.0;  //!< User input constant height offset \units{km}.
  greal userHeightAboveSurface = 0.0; //!< Used when height < -8.7 km.

  // Outputs
  greal currentOffset = 0.0;          //!< Height offset at current time and position \units{km}.

  // Internal
  MarsInterpolatorBase& marsSurface;  //!< The surface interpolator.
  MarsInterpolatorBase& marsLower;    //!< The lower atmosphere interpolator.
  MarsInterpolatorBase& marsUpper;    //!< The upper atmosphere interpolator.
  MarsDustModelBase& marsDustModel;   //!< The Mars dust model.

  greal dustOffset = 0.0;                //!< Adjustment to height offset due to dust \units{km}.
  greal upperFairingHeight = 0.0;        //!< Upper bound for smoothing between the marsLower and marsUpper models \units{km}.
  greal lowerFairingHeight = 0.0;        //!< Lower bound for smoothing between the marsLower and marsUpper models \units{km}.
  greal surfaceUpperFairingHeight = 0.0; //!< Upper bound for smoothing between the marsSurface and marsLower models \units{km}.
  greal surfaceLowerFairingHeight = 0.0; //!< Lower bound for smoothing between the marsSurface and marsLower models \units{km}.
  greal pressureScaleHeightDaily = 0.0;  //!< Mean daily pressure scale height \units{km}.
  greal specificGasConstantDaily = 0.0;  //!< Mean daily specific gas constant \units{J / (K * \text{mol})}.
  greal lowerHeight = 0.0;               //!< Height \units{km} of the lower interpolation bound of current layer.
  greal upperHeight = 0.0;               //!< Height \units{km} of the upper interpolation bound of current layer.

  AtmosphereState lower;           //!< The atmosphere state of the lower interpolation bound.
  MarsAtmosphereState lowerAtmos;  //!< Mars specific metrics added to the lower AtmosphereState.
  AtmosphereState upper;           //!< The atmosphere state of the upper interpolation bound.
  MarsAtmosphereState upperAtmos;  //!< Mars specific metrics added to the upper AtmosphereState.
  MarsAtmosphereState marsAtmos;   //!< Mars specific metrics added to the AtmosphereState.

  // Convenience References
  greal& temperatureDaily = marsAtmos.temperatureDaily;             //!< \copydoc MarsAtmosphereState::temperatureDaily
  greal& pressureDaily = marsAtmos.pressureDaily;                   //!< \copydoc MarsAtmosphereState::pressureDaily
  greal& densityDaily = marsAtmos.densityDaily;                     //!< \copydoc MarsAtmosphereState::densityDaily
  greal& ewWindDaily = marsAtmos.ewWindDaily;                       //!< \copydoc MarsAtmosphereState::ewWindDaily
  greal& nsWindDaily = marsAtmos.nsWindDaily;                       //!< \copydoc MarsAtmosphereState::nsWindDaily
  greal& densityMin = marsAtmos.densityMin;                         //!< \copydoc MarsAtmosphereState::densityMin
  greal& densityMax = marsAtmos.densityMax;                         //!< \copydoc MarsAtmosphereState::densityMax
  greal& temperatureMin = marsAtmos.temperatureMin;                 //!< \copydoc MarsAtmosphereState::temperatureMin
  greal& temperatureMax = marsAtmos.temperatureMax;                 //!< \copydoc MarsAtmosphereState::temperatureMax
  greal& heightOffset = marsAtmos.heightOffset;                     //!< \copydoc MarsAtmosphereState::heightOffset
  greal& localHeightOffset = marsAtmos.localHeightOffset;           //!< \copydoc MarsAtmosphereState::localHeightOffset
  greal& thermosphereBaseHeight = marsAtmos.thermosphereBaseHeight; //!< \copydoc MarsAtmosphereState::thermosphereBaseHeight
  greal& groundTemperature = marsAtmos.groundTemperature;           //!< \copydoc MarsAtmosphereState::groundTemperature
  greal& dustOpticalDepth = marsAtmos.dustOpticalDepth;             //!< \copydoc MarsAtmosphereState::dustOpticalDepth
  int& iceIsPresent = marsAtmos.iceIsPresent;                       //!< \copydoc MarsAtmosphereState::iceIsPresent       


#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MarsGCMBase, bltp);
  FRIEND_TEST(MarsGCMBase, updateGasConstant);
#endif // GRAM_UNIT_TEST

};

} // namespace

