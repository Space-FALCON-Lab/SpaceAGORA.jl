//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "unittest_friend.h"
#include "Atmosphere.h"
#include "MarsCommon.h"
#include "MarsAtmosphereState.h"

namespace GRAM {

//! This structure contains displacements used in interpolation.
//! It also contains dampening factors used near the poles.
struct MIBDisplacements {
  greal lat = 0.0;       //!< Latitude displacement
  greal lon = 0.0;       //!< Longitude displacement
  greal ls = 0.0;        //!< Longitude of the sun displacement
  greal f107 = 0.0;      //!< F10.7 solar flux displacement
  greal od = 0.0;        //!< Optical depth displacement
  greal wlat = 0.0;      //!< Latitude displacement for winds
  greal wlon = 0.0;      //!< Longitude displacement for winds
  greal wpolefac = 0.0;  //!< Pole dampening factor for winds
  greal tpolefac = 0.0;  //!< Pole dampening factor
};

//! This structure contains indices used in data lookup.
struct MIBIndices {
  size_t hgt = 0;   //!< Height index
  size_t lat = 0;   //!< Latitude index
  size_t wlat = 0;  //!< Latitude index for winds
  size_t lon = 0;   //!< Longitude index
  size_t wlon = 0;  //!< Longitude index for winds
  size_t ls = 0;    //!< Longitude of the sun index
  size_t od = 0;    //!< Optical depth index
  size_t tyr = 0;   //!< TES year index
  size_t f107 = 0;  //!< F10.7 solar flux index
};

//! This structure contains the five tidal parameters needed to compute a tide value.
struct MarsTideParameters {
  greal diurnalMean;  //!< Diurnal mean value
  greal amplitude1D;  //!< Amplitude of once-daily tidal variation    
  greal phase1D;      //!< Phase of once-daily tidal variation        
  greal amplitude2D;  //!< Amplitude of twice-daily tidal variation   
  greal phase2D;      //!< Phase of twice-daily tidal variation       
};


//! \brief The base class for all Mars GCM interpolation models.
//!
//! The Mars GCM models consist of tidal data and multi-dimensional interpolations.
//! This is the base class for all of the GCM models.  This class defines the interface
//! to be implemented by all of the models.  It also contains methods which implement
//! the three tidal models.
//!
//! \ingroup MarsGRAM
class MarsInterpolatorBase: public Atmosphere, public MarsCommon
{
protected:

public:
  MarsInterpolatorBase();
  MarsInterpolatorBase(const MarsInterpolatorBase& orig) = default;
  virtual ~MarsInterpolatorBase() = default;

  virtual void updateIndices(size_t heightIndexOffset) = 0;

  // Height level methods
  virtual greal getUpperFairingHeight() const { return 0.0; }
  virtual greal getLowerFairingHeight() const { return 0.0; }
  virtual greal getLowestBoundaryLayer() const { return 0.0; }
  virtual greal getThermosphereBaseHeight() const { return 0.0; }
  virtual void setHeightOffset(greal offset) { }

  // Convenience methods (for testing)
  const MIBIndices& getBaseIndex() const { return baseIndex; }
  const MIBDisplacements& getDisplacements() const { return displacements; }
  const greal getPressureScaleHeightDaily() const { return pressureScaleHeightDaily; }

protected:
  bool computeMinMax = true;             //!< Daily min/max are computed when true.
  MarsAtmosphereState marsAtmos;         //!< Mars specific metrics added to the AtmosphereState.
  MIBIndices index;                      //!< Index counters for generating interpolation cubes.
  MIBIndices baseIndex;                  //!< The indices for the "lower left" corner of the interpolation cubes.
  MIBDisplacements displacements;        //!< The normalized offsets of the current point within the interpolation cube.
  greal pressureScaleHeightDaily = 0.0;  //!< Mean daily pressure scale height \units{km}.

  //! \brief The forms available for converting tidal coefficients to values.
  enum TideType { 
    UVT,     //!< Tidal form for winds and temperature
    DP,      //!< Tidal form for density and pressure
    BIG_DP   //!< Tidal form for density and pressure with assumed large amplitudes
  };
  
  greal getTideValue(const MarsTideParameters& p, greal ltst, TideType form);
  void getSolMinMax(const MarsTideParameters& tp, TideType form, greal& minValue, greal& maxValue);

  // Convenience References
  greal& dustOpticalDepth = marsAtmos.dustOpticalDepth;              //!< \copydoc MarsAtmosphereState::dustOpticalDepth
  greal& temperatureDaily = marsAtmos.temperatureDaily;              //!< \copydoc MarsAtmosphereState::temperatureDaily
  greal& pressureDaily = marsAtmos.pressureDaily;                    //!< \copydoc MarsAtmosphereState::pressureDaily
  greal& densityDaily = marsAtmos.densityDaily;                      //!< \copydoc MarsAtmosphereState::densityDaily
  greal& ewWindDaily = marsAtmos.ewWindDaily;                        //!< \copydoc MarsAtmosphereState::ewWindDaily
  greal& nsWindDaily = marsAtmos.nsWindDaily;                        //!< \copydoc MarsAtmosphereState::nsWindDaily
  greal& densityMin = marsAtmos.densityMin;                          //!< \copydoc MarsAtmosphereState::densityMin
  greal& densityMax = marsAtmos.densityMax;                          //!< \copydoc MarsAtmosphereState::densityMax
  greal& temperatureMin = marsAtmos.temperatureMin;                  //!< \copydoc MarsAtmosphereState::temperatureMin
  greal& temperatureMax = marsAtmos.temperatureMax;                  //!< \copydoc MarsAtmosphereState::temperatureMax
  greal& thermosphereBaseHeight = marsAtmos.thermosphereBaseHeight;  //!< \copydoc MarsAtmosphereState::thermosphereBaseHeight

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MarsInterpolatorBase, getTideValue_UVT);
  FRIEND_TEST(MarsInterpolatorBase, getTideValue_DP);
  FRIEND_TEST(MarsInterpolatorBase, getTideValue_BIG_DP);
  FRIEND_TEST(MarsInterpolatorBase, getSolMinMax);
#endif // GRAM_UNIT_TEST

};

void readBinaryDataBlock(std::ifstream& binaryStream, size_t dataSize, greal* dataBlock);
void writeBinaryDataBlock(std::ofstream& binaryStream, size_t dataSize, greal* dataBlock);

} // namespace

