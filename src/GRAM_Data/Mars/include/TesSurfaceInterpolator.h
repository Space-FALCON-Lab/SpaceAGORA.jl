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
#include "TesLowerInterpolator.h"

namespace GRAM {

//! \brief The TES surface (boundary layer) model (0 to 0.03 km).
//!
//! The TesSurfaceInterpolator uses interpolation routines to evaluate the atmospheric
//! state at the current position.  Results are based on data from the 
//! NASA Ames Mars General Circulation Model (MGCM) surface data ranging from 0 to 0.03 km.
//! This class requires the data file "TES_surface_data.bin".  In order to access the surface
//! level MGCM lower atmosphere data, the class is a subclass of TesLowerInterpolator.
//!
//! \ingroup MarsGRAM
class TesSurfaceInterpolator : public TesLowerInterpolator
{
public:
  TesSurfaceInterpolator();
  TesSurfaceInterpolator(const TesSurfaceInterpolator& orig) = default;
  virtual ~TesSurfaceInterpolator() = default;

  void reset() { temperatureAtFiveMeters = 0.0; }
  void updateIndices(size_t heightIndexOffset) override;
  void update() override;

  greal getUpperFairingHeight() const override { return (oneKilometerIndex < 16) ? -6.0 + oneKilometerIndex + 1 : 10.0 + 5.0 * (oneKilometerIndex + 1 - 16.0); }
  greal getLowerFairingHeight() const override { return surfaceHeight + surfaceHeights[SURFACE_HEIGHT_SIZE - 1]; }
  greal getLowestBoundaryLayer() const override { return surfaceHeights[1]; }

protected:
  static const size_t SURFACE_HEIGHT_SIZE = 3; //!< The number of height layers in the MGCM surface data (0, 0.005, 0.03 km).
  static const size_t SURFACE_LON_SIZE = 41;   //!< Tne number of longitude layers in the MGCM surface data (0 to 360 by 9 degrees).

  MarsTideParameters getTideParameters(const greal data[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE], size_t latIndex, size_t lonIndex, greal poleFactor);
  using TesLowerInterpolator::getTideParameters;
  void initializeData();

#ifdef CONVERT_MARS_DATA
  void readLegacySurfaceData();
  void writeSurfaceData();
#endif

  size_t oneKilometerIndex = 0;       //!< Height index of surface + 1 km.
  greal temperatureAtFiveMeters = 0;  //!< Temperature at the surface + 5 meters.


  static bool isInitialized;  //!< \brief True after data initialization.
  static greal surfaceHeights[SURFACE_HEIGHT_SIZE]; //!< Surface height layers.
  static greal surfaceTemperatures[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE];  //!< Temperature tidal parameters.
  static greal surfaceEWWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE];       //!< EW wind tidal parameters.
  static greal surfaceNSWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE];       //!< NS wind tidal parameters.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TesSurfaceInterpolator, updateIndices);
  FRIEND_TEST(TesSurfaceInterpolator, initializeSurfaceData_getTideParameters);
  FRIEND_TEST(TesSurfaceInterpolator, update_0);
  FRIEND_TEST(TesSurfaceInterpolator, update_1);
#endif // GRAM_UNIT_TEST
};

} // namespace

