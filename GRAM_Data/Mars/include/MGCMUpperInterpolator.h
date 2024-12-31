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
#include "MGCMInterpolator.h"

namespace GRAM {

//! \brief The MGCM upper atmosphere model (80 to 170 km).
//!
//! The MGCMUpperInterpolator uses interpolation routines to evaluate the atmospheric
//! state at the current position.  Results are based on data from the University
//! of Michigan Mars Thermospheric General Circulation Model (MTGCM) ranging from 80 to 170 km.
//! This class requires the data file "MGCM_upper_data.bin".
//!
//! \ingroup MarsGRAM
class MGCMUpperInterpolator : public MGCMInterpolator
{
public:
  MGCMUpperInterpolator();
  MGCMUpperInterpolator(const MGCMUpperInterpolator& orig) = default;
  virtual ~MGCMUpperInterpolator() = default;

  void setHeightOffset(greal offset) override { heightOffset = offset; }
  void setF107(greal value) { F107 = value; }

  greal getHeightOffset() const { return heightOffset; }
  greal getF107() const { return F107; }

  void updateIndices(size_t heightIndexOffset) override;
  void update() override;

  greal getUpperFairingHeight() const override { return 80.0 + heightOffset + 5.0 * heightOffsetAdjustment; }
  greal getLowerFairingHeight() const override { return 75.0; }

protected:
  static const size_t MTGCM_HEIGHT_SIZE = 19; //!< The number of height layers in the MTGCM data (80 to 170 by 5 km).
  static const size_t MTGCM_LAT_SIZE = 36;    //!< The number of latitude layers in the MTGCM data (-87.7 to 87.5 by 5 degrees).
  static const size_t MTGCM_F107_SIZE = 2;    //!< The number of F10.7 layers at 1AU in the MTGCM data (70 and 130).

  MarsTideParameters getTideParameters(const greal data[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE], greal poleFactor);
  MarsTideParameters getTideParameters(const greal data[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE], greal poleFactor);
  void initializeData();

#ifdef CONVERT_MARS_DATA
  void writeUpperData();
  void readLegacyUpperData();
#endif

  greal heightOffset = 0;            //!< Height offset \units{km}.
  int heightOffsetAdjustment = 0;    //!< Adjustment for very negative offsets (was lhgtt)
  greal F107 = 68.0;                 //!< Solar activity value, F10.7 at 1 AU.

  static bool isInitialized; //!< \brief True after data initialization.
  static greal mtgcmTemperature[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE]; //!< Temperature tidal parameters.
  static greal mtgcmPressure[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE];    //!< Pressure tidal parameters.
  static greal mtgcmDensity[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE];     //!< Density tidal parameters.
  static greal mtgcmEWWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE];      //!< EW wind tidal parameters.
  static greal mtgcmNSWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE];      //!< NS wind tidal parameters.
  static greal mtgcmThermosphereBaseHeight[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE];         //!< Thermosphere base height tidal parameters.
  static greal mtgcmF107[MTGCM_F107_SIZE];  //!< F10.7 levels

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MGCMUpperInterpolator, initializeUpperData_getTideParameters);
#endif // GRAM_UNIT_TEST
};

} // namespace

