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
#include "TesInterpolator.h"

namespace GRAM {

//! \brief The TES lower atmosphere model (-5 to 80 km).
//!
//! The TesLowerInterpolator uses interpolation routines to evaluate the atmospheric
//! state at the current position.  Results are based on data from the 
//! NASA Ames Mars General Circulation Model (MGCM) ranging from -5 to 80 km.
//! This class requires the data file "TES_lower_data.bin".
//!
//! \ingroup MarsGRAM
class TesLowerInterpolator : public TesInterpolator
{
public:
  TesLowerInterpolator();
  TesLowerInterpolator(const TesLowerInterpolator& orig) = default;
  virtual ~TesLowerInterpolator() = default;

  void updateIndices(size_t heightIndexOffset) override;
  void update() override;

protected:
  static const size_t MGCM_HEIGHT_SIZE = 30; //!< The number of height layers in the MGCM data (-5 to 10 by 1km, 10 to 80 by 5 km).
  static const size_t MGCM_LAT_SIZE = 25;    //!< The number of latitude layers in the MGCM data (-90 to 90 by 7.5 degrees).

  MarsTideParameters getTideParameters(const greal data[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE], size_t latIndex, greal poleFactor);
  void initializeData();

#ifdef CONVERT_MARS_DATA
  void readLegacyLowerData();
  void writeLowerData();
#endif

  static bool isInitialized;  //!< \brief True after data initialization.
  static greal mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE]; //!< Temperature tidal parameters.
  static greal mgcmDensity[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE];     //!< Density tidal parameters.
  static greal mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE];      //!< EW wind tidal parameters.
  static greal mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE];      //!< NS wind tidal parameters.
  static greal mgcmPressure[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE];                //!< Mean pressure.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TesLowerInterpolator, initializeLowerData_getTideParameters);
#endif // GRAM_UNIT_TEST

};

} // namespace

