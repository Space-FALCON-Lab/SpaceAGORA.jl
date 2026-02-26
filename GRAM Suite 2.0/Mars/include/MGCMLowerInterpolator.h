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
#include "MGCMInterpolator.h"

namespace GRAM {

//! \brief The MGCM lower atmosphere model (0 to 80 km).
//!
//! The MGCMLowerInterpolator uses interpolation routines to evaluate the atmospheric
//! state at the current position.  Results are based on data from the 
//! NASA Ames Mars General Circulation Model (MGCM) ranging from 0 to 80 km.
//! This class requires the data file "MGCM_lower_data.bin".
//!
//! \ingroup MarsGRAM
class MGCMLowerInterpolator : public MGCMInterpolator
{
public:
  MGCMLowerInterpolator();
  MGCMLowerInterpolator(const MGCMLowerInterpolator& orig) = default;
  virtual ~MGCMLowerInterpolator() = default;

  virtual void updateIndices(size_t heightIndexOffset) override;
  virtual void update() override;

protected:
  static const size_t MGCM_HEIGHT_SIZE = 17; //!< The number of height layers in the MGCM data (0 to 80 by 5 km).
  static const size_t MGCM_LAT_SIZE = 25;    //!< The number of latitude layers in the MGCM data (-90 to 90 by 7.5 degrees).

  MarsTideParameters getTideParameters(const greal data[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE], size_t latIndex, greal poleFactor);
  void initializeData();

#ifdef CONVERT_MARS_DATA
  void writeLowerData();
  void readLegacyLowerData();
#endif

  static bool isInitialized;  //!< \brief True after data initialization.
  static greal mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE]; //!< Temperature tidal parameters.
  static greal mgcmPressure[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE];    //!< Pressure tidal parameters.
  static greal mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE];      //!< EW wind tidal parameters.
  static greal mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE];      //!< NS wind tidal parameters.
  static greal mgcmDensity[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE];                 //!< Mean density.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MGCMLowerInterpolator, initializeLowerData_getTideParameters);
#endif // GRAM_UNIT_TEST
};

} // namespace

