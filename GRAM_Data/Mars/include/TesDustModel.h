//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#include "MarsDustModelBase.h"

namespace GRAM {

//! \brief The TES dust model class.
//!
//! The TesDustModel Interpolates TES observed dust optical depth to current 
//! latitude, longitude, and Ls for the specified map year.
//!
//! \ingroup MarsGRAM
class TesDustModel : public MarsDustModelBase
{
public:
  TesDustModel();
  TesDustModel(const TesDustModel& orig) = default;
  virtual ~TesDustModel() = default;

  void setInputParameters(const MarsInputParameters& params) override;
  void setMapYear(int value) { mapYear = value; }

private:
  void updateDustOpticalDepth() override;
  void initializeData();

#ifdef CONVERT_MARS_DATA
  void writeData();
  void readLegacyData();
#endif

  // User input.
  int mapYear = 1;   //!< The selected TES map year (1 or 2).

  static const size_t MAP_YEAR_SIZE = 2;  //!< TES map year dimensions.
  static const size_t LS_SIZE = 72;       //!< Ls dimensions.
  static const size_t LAT_SIZE = 25;      //!< Latitude dimensions.
  static const size_t LON_SIZE = 40;      //!< Longitude dimensions.
  static bool isInitialized;              //!< If true, data has been initialized.
  static greal tes_tau[MAP_YEAR_SIZE][LS_SIZE][LAT_SIZE][LON_SIZE];  //!< TES dust optical depth data.

};

} // namespace
