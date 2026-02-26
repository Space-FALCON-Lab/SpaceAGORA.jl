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
#include "MarsGCMBase.h"
#include "TesSurfaceInterpolator.h"
#include "TesLowerInterpolator.h"
#include "TesUpperInterpolator.h"
#include "TesDustModel.h"

namespace GRAM {

//! \brief The TesGCM model of the Mars atmosphere.
//!
//! The TesGCM model of the Mars atmosphere is comprised of three layers of data:
//!   - Boundary layer data from NASA Ames Mars General Circulation Model (MGCM) surface data ranging from 0 to 0.03 km.
//!   - Lower atmosphere data from NASA Ames Mars General Circulation Model (MGCM) ranging from -5 to 80 km.
//!   - Upper atmosphere data from the University of Michigan Mars Thermospheric General Circulation Model (MTGCM) ranging from 80 to 240 km.
//!
//! \ingroup MarsGRAM
class TesGCM : public MarsGCMBase
{
public:
  TesGCM();
  TesGCM(const TesGCM& orig);
  virtual ~TesGCM() override = default;

  void setInputParameters(const MarsInputParameters& params) override;
  void updateDustModel() override;
  greal getMaxHeight() const override { return 240; }

private:
  greal getGlobalMeanOffset() override;
  greal getSeasonalOffset() override;

  TesSurfaceInterpolator tesSurface; //!< The Tes surface model.
  TesLowerInterpolator tesLower;     //!< The Tes lower atmosphere model.
  TesUpperInterpolator tesUpper;     //!< The Tes upper atmosphere model.
  TesDustModel dustModel;            //!< The Tes dust model.

  int mapYear = 1; //!< The selected map year (1 or 2).

  static const size_t LS_SIZE = 13;       //!< Size of the LS dimension in the height offsets.
  static const size_t MAP_YEAR_SIZE = 2;  //!< Size of the MapYear dimension in the height offsets.
  static const double mtgcmHeightOffsets[LS_SIZE][MAP_YEAR_SIZE]; //!< Height offset array.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TesGCM, bltp);
  FRIEND_TEST(TesGCM, updateHeightOffsets);
  FRIEND_TEST(TesGCM, updateFairingHeights);
  FRIEND_TEST(TesGCM, DISABLED_updateScaleHeights);
#endif // GRAM_UNIT_TEST
};

} // namespace

