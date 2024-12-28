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
#include "MGCMSurfaceInterpolator.h"
#include "MGCMLowerInterpolator.h"
#include "MGCMUpperInterpolator.h"
#include "MGCMDustModel.h"

namespace GRAM {

//! \brief The MGCM model of the Mars atmosphere.
//!
//! The MGCM model of the Mars atmosphere is comprised of three layers of data:
//!   - Boundary layer data from NASA Ames Mars General Circulation Model (MGCM) surface data ranging from 0 to 0.03 km.
//!   - Lower atmosphere data from NASA Ames Mars General Circulation Model (MGCM) ranging from 0 to 80 km.
//!   - Upper atmosphere data from the University of Michigan Mars Thermospheric General Circulation Model (MTGCM) ranging from 80 to 170 km.
//!
//! \ingroup MarsGRAM
class MGCM : public MarsGCMBase
{
public:
  MGCM();
  MGCM(const MGCM& orig);
  virtual ~MGCM() override = default;

  void setInputParameters(const MarsInputParameters& params) override;
  void updateDustModel() override;
  greal getMaxHeight() const override { return 170.0; }

  // For testing.
  void setDustParameters(greal offset, greal opticalDepth);

private:
  greal getGlobalMeanOffset() override;
  greal getSeasonalOffset() override;

  MGCMSurfaceInterpolator mgcmSurface; //!< The MGCM surface model.
  MGCMLowerInterpolator mgcmLower;     //!< The MGCM lower atmosphere model.
  MGCMUpperInterpolator mgcmUpper;     //!< The MGCM upper atmosphere model.
  MGCMDustModel dustModel;             //!< The MGCM dust model.

  static const size_t LS_SIZE = 13; //!< Size of the LS dimension in the height offsets.
  static const size_t OD_SIZE = 3;  //!< Size of the OD dimension in the height offsets.
  static const double mtgcmHeightOffsets[LS_SIZE][OD_SIZE]; //!< Height offset array.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MGCM, bltp);
  FRIEND_TEST(MGCM, updateHeightOffsets);
  FRIEND_TEST(MGCM, updateFairingHeights);
  FRIEND_TEST(MGCM, updateUpperLower_1);
  FRIEND_TEST(MGCM, updateUpperLower_2);
  FRIEND_TEST(MGCM, updateScaleHeights);
  FRIEND_TEST(MGCM, updateGasConstant);
#endif // GRAM_UNIT_TEST
};

} // namespace

