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
#include "MarsInterpolatorBase.h"

namespace GRAM {

//! \brief The base class for MGCM interpolation models.
//!
//! The MGCM models consist of tidal data and multi-dimensional interpolations.
//! This is the base class for the three MGCM models.  The MarsInterpolatorBase class
//! defines the interface to be implemented by all of the models.  This class provides
//! a common method for computing interpolation indices and displacements.
//!
//! \ingroup MarsGRAM
class MGCMInterpolator : public MarsInterpolatorBase
{
public:
  MGCMInterpolator();
  MGCMInterpolator(const MGCMInterpolator& orig) = default;
  virtual ~MGCMInterpolator() = default;

  void setDustOpticalDepth(greal depth) { dustOpticalDepth = depth; }

protected:
  void updateBaseIndices(size_t latitudeArraySize, greal latitudeOffset);

  // These sizes are common to all MGCM models
  static const size_t PARAM_SIZE = 5;  //!< The number of tidal parameters.
  static const size_t LS_SIZE = 13;    //!< Longitude of the sun layers. From 0 to 360 by 30 degrees.
  static const size_t OD_SIZE = 3;     //!< The number of optical depth layers.
  static const greal dustOD[OD_SIZE];  //!< The optical depth layers.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MGCMInterpolator, updateIndices);
  FRIEND_TEST(MGCMInterpolator, updateIndices_noOffset);
#endif // GRAM_UNIT_TEST

};

} // namespace

