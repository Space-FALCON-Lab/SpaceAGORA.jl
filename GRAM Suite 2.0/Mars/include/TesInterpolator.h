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
#include "MarsInputParameters.h"

namespace GRAM {

//! \brief The base class for TesGCM interpolation models.
//!
//! The TesGCM models consist of tidal data and multi-dimensional interpolations.
//! This is the base class for the three TesGCM models.  The MarsInterpolatorBase class
//! defines the interface to be implemented by all of the models.  This class provides
//! a common method for computing interpolation indices and displacements.
//!
//! \ingroup MarsGRAM
class TesInterpolator : public MarsInterpolatorBase
{
public:
  TesInterpolator();
  TesInterpolator(const TesInterpolator& orig) = default;
  virtual ~TesInterpolator() = default;

  void setMapYear(int year) { mapYear = year; }

protected:
  void updateBaseIndices(size_t latSize, greal offset);

  int mapYear = 1;

  // These sizes are common to all TesGCM models
  static const size_t PARAM_SIZE = 5;     //!< The number of tidal parameters.
  static const size_t LS_SIZE = 13;       //!< Longitude of the sun layers (from 0 to 360 by 30 degrees).
  static const size_t TES_YEAR_SIZE = 2;  //!< The number of TES years (1 and 2).


  static greal dust[TES_YEAR_SIZE];       //!< Optical depths for TES years 1 and 2.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TesInterpolator, updateIndices);
  FRIEND_TEST(TesInterpolator, updateIndices_noOffset);
#endif // GRAM_UNIT_TEST

};

} // namespace

