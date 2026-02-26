//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#include "gram.h"
#include "MarsDustModelBase.h"

namespace GRAM {

//! \brief The MGCM dust model class.
//!
//! The MGCMDustModel uses viking-like seasonal variation of optical depth.
//! Assumed (sinusoidal) seasonal variation for non-dust-storm optical depth.
//!
//! \ingroup MarsGRAM
class MGCMDustModel : public MarsDustModelBase
{
public:
  MGCMDustModel();
  MGCMDustModel(const MGCMDustModel& orig) = default;
  virtual ~MGCMDustModel() = default;

  void setInputParameters(const MarsInputParameters& params) override;
  void setLevels(greal constant, greal max, greal min)
  {
    mgcmConstantDustLevel = constant;
    mgcmMaxDustLevel = max;
    mgcmMinDustLevel = min;
  }
  greal getConstantDustLevel() const { return mgcmConstantDustLevel; }
  greal getMaxDustLevel() const { return mgcmMaxDustLevel; }
  greal getMinDustLevel() const { return mgcmMinDustLevel; }

private:
  void updateDustOpticalDepth() override;

  // User input parameters.
  greal mgcmConstantDustLevel = 0;  //!< If non-zero, this value is the constant optical depth.
  greal mgcmMaxDustLevel = 0;       //!< The largest optical depth value (at Ls = 270).
  greal mgcmMinDustLevel = 0;       //!< The smallest optical depth value (at Ls = 90).
};

} // namespace
