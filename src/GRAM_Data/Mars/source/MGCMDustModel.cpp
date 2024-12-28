//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include "MGCMDustModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MGCMDustModel::MGCMDustModel()
  : MarsDustModelBase()
{
}

//! \fn  MGCMDustModel::MGCMDustModel(const MGCMDustModel& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MGCMDustModel::~MGCMDustModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MGCMDustModel::setInputParameters(const MarsInputParameters& params)
{
  MarsDustModelBase::setInputParameters(params);

  mgcmConstantDustLevel = params.mgcmConstantDustLevel;
  mgcmMaxDustLevel = params.mgcmMaxDustLevel;
  mgcmMinDustLevel = params.mgcmMinDustLevel;
}

//! \brief Computes the dust optical depth (tau).
//!
//! Computes the dust optical depth  using viking-like seasonal variation.
//! Computes seasonal (function of Ls) dust optical depth.  Assumed (sinusoidal) seasonal variation 
//! for non-dust-storm optical depth  (mgcmMinDustLevel at Ls = 90; mgcmMaxDustLevel at Ls = 270)                
//!
//! \b Inputs
//! \arg #longitudeSun
//! \arg #mgcmConstantDustLevel
//! \arg #mgcmMaxDustLevel
//! \arg #mgcmMinDustLevel
//!
//! \retval #dustOpticalDepth
void MGCMDustModel::updateDustOpticalDepth()
{
  if (mgcmConstantDustLevel > 0.0) {
    // if mgcmConstantDustLevel > 0, use as OD, with minimum value of 0.1
    dustOpticalDepth = max(0.1, mgcmConstantDustLevel); 
  }
  else {
    // sinusoidal variation on Ls
    double s = sin(toRadians(longitudeSun));
    // mgcmMinDustLevel at Ls = 90; mgcmMaxDustLevel at Ls = 270
    dustOpticalDepth = (mgcmMaxDustLevel * (1.0 - s) + mgcmMinDustLevel * (1.0 + s)) / 2.0;
  }
}


} // namespace
