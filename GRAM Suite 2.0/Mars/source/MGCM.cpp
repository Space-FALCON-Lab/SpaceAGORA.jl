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
#include "MGCM.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
MGCM::MGCM()
: MarsGCMBase(mgcmSurface, mgcmLower, mgcmUpper, dustModel)
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
MGCM::MGCM(const MGCM& orig)
  : MarsGCMBase(mgcmSurface, mgcmLower, mgcmUpper, dustModel)
{
  constantHeightOffset = orig.constantHeightOffset;
  offsetModel = orig.offsetModel;
  userHeightAboveSurface = orig.userHeightAboveSurface;
  mgcmUpper.setF107(orig.mgcmUpper.getF107());
  dustModel.setLevels(orig.dustModel.getConstantDustLevel(),
    orig.dustModel.getMaxDustLevel(), orig.dustModel.getMinDustLevel());
}

//! \fn  MGCM::~MGCM()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MGCM::getMaxHeight()
//! \brief Returns the maximum height in the upper atmosphere model.

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MGCM::setInputParameters(const MarsInputParameters& params)
{
  mgcmUpper.setF107(params.F107);
  userHeightAboveSurface = params.heightAboveSurface;
  dustModel.setInputParameters(params);
  constantHeightOffset = params.constantHeightOffset;
  offsetModel = params.offsetModel;
}

//! \brief Used by testing only to set dust parameters.
//!
//! This method facilitates testing by setting dust parameters without using the dust model.
//! \param offset        The dust offset.
//! \param opticalDepth  The dust optical depth.
void MGCM::setDustParameters(greal offset, greal opticalDepth)
{
  dustOffset = offset;
  mgcmSurface.setDustOpticalDepth(opticalDepth);
  mgcmLower.setDustOpticalDepth(opticalDepth);
  mgcmUpper.setDustOpticalDepth(opticalDepth);
}

//! \brief Computes the dust optical depth and dust offset.
//!
//! This method uses the MGCMDustModel to compute the dust optical depth and dust offset.
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//!
//! \retval #dustOffset
//! \retval #dustOpticalDepth
void MGCM::updateDustModel()
{
  dustModel.setPosition(position);
  dustModel.setLongitudeSun(ephem.longitudeSun);
  dustModel.update();
  dustOffset = dustModel.getDustOffset();
  dustOpticalDepth = dustModel.getDustOpticalDepth();
  mgcmSurface.reset();
  mgcmSurface.setDustOpticalDepth(dustOpticalDepth);
  mgcmLower.setDustOpticalDepth(dustOpticalDepth);
  mgcmUpper.setDustOpticalDepth(dustOpticalDepth);
}

//! \brief Computes and returns the global mean height offset.
//!
//! This method computes and returns the global mean height offset.
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//!
//! \returns The global mean height offset.
greal MGCM::getGlobalMeanOffset()
{
  // Establis the LS and OD indices using the upper atmosphere model.
  mgcmUpper.setPosition(position);
  mgcmUpper.setEphemerisState(ephem);
  mgcmUpper.updateIndices(0);
  greal odDisp = mgcmUpper.getDisplacements().od;
  greal lsDisp = mgcmUpper.getDisplacements().ls;
  size_t odIndex = mgcmUpper.getBaseIndex().od;
  size_t lsIndex = mgcmUpper.getBaseIndex().ls;

  // Interpolate over optical depth and longitude of the sun dimensions
  Interpolator interp(lsDisp, odDisp);
  // 2-d array for interpolating offset
  greal offset[2][2];
  offset[0][0] = mtgcmHeightOffsets[lsIndex][odIndex];
  offset[0][1] = mtgcmHeightOffsets[lsIndex][odIndex + 1];
  offset[1][0] = mtgcmHeightOffsets[lsIndex + 1][odIndex];
  offset[1][1] = mtgcmHeightOffsets[lsIndex + 1][odIndex + 1];

  // linear interpolation for global offset
  return interp.linear(offset);
}

//! \brief Computes and returns a seasonal height offset.
//!
//! This method computes and returns the seasonal height offset needed for the MARS_SEASONAL offset model.
//!
//! \b Inputs
//! \arg #constantHeightOffset
//! \arg #longitudeSun
//!
//! \returns The seasonal height offset.
greal MGCM::getSeasonalOffset() {
  return constantHeightOffset + cos(toRadians(2.0 * longitudeSun));
}

const double MGCM::mtgcmHeightOffsets[LS_SIZE][OD_SIZE] =
{ // .3      1.0     3.0 Optical Depth
  { -0.002, -0.001, -0.003 }, // LS =   0
  { -0.001,  0.000,  0.003 }, // LS =  30
  {  0.001, -0.003, -0.002 }, // LS =  60
  { -0.002, -0.002, -0.000 }, // LS =  90
  {  0.001, -0.002, -0.001 }, // LS = 120
  {  0.004, -0.002,  0.011 }, // LS = 150
  { -0.003, -0.004, -0.002 }, // LS = 180
  { -0.002,  0.001, -0.001 }, // LS = 210
  { -0.002, -0.001,  0.018 }, // LS = 240
  { -0.002,  0.004, -0.003 }, // LS = 270
  { -0.001,  0.008, -0.002 }, // LS = 300
  {  0.033, -0.001, -0.003 }, // LS = 330
  { -0.002, -0.001, -0.003 }  // LS = 360
};
//{ // .3      1.0     3.0 Optical Depth
//  { 3.26, 2.88, 4.49 }, // LS =   0
//  { 2.96,   3.19,   4.26 }, // LS =  30
//  { 2.62,   3.05,   4.07 }, // LS =  60
//  { 2.35,   2.86,   3.36 }, // LS =  90
//  { 2.21,   3.06,   3.59 }, // LS = 120
//  { 2.77,   3.16,   4.07 }, // LS = 150
//  { 3.28,   3.55,   4.86 }, // LS = 180
//  { 3.35,   3.58,   4.94 }, // LS = 210
//  { 3.24,   3.52,   2.81 }, // LS = 240
//  { 3.08,   1.32,   2.81 }, // LS = 270
//  { 3.05,   1.86,   2.43 }, // LS = 300
//  { 3.21,   2.93,   4.18 }, // LS = 330
//  { 3.26,   2.88,   4.49 }  // LS = 360
//};

} // namespace