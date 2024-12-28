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
#include "TesGCM.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
TesGCM::TesGCM()
	: MarsGCMBase(tesSurface, tesLower, tesUpper, dustModel)
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
TesGCM::TesGCM(const TesGCM& orig)
  : MarsGCMBase(tesSurface, tesLower, tesUpper, dustModel)
{
  constantHeightOffset = orig.constantHeightOffset;
  offsetModel = orig.offsetModel;
  mapYear = orig.mapYear;
  userHeightAboveSurface = orig.userHeightAboveSurface;
  tesUpper.setF107(orig.tesUpper.getF107());
  tesUpper.setHeightOffset(orig.tesUpper.getHeightOffset());
  tesUpper.setMapYear(mapYear);
  tesLower.setMapYear(mapYear);
  tesSurface.setMapYear(mapYear);
  dustModel.setMapYear(mapYear);
}

//! \fn  TesGCM::~TesGCM()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  TesGCM::getMaxHeight()
//! \brief Returns the maximum height in the upper atmosphere model.

//! \copydoc PerturbedAtmosphere::setInputParameters()
void TesGCM::setInputParameters(const MarsInputParameters& params)
{
  dustModel.setInputParameters(params);
  constantHeightOffset = params.constantHeightOffset;
  offsetModel = params.offsetModel;
  mapYear = params.mapYear;
  tesUpper.setMapYear(mapYear);
  tesUpper.setF107(params.F107);
  tesLower.setMapYear(mapYear);
  tesSurface.setMapYear(mapYear);
  userHeightAboveSurface = params.heightAboveSurface;
}

//! \brief Computes the dust optical depth.
//!
//! This method uses the TesDustModel to compute the dust optical depth.  It also zeroes out the
//! dust offset which is not used in the Tes model.
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//!
//! \retval #dustOpticalDepth
void TesGCM::updateDustModel()
{
  dustModel.setPosition(position);
  dustModel.setLongitudeSun(ephem.longitudeSun);
  dustModel.update();
  dustOpticalDepth = dustModel.getDustOpticalDepth();
  dustOffset = 0.0;
  tesSurface.reset();
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
greal TesGCM::getGlobalMeanOffset()
{
  // Establis the LS index using the upper atmosphere model.
  tesUpper.setPosition(position);
  tesUpper.setEphemerisState(ephem);
  tesUpper.updateIndices(0);
  greal lsDisp = tesUpper.getDisplacements().ls;
  size_t lsIndex = tesUpper.getBaseIndex().ls;

  // Interpolate over the longitude of the sun dimension for the current map year.
  Interpolator interp(lsDisp);
  // 2-d array for interpolating offset
  greal offset[2];
  offset[0] = mtgcmHeightOffsets[lsIndex][mapYear - 1];
  offset[1] = mtgcmHeightOffsets[lsIndex + 1][mapYear - 1];

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
greal TesGCM::getSeasonalOffset() {
  return constantHeightOffset + 1.2 * sin(toRadians(2.0 * longitudeSun));
}

const double TesGCM::mtgcmHeightOffsets[LS_SIZE][MAP_YEAR_SIZE] =
{ // y1      y2      Mapyear
	{-1.22,  -1.25 }, // LS =   0
	{-1.22,   -.38 }, // LS =  30
	{-1.41,  -1.78 }, // LS =  60
	{-1.55,  -1.57 }, // LS =  90
	{-1.54,  -1.46 }, // LS = 120
	{-1.65,  -1.51 }, // LS = 150
	{-1.43,  -1.25 }, // LS = 180
	{-1.09,   2.29 }, // LS = 210
	{  .49,   1.27 }, // LS = 240
	{ -.26,   -.08 }, // LS = 270
	{-1.06,   -.80 }, // LS = 300
	{-1.12,  -1.08 }, // LS = 330
	{-1.22,  -1.25 }  // LS = 360
};

} // namespace