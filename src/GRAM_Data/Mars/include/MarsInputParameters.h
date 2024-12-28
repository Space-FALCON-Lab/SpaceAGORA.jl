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
#include "InputParameters.h"
#include "MarsCommon.h"

namespace GRAM {

//! \brief Mars input parameters.
//!
//! This class contains Mars specific input parameters as well
//! as all common input parameters. 
//! \ingroup Cpp_Mars MarsGRAM
class MarsInputParameters : public InputParameters
{
public:
  MarsInputParameters();
  MarsInputParameters(const MarsInputParameters& orig) = default;
  virtual ~MarsInputParameters() = default;

  // MGCMDustModel
  greal mgcmConstantDustLevel = 0.3;  //!< Optical depth of background dust level.   \note was dusttau
  greal mgcmMaxDustLevel = 1.0;       //!< Maximum seasonal dust level.              \note was dustmax
  greal mgcmMinDustLevel = 0.3;       //!< Minimum seasonal dust level.              \note was dustmin

  // MarsAtmosphere
  greal dustNu = 0.003;                    //!< Parameter for vertical distribution of dust density.
  greal dustDiameter = 5.0;                //!< Dust particle diameter.                         \note was dustdiam
  greal dustDensity = 3000.0;              //!< Dust particle density.                          \note was dustdens
  bool isMolaHeights = true;               //!< True for input heights relative to MOLA areoid. \note was MOLAhgts
  greal perturbationWaveLengthScale = 1.0; //!< Scale factor for perturbation wavelengths.      \note was wlscale

  greal waveDate = 0.0;        //!< Julian date for (primary) peak(s) of wave. 
  greal waveMeanOffset = 1.0;  //!< Mean term of longitude-dependent wave multiplier for density. \note was WaveA0
  greal waveAmplitude1 = 0.0;  //!< Amplitude of wave-1 component of longitude-dependent wave multiplier for density. \note was WaveA1
  greal waveAmplitude2 = 0.0;  //!< Amplitude of wave-2 component of longitude-dependent wave multiplier for density. \note was WaveA2
  greal waveAmplitude3 = 0.0;  //!< Amplitude of wave-3 component of longitude-dependent wave multiplier for density. \note was WaveA3
  greal wavePhase1 = 0.0;      //!< Phase of wave-1 component of longitude-dependent wave multiplier \note was Wavephi1
  greal wavePhase2 = 0.0;      //!< Phase of wave-2 component of longitude-dependent wave multiplier \note was Wavephi2
  greal wavePhase3 = 0.0;      //!< Phase of wave-3 component of longitude-dependent wave multiplier \note was Wavephi3
  greal wavePhase1Rate = 0.0;  //!< Rate of longitude movement for wave-1 component. \note was phi1dot
  greal wavePhase2Rate = 0.0;  //!< Rate of longitude movement for wave-2 component. \note was phi2dot
  greal wavePhase3Rate = 0.0;  //!< Rate of longitude movement for wave-3 component. \note was phi3dot
  greal waveScale = 20.0;      //!< Vertical scale of longitude-dependent wave damping at altitudes below 100 km \note was wscale
  std::string waveFile = "";   //!< File path for time-dependent wave coefficient data.

  // MarsDustModelBase
  greal stormLongitudeSun = 0.0; //!< Starting LS value for dust storm.   /note was ALS0
  greal stormDuration  = 48.0;   //!< Duration for dust storm.            /note was ALSDUR
  greal stormIntensity  = 0.0;   //!< Dust storm intensity.               /note was INTENS
  greal stormMaxRadius  = 0.0;   //!< Maximum radius of the dust storm.   /note was RADMAX
  greal stormLatitude = 0.0;     //!< Latitude for center of dust storm.  /note was DUSTLAT
  greal stormLongitude = 0.0;    //!< Longitude for center of dust storm. /note was DUSTLON

  // MarsInterpolatorBase
  greal F107 = 68.0;             //!< F10.7 cm solar flux at 1 AU.
  bool computeMinMax = true;     //!< True for daily max/min data computations.  \note was idaydata

  // StewartModel
  greal exosphericTemperatureFactor = 0.0; //!< Standard deviation for thermosphere variation. \note was STDL
  greal exosphericTemperatureOffset = 0;   //!< Adjustment for exospheric temperature.         \note was deltatx

  // MarsGCMBase
  greal constantHeightOffset = 3.25;              //!< Constant part of height offset on select offset models. \note was zoffset
  MarsOffsetModel offsetModel = MARS_GLOBAL_MEAN; //!< Model selection for height offset term.                 \note was ibougher

  // SlopeWindsModel
  greal meanWindsScale = 1.0;          //!< Scale factor for mean winds.                 \note was wmscale
  greal boundaryLayerWindsScale = 1.0; //!< Scale factor for boundary layer slope winds  \note was blwinfac

  // TesInterpolator, TesDustModel
  int mapYear = 0;                     //!< TES mapping year 1 or 2.  0 for GCM model.

  greal heightAboveSurface = 0.0;      //!< Height above surface.  \note was hgtasfcm

  greal equatorialRadius = 0.0;        //!< Equatorial radius for reference ellipsoid  \note was requa
  greal polarRadius = 0.0;             //!< Polar radius for reference ellipsoid       \note was rpole

private:

};

} // namespace

