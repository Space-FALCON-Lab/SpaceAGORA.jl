//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"

namespace GRAM {

//! \brief The augmented output state of an Earth Atmosphere model.
//!
//! The EarthAtmosphereState provides Earth model specific metrics as an
//! extension to the AtmosphereState class using AtmosphereState::planetSpecificMetrics.
//!
//! \ingroup EarthGRAM Cpp_Earth
class EarthAtmosphereState
{
public:
  EarthAtmosphereState();
  EarthAtmosphereState(const EarthAtmosphereState& orig) = default;
  virtual ~EarthAtmosphereState() = default;

  // Dynamics
  greal perturbedTemperature = 0.0;       //!< \units{ K     } \note was tmpert
  greal temperaturePerturbation = 0.0;    //!< \units{ \%    } \note was trh
  greal perturbedPressure = 0.0;          //!< \units{ N/m^2 } \note was pmpert
  greal pressurePerturbation = 0.0;       //!< \units{ \%    } \note was prh

  greal presPertSmall  = 0.0;             //!< Small scale perturbation of pressure                   \units{ \% }. \note was prhs
  greal densPertSmall  = 0.0;             //!< Small scale perturbation of density                    \units{ \% }. \note was drhs
  greal tempPertSmall  = 0.0;             //!< Small scale perturbation of temperature                \units{ \% }. \note was trhs
  greal ewWindPertSmall = 0.0;            //!< Small scale perturbation of east/west wind component   \units{ \% }. \note was urhs
  greal nsWindPertSmall = 0.0;            //!< Small scale perturbation of north/south wind component \units{ \% }. \note was vrhs

  greal presSDSmall = 0.0;                //!< Standard deviation of the small scale perturbations of pressure    \units{ \% }. \note was sphs
  greal densSDSmall = 0.0;                //!< Standard deviation of the small scale perturbations of density     \units{ \% }. \note was sdhs
  greal tempSDSmall = 0.0;                //!< Standard deviation of the small scale perturbations of temperature \units{ \% }. \note was sths
  greal ewWindSDSmall = 0.0;              //!< Standard deviation of the small scale perturbations of east/west wind component   \units{ \% }. \note was suhs
  greal nsWindSDSmall = 0.0;              //!< Standard deviation of the small scale perturbations of north/south wind component \units{ \% }. \note was svhs

  greal presPertLarge = 0.0;              //!< Large scale perturbation of pressure                   \units{ \% }. \note was prhl
  greal densPertLarge = 0.0;              //!< Large scale perturbation of density                    \units{ \% }. \note was drhl
  greal tempPertLarge = 0.0;              //!< Large scale perturbation of temperature                \units{ \% }. \note was trhl
  greal ewWindPertLarge = 0.0;            //!< Large scale perturbation of east/west wind component   \units{ \% }. \note was urhl
  greal nsWindPertLarge = 0.0;            //!< Large scale perturbation of north/south wind component \units{ \% }. \note was vrhl

  greal presSDLarge = 0.0;                //!< Large scale perturbation of pressure                   \units{ \% }. \note was sphl
  greal densSDLarge = 0.0;                //!< Large scale perturbation of density                    \units{ \% }. \note was sdhl
  greal tempSDLarge = 0.0;                //!< Large scale perturbation of temperature                \units{ \% }. \note was sthl
  greal ewWindSDLarge = 0.0;              //!< Large scale perturbation of east/west wind component   \units{ \% }. \note was suhl
  greal nsWindSDLarge = 0.0;              //!< Large scale perturbation of north/south wind component \units{ \% }. \note was svhl

  greal vaporPressure = 0.0;              //!< Mean vapor pressure of water                   \units{ Pa }. \note was eoft
  greal vaporPressureSD = 0.0;            //!< Standard deviation of water vapor pressure     \units{ \% }. \note was seoft
  greal vaporDensity = 0.0;               //!< Mean water vapor density                       \units{ kg/m^3 }. \note was rhov
  greal vaporDensitySD = 0.0;             //!< Standard deviation of the water vapor density  \units{ \% }. \note was shrov
  greal dewPoint = 0.0;                   //!< Mean dewpoint temperature                      \units{ K }.  \note was tdgh
  greal dewPointSD = 0.0;                 //!< Standard deviation of the dewpoint temperature \units{ \% }. \note was stdgh
  greal relativeHumidity = 0.0;           //!< Mean relative humidity                         \units{ \% }. \note was rhp
  greal relativeHumiditySD = 0.0;         //!< Standard deviation of the relative humidity    \units{ \% }. \note was srhp

  greal geodeticLatitude = 0.0;           //!< Latitude in the geodetic reference frame. \units{\text{degrees}}
  char rraSiteName[6] = "";               //!< Name of the active RRA site.
  greal rraWeight = 0.0;                  //!< Fairing weight between RRA and Earth model.

  // Winds
  greal windSpeed = 0.0;                  //!< Mean wind speed \units{ m/s } \note was spdgh
  greal windSpeedStandardDeviation = 0.0; //!< Standard deviation of wind speed \units{ \%  } \note was sdsph
  greal windCorrelation = 0.0;            //!< Horizontal wind correlation \note was uvt2)

  greal windSpeedAtSurface = 0.0;         //!< Mean wind speed at the surface \units{ m/s } \note was spdavsrf1
  greal windSpeedSDAtSurface = 0.0;       //!< Standard deviation of wind speed at the surface \units{ \% } \note was spdsdsrf1
  greal temperatureAtSurface = 0.0;       //!< Mean temperature at the surface \units{ K } \note was tsrf1
  greal temperatureSDAtSurface = 0.0;     //!< Standard deviation of temperature at the surface \units{ \% } \note was stsrf1
  greal pressureSDAtSurface = 0.0;        //!< Standard deviation of the pressure at the surface \units{ \% } \note was sdsrf
  greal densitySDAtSurface = 0.0;         //!< Standard deviation of the density at the surface \units{ \% } \note was spsrf
  greal densityAtSurface = 0.0;           //!< Mean density at the surface \units{ kg/m^3 } \note was dsrf
  greal ewWindAtSurface = 0.0;            //!< Mean velocity of the east/west winds at the surface \units{ m/s } \note was usrf1
  greal nsWindAtSurface = 0.0;            //!< Mean velocity of the north/south winds at the surface \units{ m/s } \note was vsrf1
  greal ewWindSDAtSurface = 0.0;          //!< Standard deviation of the east/west winds at the surface \units{ \% } \note was susrf1
  greal nsWindSDAtSurface = 0.0;          //!< Standard deviation of the north/south winds at the surface \units{ \% } \note was svsrf1
  greal windCorrelationAtSurface = 0.0;   //!< Horizontal wind correlation at the surface \note was uvtsrf1

  int severityLevel = 0;                  //!< Severe turbulence indicator. \note was isev

  // Boundary Layer
  int landCode = 99;                       //!< Land surface type code. \note was lc
  greal surfaceRoughness = 9.9999;         //!< Surface roughness length \units{ m }. \note was z0
  greal solarDays = 0.0;                   //!< Number of solar days since start \units{ days }. \note was rday
  greal solarElevation = 0.0;              //!< The angular height of the sun in the sky \units{\text{degrees}}. \note was el
  greal elevationAtMidnight = 0.0;         //!< The solar elevation at midnight \units{\text{degrees}}. \note was elmn
  greal elevationAtNoon = 0.0;             //!< The solar elevation at noon \units{\text{degrees}}. \note was elmd
  greal netRadiationIndex = 0.0;           //!< Net radiation index [unitless], for stability calculation.     \note was nri
  greal stability = 0.0;                   //!< Atmospheric stability category [unitless].                     \note was s
  greal inverseLength = 0.0;               //!< Inverse of Monin-Obukhov scale length         \units{ 1/m   }. \note was ool
  greal frictionVelocity = 0.0;            //!< Surface friction velocity                     \units{ m/s   }. \note was ustar
  greal BVFrequencySquare = 0.0;           //!< Square of Brunt-Vaisala frequency             \units{ 1/s^2 }. \note was bvfsq
  greal boundaryLayerDepth = 0.0;          //!< Current depth of boundary layer               \units{ m     }. \note was hbl
  greal neutralBoundaryLayerDepth = 0.0;   //!< Depth of neutral boundary layer               \units{ m     }. \note was hn
  greal unstableBLFactor = 1.0;            //!< Height factor for unstable BL during early daytime [unitless]. \note was blfact
  greal sigmaRatio = 0.0;                  //!< Ratio of sigmaW to friction velocity at current height.        \note was sigrat
  greal sigmaW = 0.0;                      //!< Vertical wind standard deviation at top of boundary layer          \units{ m/s }. \note was swb
  greal metersAboveSurface = 0.0;          //!< Current height above surface or height to top of boundary layer    \units{ m   }. \note was chb
  greal perturbedWindSpeedAtSurface = 0.0; //!< Surface (10m) wind speed (from mean-plus-large-scale-perturbation) \units{ m/s }. \note was spdsrf
};

} // namespace
