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

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

enum ThermosphereModelType { EG_MSIS, EG_MET, EG_JB2008 };
enum RRAYearType { RRA_1983, RRA_2006, RRA_2013, RRA_2019 };

//! \brief Earth input parameters.
//!
//! This class contains Earth specific input parameters as well
//! as all common input parameters. 
//! \ingroup Cpp_Earth EarthGRAM
class EarthInputParameters : public InputParameters
{
public:
  EarthInputParameters();
  EarthInputParameters(const EarthInputParameters& orig) = default;
  virtual ~EarthInputParameters() = default;
  EarthInputParameters& operator=(const EarthInputParameters& orig) = default;

  // Data location parameters
  std::string atmPath;   //!< Optional override of the default atmosphere data location
  std::string NCEPPath;  //!< Optional override of the default NCEP data location 
  std::string M2Path; //!< Optional override of the default MERRA-2 data location

  // Print parameters
  std::string speciesPath;       //!< Override of legacy species output file name.
  std::string boundaryLayerPath; //!< Override of legacy bltest output file name.

  ThermosphereModelType thermosphereModel = EG_MSIS;  //!< Thermosphere model selection. \note was itherm
  greal surfaceRoughness = -1.0;                      //!< Surface roughness length \units{ m }. \note was z0

  // NCEP Parameters
  bool useNCEP = false;  //!< Set to true to use NCEP data instead of MERRA-2
  int NCEPYear = 9715;   //!< NCEP year code.  \note was NCEPyr
  int NCEPHour = 5;      //!< NCEP hour code.  \note was NCEPhr

  // MERRA-2 Parameters
  int M2Hour = 9;        //!< MERRA-2 hour code.
  double minimumLatitude = -90.0;   //!< Boundary for restricted read of MERRA-2 data
  double maximumLatitude = 90.0;    //!< Boundary for restricted read of MERRA-2 data
  double minimumLongitude = 0.0;    //!< Boundary for restricted read of MERRA-2 data
  double maximumLongitude = 360.0;  //!< Boundary for restricted read of MERRA-2 data

  // RRA Parameters
  std::string rraPath;             //!< Optional override of the default RRA data location
  std::string rraSiteList = "rrasites.txt";  //!< Optional override of the default RRA sites file name
  bool useRRA = false;             //!< RRA activation flag. \note was iurra
  RRAYearType rraYear = RRA_2013;  //!< RRA data selection by year. \note was iyrrra
  greal rraInnerRadius = 0.0;      //!< Lat-lon radius (deg) from RRA site, inside which RRA data is fully weighted \note was sitenear
  greal rraOuterRadius = 0.0;      //!< Lat-lon radius (deg) from RRA site, outside which RRA data are NOT used \note was sitelim

  // JB2008, MET, MSIS model parameters
  greal dailyF10 = 230.0;    //!< Daily 10.7-cm flux \note was f10
  greal meanF10 = 230.0;     //!< Mean 10.7-cm flux  \note was f10b
  greal ap = 16.0;           //!< Geomagnetic index

  // JB2008 model parameters
  greal dailyS10   = 0.0;    //!< EUV index (26-34 nm) scaled to F10 units (0.0 -> s10=f10)       \note was s10
  greal meanS10  = 0.0;      //!< EUV 81-day center-averaged index (0.0 -> s10b = f10b)           \note was s10b
  greal dailyXM10  = 0.0;    //!< MG2 index scaled to F10 units (0.0 -> xm10 = f10)               \note was xm10
  greal meanXM10 = 0.0;      //!< MG2 81-day center-averaged index (0.0 -> xm10b = f10b)          \note was xm10b
  greal dailyY10   = 0.0;    //!< Solar X-Ray & Lya index scaled to F10 (0.0 -> y10=f10)          \note was y10
  greal meanY10  = 0.0;      //!< Solar X-Ray & Lya 81-day avg. centered index (0.0 -> y10b=f10b) \note was y10b
  greal dstdtc = 0.0;        //!< Temperature change computed from Dst index

  // Perturbation parameters
  bool initializePerturbations = false;         //!< Flags user-defined initial perturbations             \note was initpert
  greal initialDensityPerturbation = 0;         //!< Initial density perturbation value (% of mean)       \note was rdinit
  greal initialTemperaturePerturbation = 0;     //!< Initial temperature perturbation value (% of mean).  \note was rtinit
  greal initialEWWindPerturbation = 0;          //!< Initial eastward velocity perturbation (m/s)         \note was ruinit
  greal initialNSWindPerturbation = 0;          //!< Initial northward velocity perturbation (m/s)        \note was rvinit
  greal initialVerticalWindPerturbation = 0;    //!< Initial upward velocity perturbation (m/s)           \note was rwinit
  bool patchy = false;                          //!< Flags patchiness in perturbation model              
  greal randomPerturbationScale = 1.0;          //!< Random perturbation scale for density, temperature and pressure \note was rpscale
  greal horizontalWindPerturbationScale = 1.0;  //!< Random perturbation scale for horizontal winds                  \note was ruscale

  bool corrMonte = false;     //!< Flags use of EarthCorrMonte program.
  greal corrDeltaHours = 0.1; //!< Time offset for hourly dispersions
  bool corrMean = false;      //!< Flags correlation to mean values instead of perturbed.

  greal maximumHeight = 0.0;  //!< Maximum height for ballistic trajectory in EarthCorrTraj
};

} // namespace

