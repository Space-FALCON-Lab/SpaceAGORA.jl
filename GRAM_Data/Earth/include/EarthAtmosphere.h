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
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "unittest_friend.h"
#include "EarthCommon.h"
#include "PerturbedAtmosphere.h"
#include "Ephemeris.h"
#include "EarthInputParameters.h"
#include "EarthModel.h"
#include "EarthAtmosphereState.h"
#include "RRA.h"
#include "Interpolator.h"
#include "RandomNumberGenerator.h"

//! \defgroup Cpp_Earth The C++ Interface for Earth 
//! \brief Class references for the EarthGRAM interface.
//!
//! These are the classes needed for building an application that interfaces with EarthGRAM.
//! The interface for Earth is declared in the following:
//! \code #include "EarthAtmosphere.h" \endcode 
//! An example using this class can be found \ref Earth/examples/Earth.cpp "here".

//! \defgroup EarthGRAM The Full List of EarthGRAM Classes.
//! \brief Class references for all Earth models.
//!
//! These are the Earth specific classes needed for building EarthGRAM.

namespace GRAM {

//! \brief The primary interface for the Earth atmosphere model.
//!
//! The Earth atmosphere atmosphere model is a mixture of empirically based models that represent
//! different altitude ranges and the geographical and temporal variations within these altitude ranges.
//! Within this class is a perturbation model that computes variations about these means if dispersions 
//! are desired.
//! \ingroup Cpp_Earth EarthGRAM
class EarthAtmosphere : public PerturbedAtmosphere, public EarthCommon
{
//! \example Earth/examples/Earth.cpp
public:
  EarthAtmosphere();
  EarthAtmosphere(const EarthAtmosphere& orig);
  virtual ~EarthAtmosphere() override;

  void setInputParameters(const EarthInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void setThermosphereModel(ThermosphereModelType model);
  void setSurfaceRoughness(greal z0) { inputParameters.surfaceRoughness = z0; }
  void setDataPath(const std::string& path);
  void setAtmosPath(const std::string& path);
  void setUseNCEP(bool flag);
  void setNCEPPath(const std::string& path);
  void setMERRA2Path(const std::string& path);
  void setNCEPParameters(int NCEPYear, int NCEPHour);
  void setMERRA2Parameters(int M2Hour, greal latMin, greal latMax, greal lonMin, greal lonMax);
  void setRRAPath(const std::string& path);
  void setRRASiteList(const std::string& fileName);
  void setRRAParameters(RRAYearType year, greal innerRadius, greal outerRadius);
  void setUseRRA(bool useFlag);
  void setSolarParameters(greal dailyF10, greal meanF10, greal ap);
  void setJB2008Parameters(greal dailyS10, greal meanS10, greal dailyXM10, greal meanXM10, greal dailyY10, greal meanY10, greal dstdtc);
  void setInitialPerturbations(greal densPert, greal tempPert, greal ewPert, greal nsPert, greal verticalPert);
  void setPatchiness(bool usePatchiness) { inputParameters.patchy = usePatchiness; }
  void setRandomPerturbationScale(greal scaleFactor);
  void setHorizontalWindPerturbationScale(greal scaleFactor);
  void setVerticalWindPerturbationScale(greal scaleFactor) override;
  void setStartTime(const GramTime& gramTime) override;
  void setSeed(int seed) override { initializing = true; randomNumberGenerator.setSeed(seed); }

  void update() override;

  static const std::string& getVersionString();
  const EarthAtmosphereState& getEarthAtmosphereState() { return earthAtmos; }

  static void getCorrelationCoefficients(const Position& base, const Position& target, greal& baseCoefficient, greal& targetCoefficient);

private:
  // Removing these from public exposure.
  using PerturbedAtmosphere::setDensityPerturbationScale;
  using PerturbedAtmosphere::setEWWindPerturbationScale;
  using PerturbedAtmosphere::setNSWindPerturbationScale;
  using PerturbedAtmosphere::setPerturbationScales;
  using PerturbedAtmosphere::setMinRelativeStepSize;

  void updateReferenceValues() override;
  void updatePressureAtSurface() override;
  void updatePerturbations(greal meanDensity = -9999.0, greal meanEWWind = -9999.0, greal meanNSWind = -9999.0) override;
  void updateSpeedOfSound() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) {}
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) {}
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) {}

  void initializeData();
  void initializeTopographyData();
  void initializePerturbations();

  void updateStandardDeviations();
  void updatePDCorrelation();
  void updateLargeScalePerturbations(greal scaleFactor);
  void updateWindSpeed();
  void updateGases();
  void updateRRA(bool init = false);
  void updateSurfaceHeight();
  void updateBoundaryLayerModel(bool init = false);

  void getStateVariables(const AtmosphereState& em);
  void getUS76StandardAtmosphere(greal &temp, greal &pres, greal &dens);
  greal getProbabilityTail(greal x);
  void getPerturbationCoefficients(greal densSD_1, greal densSD_2, greal presSD_1, greal presSD_2,
    greal tempSD_1, greal tempSD_2, greal ewSD_1, greal ewSD_2, greal nsSD_1, greal nsSD_2,
    greal udCorr_1, greal udCorr_2, greal vdCorr_1, greal vdCorr_2, greal pdCorr,
    greal presDistCorr, greal densDistCorr, greal windDistCorr,
    greal &densPrevCoef, greal &densRandCoef, greal &presPrevCoef, greal &presDensCoef, greal &presRandCoef,
    greal &ewPrevCoef, greal &ewDensCoef, greal &ewRandCoef, greal &nsPrevCoef, greal &nsDensCoef, greal &nsRandCoef);

  static const size_t rsHeightSize = 29;  //!< Array size for variable-scale perturbation data (RS data).
  static const size_t pcHeightSize = 25;  //!< Array size for large-scale fraction data (PT,PW data) and density-velocity correlation data (CS,CL data).
  static const size_t pcLatSize = 10;     //!< Array size for large-scale fraction data (PT,PW data) and density-velocity correlation data (CS,CL data).
  static const size_t topoLonSize = 360;  //!< Array size for topographical data.
  static const size_t topoLatSize = 180;  //!< Array size for topographical data.

  static greal getCorrelation(greal x);
  static greal interpolateRSData(greal data[rsHeightSize], greal hgt);  // was intrw
  void interpolatePCData(greal xarray[pcHeightSize][pcLatSize], greal yarray[pcHeightSize][pcLatSize],
    greal h, greal phi, greal &suh, greal &svh);      // was intr25


  EarthModel* earthModelPtr = nullptr;     //!< The Earth height based data model.
  EarthInputParameters inputParameters;    //!< User supplied parameters.
  RRA* rraPtr = nullptr;                   //!< The Range Reference Atmosphere model (pointer).
  EarthAtmosphereState earthAtmos;         //!< Earth specific metrics added to the AtmosphereState.
  RandomNumberGenerator corrRandom;        //!< Provides random numbers for CorrMonte perturbations.

  static bool initialized;                 //!< Initialization flag for static data.

  // Topographical data
  static int landcd[topoLonSize][topoLatSize];   //!< Land code array for topographical data.
  static greal ztopo[topoLonSize][topoLatSize];  //!< Surface elevation array (meters above sea level).

  // RS: Variable-scale perturbation data.
  static greal xlbar[rsHeightSize];   //!< Mean value for L_h.
  static greal zlbar[rsHeightSize];   //!< Mean value for L_z.
  static greal xsigl[rsHeightSize];   //!< Standard deviation for L_h.
  static greal zsigl[rsHeightSize];   //!< Standard deviation for L_z.
  static greal xlmin[rsHeightSize];   //!< The minimum (turbulence) scale value for L_h.
  static greal zlmin[rsHeightSize];   //!< The minimum (turbulence) scale value for L_z.
  static greal xscale[rsHeightSize];  //!< The spatial scale value for L_h.
  static greal zscale[rsHeightSize];  //!< The spatial scale value for L_z.
  static greal wr[rsHeightSize];      //!< Standard deviation of small-scale vertical wind component.

  // PT: Large-scale fraction data (pdt)
  static greal plp[pcHeightSize][pcLatSize];  //!< Large-scale fraction pressure data.
  static greal dlp[pcHeightSize][pcLatSize];  //!< Large-scale fraction density data.
  static greal tlp[pcHeightSize][pcLatSize];  //!< Large-scale fraction temperature data.

  // PW: Large-scale fraction data (winds)
  static greal ulp[pcHeightSize][pcLatSize];  //!< Large-scale fraction east/west wind data.
  static greal vlp[pcHeightSize][pcLatSize];  //!< Large-scale fraction north/south wind data.

  // CS: Small-scale density-velocity correlation data
  static greal uds[pcHeightSize][pcLatSize];  //!< Small-scale east/west wind and density correlation data.
  static greal vds[pcHeightSize][pcLatSize];  //!< Small-scale north/south wind and density correlation data.

  // CL: Large-scale density-velocity correlation data
  static greal udl[pcHeightSize][pcLatSize];  //!< Large-scale east/west wind and density correlation data.
  static greal uvt[pcHeightSize][pcLatSize];  //!< Large-scale north/south wind and density correlation data.

  static const int ptData[pcHeightSize][15];  //!< Integer form of PT data: Large-scale fraction data (pdt).
  static const int pwData[pcHeightSize][10];  //!< Integer form of PW data: Large-scale fraction data (winds).
  static const int csData[pcHeightSize][10];  //!< Integer form of CS data: Small-scale density-velocity correlation data.
  static const int clData[pcHeightSize][10];  //!< Integer form of CL data: Large-scale density-velocity correlation data.

  typedef struct RsDataType {
    greal xlbar;      //!< Mean value for L_h.
    greal xsigl;      //!< Mean value for L_z.
    greal xlmin;      //!< Standard deviation for L_h.
    greal xscale;     //!< Standard deviation for L_z.
    greal zlbar;      //!< The minimum (turbulence) scale value for L_h.
    greal zsigl;      //!< The minimum (turbulence) scale value for L_z.
    greal zlmin;      //!< The spatial scale value for L_h.
    greal zscale;     //!< The spatial scale value for L_z.
    greal wr;         //!< Standard deviation of small-scale vertical wind component.
  } RsData;

  static RsData rsData[rsHeightSize];  //!< RS data: Variable-scale perturbation data.

  const greal maximumPressureVariance = 0.99;   //!< Upper limit for pressure variance.
  const greal maximumVariance = 0.85;           //!< Upper limit for density and winds variance.    \note was plmax
  const greal BuellLimit = 0.99;                //!< The Buell constraint for standard deviations.  \note was rlim
  const greal nearlyZero = 1.0e-5;              //!< A near zero value to ensure non-zero division. \note was slim
  bool initializing = true;                     //!< Set to \c true on first pass of each profile.
  greal baseDayOfYear = 0.0;
  int leapYear = 0;
  bool heightErrorThrown = false;

  greal udCorrSmall = 0.0;   //!< Small-scale correlation between east/west wind and density.        \note was uds2
  greal vdCorrSmall = 0.0;   //!< Small-scale correlation between north/south wind and density.      \note was vds2
  greal pdCorrLarge = 0.0;   //!< Large-scale correlation between pressure and density.              \note was rpdl
  greal udCorrLarge = 0.0;   //!< Large-scale correlation between east/west wind and density.        \note was udl2

  greal waveRand = 0.0;      //!< Random number used in large-scale wave perturbations.              \note was waverand
  greal phi_q = 0.0;         //!< Phase angle portion of the wave angle in equation 48.              \note was phidens
  greal A_Q = 0.0;           //!< Amplitude factor in equation 49.                                   \note was ampfact
  greal a_v = 0.0;           //!< Randomized factor in equation 51.
  greal b_v = 0.045;         //!< Constant factor in equation 51.
  greal hwn = 0.0;           //!< Horizontal wave number (m and n in equation 50).
  greal T = 0.0;             //!< A randomized wave period in equation 48.  (east/west)               \note was wper
  greal T_v = 0.0;           //!< A randomized wave period in equation 48.  (north/south)             \note was wperv
  greal ewWaveAngle = 0.0;   //!< Angle appearing as argument to cosine in equation 48. (east/west)
  greal nsWaveAngle = 0.0;   //!< Angle appearing as argument to cosine in equation 48. (north/south)
  greal ewCorrPhase = 0.0;   //!< Phase angle for large-scale e/w component.
  greal nsCorrPhase = 0.0;   //!< Phase angle for large-scale n/s component.

  greal xMeanSmall = 0.0;    //!< Mean value for L_h.                                                 \note was xbarh
  greal zMeanSmall = 0.0;    //!< Mean value for L_z.                                                 \note was xzarh
  greal xLengthSmall = 0.0;  //!< Horizontal small-scale length L_h.                                  \note was xlh
  greal zLengthSmall = 0.0;  //!< Vertical small-scale length L_z.                                    \note was zlh
  greal xSDSmall = 0.0;      //!< Horizontal small-scale standard deviation of L_h.                   \note was sxlh
  greal zSDSmall = 0.0;      //!< Vertical small-scale standard deviation of L_z.                     \note was szlh

  greal latitudePrev = 0.0;  //!< Latitude of last pertubation update.
  greal longitudePrev = 0.0; //!< Longitude of last pertubation update.
  greal heightPrev = 0.0;    //!< Height of last pertubation update.
  greal timePrev = 0.0;      //!< Elapsed time of last pertubation update.

  greal& temperatureAtSurface = earthAtmos.temperatureAtSurface;          //!< \copydoc EarthAtmosphereState::temperatureAtSurface
  greal& temperatureSDAtSurface = earthAtmos.temperatureSDAtSurface;      //!< \copydoc EarthAtmosphereState::temperatureSDAtSurface
  greal& ewWindAtSurface = earthAtmos.ewWindAtSurface;                    //!< \copydoc EarthAtmosphereState::ewWindAtSurface
  greal& nsWindAtSurface = earthAtmos.nsWindAtSurface;                    //!< \copydoc EarthAtmosphereState::nsWindAtSurface
  greal& ewWindSDAtSurface = earthAtmos.ewWindSDAtSurface;                //!< \copydoc EarthAtmosphereState::ewWindSDAtSurface
  greal& nsWindSDAtSurface = earthAtmos.nsWindSDAtSurface;                //!< \copydoc EarthAtmosphereState::nsWindSDAtSurface
  greal& windCorrelationAtSurface = earthAtmos.windCorrelationAtSurface;  //!< \copydoc EarthAtmosphereState::windCorrelationAtSurface
  greal& windSpeedAtSurface = earthAtmos.windSpeedAtSurface;              //!< \copydoc EarthAtmosphereState::windSpeedAtSurface
  greal& windSpeedSDAtSurface = earthAtmos.windSpeedSDAtSurface;          //!< \copydoc EarthAtmosphereState::windSpeedSDAtSurface

  // small and large perturbations
  greal& presPertSmall = earthAtmos.presPertSmall;     //!< \copydoc EarthAtmosphereState::presPertSmall
  greal& densPertSmall = earthAtmos.densPertSmall;     //!< \copydoc EarthAtmosphereState::densPertSmall
  greal& tempPertSmall = earthAtmos.tempPertSmall;     //!< \copydoc EarthAtmosphereState::tempPertSmall
  greal& ewWindPertSmall = earthAtmos.ewWindPertSmall; //!< \copydoc EarthAtmosphereState::ewWindPertSmall
  greal& nsWindPertSmall = earthAtmos.nsWindPertSmall; //!< \copydoc EarthAtmosphereState::nsWindPertSmall
  greal& presPertLarge = earthAtmos.presPertLarge;     //!< \copydoc EarthAtmosphereState::presPertSmall
  greal& densPertLarge = earthAtmos.densPertLarge;     //!< \copydoc EarthAtmosphereState::densPertSmall
  greal& tempPertLarge = earthAtmos.tempPertLarge;     //!< \copydoc EarthAtmosphereState::tempPertSmall
  greal& ewWindPertLarge = earthAtmos.ewWindPertLarge; //!< \copydoc EarthAtmosphereState::ewWindPertSmall
  greal& nsWindPertLarge = earthAtmos.nsWindPertLarge; //!< \copydoc EarthAtmosphereState::nsWindPertSmall

  // small and large perturbations std devs
  greal& presSDSmall = earthAtmos.presSDSmall;       //!< \copydoc EarthAtmosphereState::presSDSmall
  greal& densSDSmall = earthAtmos.densSDSmall;       //!< \copydoc EarthAtmosphereState::densSDSmall
  greal& tempSDSmall = earthAtmos.tempSDSmall;       //!< \copydoc EarthAtmosphereState::tempSDSmall
  greal& ewWindSDSmall = earthAtmos.ewWindSDSmall;   //!< \copydoc EarthAtmosphereState::ewWindSDSmall
  greal& nsWindSDSmall = earthAtmos.nsWindSDSmall;   //!< \copydoc EarthAtmosphereState::nsWindSDSmall
  greal& presSDLarge = earthAtmos.presSDLarge;       //!< \copydoc EarthAtmosphereState::presSDSmall
  greal& densSDLarge = earthAtmos.densSDLarge;       //!< \copydoc EarthAtmosphereState::densSDSmall
  greal& tempSDLarge = earthAtmos.tempSDLarge;       //!< \copydoc EarthAtmosphereState::tempSDSmall
  greal& ewWindSDLarge = earthAtmos.ewWindSDLarge;   //!< \copydoc EarthAtmosphereState::ewWindSDSmall
  greal& nsWindSDLarge = earthAtmos.nsWindSDLarge;   //!< \copydoc EarthAtmosphereState::nsWindSDSmall

  greal& pressureStandardDeviation = atmos.pressureStandardDeviation;        //!< \copydoc AtmosphereState::pressureStandardDeviation
  greal& temperatureStandardDeviation = atmos.temperatureStandardDeviation;  //!< \copydoc AtmosphereState::temperatureStandardDeviation

  greal& pressurePerturbation = earthAtmos.pressurePerturbation;         //!< \copydoc EarthAtmosphereState::pressurePerturbation
  greal& temperaturePerturbation = earthAtmos.temperaturePerturbation;   //!< \copydoc EarthAtmosphereState::temperaturePerturbation

  greal& perturbedPressure = earthAtmos.perturbedPressure;        //!< \copydoc EarthAtmosphereState::perturbedPressure
  greal& perturbedTemperature = earthAtmos.perturbedTemperature;  //!< \copydoc EarthAtmosphereState::perturbedTemperature

  int& severityLevel = earthAtmos.severityLevel;  //!< \copydoc EarthAtmosphereState::severityLevel

  greal& windCorrelation = earthAtmos.windCorrelation;                       //!< \copydoc EarthAtmosphereState::windCorrelation
  greal& windSpeed = earthAtmos.windSpeed;                                   //!< \copydoc EarthAtmosphereState::windSpeed
  greal& windSpeedStandardDeviation = earthAtmos.windSpeedStandardDeviation; //!< \copydoc EarthAtmosphereState::windSpeedStandardDeviation

  greal& vaporPressure = earthAtmos.vaporPressure;               //!< \copydoc EarthAtmosphereState::vaporPressure
  greal& vaporPressureSD = earthAtmos.vaporPressureSD;           //!< \copydoc EarthAtmosphereState::vaporPressureSD
  greal& vaporDensity = earthAtmos.vaporDensity;                 //!< \copydoc EarthAtmosphereState::vaporDensity
  greal& vaporDensitySD = earthAtmos.vaporDensitySD;             //!< \copydoc EarthAtmosphereState::vaporDensitySD
  greal& dewPoint = earthAtmos.dewPoint;                         //!< \copydoc EarthAtmosphereState::dewPoint
  greal& dewPointSD = earthAtmos.dewPointSD;                     //!< \copydoc EarthAtmosphereState::dewPointSD
  greal& relativeHumidity = earthAtmos.relativeHumidity;         //!< \copydoc EarthAtmosphereState::relativeHumidity
  greal& relativeHumiditySD = earthAtmos.relativeHumiditySD;     //!< \copydoc EarthAtmosphereState::relativeHumiditySD

  int& landCode = earthAtmos.landCode;                           //!< \copydoc EarthAtmosphereState::landCode
  greal& surfaceRoughness = earthAtmos.surfaceRoughness;         //!< \copydoc EarthAtmosphereState::surfaceRoughness
  greal& solarElevation = earthAtmos.solarElevation;             //!< \copydoc EarthAtmosphereState::solarElevation
  greal& elevationAtMidnight = earthAtmos.elevationAtMidnight;   //!< \copydoc EarthAtmosphereState::elevationAtMidnight
  greal& elevationAtNoon = earthAtmos.elevationAtNoon;           //!< \copydoc EarthAtmosphereState::elevationAtNoon
  greal& netRadiationIndex = earthAtmos.netRadiationIndex;       //!< \copydoc EarthAtmosphereState::netRadiationIndex
  greal& stability = earthAtmos.stability;                       //!< \copydoc EarthAtmosphereState::stability
  greal& inverseLength = earthAtmos.inverseLength;               //!< \copydoc EarthAtmosphereState::inverseLength
  greal& frictionVelocity = earthAtmos.frictionVelocity;         //!< \copydoc EarthAtmosphereState::frictionVelocity
  greal& BVFrequencySquare = earthAtmos.BVFrequencySquare;       //!< \copydoc EarthAtmosphereState::BVFrequencySquare
  greal& metersAboveSurface = earthAtmos.metersAboveSurface;     //!< \copydoc EarthAtmosphereState::metersAboveSurface
  greal& sigmaRatio = earthAtmos.sigmaRatio;                     //!< \copydoc EarthAtmosphereState::sigmaRatio
  greal& sigmaW = earthAtmos.sigmaW;                             //!< \copydoc EarthAtmosphereState::sigmaW
  greal& solarDays = earthAtmos.solarDays;                       //!< \copydoc EarthAtmosphereState::solarDays
  greal& boundaryLayerDepth = earthAtmos.boundaryLayerDepth;     //!< \copydoc EarthAtmosphereState::boundaryLayerDepth
  greal& unstableBLFactor = earthAtmos.unstableBLFactor;         //!< \copydoc EarthAtmosphereState::unstableBLFactor
  greal& neutralBoundaryLayerDepth = earthAtmos.neutralBoundaryLayerDepth;       //!< \copydoc EarthAtmosphereState::neutralBoundaryLayerDepth
  greal& perturbedWindSpeedAtSurface = earthAtmos.perturbedWindSpeedAtSurface;   //!< \copydoc EarthAtmosphereState::perturbedWindSpeedAtSurface

  // Input Parameters
  int& year = inputParameters.year;               //!< \copydoc EarthInputParameters::year
  int& month = inputParameters.month;             //!< \copydoc EarthInputParameters::month
  int& day = inputParameters.day;                 //!< \copydoc EarthInputParameters::day
  int& hour = inputParameters.hour;               //!< \copydoc EarthInputParameters::hour
  int& minute = inputParameters.minute;           //!< \copydoc EarthInputParameters::minute
  greal& seconds = inputParameters.seconds;       //!< \copydoc EarthInputParameters::seconds
  bool& patchy = inputParameters.patchy;          //!< \copydoc EarthInputParameters::patchy
  bool& initPerturbations = inputParameters.initializePerturbations;                        //!< \copydoc EarthInputParameters::initializePerturbations
  greal& initialDensityPerturbation = inputParameters.initialDensityPerturbation;           //!< \copydoc EarthInputParameters::initialDensityPerturbation
  greal& initialTemperaturePerturbation = inputParameters.initialTemperaturePerturbation;   //!< \copydoc EarthInputParameters::initialTemperaturePerturbation
  greal& initialEWWindPerturbation = inputParameters.initialEWWindPerturbation;             //!< \copydoc EarthInputParameters::initialEWWindPerturbation
  greal& initialNSWindPerturbation = inputParameters.initialNSWindPerturbation;             //!< \copydoc EarthInputParameters::initialNSWindPerturbation
  greal& initialVerticalWindPerturbation = inputParameters.initialVerticalWindPerturbation; //!< \copydoc EarthInputParameters::initialVerticalWindPerturbation
  greal& userSurfaceRoughness = inputParameters.surfaceRoughness;                           //!< \copydoc EarthInputParameters::surfaceRoughness \note was z0in

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(EarthAtmosphere, initializeData);
  FRIEND_TEST(EarthAtmosphere, updateSurfaceHeight);
  FRIEND_TEST(EarthAtmosphere, interpolateRSData);
  FRIEND_TEST(EarthAtmosphere, interpolatePCData);
  FRIEND_TEST(EarthAtmosphere, getUS76StandardAtmosphere);
#endif // GRAM_UNIT_TEST
};

} // namespace
