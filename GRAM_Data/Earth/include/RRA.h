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

#include <vector>
#include <string>
#include <iosfwd>
#include "unittest_friend.h"
#include "gram.h"
#include "Atmosphere.h"
#include "EarthInputParameters.h"
#include "EarthAtmosphereState.h"
#include "EarthCommon.h"


namespace GRAM {

//! \brief Reads and processes Range Reference Atmosphere (RRA) database sites.
//!
//! The RRA class determines if the current position falls within limits of one
//! the RRA sites.  If so, fairing is performed for a smooth transition from
//! the model means to the RRA data.
//! \ingroup EarthGRAM
class RRA : public Atmosphere, public EarthCommon
{
public:
	RRA();
  RRA(const RRA& orig);
  virtual ~RRA() = default;

  void setInputParameters(const EarthInputParameters& params);
  void setRRASiteList(const std::string& fileName) { rraSiteList = fileName; }
  void setMonth(int mon) { month = mon; }
  void setRRAParameters(RRAYearType year, greal inner, greal outer);
  void setUseRRA(bool useFlag);

  void setInputState(AtmosphereState& state);
  void setEpsilon(greal eps) { epsilon = eps; }
  void setPpmToND(greal ppm2nd) { ppmToND = ppm2nd; }

  void update() override;
  void updatePsychrometrics();
  void updateWindSpeed();

  void locateNearestSite();

  bool isWeighted() const { return (useRRA && siteIndex >= 0 && rraWeight != 0.0); }
  bool isActive() const { return useRRA; }

  void getFairingHeights(greal& lower, greal& upper) const;

  RRAYearType getRRAYear() const { return rraYear; }

 private:
  void readSiteList();
  void readrra1();
  void readrra2();
  void readrra3();
  void parseTopLine(std::ifstream& rraDataFile, const std::string& fileName);
  void parseToMonth(std::ifstream& rraDataFile);
  greal getHeightWeight(greal h, const std::vector<greal>& z);
  int findHeightIndex(greal h, const std::vector<greal>& z);
  greal getWeightFactor(greal delta);

  EarthAtmosphereState earthAtmos;                                            //!< \copydoc EarthAtmosphere::earthAtmos
  greal& pressureStandardDeviation = atmos.pressureStandardDeviation;         //!< \copydoc AtmosphereState::pressureStandardDeviation
  greal& temperatureStandardDeviation = atmos.temperatureStandardDeviation;   //!< \copydoc AtmosphereState::temperatureStandardDeviation
  greal& temperatureAtSurface = earthAtmos.temperatureAtSurface;              //!< \copydoc EarthAtmosphereState::temperatureAtSurface
  greal& temperatureSDAtSurface = earthAtmos.temperatureSDAtSurface;          //!< \copydoc EarthAtmosphereState::temperatureSDAtSurface
  greal& windSpeedAtSurface = earthAtmos.windSpeedAtSurface;                  //!< \copydoc EarthAtmosphereState::windSpeedAtSurface
  greal& windSpeedSDAtSurface = earthAtmos.windSpeedSDAtSurface;              //!< \copydoc EarthAtmosphereState::windSpeedSDAtSurface
  greal& ewWindAtSurface = earthAtmos.ewWindAtSurface;                        //!< \copydoc EarthAtmosphereState::ewWindAtSurface
  greal& nsWindAtSurface = earthAtmos.nsWindAtSurface;                        //!< \copydoc EarthAtmosphereState::nsWindAtSurface
  greal& ewWindSDAtSurface = earthAtmos.ewWindSDAtSurface;                    //!< \copydoc EarthAtmosphereState::ewWindSDAtSurface
  greal& nsWindSDAtSurface = earthAtmos.nsWindSDAtSurface;                    //!< \copydoc EarthAtmosphereState::nsWindSDAtSurface
  greal& windCorrelationAtSurface = earthAtmos.windCorrelationAtSurface;      //!< \copydoc EarthAtmosphereState::windCorrelationAtSurface
  greal& windCorrelation = earthAtmos.windCorrelation;                        //!< \copydoc EarthAtmosphereState::windCorrelation
  greal& densityAtSurface = earthAtmos.densityAtSurface;                      //!< \copydoc EarthAtmosphereState::densityAtSurface
  greal& densitySDAtSurface = earthAtmos.densitySDAtSurface;                  //!< \copydoc EarthAtmosphereState::densitySDAtSurface
  greal& pressureSDAtSurface = earthAtmos.pressureSDAtSurface;                //!< \copydoc EarthAtmosphereState::pressureSDAtSurface
  greal& vaporPressure = earthAtmos.vaporPressure;                            //!< \copydoc EarthAtmosphereState::vaporPressure
  greal& vaporPressureSD = earthAtmos.vaporPressureSD;                        //!< \copydoc EarthAtmosphereState::vaporPressureSD
  greal& dewPoint = earthAtmos.dewPoint;                                      //!< \copydoc EarthAtmosphereState::dewPoint
  greal& dewPointSD = earthAtmos.dewPointSD;                                  //!< \copydoc EarthAtmosphereState::dewPointSD
  greal& relativeHumidity = earthAtmos.relativeHumidity;                      //!< \copydoc EarthAtmosphereState::relativeHumidity
  greal& relativeHumiditySD = earthAtmos.relativeHumiditySD;                  //!< \copydoc EarthAtmosphereState::relativeHumiditySD
  greal& windSpeed = earthAtmos.windSpeed;                                    //!< \copydoc EarthAtmosphereState::windSpeed
  greal& windSpeedStandardDeviation = earthAtmos.windSpeedStandardDeviation;  //!< \copydoc EarthAtmosphereState::windSpeedStandardDeviation

  AtmosphereState *inState = nullptr;       //!< The input (existing) atmosphere state.
                                            
  const greal SURFACE_LIMIT = 27.0_km;      //!< Surface values are affected below this surface limit.
  const int NO_SITE = -1;                   //!< Constant to designate no site has been located.
  int siteIndex = NO_SITE;                  //!< Index of the last RRA site located.  \note was isite
  std::string siteName;                     //!< Name of the last RRA site located. \note was rra_id
  greal siteLatitude = 0.0;                 //!< Geodetic latitude of the last RRA site located \units{\text{degrees}}. \note was xlat
  greal siteLongitude = 0.0;                //!< Longitude of the last RRA site located \units{\text{degrees}}. \note was xlon
  greal siteWeight = 0.0;                   //!< Horizontal weight factor of the last RRA site located. \note was sitewgt
  int siteReadIndex = NO_SITE;              //!< Index of the last RRA site read into memory. \note was isitered
                                            
  greal presBeforeRRA = 0.0;                //!< Temporary storage for the incoming pressure.
  greal tempBeforeRRA = 0.0;                //!< Temporary storage for the incoming temperature.
                                            
  int surfaceIndex1 = 0;                    //!< Index of z1 value nearest site list surface value.
  int surfaceIndex2 = 0;                    //!< Index of z2 value nearest site list surface value.
  int surfaceIndex3 = 0;                    //!< Index of z3 value nearest site list surface value.
  greal rraWeight1 = 0.0;                   //!< Interpolation weight for wind data.
  greal rraWeight2 = 0.0;                   //!< Interpolation weight for pdt data.
  greal rraWeight3 = 0.0;                   //!< Interpolation weight for psychrometric data.
  greal rraWeight = 0.0;                    //!< For output.
                                            
  greal epsilon = 0.0;                      //!< The ratio between the gas constant for dry air and the gas constant for water vapor.
  greal ppmToND = 0.0;                      //!< Parts per million to number density conversion factor.
//  greal dewPointSurface = 0.0;   // (was tdsrf)   // Not used and not output
//  greal dewPointSDSurface = 0.0; // (was stdsrf)  // Not used and not output

  // namelist data
  bool useRRA = false;                      //!< Flag denoting activation of RRA data.   \note was iurra
  RRAYearType rraYear = RRA_2013;           //!< Denotes with year of RRA data is used.  \note was iyrrra
  greal outerRadius = 0.0;                  //!< Outer distance of RRA fairing ring.     \note was sitelim
  greal innerRadius = 0.0;                  //!< Inner distance of RRA fairing ring.     \note was sitenear
  greal ewWindPerturbationScale = 1.0;      //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;      //!< Scales North/South wind perturbations.
  greal densityPerturbationScale = 1.0;     //!< Scales density standard deviations.
  greal pressurePerturbationScale = 1.0;    //!< Scales pressure standard deviations.
  greal temperaturePerturbationScale = 1.0; //!< Scales temperature standard deviations.
  int month = 0;                            //!< Month selection for RRA data ingestion.
  std::string rraSiteList = "rrasites.txt"; //!< Name of RRA sites file.
  
  // RRA sites list (from rrasites.txt)
  std::vector<std::string> siteNameList;     //!< List of RRA site names.
  std::vector<greal> siteGeodeticLatitudeList;       //!< List of RRA site geodetic latitudes \units{\text{degrees}}.
  std::vector<greal> siteGeocentricLatitudeList;       //!< List of RRA site geocentric latitudes \units{\text{degrees}}.
  std::vector<greal> siteLongitudeList;      //!< List of RRA site longitudes \units{\text{degrees}}.
  std::vector<greal> siteSurfaceHeightList;  //!< List of RRA site surface heights \units{m}.
  std::vector<greal> siteMaxHeightList;      //!< List of RRA site maximum heights \units{km}.

  // Data for current RRA site from the T1 file
  std::vector<greal> z1;       //!< Heights for winds data \units{km}.
  std::vector<greal> u;        //!< East/west winds \units{m/s}.
  std::vector<greal> v;        //!< North/south winds \units{m/s}.
  std::vector<greal> su;       //!< East/west wind standard deviations \units{m/s}.
  std::vector<greal> sv;       //!< North/south wind standard deviations \units{m/s}.
  std::vector<greal> ruvt;     //!< Wind correlation data.  
  std::vector<greal> avspd;    //!< Wind speed data \units{m/s}.
  std::vector<greal> sdspd;    //!< Wind speed standard deviations \units{m/s}.

  // Data for current RRA site from the T2 file
  std::vector<greal> z2;       //!< Heights for pdt data \units{km}.
  std::vector<greal> p;        //!< Pressure data \units{Pa}.
  std::vector<greal> d;        //!< Density data \units{kg/m^3}.
  std::vector<greal> t;        //!< Temperature data \units{K}.
  std::vector<greal> sp;       //!< Pressure standard deviations \units{\%}.
  std::vector<greal> sd;       //!< Density standard deviations \units{\%}.
  std::vector<greal> st;       //!< Temperature standard deviations \units{\%}.

  // Data for current RRA site from the T3 file
  std::vector<greal> z3;       //!< Heights for psychrometric data \units{km}.
  std::vector<greal> vp;       //!< Vapor pressure data \units{Pa}.
  std::vector<greal> td;       //!< Dew point temperature data \units{K}.
  std::vector<greal> svp;      //!< Vapor pressure standard deviations \units{Pa}.
  std::vector<greal> std;      //!< Dew point temperature standard deviations \units{K}.

  // temperature was trra
  // density was drra
  // pressure was prra
  // temperatureStandardDeviation was strra
  // windCorrelation was uvrra
  // windSpeed was avspdrra
  // windSpeedStandardDeviation was sdspdrra

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(RRA, parseTopLine_and_readSiteList);
  FRIEND_TEST(RRA, locateNearestSite);
  FRIEND_TEST(RRA, getWeightFactor);
  FRIEND_TEST(RRA, findHeightIndex);
  FRIEND_TEST(RRA, getHeightWeight);
  FRIEND_TEST(RRA, dedt);
  FRIEND_TEST(RRA, d2edt2);
#endif // GRAM_UNIT_TEST
};

} // namespace
