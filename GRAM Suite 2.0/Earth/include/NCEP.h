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

#include <fstream>
#include <string>
#include "unittest_friend.h"
#include "gram.h"
#include "EarthCommon.h"
#include "Atmosphere.h"
#include "EarthInputParameters.h"
#include "EarthAtmosphereState.h"

namespace GRAM {

//! \brief The National Centers for Environmental Prediction (NCEP) model of a Earth lower atmosphere.
//!
//! The NCEP/NCAR Reanalysis data set is a continually updated (1948-present) globally gridded 
//! data set that represents the state of the Earth's atmosphere, incorporating observations and 
//! numerical weather prediction model output from 1948 to present.
//!
//! \ingroup EarthGRAM
class NCEP : public Atmosphere, public EarthCommon
{
public:
	NCEP();
  NCEP(const NCEP& orig);
  virtual ~NCEP() override;

  void setInputParameters(const EarthInputParameters& params);
  void setYearAndHour(int NCEPyr, int NCEPhr) { NCEPYear = NCEPyr; NCEPHour = NCEPhr; initialized = false; }

  void update() override;
  void updateSurface();

  void getHeights(greal lat, greal lon, greal &highest20mb, greal &lowest10mb);
  greal getTemperatureGradient() { return temperatureGradient; }

private:
  static const size_t presSize = 19;                               //!< Number of pressure levels.
  static const size_t latSize = 73;                                //!< Size of the latitude dimension.
  static const size_t lonSize = 145;                               //!< Size of the longitude dimension.
  static const size_t arraySize3D = presSize * latSize * lonSize;  //!< Size of the 3D arrays allocated as 1D.
  static const size_t arraySize2D = latSize * lonSize;             //!< Size of the 2D arrays allocated as 1D.

  greal getGasConstant(greal td, greal p);
  greal getGeopotentialHeight(greal g0, greal r, greal z);
  inline size_t getIndex(size_t ipres, size_t ilat, size_t ilon) { return (ipres * latSize + ilat) * lonSize + ilon; }
  inline size_t getIndex(size_t ilat, size_t ilon) { return ilat * lonSize + ilon; }
  void sortLevel(size_t lsort[2][2][presSize], greal zsort[2][2][presSize]);
  void readNCEPFile();
  void readBlock(std::ifstream& file, greal* dataBlock, size_t blockSize);

  // Input Parameters
  int NCEPYear = 0;                         //!< NCEP year code.
  int NCEPHour = 0;                         //!< NCEP hour code.
  int hour = 0;                             //!< Start time hour.
  int month = 0;                            //!< Start time month.
  greal ewWindPerturbationScale = 1.0;      //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;      //!< Scales North/South wind perturbations.
  greal densityPerturbationScale = 1.0;     //!< Scales density standard deviations.
  greal pressurePerturbationScale = 1.0;    //!< Scales pressure standard deviations.
  greal temperaturePerturbationScale = 1.0; //!< Scales temperature standard deviations.

  // This member must be of type double to read the binary files.
  double* blockOfDoubles = nullptr;  //!< Used in reading binary NCEP files.
  bool initialized = false;          //!< Set to true when NCEP data has been initialized.
  greal temperatureGradient = 0.0;   //!< Temperature gradient with respect to height. \note was dtz
  greal dtdx[presSize] = { 0.0 };    //!< Temperature gradient with respect to longitude by pressure level.
  greal dtdy[presSize] = { 0.0 };    //!< Temperature gradient with respect to latitude by pressure level.
  bool surfaceUpdate = false;

  // Pressures (mb).
  static const greal pn[presSize];   //!< Pressure levels \units{Pa}.

  // Dynamically allocated (in constructor)
  // 3D arrays (pres, lat, lon)
  greal*  temp = nullptr;  //!< Temperature.
  greal*  dens = nullptr;  //!< Density.
  greal*  dewp = nullptr;  //!< Dew point.
  greal*  uwnd = nullptr;  //!< E/W wind.
  greal*  vwnd = nullptr;  //!< N/S wind.
  greal*  geop = nullptr;  //!< Geopotential height.
  greal*  stmp = nullptr;  //!< Temperature standard deviation.
  greal*  sden = nullptr;  //!< Density standard deviation.
  greal* sdewp = nullptr;  //!< Dew point standard deviation.
  greal*  suwd = nullptr;  //!< E/W wind standard deviation.
  greal*  svwd = nullptr;  //!< N/S wind standard deviation.
  greal*  sprs = nullptr;  //!< Pressure standard deviation.
  greal*  rhum = nullptr;  //!< Relative humidity.
  greal* srhum = nullptr;  //!< Relative humidity standard deviation.
  greal*  vprs = nullptr;  //!< Vapor pressure.
  greal* svprs = nullptr;  //!< Vapor pressure standard deviation.
  greal* spdav = nullptr;  //!< Average wind speed.
  greal* spdsd = nullptr;  //!< Wind speed standard deviation.
  greal* uvcor = nullptr;  //!< Wind correlation.
  // 2D arrays (lat, lon)
  greal* slpn = nullptr;   //!< Sea-level pressure.
  greal* sfcp = nullptr;   //!< Surface-level pressure.
  greal* sslp = nullptr;   //!< Sea-level pressure standard deviation.
  greal* ssfcp = nullptr;  //!< Surface-level pressure standard deviation.

  EarthAtmosphereState earthAtmos;  //!< Earth specific metrics added to the AtmosphereState.

  greal& pressureStandardDeviation = atmos.pressureStandardDeviation;         //!< \copydoc AtmosphereState::pressureStandardDeviation        \note was spz
  greal& temperatureStandardDeviation = atmos.temperatureStandardDeviation;   //!< \copydoc AtmosphereState::temperatureStandardDeviation     \note was stz
  greal& dewPoint = earthAtmos.dewPoint;                                      //!< \copydoc EarthAtmosphereState::dewPoint                    \note was tdz
  greal& dewPointSD = earthAtmos.dewPointSD;                                  //!< \copydoc EarthAtmosphereState::dewPointSD                  \note was stdz
  greal& relativeHumidity = earthAtmos.relativeHumidity;                      //!< \copydoc EarthAtmosphereState::relativeHumidity            \note was rhn
  greal& relativeHumiditySD = earthAtmos.relativeHumiditySD;                  //!< \copydoc EarthAtmosphereState::relativeHumiditySD          \note was srhn
  greal& vaporPressure = earthAtmos.vaporPressure;                            //!< \copydoc EarthAtmosphereState::vaporPressure               \note was vpn
  greal& vaporPressureSD = earthAtmos.vaporPressureSD;                        //!< \copydoc EarthAtmosphereState::vaporPressureSD             \note was svpn
  greal& windSpeed = earthAtmos.windSpeed;                                    //!< \copydoc EarthAtmosphereState::windSpeed                   \note was spdavz
  greal& windSpeedStandardDeviation = earthAtmos.windSpeedStandardDeviation;  //!< \copydoc EarthAtmosphereState::windSpeedStandardDeviation  \note was spdsdz
  greal& windCorrelation = earthAtmos.windCorrelation;                        //!< \copydoc EarthAtmosphereState::windCorrelation             \note was uvcorrz

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(NCEP, ncepmd);
  FRIEND_TEST(NCEP, getIndex);
  FRIEND_TEST(NCEP, readNCEPFile);
  FRIEND_TEST(NCEP, wexler);
  FRIEND_TEST(NCEP, getGasConstant);
  FRIEND_TEST(NCEP, getGeopotentialHeight);
  FRIEND_TEST(NCEP, getHeights);
  FRIEND_TEST(NCEP, sortLevel);
#endif // GRAM_UNIT_TEST
};

} // namespace
