//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// MERRA-2 dataset implementation
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

//! \brief The Modern-Era Retrospective Analysis for Research and Applications, Version 2 (MERRA-2) lower atmosphere model.
//!
//! Along with the enhancements in the meteorological assimilation, MERRA-2 takes some 
//! significant steps towards GMAO's target of an Earth System reanalysis. MERRA-2 is the 
//! first long-term global reanalysis to assimilate space-based observations of aerosols 
//! and represent their interactions with other physical processes in the climate system. 
//!
//! \ingroup EarthGRAM
class MERRA2 : public Atmosphere, public EarthCommon {
public:
  MERRA2();
  MERRA2(const MERRA2& orig);
  virtual ~MERRA2() override;

  void setInputParameters(const EarthInputParameters& params);
  void setHour(int M2hr) { M2Hour = M2hr;  initialized = false; }
  void setExtents(greal latMin, greal latMax, greal lonMin, greal lonMax)
       { minimumLatitude = latMin; maximumLatitude = latMax; 
         minimumLongitude = lonMin; maximumLongitude = lonMax; }

  void update() override;
  void updateSurface();

  void getHeights(greal lat, greal lon, greal& lowest0p3mb, greal& highest0p1mb);
  greal getTemperatureGradient() { return temperatureGradient; }

private:
  static const size_t presSize = 42;  //!< Number of pressure levels.
  size_t latFullSize = 361;  //!< Size of the latitude dimension.
  size_t lonFullSize = 576;  //!< Size of the longitude dimension.
  size_t lonSurfaceFullSize = 540; //!< Size of the longitude dimension for surface variables.
  size_t arraySize3D = presSize * latFullSize * lonFullSize; //!< Size of the 3D arrays allocated as 1D.
  size_t arraySize2D = latFullSize * lonFullSize; //!< Size of the 2D arrays allocated as 1D.
  size_t arraySize2DSurface = latFullSize * lonSurfaceFullSize;

  // Used in reading data files
  size_t latStartIndex = 0;
  size_t latSize = latFullSize;
  size_t lonStartIndex = 0;
  double lonShift = 0;
  size_t lonSize = lonFullSize;
  size_t lonSurfaceStartIndex = 0;
  size_t lonSurfaceSize = lonSurfaceFullSize;
  double lonSurfaceShift = 0;
  size_t lonData[2] = { 0, 0 };
  size_t lonSkip[2] = { 0, 0 };
  size_t lonSurfaceData[2] = { 0, 0 };
  size_t lonSurfaceSkip[2] = { 0, 0 };

  // Used in interpolation update
  greal gridLat = 0.0;;
  greal gridLon = 0.0;;
  size_t latIndex = 0;
  size_t lonIndex = 0;
  size_t lonSurfaceIndex = 0;
  greal lonFactor = 0.0;
  greal latFactor = 0.0;
  greal lonSurfaceFactor = 0.0;

  std::string M2FileName;

  inline size_t getIndex(size_t ipres, size_t ilat, size_t ilon) { return (ipres * latSize + ilat) * lonSize + ilon; }
  inline size_t getIndex(size_t ilat, size_t ilon) { return ilat * lonSize + ilon; }
  inline size_t getSurfaceIndex(size_t ilat, size_t ilon) { return ilat * lonSurfaceSize + ilon; }

  greal getGasConstant(greal td, greal p);
  void updateIndices();
  void readM2File();
  void readBlock(std::ifstream& file, greal* dataBlock, size_t blockSize);
  void readInfoFile();
  void computeArraySizes();
  void allocateArrays();
  void deleteArrays();
  void buildDataFileName();

  // Input Parameters
  int M2Hour = 0;                           //!< MERRA-2 hour code.
  greal minimumLatitude = -90.0;            //!< Boundary for restricted read of MERRA-2 data.
  greal maximumLatitude = 90.0;             //!< Boundary for restricted read of MERRA-2 data.
  greal minimumLongitude = 0.0;             //!< Boundary for restricted read of MERRA-2 data.
  greal maximumLongitude = 360.0;           //!< Boundary for restricted read of MERRA-2 data.
  int hour = 0;                             //!< Start time hour.
  int month = 0;                            //!< Start time month.
  greal ewWindPerturbationScale = 1.0;      //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;      //!< Scales North/South wind perturbations.
  greal densityPerturbationScale = 1.0;     //!< Scales density standard deviations.
  greal pressurePerturbationScale = 1.0;    //!< Scales pressure standard deviations.
  greal temperaturePerturbationScale = 1.0; //!< Scales temperature standard deviations.

  // This member must be of type double to read the binary files.
  double* blockOfDoubles = nullptr; //!< Used in reading binary NCEP files.
  bool initialized = false;         //!< Set to true when NCEP data has been initialized.
  greal temperatureGradient = 0.0;  //!< Temperature gradient with respect to height. \note was dtz
  bool surfaceUpdate = false;
  greal hp1, hp3;

  // Pressures (mb).
  static const greal pn[presSize]; //!< Pressure levels \units{Pa}.

  // Dynamically allocated
  // 3D arrays (pres, lat, lon)
  greal* temp = nullptr;  //!< Temperature.
  greal* dens = nullptr;  //!< Density.
  greal* dewp = nullptr;  //!< Dew point.
  greal* uwnd = nullptr;  //!< E/W wind.
  greal* vwnd = nullptr;  //!< N/S wind.
  greal* hgt = nullptr;   //!< Edge height.
  greal* stmp = nullptr;  //!< Temperature standard deviation.
  greal* sden = nullptr;  //!< Density standard deviation.
  greal* spres = nullptr; //!< Pressure standard deviation.
  greal* sdewp = nullptr; //!< Dew point standard deviation.
  greal* suwd = nullptr;  //!< E/W wind standard deviation.
  greal* svwd = nullptr;  //!< N/S wind standard deviation.
  greal* sprs = nullptr;  //!< Pressure standard deviation.
  greal* rhum = nullptr;  //!< Relative humidity.
  greal* srhum = nullptr; //!< Relative humidity standard deviation.
  greal* vprs = nullptr;  //!< Vapor pressure.
  greal* svprs = nullptr; //!< Vapor pressure standard deviation.
  greal* spdav = nullptr; //!< Average wind speed.
  greal* spdsd = nullptr; //!< Wind speed standard deviation.
  greal* uvcor = nullptr; //!< Wind correlation.
  // 2D arrays (lat, lon)
  greal* psfc = nullptr;  //!< Surface-level pressure.
  greal* spsfc = nullptr; //!< Surface-level pressure standard deviation.
  greal* usfc = nullptr;
  greal* susfc = nullptr;
  greal* vsfc = nullptr;
  greal* svsfc = nullptr;
  greal* tsfc = nullptr;
  greal* stsfc = nullptr;
  greal* dsfc = nullptr;
  greal* sdsfc = nullptr;
  greal* tdsfc = nullptr;
  greal* stdsfc = nullptr;
  greal* vpsfc = nullptr;
  greal* svpsfc = nullptr;
  greal* ghsfc = nullptr;
  greal* spdavsfc = nullptr;
  greal* spdsdsfc = nullptr;
  greal* uvcorsfc = nullptr;
  greal* rhsfc = nullptr;
  greal* srhsfc = nullptr;

  EarthAtmosphereState earthAtmos; //!< Earth specific metrics added to the AtmosphereState.

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
  FRIEND_TEST(MERRA2, getIndex);
  FRIEND_TEST(MERRA2, readM2File);
  FRIEND_TEST(MERRA2, DISABLED_readM2File_dataPoints);
  FRIEND_TEST(MERRA2, getGasConstant);
  FRIEND_TEST(MERRA2, DISABLED_getHeights);
  FRIEND_TEST(MERRA2, updateIndices_basic);
  FRIEND_TEST(MERRA2, updateIndices_complex);
#endif // GRAM_UNIT_TEST
};

} // namespace GRAM
