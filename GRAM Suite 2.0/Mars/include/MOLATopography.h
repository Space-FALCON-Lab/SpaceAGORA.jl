//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "unittest_friend.h"
#include "Position.h"
#include "MarsCommon.h"

namespace GRAM {

//! \brief Computes topographical parameters based off of MOLA data.
//!
//! The MOLATopography class will compute the aeroid radius, the topographic height,
//! or the albedo given a latitude and longitude pair. The results are based on
//! Mars Orbiter Laser Altimeter (MOLA) data.  The file "MOLA_data.bin" is required.
//!
//! \ingroup MarsGRAM
class MOLATopography : public MarsCommon
{
public:
  MOLATopography();
  MOLATopography(const MOLATopography& orig) = default;
  virtual ~MOLATopography() = default;

  double getAreoidRadius(greal lat, greal lon);
  double getTopographicHeight(greal lat, greal lon);
  double getAlbedo(greal lat, greal lon);
  void setAreoidRadiusCallback(TopoCallback callback) { getAreoidRadiusCallback = callback; }
  void setTopographicHeightCallback(TopoCallback callback) { getTopographicHeightCallback = callback; }
  void setCallbackData(void* dataPointer) { callbackDataPointer = dataPointer; }

private:
  struct MOLAParams {
    size_t latSize = 0;    //!< Size of the latitude dimesion in the data array.
    size_t lonSize = 0;    //!< Size of the longitude dimesion in the data array.
    greal stepSize = 0.0;  //!< Angle step size in the data array (same for lat and lon assumed).
    size_t latIndex = 0;   //!< Latitude index for data interpolation.
    size_t lonIndex = 0;   //!< Longitude index for data interpolation.
    greal latDisp = 0.0;   //!< Latitude relative displacement for data interpolation.
    greal lonDisp = 0.0;   //!< Longitude relative displacement for data interpolation.
    greal previousLatitude = 999.0;   //!< The previous latitude used to compute indices.
    greal previousLongitude = 999.0;  //!< The previous longitude used to compute indices.
  };

  void initializeData();
  void updateIndices(greal lat, greal lon, MOLAParams& p);

  TopoCallback getAreoidRadiusCallback = nullptr;      //!< Override of getAreoidRadius.
  TopoCallback getTopographicHeightCallback = nullptr; //!< Override of getTopographicHeight.
  void* callbackDataPointer = nullptr;                 //!< Data pointer for overrides.

  MOLAParams topoParams;   //!< Parameters for aeroid radius and topographic height computations.
  MOLAParams albedoParams; //!< Parameters for albedo computations.

  const static size_t LAT_SIZE = 362;      //!< The size of the latitude dimension for the radius and topo height arrays.
  const static size_t LON_SIZE = 722;      //!< The size of the longitude dimension for the radius and topo height arrays.
  const static size_t ALB_LAT_SIZE = 182;  //!< The size of the latitude dimension for the albedo array.
  const static size_t ALB_LON_SIZE = 362;  //!< The size of the longitude dimension for the albedo array.

  static bool isInitialized;                        //!< True if data has been initialized.
  static double areoRadius[LAT_SIZE][LON_SIZE];     //!< The MOLA radius data.
  static double topoHeight[LAT_SIZE][LON_SIZE];     //!< The MOLA topological height data.
  static double albedo[ALB_LAT_SIZE][ALB_LON_SIZE]; //!< The MOLA albedo data.

#ifdef CONVERT_MARS_DATA
  void writeData();
  void readLegacyData();
#endif

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MOLATopography, updateIndices);
  FRIEND_TEST(MOLATopography, updateAlbedoIndices);
#endif // GRAM_UNIT_TEST
};

} // namespace

