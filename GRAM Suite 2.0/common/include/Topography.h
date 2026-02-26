//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "unittest_friend.h"
#include "Position.h"

namespace GRAM {

typedef double(*TopoCallback)(double, double, void*);

//! \brief Computes topographical height based off of PDS data.
//!
//! The Topography class will compute the topographic height,
//! given a latitude and longitude pair. 
//!
//! \ingroup CommonGRAM
class Topography 
{
public:
  Topography();
  Topography(const Topography& orig) = default;
  virtual ~Topography() = default;

  virtual double getTopographicHeight(greal lat, greal lon);
  void setTopographicHeightCallback(TopoCallback callback) { getTopographicHeightCallback = callback; }
  void setCallbackData(void* dataPointer) { callbackDataPointer = dataPointer; }

protected:
  // Data parameters
  size_t latSize = 0;          //!< Size of the latitude dimesion in the data array.
  greal latStep = 1.0;         //!< Angle step size in the data array.
  size_t lonSize = 0;          //!< Size of the longitude dimesion in the data array.
  greal lonStep = 1.0;         //!< Angle step size in the data array.
  size_t recordBytes = 0;      //!< Length of each record in PDS file.
  size_t fileRecords = 0;      //!< Total number of fixed length records in PDS file.
  size_t labelRecords = 0;     //!< Number of non-data header records in PDS file.
  size_t lines = 0;            //!< Number of latitude data records in PDS file.
  size_t lineSamples = 0;      //!< Number of longitude entries per record in PDS file.
  size_t sampleOffset = 0;     //!< Offset of first longitude entry in PDS file.
  bool eastPositive = true;    //!< Longitude flag for data in the PDS file.
  bool usePoleHeights = false; //!< Flag to use specified pole heights or compute an average.
  greal southPoleHeight = 0.0; //!< Surface height to use at the south pole.
  greal northPoleHeight = 0.0; //!< Surface height to use at the north pole.

  // Results
  size_t latIndex = 0;   //!< Latitude index for data interpolation.
  greal latDisp = 0.0;   //!< Latitude relative displacement for data interpolation.
  greal previousLatitude = 999.0;  //!< Caching helps optimize computation.

  size_t lonIndex = 0;   //!< Longitude index for data interpolation.
  greal lonDisp = 0.0;   //!< Longitude relative displacement for data interpolation.
  greal previousLongitude = 999.0;  //!< Caching helps optimize computation.

  void initializeData(const char* topoData);
  void updateIndices(greal lat, greal lon);

  TopoCallback getTopographicHeightCallback = nullptr; //!< Override of getTopographicHeight.
  void* callbackDataPointer = nullptr;                 //!< Data pointer for overrides.

  bool isInitialized = false;                   //!< True if data has been initialized.
  std::vector<std::vector<double>> topoHeight;  //!< The topological height data.


#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(Topography, updateIndices);
  FRIEND_TEST(Topography, updateIndices_oddSizes);
#endif // GRAM_UNIT_TEST
};

} // namespace

