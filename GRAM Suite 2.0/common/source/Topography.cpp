//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <algorithm>
#include <cmath>
#include "Topography.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
Topography::Topography()
{
//  initializeData();
}

//! \fn  Topography::Topography(const Topography& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  Topography::~Topography()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

//! \brief  Computes lower bounding  areoid radius data array indices and relative displacements in
//!
//! Computes lower bounding  areoid radius data array indices and relative displacements in
//! latitude and longitude for given input latitude and longitude. This method is optimized to only update
//! the indices if the latitude or longitude is different from the previous call.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//!
//! \retval p     The indices and displacements for the latitude and longitude.
void Topography::updateIndices(greal lat, greal lon)
{
  // Make sure lon is within bounds
  lon = wrapDegrees(lon);

  // Only update the longitude indices if a new longitude is requested.  Otherwise use the cached values.
  if (lon != previousLongitude) {
    previousLongitude = lon;

    greal normalizedLon = lon / lonStep;
    // Use the lower bound index
    size_t index = size_t(floor(normalizedLon));
    // ensure lon index is in proper bounds
    lonIndex = clampSize(index, lonSize);

    // Relative displacement from the index.
    lonDisp = normalizedLon - greal(lonIndex);
  }

  // Only update the latitude indices if a new latitude is requested.  Otherwise use the cached values.
  if (lat != previousLatitude) {
    previousLatitude = lat;

    // Shift so values are between 0 and 180.
    greal shiftedLat = lat + 90.0;

    greal halfLatStep = latStep * 0.5;

    // South Pole case occurs when latitude is within a half step of -90.
    if (shiftedLat < halfLatStep) {
      // The array data only takes a half step for index range 0 to 1 (lat -90 to -89.75).
      latIndex = 0;

      // Relative displacement is 0 for -90 and 1 for -90+halfLonStep.
      latDisp = (shiftedLat) / halfLatStep;
    }
    //=================================================================================================//
    // Compute latitude index and relative displacement for case near N pole                           //
    //-------------------------------------------------------------------------------------------------//
    // North Pole case occurs when latitude is within a half step of 90
    else if (lat > 90.0 - halfLatStep) {
      // The array data only takes a half step for index range 360 to 361 (lat 89.75 to 90).
      // set mola areoid radius lat index to 360 (to allow index + 1 = 361, uppermost grid index)
      latIndex = latSize - 2;

      // Relative displacement is 0 for 90-halfLonStep and 1 for 90.
      latDisp = (lat - 90.0 + halfLatStep) / halfLatStep;
    }
    //=================================================================================================//
    // Compute latitude index and relative displacement for case not near either pole                  //
    //-------------------------------------------------------------------------------------------------//
    // The array data takes whole steps for index range 1 to 360 (lat -89.75 to 89.75).
    else {
      // offset data so that -89.75 maps to 1 and 89.75 maps to 360.  Offset = 90 + MOLA_HALF_STEP
      greal offsetLat = shiftedLat + halfLatStep;
      // Normalize to the index range (1 to 360).
      greal normalizedLat = offsetLat / latStep;
      // Use the lower bound index
      size_t index = size_t(floor(normalizedLat));
      // ensure lat index is in proper bounds
      latIndex = clampSize(index, latSize - 1);

      // Relative displacement from the index.
      latDisp = normalizedLat - greal(latIndex);
    }
  }
}


//! \brief Interpolates  topographic height data to current latitude and longitude
//!
//! This method will update the appropriate lookup indices and displacements which it will use to
//! perform 2D linear interpolation on the  topographic height data.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//! \returns  topographic height \units{km}.
double Topography::getTopographicHeight(greal lat, greal lon)
{
  if (getTopographicHeightCallback) {
    return getTopographicHeightCallback(lat, lon, callbackDataPointer);
  }

  // MIBIndices must be updated prior to interpolation.
  updateIndices(lat, lon);

  // use 2D linear interpolation to get (and return) the topo height at current position
  Interpolator interp2d(latDisp, lonDisp);
  return interp2d.linear(
    topoHeight[latIndex    ][lonIndex], topoHeight[latIndex    ][lonIndex + 1],
    topoHeight[latIndex + 1][lonIndex], topoHeight[latIndex + 1][lonIndex + 1]);
}


//! \brief Reads  radius, topological height, and albedo data.
//!
//! This method reads _data.bin, a binary data file, to populate the data arrays for the
//!  radius, topological height, and albedo data.  The data is stored in static memory.
//! The data format consists three data blocks.  Each data block is preceeded by a size_t with
//! the number of doubles in the data block.
//!
//! \b Inputs
//! \arg #isInitialized
//!
//! \retval #topoHeight
void Topography::initializeData(const char* topoData)
{
  // Only perform initialization once.
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // Skip the label records
  const char* dataBlock = reinterpret_cast<const char*>(topoData + recordBytes * labelRecords);

  union {
    double value;
    char bytes[8];
  } x;
  x.bytes[7] = dataBlock[0];
  x.bytes[6] = dataBlock[1];
  x.bytes[5] = dataBlock[2];
  x.bytes[4] = dataBlock[3];
  x.bytes[3] = dataBlock[4];
  x.bytes[2] = dataBlock[5];
  x.bytes[1] = dataBlock[6];
  x.bytes[0] = dataBlock[7];

  // Convert data block to indexed 2D vector.
  topoHeight.resize(latSize);
  size_t lonSizeMinus1 = lonSize - 1;
  for (size_t line = 0; line < lines; ++line) {
    size_t row = (latSize - line - 2);
    topoHeight[row].resize(lonSize);
    auto& topoHeightRow = topoHeight[row];
    for (size_t sample = 0; sample < lonSizeMinus1; ++sample) {
      size_t col = (sample + sampleOffset) % lonSizeMinus1;
      size_t index = (sample + line * lonSizeMinus1) * 8;
      for (int i = 0; i < 8; ++i) {
        x.bytes[8 - i - 1] = dataBlock[index + i];
      }
      topoHeightRow[col] = x.value;
    }
  }

  // Compute pole heights if not supplied.  Use average of nearest lat.
  if (!usePoleHeights) {
    southPoleHeight = 0.0;
    northPoleHeight = 0.0;
    for (size_t i = 0; i < lonSizeMinus1; ++i) {
      southPoleHeight += topoHeight[1][i];
      northPoleHeight += topoHeight[latSize - 2][i];
    }
    southPoleHeight /= lonSizeMinus1;
    northPoleHeight /= lonSizeMinus1;
  }

  topoHeight[0].resize(lonSize);
  topoHeight[latSize - 1].resize(lonSize);
  // Fill in pole heights.
  for (size_t i = 0; i < lonSize; ++i) {
    topoHeight[0][i] = southPoleHeight;
    topoHeight[latSize - 1][i] = northPoleHeight;
  }

  // Copy 0 to 360.
  for (size_t i = 0; i < latSize; ++i) {
    topoHeight[i][lonSizeMinus1] = topoHeight[latSize - 1][0];
  }

}

//! \fn  Topography::setTopographicHeightCallback()
//! \brief Allow for override of  topographic height computations.
//!
//! This function sets an override callback for getTopographicHeight. The override must match the
//! signature `double(*TopoCallback)(double, double, void*)` where the arguments are the desired
//! latitude, longitude, and a data pointer.  The data pointer is provided to allow the developer
//! to pass data to the override function via setCallbackData(). The return value must be the 
//! value of the topographic height at the specified latitude and longitude.
//!
//! \param callback A pointer to the callback function or nullptr to disable.

//! \fn  Topography::setCallbackData()
//! \brief Set pointer to user data needed in overrides.
//!
//! This method allows a developer to pass data to the override functions set by 
//! setTopographicHeightCallback().
//!
//! \param dataPointer A pointer to a block of data.



} // namespace
