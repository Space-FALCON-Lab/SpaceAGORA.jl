//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <algorithm>
#include <cmath>
#include "MOLATopography.h"
#include "Interpolator.h"
#include "MarsInterpolatorBase.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool MOLATopography::isInitialized = false;
double MOLATopography::areoRadius[LAT_SIZE][LON_SIZE] = {};
double MOLATopography::topoHeight[LAT_SIZE][LON_SIZE] = {};
double MOLATopography::albedo[ALB_LAT_SIZE][ALB_LON_SIZE] = {};

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
MOLATopography::MOLATopography()
  : MarsCommon(NULL)
{
  initializeData();

  // Set the constant parameters for radius and topo height.
  topoParams.latSize = LAT_SIZE;
  topoParams.lonSize = LON_SIZE;
  topoParams.stepSize = 0.5;

  // Set the constant parameters for albedo.
  albedoParams.latSize = ALB_LAT_SIZE;
  albedoParams.lonSize = ALB_LON_SIZE;
  albedoParams.stepSize = 1.0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \fn  MOLATopography::MOLATopography(const MOLATopography& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \fn  MOLATopography::~MOLATopography()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief  Computes lower bounding mola areoid radius data array indices and relative displacements in
//!
//! Computes lower bounding mola areoid radius data array indices and relative displacements in
//! latitude and longitude for given input latitude and longitude. This method is optimized to only update
//! the indices if the latitude or longitude is different from the previous call.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//! \param[out] p A MOLAParams struct.
//!
//! \retval p     The indices and displacements for the latitude and longitude.
void MOLATopography::updateIndices(greal lat, greal lon, MOLAParams& p)
{
  const greal stepSize = p.stepSize;
  // latitude/longitude half step size for MOLA data
  const greal halfStepSize = 0.5 * stepSize;

  // Only update the longitude indices if a new longitude is requested.  Otherwise use the cached values.
  if (lon != p.previousLongitude) {
    p.previousLongitude = lon;

    // Convert to West Postive longitude
    greal lonW = (lon == 0.0) ? lon : 360.0 - lon;

    //=================================================================================================//
    // make necessary adjustments to lon value to agree with IAU 2000 epoch                            //
    // initially: double lon = pos.LON_E + 0.238; but made changes such that input lon is forced W +   //
    //-------------------------------------------------------------------------------------------------//
    lonW = wrapDegrees(lonW + 0.238);

    //=================================================================================================//
    //  compute the lower bounding mola radius data array index and relative displacement              //
    //-------------------------------------------------------------------------------------------------//
    // offset longitude a half step so that 0.25 maps to 1 and 359.75 maps to 720
    greal offsetLon = lonW + halfStepSize;
    // Normalize to the index range (0 to 721).
    greal normalizedLon = offsetLon / stepSize;
    // Use the lower bound index
    size_t index = size_t(floor(normalizedLon));
    // ensure lon index is in proper bounds
    p.lonIndex = clampSize(index, p.lonSize - 1);

    // Relative displacement from the index.
    p.lonDisp = normalizedLon - greal(p.lonIndex);
  }

  // Only update the latitude indices if a new latitude is requested.  Otherwise use the cached values.
  if (lat != p.previousLatitude) {
    p.previousLatitude = lat;
    //=================================================================================================//
    // Compute lower bounding mola radius index and relative displacement for case near S pole         //
    //-------------------------------------------------------------------------------------------------//
    // South Pole case occurs when latitude is within a half step of -90.
    if (lat < halfStepSize - 90.0) {
      // The array data only takes a half step for index range 0 to 1 (lat -90 to -89.75).
      p.latIndex = 0;

      // Relative displacement is 0 for -90 and 1 for -89.75.
      p.latDisp = (lat + 90.0) / halfStepSize;
    }
    //=================================================================================================//
    // Compute latitude index and relative displacement for case near N pole                           //
    //-------------------------------------------------------------------------------------------------//
    // North Pole case occurs when latitude is within a half step of 90
    else if (lat > 90.0 - halfStepSize) {
      // The array data only takes a half step for index range 360 to 361 (lat 89.75 to 90).
      // set mola areoid radius lat index to 360 (to allow index + 1 = 361, uppermost grid index)
      p.latIndex = p.latSize - 2;

      // Relative displacement is 0 for 89.75 and 1 for 90.
      p.latDisp = (lat - 90.0 + halfStepSize) / halfStepSize;
    }
    //=================================================================================================//
    // Compute latitude index and relative displacement for case not near either pole                  //
    //-------------------------------------------------------------------------------------------------//
    // The array data takes whole steps for index range 1 to 360 (lat -89.75 to 89.75).
    else {
      // offset data so that -89.75 maps to 1 and 89.75 maps to 360.  Offset = 90 + MOLA_HALF_STEP
      greal offsetLat = lat + 90.0 + halfStepSize;
      // Normalize to the index range (1 to 360).
      greal normalizedLat = offsetLat / stepSize;
      // Use the lower bound index
      size_t index = size_t(floor(normalizedLat));
      // ensure lat index is in proper bounds
      p.latIndex = clampSize(index, p.latSize - 1);

      // Relative displacement from the index.
      p.latDisp = normalizedLat - greal(p.latIndex);
    }
  }
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Interpolates MOLA areoid radius data to current latitude and longitude
//!
//! This method will update the appropriate lookup indices and displacements which it will 
//! then use to perform 2D linear interpolation on the MOLA areoid radius data.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//!
//! \returns Mola areoid radius \units{km}.
greal MOLATopography::getAreoidRadius(greal lat, greal lon)
{
  if (getAreoidRadiusCallback) {
    return getAreoidRadiusCallback(lat, lon, callbackDataPointer);
  }

  // MIBIndices must be updated prior to interpolation.
  updateIndices(lat, lon, topoParams);

  // use 2D linear interpolation to get (and return) the radius at current position
  Interpolator interp2d(topoParams.latDisp, topoParams.lonDisp);
  size_t ilat = topoParams.latIndex;
  size_t ilon = topoParams.lonIndex;
  return interp2d.linear(
    areoRadius[ilat    ][ilon], areoRadius[ilat    ][ilon + 1],
    areoRadius[ilat + 1][ilon], areoRadius[ilat + 1][ilon + 1]);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Interpolates MOLA topographic height data to current latitude and longitude
//!
//! This method will update the appropriate lookup indices and displacements which it will use to
//! perform 2D linear interpolation on the MOLA topographic height data.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//! \returns Mola topographic height \units{km}.
double MOLATopography::getTopographicHeight(greal lat, greal lon)
{
  if (getTopographicHeightCallback) {
    return getTopographicHeightCallback(lat, lon, callbackDataPointer);
  }

  // MIBIndices must be updated prior to interpolation.
  updateIndices(lat, lon, topoParams);

  // use 2D linear interpolation to get (and return) the topo height at current position
  Interpolator interp2d(topoParams.latDisp, topoParams.lonDisp);
  size_t ilat = topoParams.latIndex;
  size_t ilon = topoParams.lonIndex;
  return interp2d.linear(
    topoHeight[ilat    ][ilon], topoHeight[ilat    ][ilon + 1],
    topoHeight[ilat + 1][ilon], topoHeight[ilat + 1][ilon + 1]);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Interpolates albedo data to current latitude and longitude
//!
//! This method will update the appropriate lookup indices and displacements which it will 
//! then use to perform 2D linear interpolation on the albedo data.
//!
//! \param lat    planetocentric latitude in degrees
//! \param lon    planetocentric east positive longitude in degrees
//! \returns albedo
double MOLATopography::getAlbedo(greal lat, greal lon)
{
  // MIBIndices must be updated prior to interpolation.
  updateIndices(lat, lon, albedoParams);

  // use 2D linear interpolation to get (and return) albedo at current position
  Interpolator interp2d(albedoParams.latDisp, albedoParams.lonDisp);
  size_t ilat = albedoParams.latIndex;
  size_t ilon = albedoParams.lonIndex;
  return interp2d.linear(
    albedo[ilat    ][ilon], albedo[ilat    ][ilon + 1],
    albedo[ilat + 1][ilon], albedo[ilat + 1][ilon + 1]);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Reads MOLA radius, topological height, and albedo data.
//!
//! This method reads MOLA_data.bin, a binary data file, to populate the data arrays for the
//! MOLA radius, topological height, and albedo data.  The data is stored in static memory.
//! The data format consists three data blocks.  Each data block is preceeded by a size_t with
//! the number of doubles in the data block.
//!
//! \b Inputs
//! \arg #isInitialized
//!
//! \retval #areoRadius
//! \retval #topoHeight
//! \retval #albedo
void MOLATopography::initializeData()
{
  // Only perform initialization once.
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // These commented calls will convert the legacy data to the new format.
  // readLegacyData();
  // writeData();

  try {
    // Open a binary file stream.
    // Note that molaFileName is defined in MarsCommon
    ifstream binaryFile;
    binaryFile.open(dataPath + molaFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);
    }

    // Read the radius and topo height data blocks.
    size_t dataSize = LAT_SIZE * LON_SIZE;
    readBinaryDataBlock(binaryFile, dataSize, (greal*)areoRadius);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)topoHeight);

    // Read the albedo data.
    dataSize = ALB_LAT_SIZE * ALB_LON_SIZE;
    readBinaryDataBlock(binaryFile, dataSize, (greal*)albedo);

    // Close the input file.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + molaFileName;
  }
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \fn  MOLATopography::setAreoidRadiusCallback()
//! \brief Allow for override of MOLA areoid radius computations.
//!
//! This function sets an override callback for getAreiodRadius. The override must match the
//! signature `double(*TopoCallback)(double, double, void*)` where the arguments are the desired
//! latitude, longitude, and a data pointer.  The data pointer is provided to allow the developer
//! to pass data to the override function via setCallbackData(). The return value must be the 
//! value of the areoid radius at the specified latitude and longitude.
//!
//! \param callback A pointer to the callback function or nullptr to disable.
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \fn  MOLATopography::setTopographicHeightCallback()
//! \brief Allow for override of MOLA topographic height computations.
//!
//! This function sets an override callback for getTopographicHeight. The override must match the
//! signature `double(*TopoCallback)(double, double, void*)` where the arguments are the desired
//! latitude, longitude, and a data pointer.  The data pointer is provided to allow the developer
//! to pass data to the override function via setCallbackData(). The return value must be the 
//! value of the topographic height at the specified latitude and longitude.
//!
//! \param callback A pointer to the callback function or nullptr to disable.
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \fn  MOLATopography::setCallbackData()
//! \brief Set pointer to user data needed in overrides.
//!
//! This method allows a developer to pass data to the override functions set by 
//! setAreoidRadiusCallback() and setTopographicHeightCallback().
//!
//! \param dataPointer A pointer to a block of data.
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#ifdef CONVERT_MARS_DATA
//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Write the new binary data file format.
//!
//! The new file format consists of a repeated sequence of (block size, data block) pairs.
//! Output is a binary file consisting of greals.  If the definition of greal is modified,
//! then this file will need to be recreated.
void MOLATopography::writeData()
{
  ofstream binaryFile;
  binaryFile.open(dataPath + "MOLA_data.bin", ios::binary | ios::out);

  size_t dataSize = LAT_SIZE *LON_SIZE * sizeof(double);
  binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
  binaryFile.write(reinterpret_cast<char*>(&areoRadius), dataSize);

  binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
  binaryFile.write(reinterpret_cast<char*>(&topoHeight), dataSize);

  dataSize = ALB_LAT_SIZE * ALB_LON_SIZE * sizeof(double);
  binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
  binaryFile.write(reinterpret_cast<char*>(&albedo), dataSize);

  binaryFile.close();
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void MOLATopography::readLegacyData()
{
  int bytes;
  ifstream binaryFile[3];
  double input[LAT_SIZE];

  binaryFile[0].open("molatoph.bin", ios::binary | ios::in);
  binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
  for (size_t ilon = 0; ilon < LON_SIZE; ++ilon) {
    binaryFile[0].read(reinterpret_cast<char*>(input), sizeof(double) * LAT_SIZE);
    for (size_t ilat = 0; ilat < LAT_SIZE; ++ilat) {
      areoRadius[ilat][ilon] = input[ilat];
    }
  }
  for (size_t ilon = 0; ilon < LON_SIZE; ++ilon) {
    binaryFile[0].read(reinterpret_cast<char*>(input), sizeof(double) * LAT_SIZE);
    for (size_t ilat = 0; ilat < LAT_SIZE; ++ilat) {
      topoHeight[ilat][ilon] = input[ilat];
    }
  }
  binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
  binaryFile[0].close();

  binaryFile[1].open("albedo1.bin", ios::binary | ios::in);
  binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
  for (size_t ilon = 0; ilon < ALB_LON_SIZE; ++ilon) {
    binaryFile[1].read(reinterpret_cast<char*>(input), sizeof(double) * ALB_LAT_SIZE);
    for (size_t ilat = 0; ilat < ALB_LAT_SIZE; ++ilat) {
      albedo[ilat][ilon] = input[ilat];
    }
  }
  binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
  binaryFile[1].close();
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#endif

} // namespace
