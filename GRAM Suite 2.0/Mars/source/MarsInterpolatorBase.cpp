//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "MarsInterpolatorBase.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MarsInterpolatorBase::MarsInterpolatorBase()
  : Atmosphere(), MarsCommon(this)
{
  atmos.setPlanetSpecificMetrics(marsAtmos);
}

//! \fn  MarsInterpolatorBase::MarsInterpolatorBase(const MarsInterpolatorBase& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MarsInterpolatorBase::~MarsInterpolatorBase()
//! \copydoc Atmosphere::~Atmosphere()

//!  \brief Computes net values based on a set of tidal parameters.
//!
//! Computes net values for a given local true solar time based on a set of tidal parameters
//!
//! \param tp   Tidal parameters.
//! \param ltst Solar time.
//! \param form The form of the equation.
//!
//! \returns The net value.
greal MarsInterpolatorBase::getTideValue(const MarsTideParameters& tp, greal ltst, TideType form)
{
  greal tideValue = 0.0;  // The return value to be computed.
  greal cos1D = cos(PI * (ltst - tp.phase1D) / 12.0); // cosine term for once/day tides
  greal cos2D = cos(PI * (ltst - tp.phase2D) / 6.0);  // cosine term for twice/day tides

  switch (form) {

  // form=UVT: for MGCM/TesGCM u, v, and t
  case UVT:
    tideValue = tp.diurnalMean + tp.amplitude1D * cos1D + tp.amplitude2D * cos2D;
    break;

  // form=DP: for MGCM/TesGCM d, p (TesGCMSurface and TesGCMLower)
  case DP:
    tideValue = max(0.1 * tp.diurnalMean,
                    tp.diurnalMean * (1.0 + (tp.amplitude1D * cos1D + tp.amplitude2D * cos2D) / 100.0));
    break;

  // form=BIG_DP: for TesGCMUpperInterpolator data where amplitude1D, amplitude2D are assumed large
  case BIG_DP:
    tideValue = tp.diurnalMean * pow(1.0 + 0.01 * tp.amplitude1D, cos1D) 
                               * pow(1.0 + 0.01 * tp.amplitude2D, cos2D);
    break;

  default:
    break;
  }

  return tideValue;
}

//! \brief  Computes tide extrema over one sol.
//!
//! This method computes the minimum and maximum values over the course of one sol given
//! a set of tidal parameters.  The sampling rate is once per local hour.
//!
//! \param tp   Tidal parameters.
//! \param form The form of the equation.
//! \param[out] minValue A greal.
//! \param[out] maxValue A greal.
//! \retval minValue
//! \retval maxValue
void MarsInterpolatorBase::getSolMinMax(const MarsTideParameters& tp, TideType form, greal& minValue, greal& maxValue)
{
  // Initialize extrema
  maxValue = -9999.0;
  minValue = 9999.0;

  // loop over 24 hours in a sol
  for (int ihour = 0; ihour < 24; ihour++) {

    // get current value at loop hour
    greal hour = ihour;
    greal value = getTideValue(tp, hour, form);

    // update extrema
    maxValue = max(value, maxValue);
    minValue = min(value, minValue);
  }
}

//! \fn  MarsInterpolatorBase::updateIndices(size_t heightIndexOffset)
//! \brief Updates the base indices and displacements.
//!
//! This pure virtual function must be implemented by each derived class. The implementation
//! should populate the baseIndex and displacements members in preparation for interpolation.
//! The heightIndexOffset affects the base height index only for the computations of density
//! and pressure scale heights.
//!
//! \param heightIndexOffset Value added to the height index (0 or 1).
//!
//! \retval #baseIndex
//! \retval #displacements

//! \fn  MarsInterpolatorBase::getUpperFairingHeight()
//! \brief Get the upper height where fairing ends
//!
//! Fairing (or smoothing) between model layers occurs between the lower and upper fairing heights.
//!
//! \returns The upper height \units{km} where fairing ends.

//! \fn  MarsInterpolatorBase::getLowerFairingHeight()
//! \brief Get the lower height where fairing begins.
//!
//! Fairing (or smoothing) between model layers occurs between the lower and upper fairing heights.
//!
//! \returns The lower height \units{km} where fairing begins.

//! \fn  MarsInterpolatorBase::getLowestBoundaryLayer()
//! \brief Get the lowest boundary layer level.
//!
//! This interface method is implemented in the surface layer models.  It returns the lowest
//! boundary layer level (above 0).
//!
//! \returns The lowest boundary layer level \units{km}.

//! \fn  MarsInterpolatorBase::getThermosphereBaseHeight()
//! \brief Get the thermosphere base height.
//!
//! This interface method is implemented in the upper layer models.  It returns the
//! base height of the thermosphere.
//!
//! \returns The thermosphere base height \units{km}.

//! \fn  MarsInterpolatorBase::setHeightOffset(greal offset)
//! \brief Set the height offset.
//!
//! This interface method is implemented in the upper layer models.  The height is 
//! decremented by this offset before indices are updated.
//!
//! param offset The height offset \units{km}.

//! \fn  MarsInterpolatorBase::getBaseIndex()
//! \brief Get the current base index structure.
//!
//! This method is primarily used by unit tests. It was made public because testing was
//! too difficult with a protected method.
//!
//! \returns #baseIndex

//! \fn  MarsInterpolatorBase::getDisplacements()
//! \brief Get the current displacements structure.
//!
//! This method is primarily used by unit tests.  It was made public because testing was
//! too difficult with a protected method.
//!
//! \returns #displacements

//! \fn  MarsInterpolatorBase::getPressureScaleHeightDaily()
//! \brief Get the mean daily pressure scale height.
//!
//! \returns #pressureScaleHeightDaily

//! \brief Utility method for reading a block of binary data.
//!
//! This method reads a block of binary data from the specified stream. The number
//! of greals in the block is specified by \p dataSize.  This value is verified against
//! the size indicator preceeding the data block within the binary stream.  The data
//! is stored in the \p dataBlock parameter which must be pre-allocated and large enough
//! to hold the entire data block.
//!
//! \param binaryStream The binary data stream to be read.
//! \param dataSize The number of greals to be read.
//! \param[out] dataBlock A pre-allocated array of greals of size dataSize.
//!
//! \retval dataBlock Populated with data from the binary stream.
void readBinaryDataBlock(ifstream& binaryStream, size_t dataSize, greal* dataBlock)
{
  // Convert dataSize to bytes.
  dataSize *= sizeof(greal);

  // Get the size of the data block from the binary file.
  size_t verifySize;
  binaryStream.read(reinterpret_cast<char*>(&verifySize), sizeof(size_t));

  // Verify the data block size.
  if (!binaryStream || dataSize != verifySize) {
    throw string("Error: File read error.\n"
                 "       Data block size error in binary data format.\n"
                 "       This is an unrecoverable error.\n"
                 "       File: ");  // catch and append the file name
  }

  // Read the data block.
  binaryStream.read(reinterpret_cast<char*>(dataBlock), dataSize);

  // Check the stream for errors.
  if (!binaryStream) {
    throw string("Error: File read error.\n"
                 "       An unknown error occurred while reading a binary file.\n"
                 "       This is an unrecoverable error.\n"
                 "       File: ");  // catch and append the file name
  }
}

//! \brief Utility method for writing a block of binary data.
//!
//! blah
//!
//! \param binaryStream The binary data stream for output.
//! \param dataSize The number of greals to be written.
//! \param dataBlock An array of greals of size dataSize.
void writeBinaryDataBlock(ofstream& binaryStream, size_t dataSize, greal* dataBlock)
{
  // Write the block size
  binaryStream.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));

  // Write the data block.
  binaryStream.write(reinterpret_cast<char*>(dataBlock), dataSize * sizeof(greal));
}

} // namespace