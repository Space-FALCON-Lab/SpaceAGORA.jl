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
#include "TesLowerInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool TesLowerInterpolator::isInitialized = false;
greal TesLowerInterpolator::mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesLowerInterpolator::mgcmDensity[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesLowerInterpolator::mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesLowerInterpolator::mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesLowerInterpolator::mgcmPressure[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};

//! \copydoc Atmosphere::Atmosphere()
TesLowerInterpolator::TesLowerInterpolator()
: TesInterpolator()
{
  initializeData();
}

//! \fn  TesLowerInterpolator::TesLowerInterpolator(const TesLowerInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TesLowerInterpolator::~TesLowerInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Updates the base indices and displacements.
//!
//! Populates the baseIndex and displacements members in preparation for interpolation.
//! The heightIndexOffset affects the base height index only for the computations of density
//! and pressure scale heights.
//!
//! \param heightIndexOffset Value added to the height index (0 or 1).
//!
//! \retval #baseIndex
//! \retval #displacements
void TesLowerInterpolator::updateIndices(size_t heightIndexOffset)
{
  // Height indices start at -5 and increase by 1 for the first 10 layers, then by 5.

  // Find the height index of the first MGCM level above topographic surface height + 1 km.                                  
  // Heights might be negative here, so use int instead of size_t
  int idx = int(6.0 + floor(surfaceHeight + 0.3));
  size_t surfaceHeightIndex = size_t(max(0, idx));
  greal surfaceLayerHeight = greal(surfaceHeightIndex) - 5.0;
  if (surfaceHeightIndex > 15) {
    surfaceLayerHeight = 10.0 + 5.0 * (greal(surfaceHeightIndex) - 15.0);
  }

  // Compute the height index if it is above the surface layer
  if (height >= surfaceLayerHeight) {
    idx = int(floor(height + 5.0));
    if (idx > 15) {
     idx = 15 + int(floor((height - 10.0) / 5.0));
    }
    baseIndex.hgt = size_t(clampSize(idx, int(MGCM_HEIGHT_SIZE - 1)));
  }
  // Below the surface layer, just use the surface height index.
  else {
    baseIndex.hgt = surfaceHeightIndex;
  }

  // Apply the height index offset
  index.hgt = baseIndex.hgt + heightIndexOffset;

  // Update lat and ls indices.
  updateBaseIndices(MGCM_LAT_SIZE, 0.0);
}

//! \brief Gets tide parameters corresponding to the current #index.
//!
//! This method looks up tidal parameters in the supplied MGCM \p data array.  The current
//! #index is used for lookup within the array. The parameter \p latIndex is used as the
//! latitide index (typically index.lat or index.wlat).  The \p poleFactor is a dampening
//! factor that is applied to the amplitudes to smooth data near the poles.
//!
//! \param data Multi-dimensional array of tide parameters.
//! \param latIndex Latitude index is either index.lat or index.wlat.
//! \param poleFactor An amplitude dampening factor (0 to 1).
//!
//! \b Inputs
//! \arg #index
//!
//! \returns Tidal parameters.
MarsTideParameters TesLowerInterpolator::getTideParameters(const greal data[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE], size_t latIndex, greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][latIndex][index.ls][index.tyr];
  tp.amplitude1D = data[1][index.hgt][latIndex][index.ls][index.tyr] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][latIndex][index.ls][index.tyr] * poleFactor;
  tp.phase1D = data[2][index.hgt][latIndex][index.ls][index.tyr];
  tp.phase2D = data[4][index.hgt][latIndex][index.ls][index.tyr];
  return tp;
}

//! \brief Computes TPD by interpolating MGCM data in 2-D (latitude, Ls).
//!
//! Interpolates MGCM data in 2-D (latitude, Ls) for a given TES year, MGCM height index, and
//! local solar time.
//!
//! \b Inputs
//! \arg #baseIndex
//! \arg #solarTime
//!
//! \retval #temperature
//! \retval #pressure
//! \retval #density
//! \retval #ewWind
//! \retval #nsWind
//! \retval #specificGasConstant
//! \retval #temperatureDaily
//! \retval #pressureDaily
//! \retval #densityDaily
//! \retval #ewWindDaily
//! \retval #nsWindDaily
//! \retval #temperatureMin
//! \retval #temperatureMax
//! \retval #height
void TesLowerInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 2-D interpolation
  //-------------------------------------------------------------------------------------------------//
  greal tday[2][2] = {};  // 2-D 'box' array for interpolating daily mean T
  greal tmax[2][2] = {};  // 2-D 'box' array for interpolating daily maximum T
  greal tmin[2][2] = {};  // 2-D 'box' array for interpolating daily minimum T
  greal dmax[2][2] = {};  // 2-D 'box' array for interpolating daily maximum D
  greal dmin[2][2] = {};  // 2-D 'box' array for interpolating daily minimum D
  greal dday[2][2] = {};  // 2-D 'box' array for interpolating daily mean P
  greal uday[2][2] = {};  // 2-D 'box' array for interpolating daily mean U
  greal vday[2][2] = {};  // 2-D 'box' array for interpolating daily mean V
  greal r0[2][2] = {};    // 2-D 'box' array for interpolating current gas constant
  greal dm[2][2] = {};    // 2-D 'box' array for interpolating current P
  greal tm[2][2] = {};    // 2-D 'box' array for interpolating current T
  greal um[2][2] = {};    // 2-D 'box' array for interpolating current U
  greal vm[2][2] = {};    // 2-D 'box' array for interpolating current V

  //=================================================================================================//
  // Populate the interpolation cubes.  Loop over lat and ls.
  //-------------------------------------------------------------------------------------------------//
  // no interpolation over the TES year
  index.tyr = baseIndex.tyr;
  // Start LAT loop
  for (int ilat = 0; ilat < 2; ilat++) {
    // Pole factors are applied to the tidal amplitudes. 
    // For TPD: Away from the poles, the factor is 1.
    greal poleFactor = 1.0;
    // At the poles, set the factor to 0.
    // That's south pole at the bottom of the layer
    // Or the north pole at the top of the layer
    if ((baseIndex.lat == 0 && ilat == 0) || (baseIndex.lat == MGCM_LAT_SIZE - 2 && ilat == 1)) {
      poleFactor = 0.0;
    }

    // A different pole factor is used for winds.
    greal windPoleFactor = 1.0;
    if ((baseIndex.wlat == 0 && ilat == 0) || (baseIndex.wlat == MGCM_LAT_SIZE - 2 && ilat == 1)) {
      // This pole factor is computed in MGCMInterpolator::updateBaseIndices()
      windPoleFactor = displacements.wpolefac;
    }

    // Start LS loop
    for (int ils2 = 0; ils2 < 2; ils2++) {
      // Update index values (base + counter).  Counter = 0 means the bottom of the 
      // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
      index.ls = baseIndex.ls + ils2;
      index.lat = baseIndex.lat + ilat;
      index.wlat = baseIndex.wlat + ilat;

      // Cubes are populated by getting the tide parameters for the current index.
      // Then use the tide parameters to get the net value for the current solar time.
      // Temperature
      MarsTideParameters t = getTideParameters(mgcmTemperature, index.lat, poleFactor);
      tm[ilat][ils2] = getTideValue(t, solarTime, UVT);

      // Density
      MarsTideParameters d = getTideParameters(mgcmDensity, index.lat, poleFactor);
      dm[ilat][ils2] = getTideValue(d, solarTime, DP);

      // Diurnal mean for Pressure
      double pressureDiurnalMean = mgcmPressure[index.hgt][index.lat][index.ls][index.tyr];

      // Use the ideal gas law to compute R.
      if (abs(t.diurnalMean) > 0.0 && abs(d.diurnalMean) > 0.0) {
        r0[ilat][ils2] = pressureDiurnalMean / (t.diurnalMean * d.diurnalMean);
      }
      else {
        r0[ilat][ils2] = 190.0;
      }

      // EW Wind (note that we are using LAT and the wind pole factor)
      MarsTideParameters u = getTideParameters(mgcmEWWind, index.lat, windPoleFactor);
      um[ilat][ils2] = getTideValue(u, solarTime, UVT);

      // NS Wind (note that we are using WLAT and the wind pole factor)
      MarsTideParameters v = getTideParameters(mgcmNSWind, index.wlat, windPoleFactor);
      vm[ilat][ils2] = getTideValue(v, solarTime, UVT);

      // Save the daily means for Temperature, Density, and winds.
      tday[ilat][ils2] = t.diurnalMean;
      dday[ilat][ils2] = d.diurnalMean;
      uday[ilat][ils2] = u.diurnalMean;
      vday[ilat][ils2] = v.diurnalMean;

      // If directed, compute daily extrema for Temperature and Density
      if (computeMinMax) {
        // Initialize to extreme nonsense.
        tmax[ilat][ils2] = -9999.0;
        tmin[ilat][ils2] = 9999.0;
        dmax[ilat][ils2] = -9999.0;
        dmin[ilat][ils2] = 9999.0;

        // Sample once per local hour over one sol.
        for (int ihour = 0; ihour < 24; ihour++) {
          // Use current tide parameters and sample time to get
          // Temperature and Pressure
          greal hourlyTemp = getTideValue(t, greal(ihour), UVT);
          greal hourlyDensity = getTideValue(d, greal(ihour), DP);

          // Save extrema.
          tmax[ilat][ils2] = max(hourlyTemp, tmax[ilat][ils2]);
          tmin[ilat][ils2] = min(hourlyTemp, tmin[ilat][ils2]);
          dmax[ilat][ils2] = max(hourlyDensity, dmax[ilat][ils2]);
          dmin[ilat][ils2] = min(hourlyDensity, dmin[ilat][ils2]);
        }
      }
    }  // end ls loop
  }  // end lat loop

  //=================================================================================================//
  //  use 2-D interpolation (in LAT and LS) to compute final values                                  //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp(displacements.lat, displacements.ls);

  // Temperature (current and daily average)
  temperature = interp.linear(tm);
  temperatureDaily = interp.linear(tday);

  // Density (current and daily average)
  //density = interp.linear(dm);
  //densityDaily = interp.linear(dday);
  density = interp.log(dm);
  densityDaily = interp.log(dday);

  // Gas Constant (mean)
  specificGasConstant = interp.linear(r0);

  // Pressure (current and daily average from gas law)
  pressure = density * specificGasConstant * temperature;
  pressureDaily = densityDaily * specificGasConstant * temperatureDaily;

  // If directed, Temperature and Density extrema
  if (computeMinMax) {
    temperatureMax = interp.linear(tmax);
    temperatureMin = interp.linear(tmin);
    //densityMax = interp.linear(dmax);
    //densityMin = interp.linear(dmin);
    densityMax = interp.log(dmax);
    densityMin = interp.log(dmin);
  }

  // EW Wind (current and daily average)
  ewWind = interp.linear(um);
  ewWindDaily = interp.linear(uday);

  // Create an interpolator for NS winds using wlat.
  Interpolator interpW(displacements.wlat, displacements.ls);

  // NS Wind (current and daily average)
  nsWind = interpW.linear(vm);
  nsWindDaily = interpW.linear(vday);

  // Set height to the value of the bottom of the layer.
  if (index.hgt > 15) {
		height = 10.0 + 5.0 * (index.hgt - 15.0);
	}
	else {
		height = index.hgt - 5.0;
	}
}

//! \brief Reads TES lower atmosphere tidal parameter data.
//!
//! This method reads TES_lower_data.bin, a binary data file, to populate the data arrays for the
//! MGCM temperature, pressure, density, and winds data.  The data is stored in static memory.
//! The data format consists five data blocks.  Each data block is preceeded by a size_t with
//! the number of doubles in the data block.
//!
//! \b Inputs
//! \arg #isInitialized
//!
//! \retval #mgcmDensity
//! \retval #mgcmTemperature
//! \retval #mgcmPressure
//! \retval #mgcmEWWind
//! \retval #mgcmNSWind
void TesLowerInterpolator::initializeData() {
  // Only initialize the date once
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // These lines facilitate conversion of file formats.
  //readLegacyLowerData();
  //writeLowerData();

  try {
    // Open a binary file stream.
    ifstream binaryFile;
    binaryFile.open(dataPath + tesLowerFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the Pressure data block (means only)
    size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE;
    // Read the Pressure data
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmPressure);

    // Compute the size of the remaining data blocks (augment with param dimension).
    dataSize *= PARAM_SIZE;
    // Read the remaining data blocks.
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmTemperature);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmDensity);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmEWWind);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmNSWind);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + tesLowerFileName;
  }
}

#ifdef CONVERT_MARS_DATA
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void TesLowerInterpolator::readLegacyLowerData()
{
  int bytes, ils, ihgt;
  double ylat;


  MIBIndices idx;
  for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
    string fileName = "";
    if (idx.tyr == 0) fileName = "tpdloy11.bin";
    if (idx.tyr == 1) fileName = "tpdloy21.bin";
    ifstream binaryFile[3];
    binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);

    if (idx.tyr == 0) fileName = "uvloy11.bin";
    if (idx.tyr == 1) fileName = "uvloy21.bin";
    binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);

    for (idx.ls = 1; idx.ls < LS_SIZE; ++idx.ls) {
      for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
        for (size_t ix = 0; ix < MGCM_HEIGHT_SIZE; ++ix) {
          idx.hgt = MGCM_HEIGHT_SIZE - ix - 1;
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ils), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ylat), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[0][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[1][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[2][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[3][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[4][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[0][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[1][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[2][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[3][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[4][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));


          binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ils), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ylat), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[0][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[1][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[2][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[3][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[4][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[0][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[1][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[2][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[3][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[4][idx.hgt][idx.lat][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
        }
      }
    }
    binaryFile[0].close();
    binaryFile[1].close();
  }

  for (idx.hgt = 0; idx.hgt < MGCM_HEIGHT_SIZE; ++idx.hgt) {
    for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
      for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          mgcmTemperature[p][idx.hgt][idx.lat][0][idx.tyr] = mgcmTemperature[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr];
          mgcmDensity[p][idx.hgt][idx.lat][0][idx.tyr] = mgcmDensity[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr];
          mgcmEWWind[p][idx.hgt][idx.lat][0][idx.tyr] = mgcmEWWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr];
          mgcmNSWind[p][idx.hgt][idx.lat][0][idx.tyr] = mgcmNSWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr];
        }
        mgcmPressure[idx.hgt][idx.lat][0][idx.tyr] = mgcmPressure[idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr];
      }
    }
  }
}

//! \brief Write the new binary data file format.
//!
//! The new file format consists of a repeated sequence of (block size, data block) pairs.
//! Output is a binary file consisting of greals.  If the definition of greal is modified,
//! then this file will need to be recreated.
void TesLowerInterpolator::writeLowerData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "TES_lower_data.bin", ios::binary | ios::out);

  size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmPressure);

  dataSize *= PARAM_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmTemperature);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmDensity);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmEWWind);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmNSWind);

  binaryFile.close();
}
#endif

} // namespace