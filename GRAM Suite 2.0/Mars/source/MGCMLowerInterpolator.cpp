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
#include "MGCMLowerInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool MGCMLowerInterpolator::isInitialized = false;
greal MGCMLowerInterpolator::mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMLowerInterpolator::mgcmPressure[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMLowerInterpolator::mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMLowerInterpolator::mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMLowerInterpolator::mgcmDensity[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][3] = {};

//! \copydoc Atmosphere::Atmosphere()
MGCMLowerInterpolator::MGCMLowerInterpolator()
: MGCMInterpolator()
{
  initializeData();
}

//! \fn  MGCMLowerInterpolator::MGCMLowerInterpolator(const MGCMLowerInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MGCMLowerInterpolator::~MGCMLowerInterpolator()
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
void MGCMLowerInterpolator::updateIndices(size_t heightIndexOffset)
{
  constexpr greal heightStepSize = 5.0;

  // Find the height index of the first MGCM level above topographic surface height + 1 km.                                  
  // Heights might be negative here, so use int instead of size_t
  int idx = int(1.0 + floor((surfaceHeight + 1.0) / heightStepSize));
  size_t surfaceHeightIndex = size_t(max(0, idx));
  greal surfaceLayerHeight = heightStepSize * greal(surfaceHeightIndex);

  // Compute the height index if it is above the surface layer
  if (height >= surfaceLayerHeight) {
    idx = int(floor(height / heightStepSize));
    baseIndex.hgt = size_t(clampSize(idx, int(MGCM_HEIGHT_SIZE - 1)));
  }
  // Below the surface layer, just use the surface height index.
  else {
    baseIndex.hgt = surfaceHeightIndex;
  }
  // Apply the height index offset
  index.hgt = baseIndex.hgt + heightIndexOffset;

  // Update lat, ls, and od indices.
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
MarsTideParameters MGCMLowerInterpolator::getTideParameters(const greal data[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE], size_t latIndex, greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][latIndex][index.ls][index.od];
  tp.amplitude1D = data[1][index.hgt][latIndex][index.ls][index.od] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][latIndex][index.ls][index.od] * poleFactor;
  tp.phase1D = data[2][index.hgt][latIndex][index.ls][index.od];
  tp.phase2D = data[4][index.hgt][latIndex][index.ls][index.od];
  return tp;
}

//! \brief Computes TPD by interpolating MGCM data in 3-D (latitude, Ls, and dust OD).
//!
//! Interpolates MGCM data in 3-D (latitude, Ls, and dust OD) for a given MGCM height index and
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
//! \retval #densityMin
//! \retval #densityMax
//! \retval #height
void MGCMLowerInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 3-D interpolation
  //-------------------------------------------------------------------------------------------------//
  greal tday[2][2][2] = {};  // 3-D 'box' array for interpolating daily mean T
  greal tmax[2][2][2] = {};  // 3-D 'box' array for interpolating daily maximum T
  greal tmin[2][2][2] = {};  // 3-D 'box' array for interpolating daily minimum T
  greal dmax[2][2][2] = {};  // 3-D 'box' array for interpolating daily maximum D
  greal dmin[2][2][2] = {};  // 3-D 'box' array for interpolating daily minimum D
  greal pday[2][2][2] = {};  // 3-D 'box' array for interpolating daily mean P
  greal uday[2][2][2] = {};  // 3-D 'box' array for interpolating daily mean U
  greal vday[2][2][2] = {};  // 3-D 'box' array for interpolating daily mean V
  greal r0[2][2][2] = {};    // 3-D 'box' array for interpolating current gas constant
  greal pm[2][2][2] = {};    // 3-D 'box' array for interpolating current P
  greal tm[2][2][2] = {};    // 3-D 'box' array for interpolating current T
  greal um[2][2][2] = {};    // 3-D 'box' array for interpolating current U
  greal vm[2][2][2] = {};    // 3-D 'box' array for interpolating current V

  //=================================================================================================//
  // Populate the interpolation cubes.  Loop over lat, ls, and od.
  //-------------------------------------------------------------------------------------------------//
  // Start LAT loop
  for (size_t ilat = 0; ilat < 2; ilat++) {
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
    if ((baseIndex.wlat == 0 && ilat == 0) || (baseIndex.wlat == MGCM_LAT_SIZE - 2 && ilat == 1)){
      // This pole factor is computed in MGCMInterpolator::updateBaseIndices()
      windPoleFactor = displacements.wpolefac;
    }
    
    // Start LS loop
    for (size_t ils = 0; ils < 2; ils++) {
      // Start OD loop
      for (size_t iod = 0; iod < 2; iod++) {
        // Update index values (base + counter).  Counter = 0 means the bottom of the 
        // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
        index.ls = baseIndex.ls + ils;
        index.od = baseIndex.od + iod;
        index.lat = baseIndex.lat + ilat;
        index.wlat = baseIndex.wlat + ilat;

        // Cubes are populated by getting the tide parameters for the current index.
        // Then use the tide parameters to get the net value for the current solar time.
        // Temperature
        MarsTideParameters t = getTideParameters(mgcmTemperature, index.lat, poleFactor);
        tm[ilat][ils][iod] = getTideValue(t, solarTime, UVT);

        // Pressure
        MarsTideParameters p = getTideParameters(mgcmPressure, index.lat, poleFactor);
        pm[ilat][ils][iod] = getTideValue(p, solarTime, DP);

        // For Density, get the mean value.
        double densityDiurnalMean = mgcmDensity[index.hgt][index.lat][index.ls][index.od];

        // Use the ideal gas law to compute R.
        if (abs(t.diurnalMean) > 0.0 && abs(densityDiurnalMean) > 0.0) {
          r0[ilat][ils][iod] = p.diurnalMean / (t.diurnalMean * densityDiurnalMean);
        }
        else {
          r0[ilat][ils][iod] = 190.0;
        }

        // EW Wind (note that we are using wlat and the wind pole factor)
        MarsTideParameters u = getTideParameters(mgcmEWWind, index.wlat, windPoleFactor);
        um[ilat][ils][iod] = getTideValue(u, solarTime, UVT);

        // NS Wind (note that we are using wlat and the wind pole factor)
        MarsTideParameters v = getTideParameters(mgcmNSWind, index.wlat, windPoleFactor);
        vm[ilat][ils][iod] = getTideValue(v, solarTime, UVT);

        // Save the daily means for Temperature, Pressure, and winds.
        tday[ilat][ils][iod] = t.diurnalMean;
        pday[ilat][ils][iod] = p.diurnalMean;
        uday[ilat][ils][iod] = u.diurnalMean;
        vday[ilat][ils][iod] = v.diurnalMean;

        // If directed, compute daily extrema for Temperature and Density
        if (computeMinMax) {
          // Initialize to extreme nonsense.
          tmax[ilat][ils][iod] = -9999.0;
          tmin[ilat][ils][iod] = 9999.0;
          dmax[ilat][ils][iod] = -9999.0;
          dmin[ilat][ils][iod] = 9999.0;

          // Sample once per local hour over one sol.
          for (int ihour = 0; ihour < 24; ihour++) {
            // Use current tide parameters and sample time to get
            // Temperature and Pressure
            greal hourlyTemp = getTideValue(t, greal(ihour), UVT);
            greal hourlyPres = getTideValue(p, greal(ihour), DP);

            // Compute Density using ideal gas law.
            greal hourlyDens = hourlyPres / (hourlyTemp * r0[ilat][ils][iod]);

            // Save extrema.
            tmax[ilat][ils][iod] = max(hourlyTemp, tmax[ilat][ils][iod]);
            tmin[ilat][ils][iod] = min(hourlyTemp, tmin[ilat][ils][iod]);
            dmax[ilat][ils][iod] = max(hourlyDens, dmax[ilat][ils][iod]);
            dmin[ilat][ils][iod] = min(hourlyDens, dmin[ilat][ils][iod]);
          }
        }
      }  // end od loop
    }  // end ls loop
  }  // end lat loop

  //=================================================================================================//
  //  use 3-D interpolation (in LAT, LS, and OD) to compute final values                             //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp(displacements.lat, displacements.ls, displacements.od);

  // Temperature (current and daily average)
  temperature = interp.linear(tm);
  temperatureDaily = interp.linear(tday);

  // Pressure (current and daily average)
  pressure = interp.log(pm);
  pressureDaily = interp.log(pday);

  // Gas Constant (mean)
  specificGasConstant = interp.linear(r0);

  // Density (current and daily average from gas law)
  density = pressure / (specificGasConstant * temperature);
  densityDaily = pressureDaily / (specificGasConstant * temperatureDaily);

  // If directed, Temperature and Density extrema
  if (computeMinMax) {
    temperatureMax = interp.linear(tmax);
    temperatureMin = interp.linear(tmin);
    densityMax = interp.log(dmax);
    densityMin = interp.log(dmin);
  }

  // Create an interpolator for winds using wlat.
  Interpolator interpW(displacements.wlat, displacements.ls, displacements.od);

  // EW Wind (current and daily average)
  ewWind = interpW.linear(um);
  ewWindDaily = interpW.linear(uday);

  // NS Wind (current and daily average)
  nsWind = interpW.linear(vm);
  nsWindDaily = interpW.linear(vday);

  // Set height to the value of the bottom of the layer.
  height = 5.0 * index.hgt;
}

//! \brief Reads MGCM lower atmosphere tidal parameter data.
//!
//! This method reads MGCM_lower_data.bin, a binary data file, to populate the data arrays for the
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
void MGCMLowerInterpolator::initializeData() {
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
    binaryFile.open(dataPath + mgcmLowerFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * OD_SIZE;
    // Read the Density data
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmDensity);

    // Compute the size of the remaining data blocks (augment with param dimension).
    dataSize *= PARAM_SIZE;
    // Read the remaining data blocks.
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmTemperature);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmPressure);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmEWWind);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmNSWind);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + mgcmLowerFileName;
  }
}


#ifdef CONVERT_MARS_DATA
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void MGCMLowerInterpolator::readLegacyLowerData()
{
  static const greal DFA[8][OD_SIZE] = {
    { 1.03050e+00,   1.01386e+00,  1.2029e+0},
    {-6.81901e-04,  -1.52668e-04,  0.0000e+0},
    { 6.67896e-06,   5.78641e-06,  4.5266e-5},
    { 1.69052e-02,   5.57000e-03,  0.0000e+0},
    { 2.83632e-04,  -8.23194e-05,  0.0000e+0},
    { 4.11891e-06,  -7.11778e-06,  0.0000e+0},
    {-1.05497e-01,   0.00000e+0,   0.0000e+0},
    { 9.17778e-05,   0.00000e+0,   0.0000e+0}
  };
  static const greal DFB[8][OD_SIZE] = {
    {-3.91632e+00,  -1.63780e+00,  18.1260e+0},
    { 1.02837e-01,   6.76088e-02,   0.0000e+0},
    {-2.73027e-04,  -3.23646e-05,  -9.5731e-3},
    { 1.02129e+00,   8.20808e-02,   0.0000e+0},
    {-3.11446e-02,  -4.21906e-02,   0.0000e+0},
    {-5.95152e-04,   1.97186e-04,   0.0000e+0},
    { 8.26562e+00,   0.00000e+0,    0.0000e+0},
    {-3.37582e-02,   0.00000e+0,    0.0000e+0}
  };
  static const greal DFC[3][OD_SIZE] = {
    {-0.1307e0,  -0.1657e0,  -0.3931e0},
    { 1.6074e0,   1.6981e0,   2.9640e0},
    {-3.7984e0,  -4.1688e0,  -5.6292e0}
  };

  int bytes, ils, ihgt;
  double ylat;


  MIBIndices idx;
  for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
    string fileName = "";
    if (idx.od == 0) fileName = "tpdlo031.bin";
    if (idx.od == 1) fileName = "tpdlo101.bin";
    if (idx.od == 2) fileName = "tpdlo301.bin";
    ifstream binaryFile[3];
    binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);

    if (idx.od == 0) fileName = "uvlo031.bin";
    if (idx.od == 1) fileName = "uvlo101.bin";
    if (idx.od == 2) fileName = "uvlo301.bin";
    binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);

    for (idx.ls = 1; idx.ls < LS_SIZE; ++idx.ls) {
      greal SINLS = sin(toRadians(30.0_deg * idx.ls));
      for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
        greal XLAT = 7.5 * idx.lat - 90.0_deg;
        for (size_t ix = 0; ix < MGCM_HEIGHT_SIZE; ++ix) {
          idx.hgt = MGCM_HEIGHT_SIZE - ix - 1;
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ils), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ylat), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[0][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[1][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[2][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[3][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmTemperature[4][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[1][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[2][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[3][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmPressure[4][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&mgcmDensity[idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));

          // Adjust pressure and density to agree with Map Year 2
          greal dfa = DFA[0][idx.od] + DFA[1][idx.od] * XLAT + DFA[2][idx.od] * XLAT * XLAT
            + (DFA[3][idx.od] + DFA[4][idx.od] * XLAT + DFA[5][idx.od] * XLAT * XLAT) * SINLS
            + (DFA[6][idx.od] + DFA[7][idx.od] * XLAT) * SINLS * SINLS;
          greal dfb = DFB[0][idx.od] + DFB[1][idx.od] * XLAT + DFB[2][idx.od] * XLAT * XLAT
            + (DFB[3][idx.od] + DFB[4][idx.od] * XLAT + DFB[5][idx.od] * XLAT * XLAT) * SINLS
            + (DFB[6][idx.od] + DFB[7][idx.od] * XLAT) * SINLS * SINLS;
          greal dfm = dfa + dfb * mgcmDensity[idx.hgt][idx.lat][idx.ls][idx.od];
          greal ZZ = greal(idx.hgt) / 20.0;
          greal dfc = 1.0 + DFC[0][idx.od] * ZZ + DFC[1][idx.od] * ZZ * ZZ + DFC[2][idx.od] * ZZ * ZZ * ZZ;
          dfc = max(0.5, dfc);
          mgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.od] = mgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.od] * dfm / dfc;
          mgcmDensity[idx.hgt][idx.lat][idx.ls][idx.od] = mgcmDensity[idx.hgt][idx.lat][idx.ls][idx.od] * dfm / dfc;

          binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ils), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
          binaryFile[1].read(reinterpret_cast<char*>(&ylat), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[0][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[1][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[2][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[3][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmEWWind[4][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[0][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[1][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[2][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[3][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&mgcmNSWind[4][idx.hgt][idx.lat][idx.ls][idx.od]), sizeof(double));
          binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
        }
      }
    }
    binaryFile[0].close();
    binaryFile[1].close();
  }

  for (idx.hgt = 0; idx.hgt < MGCM_HEIGHT_SIZE; ++idx.hgt) {
    for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
      for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          mgcmTemperature[p][idx.hgt][idx.lat][0][idx.od] = mgcmTemperature[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od];
          mgcmPressure[p][idx.hgt][idx.lat][0][idx.od] = mgcmPressure[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od];
          mgcmEWWind[p][idx.hgt][idx.lat][0][idx.od] = mgcmEWWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od];
          mgcmNSWind[p][idx.hgt][idx.lat][0][idx.od] = mgcmNSWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od];
        }
        mgcmDensity[idx.hgt][idx.lat][0][idx.od] = mgcmDensity[idx.hgt][idx.lat][LS_SIZE - 1][idx.od];
      }
    }
  }
}

//! \brief Write the new binary data file format.
//!
//! The new file format consists of a repeated sequence of (block size, data block) pairs.
//! Output is a binary file consisting of greals.  If the definition of greal is modified,
//! then this file will need to be recreated.
void MGCMLowerInterpolator::writeLowerData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "MGCM_lower_data.bin", ios::binary | ios::out);

  size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * OD_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmDensity);

  dataSize *= PARAM_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmTemperature);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmPressure);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmEWWind);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mgcmNSWind);

  binaryFile.close();
}
#endif

} // namespace