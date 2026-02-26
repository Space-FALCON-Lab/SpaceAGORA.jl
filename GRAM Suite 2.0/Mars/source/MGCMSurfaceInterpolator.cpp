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
#include "MGCMSurfaceInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool MGCMSurfaceInterpolator::isInitialized = false;
greal MGCMSurfaceInterpolator::surfaceTemperatures[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMSurfaceInterpolator::surfaceEWWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMSurfaceInterpolator::surfaceNSWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
greal MGCMSurfaceInterpolator::surfaceHeights[SURFACE_HEIGHT_SIZE] = { greal(0.0),  greal(0.005),  greal(0.03) };

//! \copydoc Atmosphere::Atmosphere()
MGCMSurfaceInterpolator::MGCMSurfaceInterpolator()
: MGCMLowerInterpolator()
{
  initializeData();
}

//! \fn  MGCMSurfaceInterpolator::MGCMSurfaceInterpolator(const MGCMSurfaceInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MGCMSurfaceInterpolator::~MGCMSurfaceInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MGCMSurfaceInterpolator::getUpperFairingHeight()
//! \brief Gets the height of the upper bound of the fairing zone.
//! \returns The height \units{km} of the upper bound of the fairing zone.

//! \fn  MGCMSurfaceInterpolator::getLowerFairingHeight()
//! \brief Gets the height of the lower bound of the fairing zone.
//! \returns The height \units{km} of the lower bound of the fairing zone.

//! \fn  MGCMSurfaceInterpolator::getLowestBoundaryLayer()
//! \brief Gets the height of the upper bound of lowest boundary layer.
//! \returns The height \units{km} of the upper bound of lowest boundary layer.

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
void MGCMSurfaceInterpolator::updateIndices(size_t heightIndexOffset)
{
  // Find the height index of the first MGCM level above topographic surface height + 1 km.                                  
  // Heights might be negative here, so use int instead of size_t
  int sidx = int(1.0 + floor((surfaceHeight + 1.0) / 5.0));
  oneKilometerIndex = max(0, sidx);  // formerly k1st

  // Find lower boundary height index
  greal heightAboveSurface = height - surfaceHeight;
  if (heightAboveSurface < surfaceHeights[1]) {
    baseIndex.hgt = 0;
  }
  else if (heightAboveSurface < surfaceHeights[2]) {
    baseIndex.hgt = 1;
  }
  else {
    baseIndex.hgt = 2;
  }
  // Apply the height index offset
  index.hgt = baseIndex.hgt + heightIndexOffset;

  // Compute the longitude index and displacement for MGCM data
  // longitude step size for surface data; evaluates to 9.0
  constexpr greal longitudeStepSize = 360.0_deg / greal(SURFACE_LON_SIZE - 1);
  size_t idx = size_t(floor(position.getLongitude(WEST_POSITIVE) / longitudeStepSize));
  baseIndex.lon = clampSize(idx, SURFACE_LON_SIZE - 1);
  displacements.lon = (position.getLongitude(WEST_POSITIVE) - longitudeStepSize * greal(baseIndex.lon)) / longitudeStepSize;

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
MarsTideParameters MGCMSurfaceInterpolator::getTideParameters(const greal data[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE], size_t latIndex, greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][latIndex][index.lon][index.ls][index.od];
  tp.amplitude1D = data[1][index.hgt][latIndex][index.lon][index.ls][index.od] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][latIndex][index.lon][index.ls][index.od] * poleFactor;
  tp.phase1D = data[2][index.hgt][latIndex][index.lon][index.ls][index.od];
  tp.phase2D = data[4][index.hgt][latIndex][index.lon][index.ls][index.od];
  return tp;
}


//! \brief Computes TPD by interpolating MGCM data in 4-D (latitude, longitude, Ls, and dust OD).
//!
//! Interpolates MGCM data in 4-D (latitude, longitude, Ls, and dust OD) for a given MGCM height index and
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
//! \retval #pressureScaleHeightDaily
//! \retval #ewWindDaily
//! \retval #nsWindDaily
//! \retval #temperatureMin
//! \retval #temperatureMax
//! \retval #densityMin
//! \retval #densityMax
//! \retval #height
void MGCMSurfaceInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 4-D interpolation for surface data or 3-D data from the Lower MGCM domain.
  //-------------------------------------------------------------------------------------------------//
  greal   tm[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of current T at current height index
  greal   um[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of current U at current height index
  greal   vm[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of current V at current height index
  greal tday[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of daily mean T at current height index
  greal tmax[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of daily max T at current height index
  greal tmin[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of daily min T at current height index
  greal uday[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of daily mean U at current height index
  greal vday[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of daily mean V at current height index
  greal tsrf[2][2][2][2];  // 4-d (lat, lon, Ls, OD) array 'box' of current T at ground surface index

  greal    p1km[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of current P at the 1 km index
  greal    t1km[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of current T at the 1 km index
  greal p1kmday[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean P at the 1 km index
  greal t1kmday[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean T at the 1 km index
  greal r1kmday[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean R at the 1 km index
  greal t1kmmax[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily max T at the 1 km index
  greal t1kmmin[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily min T at the 1 km index
  greal p1kmmax[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily max P at the 1 km index
  greal p1kmmin[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily min P at the 1 km index
  greal  plower[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean P at bottom of the surface layer
  greal  pupper[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean P at top of the surface layer
  greal  dlower[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean D at bottom of the surface layer
  greal  dupper[2][2][2];  // 3-d (lat, Ls, OD) array 'box' of daily mean D at top of the surface layer

  //=================================================================================================//
  // Populate the 4-D interpolation cubes.  Loop over lat, lon, ls, and od.
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
    if ((baseIndex.wlat == 0 && ilat == 0) || (baseIndex.wlat == MGCM_LAT_SIZE - 2 && ilat == 1)) {
      windPoleFactor = displacements.wpolefac;
    }

    // Start LON loop
    for (size_t ilon = 0; ilon < 2; ilon++) {  
      // Start LS loop
      for (size_t ils = 0; ils < 2; ils++) {   
        // Start OD loop
        for (size_t iod = 0; iod < 2; iod++) { 
          // Update index values (base + counter).  Counter = 0 means the bottom of the 
          // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
          index.lat = baseIndex.lat + ilat;  
          index.wlat = baseIndex.wlat + ilat;
          index.ls = baseIndex.ls + ils;     
          index.lon = baseIndex.lon + ilon;  
          index.od = baseIndex.od + iod;     

          // Cubes are populated by getting the tide parameters for the current index.
          // Then use the tide parameters to get the net value for the current solar time.
          // Temperature
          MarsTideParameters t = getTideParameters(surfaceTemperatures, index.lat, poleFactor);
          tm[ilat][ilon][ils][iod] = getTideValue(t, solarTime, UVT);

          // EW Wind (note that we are using wlat and the wind pole factor)
          MarsTideParameters u = getTideParameters(surfaceEWWinds, index.wlat, windPoleFactor);
          um[ilat][ilon][ils][iod] = getTideValue(u, solarTime, UVT);

          // NS Wind (note that we are using wlat and the wind pole factor)
          MarsTideParameters v = getTideParameters(surfaceNSWinds, index.wlat, windPoleFactor);
          vm[ilat][ilon][ils][iod] = getTideValue(v, solarTime, UVT);

          // Mean Surface (Ground) Temperature at the current index (height index = 0)
          tsrf[ilat][ilon][ils][iod] = surfaceTemperatures[0][0][index.lat][index.lon][index.ls][index.od];

          // Save the daily means for Temperature and winds.
          tday[ilat][ilon][ils][iod] = t.diurnalMean;
          uday[ilat][ilon][ils][iod] = u.diurnalMean;
          vday[ilat][ilon][ils][iod] = v.diurnalMean;

          // If directed, compute daily extrema for Temperature
          if (computeMinMax) {
            getSolMinMax(t, UVT, tmin[ilat][ilon][ils][iod], tmax[ilat][ilon][ils][iod]);
          }
        } // end loop dust OD
      } // end loop Ls
    } // end loop lon
  } // end loop lat

  //=================================================================================================//
  //	use 4-d interpolation to get cur and day T, U, and V; max and min T; and surface T             //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp4d(displacements.lat, displacements.lon, displacements.ls, displacements.od);

  // Temperature (current and daily average)
  temperature = interp4d.linear(tm);
  temperatureDaily = interp4d.linear(tday);

  // Surface Temperature
  greal surfaceTemperature = interp4d.linear(tsrf);

  // If directed, Temperature extrema
  if (computeMinMax) {
    temperatureMax = interp4d.linear(tmax);
    temperatureMin = interp4d.linear(tmin);
  }

  // Create an interpolator for winds using wlat.
  Interpolator interp4dw(displacements.wlat, displacements.lon, displacements.ls, displacements.od);

  // EW Wind (current and daily average)
  ewWind = interp4dw.linear(um);
  ewWindDaily = interp4dw.linear(uday);

  // NS Wind (current and daily average)
  nsWind = interp4dw.linear(vm);
  nsWindDaily = interp4dw.linear(vday);

  //=================================================================================================//
  // Now prepare to use MGCM data to compute scale heights to get P and D in the surface layer
  // Populate the 3-D interpolation cubes.  Loop over lat, ls, and od.
  //-------------------------------------------------------------------------------------------------//
  // MGCM height index just below surface height (or 0 if height is below MOLA areoid)
  size_t surfaceLayerIndex = (oneKilometerIndex == 0) ? oneKilometerIndex : oneKilometerIndex - 1;

  // Save the current height index and substitute the surface height index.
  size_t saveIndex = index.hgt;
  index.hgt = oneKilometerIndex;

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

    // Start LS loop
    for (size_t ils = 0; ils < 2; ils++) {
      // Start OD loop
      for (size_t iod = 0; iod < 2; iod++) {
        // Update index values (base + counter).  Counter = 0 means the bottom of the 
        // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
        index.lat = baseIndex.lat + ilat;
        index.ls = baseIndex.ls + ils;       
        index.od = baseIndex.od + iod;       

        // Cubes are populated by getting the tide parameters for the current index.
        // Then use the tide parameters to get the net value for the current solar time.
        // Temperature
        MarsTideParameters t = getTideParameters(mgcmTemperature, index.lat, poleFactor);
        t1km[ilat][ils][iod] = getTideValue(t, solarTime, UVT);

        // Pressure
        MarsTideParameters p = getTideParameters(mgcmPressure, index.lat, poleFactor);
        p1km[ilat][ils][iod] = getTideValue(p, solarTime, DP);

        // For Density, get the mean value at the 1km level.
        greal densityDiurnalMean = mgcmDensity[oneKilometerIndex][index.lat][index.ls][index.od];

        // Save the daily means for Temperature and Pressure.
        p1kmday[ilat][ils][iod] = p.diurnalMean;
        t1kmday[ilat][ils][iod] = t.diurnalMean;

        // Use the ideal gas law to compute R.
        if (abs(t.diurnalMean) > 0.0 && abs(densityDiurnalMean) > 0.0) {
          r1kmday[ilat][ils][iod] = p.diurnalMean / (t.diurnalMean*densityDiurnalMean);
        }
        else {
          r1kmday[ilat][ils][iod] = 190.0;
        }

        // Get mean pressure and density at the Surface
        plower[ilat][ils][iod] = mgcmPressure[0][surfaceLayerIndex][index.lat][index.ls][index.od];                              // daily mean P at k1h level
        pupper[ilat][ils][iod] = mgcmPressure[0][surfaceLayerIndex + 1][index.lat][index.ls][index.od];                           // daily mean P at k1h+1 level
        dlower[ilat][ils][iod] = mgcmDensity[surfaceLayerIndex][index.lat][index.ls][index.od];                               // daily mean D at k1h level
        dupper[ilat][ils][iod] = mgcmDensity[surfaceLayerIndex + 1][index.lat][index.ls][index.od];                            // daily mean D at k1h+1 level

        // If directed, compute daily extrema for Temperature and Pressure
        if (computeMinMax) {
          getSolMinMax(p, DP, p1kmmin[ilat][ils][iod], p1kmmax[ilat][ils][iod]);
          getSolMinMax(t, UVT, t1kmmin[ilat][ils][iod], t1kmmax[ilat][ils][iod]);
        }
      } // end loop dust OD
    } // end loop Ls
  } // end loop latitude

  // Restore the current height index
  index.hgt = saveIndex;

  //=================================================================================================//
  //  use 3-D interpolation (in LAT, LS, and OD) to compute final values                             //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp3d(displacements.lat, displacements.ls, displacements.od);

  // Pressure scale height at surface latyer
  greal plow =  interp3d.log(plower);        // daily mean P at bottom of the surface layer
  greal phigh = interp3d.log(pupper);        // daily mean P at top of the surface layer
  pressureScaleHeightDaily = 5.0 / log(plow / phigh);

  // Density scale height at surface layer
  greal dlow = interp3d.log(dlower);         // daily mean D at bottom of the surface layer
  greal dhigh = interp3d.log(dupper);        // daily mean D at top of the surface layer
  greal densityScaleHeightDaily = 5.0 / log(dlow / dhigh);

  // Current T and P at the 1 km height
  greal t1k = interp3d.linear(t1km);         // current T at the 1 km height
  greal p1k = interp3d.log(p1km);            // current P at the 1 km height

  // Daily means at the 1 km height
  greal t1kday = interp3d.linear(t1kmday);        // daily mean T at the 1 km height
  greal p1kday = interp3d.log(p1kmday);           // daily mean P at the 1 km height
  specificGasConstant = interp3d.linear(r1kmday); // daily average gas constant

  // layer average daily mean T between 1 km height and the surface height index
  greal tbarsurf = (t1kday + surfaceTemperature) / 2.0;

  //=================================================================================================//
  //  compute final P and D values                                                                   //
  //-------------------------------------------------------------------------------------------------//
  // if hgtIndex is 1 (5 m level) set 5m T to current T
  // Note that this must be forced to occur in the first call to function.
  if (index.hgt == 1 && temperatureAtFiveMeters <= 0.0) {
    temperatureAtFiveMeters = temperature; 
  }

  // layer average current T between the 5m level and surface + 1km level
  greal tbar = (t1k + temperatureAtFiveMeters) / 2.0;
  // adjust scale heights for differences in current, daily mean layer average T
  pressureScaleHeight = pressureScaleHeightDaily * tbar / tbarsurf;
  densityScaleHeight = densityScaleHeightDaily * tbar / tbarsurf;

  // recompute height (wrt MOLA areoid) of hgtIndex surface level
  height = surfaceHeight + surfaceHeights[index.hgt];    
  // height (wrt MOLA areoid) of the 1 km index level
  greal z1k = 5.0 * greal(oneKilometerIndex);           
  greal deltaHeight = z1k - height;

  // extrapolate downward from 1 km to current hgtIndex to get current P
  pressure = p1k * exp(deltaHeight / pressureScaleHeight);
  // get current D at hgtIndex from gas law
  density = pressure / (specificGasConstant * temperature);     

  // extrapolate downward from 1 km to current hgtIndex to get daily mean P
  pressureDaily = p1kday * exp(deltaHeight / pressureScaleHeightDaily);
  // get daily mean D at hgtIndex from gas law
  densityDaily = pressureDaily / (specificGasConstant * temperatureDaily);                                           // get daily mean D at hgtIndex from gas law

  // If directed, Temperature and Density extrema
  if (computeMinMax) {
    // daily max/min T at the 1 km height
    greal t1kmax = interp3d.linear(t1kmmax);

    // get max/min P scale height
    greal hp1kmax = pressureScaleHeightDaily * 0.5 * (t1kmax + temperatureMin) / tbarsurf;
    greal hp1kmin = pressureScaleHeightDaily * 0.5 * (t1kmax + temperatureMax) / tbarsurf;

    // extrapolate max/min P downward from i km to hgtIndex to get daily max/min P
    greal p1kmax = interp3d.log(p1kmmax) * exp(deltaHeight / hp1kmax);
    greal p1kmin = interp3d.log(p1kmmin) * exp(deltaHeight / hp1kmin);

    // get density max/min from gas law
    densityMax = densityDaily * (p1kmax / pressureDaily) * (temperatureDaily / temperatureMin);
    densityMin = densityDaily * (p1kmin / pressureDaily) * (temperatureDaily / temperatureMax);
  }
}

//! \brief Reads MGCM surface atmosphere tidal parameter data.
//!
//! This method reads MGCM_surface_data.bin, a binary data file, to populate the data arrays for the
//! MGCM temperature and winds data.  The data is stored in static memory.
//! The data format consists three data blocks.  Each data block is preceeded by a size_t with
//! the number of doubles in the data block.
//!
//! \b Inputs
//! \arg #isInitialized
//!
//! \retval #surfaceTemperatures
//! \retval #surfaceEWWinds
//! \retval #surfaceNSWinds
void MGCMSurfaceInterpolator::initializeData() {
  // Only initialize the date once
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // These lines facilitate conversion of file formats.
  //readLegacySurfaceData();
  //writeSurfaceData();

  try {
    // Open a binary file stream.
    ifstream binaryFile;
    binaryFile.open(dataPath + mgcmSurfaceFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * OD_SIZE;

    // Read the data blocks
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceTemperatures);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceEWWinds);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceNSWinds);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + mgcmSurfaceFileName;
  }
}

#ifdef CONVERT_MARS_DATA
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void MGCMSurfaceInterpolator::readLegacySurfaceData()
{

  int bytes, ils, jlon;
  double ylat;


  MIBIndices idx;
  for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
    string fileName = "";
    if (idx.od == 0) fileName = "sfc00031.bin";
    if (idx.od == 1) fileName = "sfc00101.bin";
    if (idx.od == 2) fileName = "sfc00301.bin";
    ifstream binaryFile[3];
    binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);
    if (!binaryFile[0].is_open()) exit(1);

    if (idx.od == 0) fileName = "sfc05031.bin";
    if (idx.od == 1) fileName = "sfc05101.bin";
    if (idx.od == 2) fileName = "sfc05301.bin";
    binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);
    if (!binaryFile[1].is_open()) exit(1);

    if (idx.od == 0) fileName = "sfc30031.bin";
    if (idx.od == 1) fileName = "sfc30101.bin";
    if (idx.od == 2) fileName = "sfc30301.bin";
    binaryFile[2].open(fileName.c_str(), ios::binary | ios::in);
    if (!binaryFile[2].is_open()) exit(1);

    for (idx.ls = 1; idx.ls < LS_SIZE; ++idx.ls) {
      for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
        for (idx.lon = SURFACE_LON_SIZE - 1; idx.lon > 0; --idx.lon) {
          idx.hgt = 0;
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ils), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&ylat), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&jlon), sizeof(int));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od] = 0.0;
            surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od] = 0.0;
          }

          for (idx.hgt = 1; idx.hgt < SURFACE_HEIGHT_SIZE; ++idx.hgt) {
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&jlon), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          }
        }
      }
    }
    binaryFile[0].close();
    binaryFile[1].close();
    binaryFile[2].close();
  }

  // Copy LS = 360 (at index LS_SIZE - 1) into LS = 0 (at index 0)
  // Copy LON = 360 (at index SURFACE_LON_SIZE - 1) into LON = 0 (at index 0)
  for (size_t p = 0; p < PARAM_SIZE; ++p) {
    for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
      for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
        for (idx.hgt = 0; idx.hgt < SURFACE_HEIGHT_SIZE; ++idx.hgt) {
          for (idx.lon = 1; idx.lon < SURFACE_LON_SIZE; ++idx.lon) {
            surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][0][idx.od] = surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.od];
            surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][0][idx.od] = surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.od];
            surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][0][idx.od] = surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.od];
            if (idx.lon == SURFACE_LON_SIZE - 1) {
              for (idx.ls = 0; idx.ls < LS_SIZE; ++idx.ls) {
                surfaceTemperatures[p][idx.hgt][idx.lat][0][idx.ls][idx.od] = surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od];
                surfaceEWWinds[p][idx.hgt][idx.lat][0][idx.ls][idx.od] = surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od];
                surfaceNSWinds[p][idx.hgt][idx.lat][0][idx.ls][idx.od] = surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.od];
              }
            }
          }
        }
      }
    }
  }

}

//! \brief Write the new binary data file format.
//!
//! The new file format consists of a repeated sequence of (block size, data block) pairs.
//! Output is a binary file consisting of greals.  If the definition of greal is modified,
//! then this file will need to be recreated.
void MGCMSurfaceInterpolator::writeSurfaceData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "MGCM_surface_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * OD_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceTemperatures);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceEWWinds);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceNSWinds);

  binaryFile.close();
}
#endif

} // namespace