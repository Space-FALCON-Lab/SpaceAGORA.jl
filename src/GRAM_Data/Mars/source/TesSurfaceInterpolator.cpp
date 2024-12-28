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
#include "TesSurfaceInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool TesSurfaceInterpolator::isInitialized = false;
greal TesSurfaceInterpolator::surfaceHeights[SURFACE_HEIGHT_SIZE] = { greal(0.0),  greal(0.005),  greal(0.03) };
greal TesSurfaceInterpolator::surfaceTemperatures[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesSurfaceInterpolator::surfaceEWWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
greal TesSurfaceInterpolator::surfaceNSWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};

//! \copydoc Atmosphere::Atmosphere()
TesSurfaceInterpolator::TesSurfaceInterpolator()
: TesLowerInterpolator()
{
  initializeData();
}

//! \fn  TesSurfaceInterpolator::TesSurfaceInterpolator(const TesSurfaceInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TesSurfaceInterpolator::~TesSurfaceInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  TesSurfaceInterpolator::getUpperFairingHeight()
//! \brief Gets the height of the upper bound of the fairing zone.
//! \returns The height \units{km} of the upper bound of the fairing zone.

//! \fn  TesSurfaceInterpolator::getLowerFairingHeight()
//! \brief Gets the height of the lower bound of the fairing zone.
//! \returns The height \units{km} of the lower bound of the fairing zone.

//! \fn  TesSurfaceInterpolator::getLowestBoundaryLayer()
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
void TesSurfaceInterpolator::updateIndices(size_t heightIndexOffset)
{
  // Find the height index of the first MGCM level above topographic surface height + 1 km.                                  
  // Heights might be negative here, so use int instead of size_t
  int kidx = int(floor(7.0 + surfaceHeight + 0.3));
  if (kidx >= 16) {
    kidx = int(floor(15.0 + (surfaceHeight + 1.0) / 5.0));
  }
  oneKilometerIndex = max(0, kidx - 1);  // formerly k1st

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

  // Compute the longitude WIND index and displacement for MGCM data
  size_t widx = size_t(floor(position.getLongitude(WEST_POSITIVE) / longitudeStepSize + 0.5));
  baseIndex.wlon = clampSize(widx, SURFACE_LON_SIZE - 1);
  displacements.wlon = (position.getLongitude(WEST_POSITIVE) - longitudeStepSize * (greal(baseIndex.wlon) - 0.5)) / longitudeStepSize;
  if (displacements.wlon > 40.0) {
    displacements.wlon -= 40.0;
  }

  // Update lat, ls, and tyr indices.
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
//! \param lonIndex Longitude index is either index.lon or index.wlon.
//! \param poleFactor An amplitude dampening factor (0 to 1).
//!
//! \b Inputs
//! \arg #index
//!
//! \returns Tidal parameters.
MarsTideParameters TesSurfaceInterpolator::getTideParameters(const greal data[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE], size_t latIndex, size_t lonIndex, greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][latIndex][lonIndex][index.ls][index.tyr];
  tp.amplitude1D = data[1][index.hgt][latIndex][lonIndex][index.ls][index.tyr] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][latIndex][lonIndex][index.ls][index.tyr] * poleFactor;
  tp.phase1D = data[2][index.hgt][latIndex][lonIndex][index.ls][index.tyr];
  tp.phase2D = data[4][index.hgt][latIndex][lonIndex][index.ls][index.tyr];
  return tp;
}


//! \brief Computes TPD by interpolating MGCM data in 3-D (latitude, longitude, and Ls).
//!
//! Interpolates MGCM data in 3-D (latitude, longitude, and Ls) for a given TES year and MGCM height index and
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
//! \retval #pressureScaleHeightDaily
//! \retval #temperatureMin
//! \retval #temperatureMax
//! \retval #height
void TesSurfaceInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 3-D interpolation for surface data or 2-D data from the Lower MGCM domain.
  //-------------------------------------------------------------------------------------------------//
  greal   tm[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of current T near surface
  greal   um[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of current U near surface
  greal   vm[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of current V near surface
  greal tday[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of daily mean T near surface
  greal tmax[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of daily max T near surface
  greal tmin[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of daily min T near surface
  greal uday[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of daily mean U near surface
  greal vday[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of daily mean V near surface
  greal tsrf[2][2][2];  // 3-d (lat, lon, Ls) array 'box' of current T at ground surface

  greal    t1km[2][2];  // 2-d (lat, Ls) array 'box' of current T at the 1 km index
  greal    d1km[2][2];  // 2-d (lat, Ls) array 'box' of current D at the 1 km index
  greal t1kmday[2][2];  // 2-d (lat, Ls) array 'box' of daily mean T at the 1 km index
  greal p1kmday[2][2];  // 2-d (lat, Ls) array 'box' of daily mean P at the 1 km index
  greal r1kmday[2][2];  // 2-d (lat, Ls) array 'box' of daily mean R at the 1 km index
  greal t1kmmax[2][2];  // 2-d (lat, Ls) array 'box' of daily max T at the 1 km index
  greal t1kmmin[2][2];  // 2-d (lat, Ls) array 'box' of daily min T at the 1 km index
  greal d1kmmax[2][2];  // 2-d (lat, Ls) array 'box' of daily max P at the 1 km index
  greal d1kmmin[2][2];  // 2-d (lat, Ls) array 'box' of daily min P at the 1 km index
  greal  plower[2][2];  // 2-d (lat, Ls) array 'box' of daily mean P at bottom of the surface layer
  greal  pupper[2][2];  // 2-d (lat, Ls) array 'box' of daily mean P at top of the surface layer
  greal  dlower[2][2];  // 2-d (lat, Ls) array 'box' of daily mean D at bottom of the surface layer
  greal  dupper[2][2];  // 2-d (lat, Ls) array 'box' of daily mean D at top of the surface layer

  //=================================================================================================//
  // Populate the 3-D interpolation cubes.  Loop over lat, lon, and ls.
  //-------------------------------------------------------------------------------------------------//
  index.tyr = baseIndex.tyr;
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
        // Update index values (base + counter).  Counter = 0 means the bottom of the 
        // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
        index.lat = baseIndex.lat + ilat;
        index.wlat = baseIndex.wlat + ilat;
        index.ls = baseIndex.ls + ils;
        index.lon = baseIndex.lon + ilon;
        index.wlon = baseIndex.wlon + ilon;

        // Cubes are populated by getting the tide parameters for the current index.
        // Then use the tide parameters to get the net value for the current solar time.
        // Temperature
        MarsTideParameters t = getTideParameters(surfaceTemperatures, index.lat, index.lon, poleFactor);
        tm[ilat][ilon][ils] = getTideValue(t, solarTime, UVT);

        // Mean Surface (Ground) Temperature at the current index (height index = 0)
        tsrf[ilat][ilon][ils] = surfaceTemperatures[0][0][index.lat][index.lon][index.ls][index.tyr];

        // EW Wind (note that we are using LAT, WLON, and the wind pole factor)
        MarsTideParameters u = getTideParameters(surfaceEWWinds, index.lat, index.wlon, windPoleFactor);
        um[ilat][ilon][ils] = getTideValue(u, solarTime, UVT);

        // NS Wind (note that we are using WLAT, WLON, and the wind pole factor)
        MarsTideParameters v = getTideParameters(surfaceNSWinds, index.wlat, index.lon, windPoleFactor);
        vm[ilat][ilon][ils] = getTideValue(v, solarTime, UVT);

        // Save the daily means for Temperature and winds.
        tday[ilat][ilon][ils] = t.diurnalMean;
        uday[ilat][ilon][ils] = u.diurnalMean;
        vday[ilat][ilon][ils] = v.diurnalMean;

        // If directed, compute daily extrema for Temperature
        if (computeMinMax) {
          getSolMinMax(t, UVT, tmin[ilat][ilon][ils], tmax[ilat][ilon][ils]);
        }
      } // end loop for Ls
    } // end loop for lon
  } // end loop for lat

  //=================================================================================================//
  //	use 3-d interpolation to get cur and day T, U, and V; max and min T; and surface T             //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp3d(displacements.lat, displacements.lon, displacements.ls);
  // Temperature (current and daily average)
  temperature = interp3d.linear(tm);
  temperatureDaily = interp3d.linear(tday);

  // Surface Temperature
  greal surfaceTemperature = interp3d.linear(tsrf);

  // If directed, Temperature extrema
  if (computeMinMax) {
    temperatureMax = interp3d.linear(tmax);
    temperatureMin = interp3d.linear(tmin);
  }

  // Create an interpolator for EW winds using wlon.
  Interpolator interp3dew(displacements.lat, displacements.wlon, displacements.ls);
  // EW Wind (current and daily average)
  ewWind = interp3dew.linear(um);
  ewWindDaily = interp3dew.linear(uday);

  // Create an interpolator for NW winds using wlat and wlon.
  Interpolator interp3dns(displacements.wlat, displacements.lon, displacements.ls);
  // NS Wind (current and daily average)
  nsWind = interp3dns.linear(vm);
  nsWindDaily = interp3dns.linear(vday);

  //=================================================================================================//
  // Now prepare to read MGCM data to compute scale heights to get P and D in surface layer          //
  // Populate the 2-D interpolation cubes.  Loop over lat and ls.
  //-------------------------------------------------------------------------------------------------//
  // MGCM height index just below surface height (or 0 if height is below MOLA areoid)
  size_t surfaceLayerIndex = size_t(floor(max(0.0, double(oneKilometerIndex) - 1.0)));

  // Save the current height index and substitute the surface height index.
  size_t saveIndex = index.hgt;
  index.hgt = oneKilometerIndex;

  // no interpolation over the TES year
  index.tyr = baseIndex.tyr;
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
      // Update index values (base + counter).  Counter = 0 means the bottom of the 
      // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
      index.lat = baseIndex.lat + ilat; 
      index.ls = baseIndex.ls + ils;

      // Cubes are populated by getting the tide parameters for the current index.
      // Then use the tide parameters to get the net value for the current solar time.
      // Temperature
      MarsTideParameters t = getTideParameters(mgcmTemperature, index.lat, poleFactor);
      t1km[ilat][ils] = getTideValue(t, solarTime, UVT);

      // Density
      MarsTideParameters d = getTideParameters(mgcmDensity, index.lat, poleFactor);
      d1km[ilat][ils] = getTideValue(d, solarTime, DP);

      // Get the daily mean value at the 1 km level
      p1kmday[ilat][ils] = mgcmPressure[oneKilometerIndex][index.lat][index.ls][index.tyr];
      t1kmday[ilat][ils] = t.diurnalMean;  

      // Use the ideal gas law to compute R.
      if (abs(t.diurnalMean) > 0.0 && abs(d.diurnalMean) > 0.0) {
        r1kmday[ilat][ils] = p1kmday[ilat][ils] / (t.diurnalMean * d.diurnalMean);
      }
      else {
        r1kmday[ilat][ils] = 190.0;
      }

      // Get mean pressure and density at the Surface
      dlower[ilat][ils] = mgcmDensity[0][surfaceLayerIndex][index.lat][index.ls][index.tyr];     
      plower[ilat][ils] = mgcmPressure[surfaceLayerIndex][index.lat][index.ls][index.tyr];       
      dupper[ilat][ils] = mgcmDensity[0][surfaceLayerIndex + 1][index.lat][index.ls][index.tyr]; 
      pupper[ilat][ils] = mgcmPressure[surfaceLayerIndex + 1][index.lat][index.ls][index.tyr];   

      // If directed, compute daily extrema for Temperature and Pressure
      if (computeMinMax) {
        getSolMinMax(d, DP, d1kmmin[ilat][ils], d1kmmax[ilat][ils]);
        getSolMinMax(t, UVT, t1kmmin[ilat][ils], t1kmmax[ilat][ils]);
      }
    } // end loop for Ls
  } // end loop for latitude

  // Restore the current height index
  index.hgt = saveIndex;

  //=================================================================================================//
  //  use 2-D interpolation (in LAT and LS) to compute final values                                  //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp2d(displacements.lat, displacements.ls);

  // Change in height over the surface layer.
  greal dHeight = 1.0;
  if (oneKilometerIndex >= 16) {
    dHeight = 5.0;
  }

  // Pressure scale height at surface latyer
  //greal plow = interp2d.linear(plower);           // daily mean P at bottom of the surface layer
  //greal phigh = interp2d.linear(pupper);           // daily mean P at top of the surface layer
  greal plow = interp2d.log(plower);           // daily mean P at bottom of the surface layer
  greal phigh = interp2d.log(pupper);           // daily mean P at top of the surface layer
  pressureScaleHeightDaily = dHeight / log(plow / phigh);

  // Density scale height at surface layer
  //greal dlow = interp2d.linear(dlower);            // daily mean D at bottom of the surface layer
  //greal dhigh = interp2d.linear(dupper);           // daily mean D at top of the surface layer
  greal dlow = interp2d.log(dlower);            // daily mean D at bottom of the surface layer
  greal dhigh = interp2d.log(dupper);           // daily mean D at top of the surface layer
  greal densityScaleHeightDaily = dHeight / log(dlow / dhigh);

  // daily average gas constant
  specificGasConstant = interp2d.linear(r1kmday);

  // Current T, D, and P at the 1 km height
  greal t1k = interp2d.linear(t1km);             // current T at the 1 km height
  //greal d1k = interp2d.linear(d1km);             // current D at the 1 km height
  greal d1k = interp2d.log(d1km);             // current D at the 1 km height
  greal p1k = specificGasConstant * d1k * t1k; // current P at the 1 km height

  // Daily means at the 1 km height
  greal t1kday = interp2d.linear(t1kmday);             // daily mean T at the 1 km height
  //greal p1kday = interp2d.linear(p1kmday);             // daily mean P at the 1 km height
  greal p1kday = interp2d.log(p1kmday);             // daily mean P at the 1 km height

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
  greal z1k = -5.0 + greal(oneKilometerIndex);
  if (oneKilometerIndex >= 16) {
	  z1k = 10.0 + 5.0 * (oneKilometerIndex - 16.0);
  }
  greal deltaHeight = z1k - height;

  // extrapolate downward from 1 km to current hgtIndex to get current P
  pressure = p1k * exp(deltaHeight / pressureScaleHeight);
  // get current D at hgtIndex from gas law
  density = pressure / (specificGasConstant * temperature);

  // extrapolate downward from 1 km to current hgtIndex to get daily mean P
  pressureDaily = p1kday * exp(deltaHeight / pressureScaleHeightDaily);                                                 // extrapolate downward from k1st to hgtIndex to get daily mean P
  // get daily mean D at hgtIndex from gas law
  densityDaily = pressureDaily / (specificGasConstant * temperatureDaily);                                           // get daily mean D at hgtIndex from gas law

  // If directed, Temperature and Density extrema
  if (computeMinMax) {
    // daily max/min T at the 1 km height
    greal tmax1 = interp2d.linear(t1kmmax); 

    // daily max/min D at the 1 km height
    //greal d1max = interp2d.linear(d1kmmax);
    //greal d1min = interp2d.linear(d1kmmin);
    greal d1max = interp2d.log(d1kmmax);
    greal d1min = interp2d.log(d1kmmin);

    // get max/min D scale height
    greal hdmax = densityScaleHeightDaily * 0.5 * (tmax1 + temperatureMin) / tbarsurf;                                        // get maximum P scale height
    greal hdmin = densityScaleHeightDaily * 0.5 * (tmax1 + temperatureMax) / tbarsurf;                                        // get minimum P scale height

    // extrapolate max/min D downward from i km to hgtIndex to get daily max/min D
    densityMax = d1max * exp(deltaHeight / hdmax);
    densityMin = d1min * exp(deltaHeight / hdmin);
  }
}

//! \brief Reads TES surface atmosphere tidal parameter data.
//!
//! This method reads TES_surface_data.bin, a binary data file, to populate the data arrays for the
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
void TesSurfaceInterpolator::initializeData() {
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  //readLegacySurfaceData();
  //writeSurfaceData();

  try {
    // Open a binary file stream.
    ifstream binaryFile;
    binaryFile.open(dataPath + tesSurfaceFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * TES_YEAR_SIZE;
    // Read the data blocks.
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceTemperatures);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceEWWinds);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceNSWinds);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + tesSurfaceFileName;
  }
}

#ifdef CONVERT_MARS_DATA
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void TesSurfaceInterpolator::readLegacySurfaceData()
{

  int bytes, ils, jlon;
  double ylat;


  MIBIndices idx;
  for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
    string fileName = "";
    if (idx.tyr == 0) fileName = "sfc00y11.bin";
    if (idx.tyr == 1) fileName = "sfc00y21.bin";
    ifstream binaryFile[3];
    binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);
    if (!binaryFile[0].is_open()) exit(1);

    if (idx.tyr == 0) fileName = "sfc05y11.bin";
    if (idx.tyr == 1) fileName = "sfc05y21.bin";
    binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);
    if (!binaryFile[1].is_open()) exit(1);

    if (idx.tyr == 0) fileName = "sfc30y11.bin";
    if (idx.tyr == 1) fileName = "sfc30y21.bin";
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
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&surfaceTemperatures[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
          binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr] = 0.0;
            surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr] = 0.0;
          }

          for (idx.hgt = 1; idx.hgt < SURFACE_HEIGHT_SIZE; ++idx.hgt) {
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&jlon), sizeof(int));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceTemperatures[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceEWWinds[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[0][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[1][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[2][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[3][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&surfaceNSWinds[4][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr]), sizeof(double));
            binaryFile[idx.hgt].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          }
        }
      }
    }
    binaryFile[0].close();
    binaryFile[1].close();
    binaryFile[2].close();
  }
  for (size_t p = 0; p < PARAM_SIZE; ++p) {
    for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
      for (idx.lat = 0; idx.lat < MGCM_LAT_SIZE; ++idx.lat) {
        for (idx.hgt = 0; idx.hgt < SURFACE_HEIGHT_SIZE; ++idx.hgt) {
          for (idx.lon = 1; idx.lon < SURFACE_LON_SIZE; ++idx.lon) {
            surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][0][idx.tyr] = surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.tyr];
            surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][0][idx.tyr] = surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.tyr];
            surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][0][idx.tyr] = surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][LS_SIZE - 1][idx.tyr];
            if (idx.lon == SURFACE_LON_SIZE - 1) {
              for (idx.ls = 0; idx.ls < LS_SIZE; ++idx.ls) {
                surfaceTemperatures[p][idx.hgt][idx.lat][0][idx.ls][idx.tyr] = surfaceTemperatures[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr];
                surfaceEWWinds[p][idx.hgt][idx.lat][0][idx.ls][idx.tyr] = surfaceEWWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr];
                surfaceNSWinds[p][idx.hgt][idx.lat][0][idx.ls][idx.tyr] = surfaceNSWinds[p][idx.hgt][idx.lat][idx.lon][idx.ls][idx.tyr];
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
void TesSurfaceInterpolator::writeSurfaceData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "TES_surface_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * TES_YEAR_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceTemperatures);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceEWWinds);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)surfaceNSWinds);

  binaryFile.close();
}
#endif

} // namespace