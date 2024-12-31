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
#include "TesUpperInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool TesUpperInterpolator::isInitialized = false;
greal TesUpperInterpolator::mtgcmTemperature[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmPressure[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmDensity[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmEWWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmNSWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmThermosphereBaseHeight[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
greal TesUpperInterpolator::mtgcmF107[MTGCM_F107_SIZE] = { 70.0, 130.0, 200.0 };

//! \copydoc Atmosphere::Atmosphere()
TesUpperInterpolator::TesUpperInterpolator()
: TesInterpolator()
{
  initializeData();
}

//! \fn  TesUpperInterpolator::TesUpperInterpolator(const TesUpperInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TesUpperInterpolator::~TesUpperInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  TesUpperInterpolator::setHeightOffset(greal offset)
//! \brief Set the height offset.
//! \param offset The height offset \units{km}.

//! \fn  TesUpperInterpolator::setF107(greal value)
//! \brief Set the solar flux F10.7 at 1 AU value.
//! \param value F10.7 at 1 AU.

//! \fn  TesUpperInterpolator::getUpperFairingHeight()
//! \brief Gets the height of the upper bound of the fairing zone.
//! \returns The height \units{km} of the upper bound of the fairing zone.

//! \fn  TesUpperInterpolator::getLowerFairingHeight()
//! \brief Gets the height of the lower bound of the fairing zone.
//! \returns The height \units{km} of the lower bound of the fairing zone.

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
void TesUpperInterpolator::updateIndices(size_t heightIndexOffset)
{
  // Compute the base height index for interpolating MTGCM data
  int idx = int(floor(((height - heightOffset) - 80.0) / 5.0));
  idx = max(0, idx);
  // If the heightOffset is very negative, set the adjustment factor.
  if (idx == 0 && heightOffset < -4.0) {
    idx = 1;
    heightOffsetAdjustment = 1;
  }
  else {
    heightOffsetAdjustment = 0;
  }
  baseIndex.hgt = min(MTGCM_HEIGHT_SIZE - 2, size_t(idx));
  // Apply the height index offset
  index.hgt = baseIndex.hgt + heightIndexOffset;

  // array index and displacement for solar activity (only 3 data levels for MTGCM data)
  baseIndex.f107 = 0;
  if (F107 > mtgcmF107[1]) {
    baseIndex.f107 = 1;
  }
  displacements.f107 = log(F107 / mtgcmF107[baseIndex.f107]) / log(mtgcmF107[baseIndex.f107 + 1] / mtgcmF107[baseIndex.f107]);  // relative displacement in solar activity, F10.7
// NOTE: Legacy MarsGRAM does not clamp to bounds (that is, it extrapolates)
//  displacements.f107 = clamp(displacements.f107, 0.0, 1.0);

  updateBaseIndices(MTGCM_LAT_SIZE, 2.5);
}

//! \brief Gets tide parameters corresponding to the current #index.
//!
//! This method looks up tidal parameters in the supplied MTGCM \p data array.  The current
//! #index is used for lookup within the array. The \p poleFactor is a dampening
//! factor that is applied to the amplitudes to smooth data near the poles.
//!
//! \param data Multi-dimensional array of tide parameters.
//! \param poleFactor An amplitude dampening factor (0 to 1).
//!
//! \b Inputs
//! \arg #index
//!
//! \returns Tidal parameters.
MarsTideParameters TesUpperInterpolator::getTideParameters(const greal data[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE], greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][index.lat][index.ls][index.tyr][index.f107];
  tp.amplitude1D = data[1][index.hgt][index.lat][index.ls][index.tyr][index.f107] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][index.lat][index.ls][index.tyr][index.f107] * poleFactor;
  tp.phase1D = data[2][index.hgt][index.lat][index.ls][index.tyr][index.f107];
  tp.phase2D = data[4][index.hgt][index.lat][index.ls][index.tyr][index.f107];
  return tp;
}

//! \brief Gets tide parameters corresponding to the current #index.
//!
//! This method looks up tidal parameters in the supplied MTGCM \p mtgcmThermosphereBaseHeight array.  The current
//! #index is used for lookup within the array. The \p poleFactor is a dampening
//! factor that is applied to the amplitudes to smooth data near the poles.
//!
//! \param data Multi-dimensional array of tide parameters.
//! \param poleFactor An amplitude dampening factor (0 to 1).
//!
//! \b Inputs
//! \arg #index
//!
//! \returns Tidal parameters.
MarsTideParameters TesUpperInterpolator::getTideParameters(const greal data[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE], greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.lat][index.ls][index.tyr][index.f107];
  tp.amplitude1D = data[1][index.lat][index.ls][index.tyr][index.f107] * poleFactor;
  tp.amplitude2D = data[3][index.lat][index.ls][index.tyr][index.f107] * poleFactor;
  tp.phase1D = data[2][index.lat][index.ls][index.tyr][index.f107];
  tp.phase2D = data[4][index.lat][index.ls][index.tyr][index.f107];
  return tp;
}


//! \brief Computes TPD by interpolating MTGCM data in 3-D (latitude, Ls, and F10.7).
//!
//! Interpolates MTGCM data in 4-D (latitude, Ls,and F10.7) for a given TES year, MTGCM height index and
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
void TesUpperInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 3-D interpolation
  //-------------------------------------------------------------------------------------------------//
  greal tday[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily mean T values
	greal tmax[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily maximum T values
	greal tmin[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily minimum T values
	greal dmax[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily maximum D values
	greal dmin[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily minimum D values
	greal pday[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily mean P values
	greal uday[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily mean U values
	greal vday[2][2][2] = {};    // 3-D 'box' array (lat, Ls, F10.7) of daily mean V values
	greal r0[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of gas constant values
	greal pm[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current P values
	greal dm[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current D values
	greal tm[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current T values
	greal um[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current U values
  greal vm[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current V values
  greal zt[2][2][2] = {};      // 3-D 'box' array (lat, Ls, F10.7) of current Zf values

  //=================================================================================================//
  // Populate the interpolation cubes.  Loop over lat, ls, and f107.
  //-------------------------------------------------------------------------------------------------//
  // no interpolation over the TES year
  index.tyr = baseIndex.tyr;
  // Start LAT loop
  for (int ilat = 0; ilat < 2; ilat++) {                                                      // begin for loop over 2 latitude data array indices
    // Pole factors are applied to the tidal amplitudes. 
    // Away from the poles, the factor is 1.
    greal poleFactor = 1.0;
    // At the poles, set the factor to tpolefac.
    // That's south pole at the bottom of the layer
    // Or the north pole at the top of the layer
    if ((baseIndex.lat == 0 && ilat == 0) || (baseIndex.lat == MTGCM_LAT_SIZE - 2 && ilat == 1)) {
      poleFactor = displacements.tpolefac;
    }

    // Start LS loop
    for (int ils = 0; ils < 2; ils++) {
      // Start F107 loop
      for (int if10 = 0; if10 < 2; if10++) {
        // Update index values (base + counter).  Counter = 0 means the bottom of the 
        // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
        index.lat = baseIndex.lat + ilat;     
					index.ls = baseIndex.ls + ils;      
					index.f107 = baseIndex.f107 + if10; 

          // Cubes are populated by getting the tide parameters for the current index.
          // Then use the tide parameters to get the net value for the current solar time.
          // Temperature
          MarsTideParameters t = getTideParameters(mtgcmTemperature, poleFactor);
					tm[ilat][ils][if10] = getTideValue(t, solarTime, UVT);               

          // Pressure
          MarsTideParameters p = getTideParameters(mtgcmPressure, poleFactor);
          pm[ilat][ils][if10] = getTideValue(p, solarTime, BIG_DP);            

          // Density
          MarsTideParameters d = getTideParameters(mtgcmDensity, poleFactor);
					dm[ilat][ils][if10] = getTideValue(d, solarTime, BIG_DP);            

          // compute gas constant from gas law 
          // if both T and D mean values are GT 0...
          if (abs(t.diurnalMean) > 0.0 && abs(d.diurnalMean) > 0.0) {
            r0[ilat][ils][if10] = p.diurnalMean / (t.diurnalMean*d.diurnalMean);
          }
          else {
            // else set gas constant to global mean value
            r0[ilat][ils][if10] = 190.0;
          }

          // EW Wind
          MarsTideParameters u = getTideParameters(mtgcmEWWind, poleFactor);
					um[ilat][ils][if10] = getTideValue(u, solarTime, UVT);

          // NS Wind
          MarsTideParameters v = getTideParameters(mtgcmNSWind, poleFactor);
          vm[ilat][ils][if10] = getTideValue(v, solarTime, UVT); 

          // Process tide coefficients for Thermosphere Base Height at current indices
          MarsTideParameters z = getTideParameters(mtgcmThermosphereBaseHeight, poleFactor);
          zt[ilat][ils][if10] = getTideValue(z, solarTime, UVT); 

          // Save the daily means for Temperature, Pressure, Density, and winds.
          tday[ilat][ils][if10] = t.diurnalMean;
          pday[ilat][ils][if10] = p.diurnalMean; 
          uday[ilat][ils][if10] = u.diurnalMean; 
          vday[ilat][ils][if10] = v.diurnalMean; 

          // If directed, compute daily extrema for Temperature and Density
          if (computeMinMax) {
            // Max and Min temperatures at corners of the cube
            getSolMinMax(t, UVT, tmin[ilat][ils][if10], tmax[ilat][ils][if10]);
            // Max and Min at corners of the cube
            getSolMinMax(d, BIG_DP, dmin[ilat][ils][if10], dmax[ilat][ils][if10]);
          }
        } // end loop over F10.7 indices
		} // end loop over Ls indices
	} // end loop over latitude indices

  //=================================================================================================//
  //  use 3-D interpolation (in LAT, LS,  F107) routine to get final interpolated values             //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp(displacements.lat, displacements.ls, displacements.f107);

  // Temperature (current and daily average)
  temperature = interp.linear(tm);
	temperatureDaily = interp.linear(tday);

  // Pressure (current and daily average)
  pressure = interp.log(pm);
  pressureDaily = interp.log(pday);

  // Gas Constant (mean)
  specificGasConstant = interp.linear(r0);

  // Density (current from interpolation and daily average from gas law)
  density = interp.log(dm);
  densityDaily = pressureDaily / (specificGasConstant * temperatureDaily);

  // If directed, Temperature and Density extrema
  if (computeMinMax) {
    temperatureMax = interp.linear(tmax); 
    temperatureMin = interp.linear(tmin); 
    densityMax = interp.log(dmax);
    densityMin = interp.log(dmin);
  }

  // EW Wind (current and daily average)
  ewWind = interp.linear(um);   // current U
  ewWindDaily = interp.linear(uday); // daily mean U

  // NS Wind (current and daily average)
  nsWind = interp.linear(vm);   // current V
  nsWindDaily = interp.linear(vday); // daily mean V

  // Thermosphere base height
  thermosphereBaseHeight    = interp.linear(zt);

  // Set height to the value of the bottom of the layer.
  height = 80.0 + 5.0 * index.hgt;
}

//! \brief Reads TES upper atmosphere tidal parameter data.
//! 
//! This method reads TES_upper_data.bin, a binary data file, to populate the data arrays for the 
//! upper atmosphere data.  The data is stored in static memory.
//! The data format consists three data blocks.  Each data block is preceeded by a size_t with 
//! the number of doubles in the data block.
//!
//! \b Inputs
//! \arg #isInitialized          
//!
//! \retval #mtgcmThermosphereBaseHeight     
//! \retval #mtgcmTemperature     
//! \retval #mtgcmPressure     
//! \retval #mtgcmDensity     
//! \retval #mtgcmEWWind     
//! \retval #mtgcmNSWind     
void TesUpperInterpolator::initializeData() {
  // Only initialize the date once
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // These lines facilitate conversion of file formats.
  //readLegacyUpperData();
  //writeUpperData();

  try {
    // Open a binary file stream.
    ifstream binaryFile;
    binaryFile.open(dataPath + tesUpperFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE * MTGCM_F107_SIZE;

    // Read the thermosphere base height data
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmThermosphereBaseHeight);

    // Compute the size of the remaining data blocks
    dataSize *= MTGCM_HEIGHT_SIZE;

    // Read the remaining data blocks.
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmTemperature);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmPressure);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmDensity);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmEWWind);
    readBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmNSWind);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg +  dataPath + tesUpperFileName;
  }
}

#ifdef CONVERT_MARS_DATA

//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void TesUpperInterpolator::readLegacyUpperData()
{
  int bytes, ils, ihgt;
  double ylat;


  MIBIndices idx;
  for (idx.f107 = 0; idx.f107 < MTGCM_F107_SIZE; ++idx.f107) {
    string fileName = "";
    if (idx.f107 == 0) fileName = "zfTESls1.txt";
    if (idx.f107 == 1) fileName = "zfTESms1.txt";
    if (idx.f107 == 2) fileName = "zfTEShs1.txt";
    ifstream textFile(fileName.c_str());
    for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
      if (idx.f107 == 0) {
        if (idx.tyr == 0) fileName = "tpdlsy11.bin";
        if (idx.tyr == 1) fileName = "tpdlsy21.bin";
      }
      if (idx.f107 == 1) {
        if (idx.tyr == 0) fileName = "tpdmsy11.bin";
        if (idx.tyr == 1) fileName = "tpdmsy21.bin";
      }
      if (idx.f107 == 2) {
        if (idx.tyr == 0) fileName = "tpdhsy11.bin";
        if (idx.tyr == 1) fileName = "tpdhsy21.bin";
      }
      ifstream binaryFile[3];
      binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);

      if (idx.f107 == 0) {
        if (idx.tyr == 0) fileName = "uvlsy11.bin";
        if (idx.tyr == 1) fileName = "uvlsy21.bin";
      }
      if (idx.f107 == 1) {
        if (idx.tyr == 0) fileName = "uvmsy11.bin";
        if (idx.tyr == 1) fileName = "uvmsy21.bin";
      }
      if (idx.f107 == 2) {
        if (idx.tyr == 0) fileName = "uvhsy11.bin";
        if (idx.tyr == 1) fileName = "uvhsy21.bin";
      }
      binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);

      for (idx.ls = 1; idx.ls < LS_SIZE; ++idx.ls) {
        for (idx.lat = 0; idx.lat < MTGCM_LAT_SIZE; ++idx.lat) {
          for (idx.hgt = 0; idx.hgt < MTGCM_HEIGHT_SIZE; ++idx.hgt) {
            binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[0][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[1][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[2][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[3][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[4][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[1][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[2][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[3][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[4][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[0][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[1][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[2][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[3][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[4][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));


            binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[0][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[1][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[2][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[3][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[4][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[0][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[1][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[2][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[3][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[4][idx.hgt][idx.lat][idx.ls][idx.tyr][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          }
          greal iyr;
          int nhgttop;
          textFile >> iyr >> ils >> ylat
            >> mtgcmThermosphereBaseHeight[0][idx.lat][idx.ls][idx.tyr][idx.f107]
            >> mtgcmThermosphereBaseHeight[1][idx.lat][idx.ls][idx.tyr][idx.f107]
            >> mtgcmThermosphereBaseHeight[2][idx.lat][idx.ls][idx.tyr][idx.f107]
            >> mtgcmThermosphereBaseHeight[3][idx.lat][idx.ls][idx.tyr][idx.f107]
            >> mtgcmThermosphereBaseHeight[4][idx.lat][idx.ls][idx.tyr][idx.f107]
            >> nhgttop;

        }
      }
      binaryFile[0].close();
      binaryFile[1].close();
    }
    textFile.close();
  }

  for (idx.f107 = 0; idx.f107 < MTGCM_F107_SIZE; ++idx.f107) {
    for (idx.tyr = 0; idx.tyr < TES_YEAR_SIZE; ++idx.tyr) {
      for (idx.lat = 0; idx.lat < MTGCM_LAT_SIZE; ++idx.lat) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          for (idx.hgt = 0; idx.hgt < MTGCM_HEIGHT_SIZE; ++idx.hgt) {
            mtgcmTemperature[p][idx.hgt][idx.lat][0][idx.tyr][idx.f107] = mtgcmTemperature[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
            mtgcmPressure[p][idx.hgt][idx.lat][0][idx.tyr][idx.f107] = mtgcmPressure[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
            mtgcmDensity[p][idx.hgt][idx.lat][0][idx.tyr][idx.f107] = mtgcmDensity[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
            mtgcmEWWind[p][idx.hgt][idx.lat][0][idx.tyr][idx.f107] = mtgcmEWWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
            mtgcmNSWind[p][idx.hgt][idx.lat][0][idx.tyr][idx.f107] = mtgcmNSWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
          }
          mtgcmThermosphereBaseHeight[p][idx.lat][0][idx.tyr][idx.f107] = mtgcmThermosphereBaseHeight[p][idx.lat][LS_SIZE - 1][idx.tyr][idx.f107];
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
void TesUpperInterpolator::writeUpperData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "TES_upper_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE * MTGCM_F107_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmThermosphereBaseHeight);

  dataSize *= MTGCM_HEIGHT_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmTemperature);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmPressure);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmDensity);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmEWWind);
  writeBinaryDataBlock(binaryFile, dataSize, (greal*)mtgcmNSWind);

  binaryFile.close();
}
#endif

} // namespace