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
#include "MGCMUpperInterpolator.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool MGCMUpperInterpolator::isInitialized = false;
greal MGCMUpperInterpolator::mtgcmTemperature[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmPressure[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmDensity[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmEWWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmNSWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmThermosphereBaseHeight[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
greal MGCMUpperInterpolator::mtgcmF107[MTGCM_F107_SIZE] = { 70.0, 130.0 };

//! \copydoc Atmosphere::Atmosphere()
MGCMUpperInterpolator::MGCMUpperInterpolator()
: MGCMInterpolator()
{
  initializeData();
}

//! \fn  MGCMUpperInterpolator::MGCMUpperInterpolator(const MGCMUpperInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MGCMUpperInterpolator::~MGCMUpperInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MGCMUpperInterpolator::setHeightOffset(greal offset)
//! \brief Set the height offset.
//! \param offset The height offset \units{km}.

//! \fn  MGCMUpperInterpolator::setF107(greal value)
//! \brief Set the solar flux F10.7 at 1 AU value.
//! \param value F10.7 at 1 AU.

//! \fn  MGCMUpperInterpolator::getUpperFairingHeight()
//! \brief Gets the height of the upper bound of the fairing zone.
//! \returns The height \units{km} of the upper bound of the fairing zone.

//! \fn  MGCMUpperInterpolator::getLowerFairingHeight()
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
void MGCMUpperInterpolator::updateIndices(size_t heightIndexOffset)
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

  // array index and displacement for solar activity (only 2 data levels for MTGCM data)
  baseIndex.f107 = 0;
  displacements.f107 = log(F107 / mtgcmF107[baseIndex.f107]) / log(mtgcmF107[baseIndex.f107 + 1] / mtgcmF107[baseIndex.f107]);  // relative displacement in solar activity, F10.7
// NOTE: Legacy MarsGRAM does not clamp to bounds (that is, it extrapolates)
//  displacements.f107 = clamp(displacements.f107, 0.0, 1.0);

  // Update lat, ls, and od indices.
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
MarsTideParameters MGCMUpperInterpolator::getTideParameters(const greal data[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE], greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.hgt][index.lat][index.ls][index.od][index.f107];
  tp.amplitude1D = data[1][index.hgt][index.lat][index.ls][index.od][index.f107] * poleFactor;
  tp.amplitude2D = data[3][index.hgt][index.lat][index.ls][index.od][index.f107] * poleFactor;
  tp.phase1D = data[2][index.hgt][index.lat][index.ls][index.od][index.f107];
  tp.phase2D = data[4][index.hgt][index.lat][index.ls][index.od][index.f107];
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
MarsTideParameters MGCMUpperInterpolator::getTideParameters(const greal data[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE], greal poleFactor)
{
  MarsTideParameters tp;
  tp.diurnalMean = data[0][index.lat][index.ls][index.od][index.f107];
  tp.amplitude1D = data[1][index.lat][index.ls][index.od][index.f107] * poleFactor;
  tp.amplitude2D = data[3][index.lat][index.ls][index.od][index.f107] * poleFactor;
  tp.phase1D = data[2][index.lat][index.ls][index.od][index.f107];
  tp.phase2D = data[4][index.lat][index.ls][index.od][index.f107];
  return tp;
}

//! \brief Computes TPD by interpolating MTGCM data in 4-D (latitude, Ls, dust OD, and F10.7).
//!
//! Interpolates MTGCM data in 4-D (latitude, Ls, dust OD, and F10.7) for a given MTGCM height index and
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
void MGCMUpperInterpolator::update()
{
  //=================================================================================================//
  // Declare and initialize interpolation arrays.  These will hold the corner points
  // of the 4-D interpolation
  //-------------------------------------------------------------------------------------------------//
  greal tday[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily mean T values
	greal tmax[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily maximum T values
	greal tmin[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily minimum T values
	greal dmax[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily maximum D values
	greal dmin[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily minimum D values
	greal pday[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily mean P values
	greal uday[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily mean U values
	greal vday[2][2][2][2] = {};    // 4-D 'box' array (lat, Ls, OD, F10.7) of daily mean V values
	greal r0[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of gas constant values
	greal pm[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current P values
	greal dm[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current D values
	greal tm[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current T values
	greal um[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current U values
  greal vm[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current V values
  greal zt[2][2][2][2] = {};      // 4-D 'box' array (lat, Ls, OD, F10.7) of current Zf values

  //=================================================================================================//
  // Populate the interpolation cubes.  Loop over lat, ls, od, and f107.
  //-------------------------------------------------------------------------------------------------//
  // Start LAT loop
  for (int ilat = 0; ilat < 2; ilat++) {
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
      // Start OD loop
      for (int iod = 0; iod < 2; iod++) {
        // Start F107 loop
				for (int if10 = 0; if10 < 2; if10++) {
          // Update index values (base + counter).  Counter = 0 means the bottom of the 
          // baseIndex layer. Counter = 1 means the top of the baseIndex layer.
          index.lat = baseIndex.lat + ilat;
					index.ls = baseIndex.ls + ils;      
					index.od = baseIndex.od + iod;      
					index.f107 = baseIndex.f107 + if10; 

          // Cubes are populated by getting the tide parameters for the current index.
          // Then use the tide parameters to get the net value for the current solar time.
          // Temperature
          MarsTideParameters t = getTideParameters(mtgcmTemperature, poleFactor);
					tm[ilat][ils][iod][if10] = getTideValue(t, solarTime, UVT);

          // Pressure
          MarsTideParameters p = getTideParameters(mtgcmPressure, poleFactor);
          pm[ilat][ils][iod][if10] = getTideValue(p, solarTime, DP);

          // Density
          MarsTideParameters d = getTideParameters(mtgcmDensity, poleFactor);
					dm[ilat][ils][iod][if10] = getTideValue(d, solarTime, DP);

          // compute gas constant from gas law 
          // if both T and D mean values are GT 0...
          if (abs(t.diurnalMean) > 0.0 && abs(d.diurnalMean) > 0.0) {
            r0[ilat][ils][iod][if10] = p.diurnalMean / (t.diurnalMean * d.diurnalMean);
          }
          else {
            // else set gas constant to global mean value
            r0[ilat][ils][iod][if10] = 190.0; 
          }

          // EW Wind
          MarsTideParameters u = getTideParameters(mtgcmEWWind, poleFactor);
					um[ilat][ils][iod][if10] = getTideValue(u, solarTime, UVT);

          // NS Wind
          MarsTideParameters v = getTideParameters(mtgcmNSWind, poleFactor);
          vm[ilat][ils][iod][if10] = getTideValue(v, solarTime, UVT);

          // Process tide coefficients for Thermosphere Base Height at current indices
          MarsTideParameters z = getTideParameters(mtgcmThermosphereBaseHeight, poleFactor);
          zt[ilat][ils][iod][if10] = getTideValue(z, solarTime, UVT);

          // Save the daily means for Temperature, Pressure, Density, and winds.
          tday[ilat][ils][iod][if10] = t.diurnalMean;
          pday[ilat][ils][iod][if10] = p.diurnalMean;
          uday[ilat][ils][iod][if10] = u.diurnalMean;
          vday[ilat][ils][iod][if10] = v.diurnalMean;

          // If directed, compute daily extrema for Temperature and Density
          if (computeMinMax) {
            // Max and Min temperatures at corners of the cube
            getSolMinMax(t, UVT, tmin[ilat][ils][iod][if10], tmax[ilat][ils][iod][if10]);
            // Max and Min at corners of the cube
            getSolMinMax(d, DP, dmin[ilat][ils][iod][if10], dmax[ilat][ils][iod][if10]);
          }
        } // end loop over F10.7 indices
			} // end loop over dust OD indices
		} // end loop over Ls indices
	} // end loop over latitude indices

  //=================================================================================================//
  //  use 4-D interpolation (in LAT, LS, OD, F107) to get final interpolated values                  //
  //-------------------------------------------------------------------------------------------------//
  // Create an interpolator using displacements.
  Interpolator interp(displacements.lat, displacements.ls, displacements.od, displacements.f107);

  // Temperature (current and daily average)
  temperature = interp.linear(tm);
	temperatureDaily = interp.linear(tday);

  // Pressure (current and daily average)
  pressure = interp.log(pm); 
  pressureDaily = interp.log(pday);

  // Gas Constant (mean)
  specificGasConstant = interp.linear(r0);   // current gas constant

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
  ewWind = interp.linear(um);
  ewWindDaily = interp.linear(uday);

  // NS Wind (current and daily average)
  nsWind = interp.linear(vm);
  nsWindDaily = interp.linear(vday);

  // Thermosphere base height
  thermosphereBaseHeight    = interp.linear(zt);

  // Set height to the value of the bottom of the layer.
  height = 80.0 + 5.0 * index.hgt;
}

//! \brief Reads MGCM upper atmosphere tidal parameter data.
//! 
//! This method reads MGCM_upper_data.bin, a binary data file, to populate the data arrays for the 
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
void MGCMUpperInterpolator::initializeData() {
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
    binaryFile.open(dataPath + mgcmUpperFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * OD_SIZE * MTGCM_F107_SIZE;

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
    throw msg + dataPath + mgcmUpperFileName;
  }
}

#ifdef CONVERT_MARS_DATA
//! \brief Read the legacy data format.
//!
//! This method is provided to help convert the legacy text data into a binary format.
//! Use the writeLowerData() method to create the binary file.
void MGCMUpperInterpolator::readLegacyUpperData()
{
  static const greal Hf80[OD_SIZE][LS_SIZE][7] = {
    {
      {1.3100e+00, -1.7173e-03, 2.2086e-04, 1.1335e-06, -4.7842e-08, -1.4020e-10, 1.9581e-12   },
      {1.2855e+00, 2.5407e-03, 1.1140e-08, -3.3227e-06, 5.6398e-08, 4.1167e-10, -7.4666e-12    },
      {1.2491e+00, 3.5334e-03, -8.0617e-05, -3.5076e-06, 8.6122e-08, 4.0803e-10, -9.8090e-12   },
      {1.1595e+00, 1.9913e-03, 2.8612e-05, -2.5338e-06, 4.4030e-08, 3.1054e-10, -5.8952e-12    },
      {1.1046e+00, 3.3547e-03, -1.4369e-06, -3.5880e-06, 6.7749e-08, 3.3590e-10, -7.5981e-12   },
      {1.2877e+00, 4.6463e-03, -1.7076e-04, -4.6047e-06, 1.2693e-07, 4.6269e-10, -1.3291e-11   },
      {1.2717e+00, -2.5201e-04, 1.4159e-04, -2.8328e-07, 8.3218e-09, 5.5739e-12, -3.3071e-12   },
      {1.3287e+00, -1.1228e-03, 1.4243e-04, 7.3074e-07, -2.0977e-08, -1.4763e-10, -1.3986e-13  },
      {1.3602e+00, -1.9213e-03, 4.4605e-05, -2.4831e-07, -2.0548e-09, -9.2112e-12, -4.9610e-13 },
      {1.3292e+00, -2.0575e-03, 3.0917e-05, -6.1131e-10, 8.5858e-09, -6.7420e-11, -1.9742e-12  },
      {1.2491e+00, -2.6908e-03, 6.5230e-05, 9.6963e-07, 2.3130e-08, -1.7627e-10, -4.1092e-12   },
      {1.3873e+00, -1.8219e-03, 6.3232e-05, -1.9586e-08, -1.2456e-08, -4.4703e-11, -1.4046e-13 },
      {1.3100e+00, -1.7173e-03, 2.2086e-04, 1.1335e-06, -4.7842e-08, -1.4020e-10, 1.9581e-12   }
    },
    {
      {1.0691e+00, -2.2073e-03, 2.0158e-04, 2.4884e-06, 3.8384e-08, -2.0940e-10, -8.8724e-12   },
      {1.2242e+00, 4.0856e-03, 9.6936e-05, -5.5182e-06, 6.5142e-08, 6.7364e-10, -1.0937e-11    },
      {1.2203e+00, 2.5951e-03, 1.0236e-04, -2.4933e-06, 3.2277e-08, 2.9808e-10, -5.2909e-12    },
      {1.2738e+00, 3.5809e-04, 5.7095e-05, -5.9205e-07, 1.1426e-08, -3.6894e-11, -6.3835e-13   },
      {1.2828e+00, 1.7465e-03, 6.3646e-05, -2.0767e-06, 2.4885e-08, 1.4811e-10, -2.5337e-12    },
      {1.2416e+00, 3.0057e-03, 1.4008e-04, -3.7762e-06, 2.4889e-08, 5.1451e-10, -6.1947e-12    },
      {1.2061e+00, -3.3591e-04, 2.8237e-04, -6.3895e-07, -5.5637e-09, 3.0563e-11, -4.7434e-12  },
      {1.2602e+00, 7.0066e-04, 2.6984e-04, 1.3802e-06, -3.9678e-08, -2.8090e-10, 2.5941e-13    },
      {1.2806e+00, -6.8090e-04, 1.7769e-04, 5.8416e-07, -1.4360e-08, -1.4034e-10, -7.3100e-13  },
      {8.4267e-01, -1.5830e-03, 7.8998e-05, 1.5957e-06, 4.1777e-08, 2.7419e-10, -9.1151e-13    },
      {9.1483e-01, -4.4468e-03, 7.3244e-05, 4.7337e-06, 7.2785e-08, -3.7216e-10, -7.8868e-12   },
      {1.1589e+00, -1.0558e-03, 2.2047e-04, 1.6263e-06, -2.2203e-08, -3.0259e-10, -1.0585e-12  },
      {1.0691e+00, -2.2073e-03, 2.0158e-04, 2.4884e-06, 3.8384e-08, -2.0940e-10, -8.8724e-12   }
    },
    {
      {1.3214e+00, -4.1666e-04, 3.7428e-04, -1.1720e-07, -2.3338e-08, 3.5305e-11, -3.2754e-12  },
      {1.1866e+00, 7.2828e-03, 3.0029e-04, -9.4224e-06, 8.9523e-08, 1.1680e-09, -1.6491e-11    },
      {1.3759e+00, 3.0491e-03, 1.7077e-04, -2.5388e-06, 2.5061e-08, 2.4012e-10, -3.7566e-12    },
      {1.2639e+00, 2.2375e-03, 7.6136e-05, -1.7501e-06, 3.9701e-08, -1.3015e-10, -1.3353e-12   },
      {1.3175e+00, 1.0165e-03, 6.5074e-05, -1.3655e-06, 4.5156e-08, -2.1937e-10, -1.2701e-12   },
      {1.2192e+00, 7.6473e-03, 1.2184e-04, -8.6728e-06, 1.3123e-07, 7.5499e-10, -1.4887e-11    },
      {1.4400e+00, 5.3844e-04, 2.5964e-04, -1.5433e-06, 1.5967e-08, 6.5136e-11, -5.9000e-12    },
      {1.4992e+00, 3.2545e-04, 4.0814e-04, -3.1428e-07, -9.1290e-08, -1.1960e-10, 5.6280e-12   },
      {9.3750e-01, -9.1505e-03, 1.0130e-04, 8.7482e-06, 1.3768e-07, -8.4349e-10, -1.6221e-11   },
      {1.0283e+00, 6.4878e-04, 1.8831e-04, -1.9329e-06, -2.0171e-08, 8.8323e-10, 8.6919e-12    },
      {9.6642e-01, -1.9405e-03, 1.2465e-04, 5.3758e-07, 3.1163e-08, 4.6631e-10, 1.3538e-12     },
      {1.1660e+00, -3.2232e-03, 4.0286e-04, 3.7747e-06, -9.7224e-09, -5.6825e-10, -4.0148e-12  },
      {1.3214e+00, -4.1666e-04, 3.7428e-04, -1.1720e-07, -2.3338e-08, 3.5305e-11, -3.2754e-12  }
    }
  };
  static const greal AC[9] = {
      -2.8476225,  5.1614397,  11.808264,  1.5427184,  8.2490367,
      -6.7993216,  -7.6310943,  -7.1398243,  -11.384270
  };
  static const greal BC[9] = {
      4.7110391,  -6.9787011,  -19.723706,  -2.5946215,  -13.403260,
      6.8682539,  15.810774,  9.5193664,  19.050826
  };
  int bytes, ils, ihgt;
  double ylat;


  MIBIndices idx;
  for (idx.f107 = 0; idx.f107 < MTGCM_F107_SIZE; ++idx.f107) {
    string fileName = "";
    if (idx.f107 == 0) fileName = "zfhtls1.txt";
    if (idx.f107 == 1) fileName = "zfhtms1.txt";
    ifstream textFile(fileName.c_str());
    for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
      if (idx.f107 == 0) {
        if (idx.od == 0) fileName = "tpdls031.bin";
        if (idx.od == 1) fileName = "tpdls101.bin";
        if (idx.od == 2) fileName = "tpdls301.bin";
      }
      if (idx.f107 == 1) {
        if (idx.od == 0) fileName = "tpdms031.bin";
        if (idx.od == 1) fileName = "tpdms101.bin";
        if (idx.od == 2) fileName = "tpdms301.bin";
      }
      ifstream binaryFile[3];
      binaryFile[0].open(fileName.c_str(), ios::binary | ios::in);

      if (idx.f107 == 0) {
        if (idx.od == 0) fileName = "uvls031.bin";
        if (idx.od == 1) fileName = "uvls101.bin";
        if (idx.od == 2) fileName = "uvls301.bin";
      }
      if (idx.f107 == 1) {
        if (idx.od == 0) fileName = "uvms031.bin";
        if (idx.od == 1) fileName = "uvms101.bin";
        if (idx.od == 2) fileName = "uvms301.bin";
      }
      binaryFile[1].open(fileName.c_str(), ios::binary | ios::in);

      for (idx.ls = 1; idx.ls < LS_SIZE; ++idx.ls) {
        greal sls = sin(toRadians(30.0_deg * idx.ls));
        greal cls = cos(toRadians(30.0_deg * idx.ls));
        for (idx.lat = 0; idx.lat < MTGCM_LAT_SIZE; ++idx.lat) {
          greal xlat = 5.0 * idx.lat - 87.5_deg;
          greal zlat = clamp(xlat, -82.5_deg, 82.5_deg);
          greal FitLat = Hf80[idx.od][idx.ls][0]
            + zlat * (Hf80[idx.od][idx.ls][1] + zlat * (Hf80[idx.od][idx.ls][2]
              + zlat * (Hf80[idx.od][idx.ls][3] + zlat * (Hf80[idx.od][idx.ls][4]
                + zlat * (Hf80[idx.od][idx.ls][5] + zlat * Hf80[idx.od][idx.ls][6])))));
          if (FitLat < 0.5) FitLat = 0.5;

          for (idx.hgt = 0; idx.hgt < MTGCM_HEIGHT_SIZE; ++idx.hgt) {
            binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
            binaryFile[0].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[1][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[2][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[3][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmTemperature[4][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[1][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[2][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[3][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmPressure[4][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[1][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[2][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[3][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&mtgcmDensity[4][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[0].read(reinterpret_cast<char*>(&bytes), sizeof(int));

            greal zm80 = 5.0 * idx.hgt / 100.0;
            if (zm80 > 0.55) zm80 = 0.55;
            greal phi = zlat / 100.0;
            greal phi2 = phi * phi;
            greal Acalc = AC[0] + AC[1] * phi + AC[2] * phi2 + AC[3] * sls +
              AC[4] * cls + AC[5] * phi * sls + AC[6] * phi * cls +
              AC[7] * phi2 * sls + AC[8] * phi2 * cls;
            greal Bcalc = BC[0] + BC[1] * phi + BC[2] * phi2 + BC[3] * sls +
              BC[4] * cls + BC[5] * phi * sls + BC[6] * phi * cls +
              BC[7] * phi2 * sls + BC[8] * phi2 * cls;
            greal FitHgt = (1.0 + Acalc * zm80 + Bcalc * zm80 * zm80) * 0.9883;
            if (FitHgt < 0.5) FitHgt = 0.5;
            mtgcmDensity[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107] *= FitLat * FitHgt;
            mtgcmPressure[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107] *= FitLat * FitHgt;

            binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ils), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ihgt), sizeof(int));
            binaryFile[1].read(reinterpret_cast<char*>(&ylat), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[1][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[2][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[3][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmEWWind[4][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[0][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[1][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[2][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[3][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&mtgcmNSWind[4][idx.hgt][idx.lat][idx.ls][idx.od][idx.f107]), sizeof(double));
            binaryFile[1].read(reinterpret_cast<char*>(&bytes), sizeof(int));
          }
          greal zfdust;
          textFile >> zfdust >> ils >> ylat
            >> mtgcmThermosphereBaseHeight[0][idx.lat][idx.ls][idx.od][idx.f107]
            >> mtgcmThermosphereBaseHeight[1][idx.lat][idx.ls][idx.od][idx.f107]
            >> mtgcmThermosphereBaseHeight[2][idx.lat][idx.ls][idx.od][idx.f107]
            >> mtgcmThermosphereBaseHeight[3][idx.lat][idx.ls][idx.od][idx.f107]
            >> mtgcmThermosphereBaseHeight[4][idx.lat][idx.ls][idx.od][idx.f107];

          int k1 = int((mtgcmThermosphereBaseHeight[0][idx.lat][idx.ls][idx.od][idx.f107] - 75.0) / 5.0) - 1;
          greal ZFX = mtgcmThermosphereBaseHeight[0][idx.lat][idx.ls][idx.od][idx.f107];
          mtgcmThermosphereBaseHeight[0][idx.lat][idx.ls][idx.od][idx.f107] = ZFX + 5.0 * log(FitLat) /
            log(mtgcmDensity[0][k1][idx.lat][idx.ls][idx.od][idx.f107] / mtgcmDensity[0][k1 + 1][idx.lat][idx.ls][idx.od][idx.f107]);

        }
      }
      binaryFile[0].close();
      binaryFile[1].close();
    }
    textFile.close();
  }

  for (idx.f107 = 0; idx.f107 < MTGCM_F107_SIZE; ++idx.f107) {
    for (idx.od = 0; idx.od < OD_SIZE; ++idx.od) {
      for (idx.lat = 0; idx.lat < MTGCM_LAT_SIZE; ++idx.lat) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          for (idx.hgt = 0; idx.hgt < MTGCM_HEIGHT_SIZE; ++idx.hgt) {
            mtgcmTemperature[p][idx.hgt][idx.lat][0][idx.od][idx.f107] = mtgcmTemperature[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
            mtgcmPressure[p][idx.hgt][idx.lat][0][idx.od][idx.f107] = mtgcmPressure[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
            mtgcmDensity[p][idx.hgt][idx.lat][0][idx.od][idx.f107] = mtgcmDensity[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
            mtgcmEWWind[p][idx.hgt][idx.lat][0][idx.od][idx.f107] = mtgcmEWWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
            mtgcmNSWind[p][idx.hgt][idx.lat][0][idx.od][idx.f107] = mtgcmNSWind[p][idx.hgt][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
          }
          mtgcmThermosphereBaseHeight[p][idx.lat][0][idx.od][idx.f107] = mtgcmThermosphereBaseHeight[p][idx.lat][LS_SIZE - 1][idx.od][idx.f107];
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
void MGCMUpperInterpolator::writeUpperData() {
  ofstream binaryFile;
   binaryFile.open(dataPath + "MGCM_upper_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * OD_SIZE * MTGCM_F107_SIZE;
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