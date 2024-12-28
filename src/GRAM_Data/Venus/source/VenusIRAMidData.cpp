//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//
// Adapted from Venus-GRAM 2005 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include "VenusIRAMid.h"

using namespace std;

namespace GRAM {

bool VenusIRAMid::initialized = false;
const size_t VenusIRAMid::heightSize;
const size_t VenusIRAMid::lstSize;
std::vector<greal> VenusIRAMid::zmd;
std::vector<greal> VenusIRAMid::lstmd;
std::vector<std::vector<greal>> VenusIRAMid::dmd;
std::vector<std::vector<greal>> VenusIRAMid::pmd;
std::vector<std::vector<greal>> VenusIRAMid::tmd;
std::vector<std::vector<greal>> VenusIRAMid::yco2md;
std::vector<std::vector<greal>> VenusIRAMid::yn2md;
std::vector<std::vector<greal>> VenusIRAMid::yomd;
std::vector<std::vector<greal>> VenusIRAMid::ycomd;
std::vector<std::vector<greal>> VenusIRAMid::yhemd;
std::vector<std::vector<greal>> VenusIRAMid::ynmd;

typedef struct VMidData {
  double z;
  double lst;
  double ndco2;
  double ndn2;
  double ndo;
  double ndco;
  double ndhe;
  double ndn;
  double dens;
  double pres;
  double temp;
  double mmwt;
} VIRAMidData;

const VIRAMidData VenusIRAMid::data[heightSize * lstSize] = {
//    Zkm     LST    ndCO2m3    ndN2m3     ndOm3      ndCOm3     ndHem3      ndNm3     rhokgm3     PNm2        TK     MMWT
  { 100.0,  0.0,  1.07E+21,  3.96E+19,  2.18E+17,  1.08E+18,  1.00E+16,  1.72E+12,  8.03E-05,  2.615E+00,  170.0,  43.41 },
  { 105.0,  0.0,  2.94E+20,  1.14E+19,  2.55E+17,  8.06E+17,  4.23E+15,  1.56E+13,  2.21E-05,  6.614E-01,  156.0,  43.34 },
  { 110.0,  0.0,  7.22E+19,  3.04E+18,  2.79E+17,  4.27E+17,  1.64E+15,  1.27E+14,  5.44E-06,  1.446E-01,  138.0,  43.17 },
  { 115.0,  0.0,  1.44E+19,  6.40E+17,  2.56E+17,  1.67E+17,  5.55E+14,  6.99E+14,  1.09E-06,  2.780E-02,  131.0,  42.70 },
  { 120.0,  0.0,  2.60E+18,  1.35E+17,  1.33E+17,  5.69E+16,  1.87E+14,  6.65E+14,  2.02E-07,  5.200E-03,  129.0,  41.67 },
  { 125.0,  0.0,  4.57E+17,  2.98E+16,  4.67E+16,  1.90E+16,  7.61E+13,  2.53E+14,  3.69E-08,  9.768E-04,  128.0,  40.20 },
  { 130.0,  0.0,  7.94E+16,  7.34E+15,  1.55E+16,  6.25E+15,  3.77E+13,  8.61E+13,  6.85E-09,  1.914E-04,  127.6,  37.96 },
  { 135.0,  0.0,  1.37E+16,  1.96E+15,  6.30E+15,  2.06E+15,  2.50E+13,  3.18E+13,  1.36E-09,  4.250E-05,  127.5,  33.92 },
  { 140.0,  0.0,  2.38E+15,  5.64E+14,  3.04E+15,  6.76E+14,  2.23E+13,  1.77E+13,  3.13E-10,  1.179E-05,  127.4,  28.12 },
  { 145.0,  0.0,  4.11E+14,  1.80E+14,  1.61E+15,  2.21E+14,  2.21E+13,  1.01E+13,  9.19E-11,  4.320E-06,  127.4,  22.53 },
  { 150.0,  0.0,  7.10E+13,  5.91E+13,  8.51E+14,  7.24E+13,  1.89E+13,  5.80E+12,  3.42E-11,  1.898E-06,  127.4,  19.09 },
  { 100.0, 12.0,  1.08E+21,  4.61E+19,  1.27E+17,  5.42E+17,  7.81E+15,  4.58E+09,  8.14E-05,  2.733E+00,  175.0,  43.33 },
  { 105.0, 12.0,  2.99E+20,  1.29E+19,  3.17E+17,  6.98E+17,  2.29E+15,  2.51E+11,  2.25E-05,  7.565E-01,  175.0,  43.28 },
  { 110.0, 12.0,  8.21E+19,  3.57E+18,  3.45E+17,  6.67E+17,  6.81E+14,  6.86E+11,  6.21E-06,  2.132E-01,  178.0,  43.11 },
  { 115.0, 12.0,  2.31E+19,  1.06E+18,  2.32E+17,  3.29E+17,  2.20E+14,  4.49E+13,  1.76E-06,  6.217E-02,  182.0,  42.84 },
  { 120.0, 12.0,  6.62E+18,  3.23E+17,  1.20E+17,  1.43E+17,  8.19E+13,  1.30E+14,  5.08E-07,  1.878E-02,  189.0,  42.50 },
  { 125.0, 12.0,  1.98E+18,  1.04E+17,  4.42E+16,  6.68E+16,  3.65E+13,  2.00E+14,  1.54E-07,  6.009E-03,  198.0,  42.19 },
  { 130.0, 12.0,  6.26E+17,  3.55E+16,  2.76E+16,  3.26E+16,  1.85E+13,  1.72E+14,  4.97E-08,  2.056E-03,  206.0,  41.41 },
  { 135.0, 12.0,  2.07E+17,  1.34E+16,  1.54E+16,  1.64E+16,  1.06E+13,  1.24E+14,  1.69E-08,  7.481E-04,  215.0,  40.39 },
  { 140.0, 12.0,  7.17E+16,  5.54E+15,  9.13E+15,  8.42E+15,  7.03E+12,  9.02E+13,  6.13E-09,  2.946E-04,  225.0,  38.92 },
  { 145.0, 12.0,  2.59E+16,  2.48E+15,  5.83E+15,  4.39E+15,  5.65E+12,  6.49E+13,  2.37E-09,  1.261E-04,  236.0,  36.88 },
  { 150.0, 12.0,  9.81E+15,  1.32E+15,  4.01E+15,  2.34E+15,  5.02E+12,  4.65E+13,  9.94E-10,  5.963E-05,  246.5,  34.16 }

};

//! \brief Populate the middle atmosphere tables.
//!
//! The middle atmosphere data is stored in static structures (shared by all instances).
//! Initially, all data resides in a vector of VIRAMidData structs.  This data is relocated
//! and possibly modified into vectors dimensioned by height and solar time.
//!
//! \b Inputs
//! \arg #data The initial structure of the data.
//! \arg #initialized Flags when data has been initialized.
//!
//! \returns The static data tables are populated.
void VenusIRAMid::initializeData() {
  if (initialized) {
    return;
  }
  initialized = true;

  // Reorganize the data table into indexable containers
  // vectors are resizable, resize the first dimension
  lstmd.resize(lstSize);
  zmd.resize(heightSize);
  pmd.resize(lstSize);
  tmd.resize(lstSize);
  dmd.resize(lstSize);
  yco2md.resize(lstSize);
  yn2md.resize(lstSize);
  yomd.resize(lstSize);
  ycomd.resize(lstSize);
  yhemd.resize(lstSize);
  ynmd.resize(lstSize);

  size_t dataIndex = 0;
  // loop over lst
  for (size_t lstIndex = 0; lstIndex < lstSize; ++lstIndex) {
    // store the lst values
    lstmd[lstIndex] = data[dataIndex].lst;
    // resize the second dimension
    dmd[lstIndex].resize(heightSize);
    pmd[lstIndex].resize(heightSize);
    tmd[lstIndex].resize(heightSize);
    yco2md[lstIndex].resize(heightSize);
    yn2md[lstIndex].resize(heightSize);
    yomd[lstIndex].resize(heightSize);
    ycomd[lstIndex].resize(heightSize);
    yhemd[lstIndex].resize(heightSize);
    ynmd[lstIndex].resize(heightSize);
    // loop over heights
    for (size_t zIndex = 0; zIndex < heightSize; ++zIndex) {
      if (lstIndex == 0) {
        // only store the height data once
        zmd[zIndex] = data[dataIndex].z;
      }
      // store the data
      dmd[lstIndex][zIndex] = data[dataIndex].dens;
      pmd[lstIndex][zIndex] = data[dataIndex].pres;
      tmd[lstIndex][zIndex] = data[dataIndex].temp;
      yco2md[lstIndex][zIndex] = data[dataIndex].ndco2;
      yn2md[lstIndex][zIndex] = data[dataIndex].ndn2;
      yomd[lstIndex][zIndex] = data[dataIndex].ndo;
      ycomd[lstIndex][zIndex] = data[dataIndex].ndco;
      yhemd[lstIndex][zIndex] = data[dataIndex].ndhe;
      ynmd[lstIndex][zIndex] = data[dataIndex].ndn;
      // go to the next row of data
      ++dataIndex;
    }
  }
}

} // namespace