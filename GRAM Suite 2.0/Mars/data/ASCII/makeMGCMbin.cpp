#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

const double TO_RADIANS = (0.0174532925199432950);  // PI / 180

const size_t SURFACE_HEIGHT_SIZE = 3; // The number of height layers in the MGCM surface data (0, 0.005, 0.03 km).
const size_t SURFACE_LON_SIZE = 41;   // Tne number of longitude layers in the MGCM surface data (0 to 360 by 9 degrees).
const size_t MGCM_HEIGHT_SIZE = 17;   // The number of height layers in the MGCM data (0 to 80 by 5 km).
const size_t MGCM_LAT_SIZE = 25;      // The number of latitude layers in the MGCM data (-90 to 90 by 7.5 degrees).
const size_t MTGCM_HEIGHT_SIZE = 19;  // The number of height layers in the MTGCM data (80 to 170 by 5 km).
const size_t MTGCM_LAT_SIZE = 36;     // The number of latitude layers in the MTGCM data (-87.7 to 87.5 by 5 degrees).
const size_t MTGCM_F107_SIZE = 2;     // The number of F10.7 layers at 1AU in the MTGCM data (70 and 130).
const size_t PARAM_SIZE = 5;          // The number of tidal parameters.
const size_t LS_SIZE = 13;            // Longitude of the sun layers. From 0 to 360 by 30 degrees.
const size_t OD_SIZE = 3;             // The number of optical depth layers.
const static size_t MOLA_LAT_SIZE = 362;      //!< The size of the latitude dimension for the radius and topo height arrays.
const static size_t MOLA_LON_SIZE = 722;      //!< The size of the longitude dimension for the radius and topo height arrays.
const static size_t ALB_LAT_SIZE = 182;  //!< The size of the latitude dimension for the albedo array.
const static size_t ALB_LON_SIZE = 362;  //!< The size of the longitude dimension for the albedo array.
double areoRadius[MOLA_LAT_SIZE][MOLA_LON_SIZE] = {};
double topoHeight[MOLA_LAT_SIZE][MOLA_LON_SIZE] = {};
double albedo[ALB_LAT_SIZE][ALB_LON_SIZE] = {};

double surfaceTemperatures[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
double surfaceEWWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
double surfaceNSWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][OD_SIZE] = {};
double surfaceHeights[SURFACE_HEIGHT_SIZE] = { 0.0,  0.005,  0.03 };
double mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
double mgcmPressure[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
double mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
double mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][OD_SIZE] = {};
double mgcmDensity[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][3] = {};
double mtgcmTemperature[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmPressure[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmDensity[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmEWWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmNSWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmThermosphereBaseHeight[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][OD_SIZE][MTGCM_F107_SIZE] = {};
const string dust[OD_SIZE] = { "03", "10", "30" };   // "y1", "y2"};

void surfbin(char version);
void mgcmbin(char version);
void mtgcmbin(char version);
void molabin(char version);


// Makes binary files from ASCII versions of NASA Ames Mars Global
// Circulation Model (MGCM) and University of Michigan Mars       
// Thermospheric Global Circulation Model (MTGCM) data, for fast  
// reading on user system                                         
int main(int argc, char** argv)
{
  char version = '1';

  cout << " Reading MGCM surface data" << endl;
  surfbin(version);

  cout << " Reading MGCM -5 to 80 km data" << endl;
  mgcmbin(version);

  cout << " Reading MTGCM 80 to 240 km data" << endl;
  mtgcmbin(version);

  cout << " Reading MOLA data" << endl;
  molabin(version);

  return 0;
}

// Utility method for writing a block of binary data.
//
// binaryStream The binary data stream for output.
// dataSize     The number of doubles to be written.
// dataBlock    An array of doubles of size dataSize.
void writeBinaryDataBlock(ofstream& binaryStream, size_t dataSize, double* dataBlock)
{
  // Write the block size
  size_t byteSize = dataSize * sizeof(double);
  binaryStream.write(reinterpret_cast<char*>(&byteSize), sizeof(size_t));

  // Write the data block.
  binaryStream.write(reinterpret_cast<char*>(dataBlock), byteSize);
}

// Reads ASCII version MGCM surface data and writes binary version
void surfbin(char version) {
  // Step through all dust optical depths
  for (size_t od = 0; od < OD_SIZE; ++od) {
    // Open ASCII version and binary version data files at 5m and 30m levels
    ifstream level_0("sfc00" + dust[od] + version + ".txt");
    ifstream level_5("sfc05" + dust[od] + version + ".txt");
    ifstream level_30("sfc30" + dust[od] + version + ".txt");

    // Read (and ignore) the header line from the ASCII files
    string dummy;
    getline(level_0, dummy);
    getline(level_5, dummy);
    getline(level_30, dummy);

    // Step through all Ls values  (360 will be copied to 0 later)
    for (size_t ls = 1; ls < LS_SIZE; ls++) {
      int xls = ls * 30;
      // Step through all latitudes
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        double xlat = -90.0 + lat * 7.5;
        // Step through all longitudes
        for (size_t lon = SURFACE_LON_SIZE - 1; lon > 0; --lon) {
          int xlon = lon * 9;
          // Read ASCII version tide coefficient amplitudes and phases
          // for surface temperature
          int ils, ilon;
          double ilat, tp[PARAM_SIZE], up[PARAM_SIZE], vp[PARAM_SIZE];
          level_0 >> ils >> ilat >> ilon >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4];
          if (ils != xls) cout << " Bad surface Ls" << endl;
          if (ilat != xlat) cout << " Bad surface Latitude" << endl;
          if (ilon != xlon) cout << " Bad surface Longitude" << endl;
          // Store tide amplitudes and phases for
          // temperature, EW wind, and NS wind
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            surfaceTemperatures[p][0][lat][lon][ls][od] = tp[p];
            surfaceEWWinds[p][0][lat][lon][ls][od] = 0.0;
            surfaceNSWinds[p][0][lat][lon][ls][od] = 0.0;
          }
          // Read ASCII version tide coefficient amplitudes and phases
          // for temperature, EW wind, and NS wind
          level_5 >> ils >> ilat >> ilon
            >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4]
            >> up[0] >> up[1] >> up[2] >> up[3] >> up[4]
            >> vp[0] >> vp[1] >> vp[2] >> vp[3] >> vp[4];
          if (ils != xls) cout << " Bad surface Ls" << endl;
          if (ilat != xlat) cout << " Bad surface Latitude" << endl;
          if (ilon != xlon) cout << " Bad surface Longitude" << endl;
          // Store tide amplitudes and phases for
          // temperature, EW wind, and NS wind
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            surfaceTemperatures[p][1][lat][lon][ls][od] = tp[p];
            surfaceEWWinds[p][1][lat][lon][ls][od] = up[p];
            surfaceNSWinds[p][1][lat][lon][ls][od] = vp[p];
          }
          // Read ASCII version tide coefficient amplitudes and phases
          // for temperature, EW wind, and NS wind
          level_30 >> ils >> ilat >> ilon
            >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4]
            >> up[0] >> up[1] >> up[2] >> up[3] >> up[4]
            >> vp[0] >> vp[1] >> vp[2] >> vp[3] >> vp[4];
          if (ils != xls) cout << " Bad surface Ls" << endl;
          if (ilat != xlat) cout << " Bad surface Latitude" << endl;
          if (ilon != xlon) cout << " Bad surface Longitude" << endl;
          // Store tide amplitudes and phases for
          // temperature, EW wind, and NS wind
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            surfaceTemperatures[p][2][lat][lon][ls][od] = tp[p];
            surfaceEWWinds[p][2][lat][lon][ls][od] = up[p];
            surfaceNSWinds[p][2][lat][lon][ls][od] = vp[p];
          }
        }
      }
    }
    // Close all MGCM surface data input files
    level_0.close();
    level_5.close();
    level_30.close();
  }

  // Copy LS = 360 (at index LS_SIZE - 1) into LS = 0 (at index 0)
  // Copy LON = 360 (at index SURFACE_LON_SIZE - 1) into LON = 0 (at index 0)
  for (size_t p = 0; p < PARAM_SIZE; ++p) {
    for (size_t od = 0; od < OD_SIZE; ++od) {
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        for (size_t hgt = 0; hgt < SURFACE_HEIGHT_SIZE; ++hgt) {
          for (size_t lon = 1; lon < SURFACE_LON_SIZE; ++lon) {
            surfaceTemperatures[p][hgt][lat][lon][0][od] = surfaceTemperatures[p][hgt][lat][lon][LS_SIZE - 1][od];
            surfaceEWWinds[p][hgt][lat][lon][0][od] = surfaceEWWinds[p][hgt][lat][lon][LS_SIZE - 1][od];
            surfaceNSWinds[p][hgt][lat][lon][0][od] = surfaceNSWinds[p][hgt][lat][lon][LS_SIZE - 1][od];
            if (lon == SURFACE_LON_SIZE - 1) {
              for (size_t ls = 0; ls < LS_SIZE; ++ls) {
                surfaceTemperatures[p][hgt][lat][0][ls][od] = surfaceTemperatures[p][hgt][lat][lon][ls][od];
                surfaceEWWinds[p][hgt][lat][0][ls][od] = surfaceEWWinds[p][hgt][lat][lon][ls][od];
                surfaceNSWinds[p][hgt][lat][0][ls][od] = surfaceNSWinds[p][hgt][lat][lon][ls][od];
              }
            }
          }
        }
      }
    }
  }

  cout << " Writing binary file MGCM_surface_data.bin." << endl;
  ofstream binaryFile;
  binaryFile.open("MGCM_surface_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * OD_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceTemperatures);
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceEWWinds);
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceNSWinds);

  binaryFile.close();
}

// Reads ASCII version MGCM 0 - 80 km data and writes binary version
void mgcmbin(char version) {
  static const double DFA[8][OD_SIZE] = {
    { 1.03050e+00,   1.01386e+00,  1.2029e+0},
    {-6.81901e-04,  -1.52668e-04,  0.0000e+0},
    { 6.67896e-06,   5.78641e-06,  4.5266e-5},
    { 1.69052e-02,   5.57000e-03,  0.0000e+0},
    { 2.83632e-04,  -8.23194e-05,  0.0000e+0},
    { 4.11891e-06,  -7.11778e-06,  0.0000e+0},
    {-1.05497e-01,   0.00000e+0,   0.0000e+0},
    { 9.17778e-05,   0.00000e+0,   0.0000e+0}
  };
  static const double DFB[8][OD_SIZE] = {
    {-3.91632e+00,  -1.63780e+00,  18.1260e+0},
    { 1.02837e-01,   6.76088e-02,   0.0000e+0},
    {-2.73027e-04,  -3.23646e-05,  -9.5731e-3},
    { 1.02129e+00,   8.20808e-02,   0.0000e+0},
    {-3.11446e-02,  -4.21906e-02,   0.0000e+0},
    {-5.95152e-04,   1.97186e-04,   0.0000e+0},
    { 8.26562e+00,   0.00000e+0,    0.0000e+0},
    {-3.37582e-02,   0.00000e+0,    0.0000e+0}
  };
  static const double DFC[3][OD_SIZE] = {
    {-0.1307e0,  -0.1657e0,  -0.3931e0},
    { 1.6074e0,   1.6981e0,   2.9640e0},
    {-3.7984e0,  -4.1688e0,  -5.6292e0}
  };

  // Step through all dust optical depths
  for (size_t od = 0; od < OD_SIZE; ++od) {
    // Open ASCII version and binary version data files at 5m and 30m levels
    ifstream tpdlo("tpdlo" + dust[od] + version + ".txt");
    ifstream uvlo("uvlo" + dust[od] + version + ".txt");

    // Read (and ignore) the header line from the ASCII files
    string dummy;
    getline(tpdlo, dummy);
    getline(uvlo, dummy);

    // Step through all Ls values  (360 will be copied to 0 later)
    for (size_t ls = 1; ls < LS_SIZE; ls++) {
      int xls = ls * 30;
      double sinLs = sin(TO_RADIANS * xls);
      // Step through all latitudes
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        double xlat = -90.0 + lat * 7.5;
        // Step through all heights
        for (int hgt = MGCM_HEIGHT_SIZE - 1; hgt >= 0; --hgt) {
          int xhgt = hgt * 5;
          // Read ASCII version tide coefficient amplitudes and phases
          // for temperature, pressure, and density.
          int ils, ihgt;
          double ilat, tp[PARAM_SIZE], pp[PARAM_SIZE], da0;
          tpdlo >> ils >> ihgt >> ilat
            >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4]
            >> pp[0] >> pp[1] >> pp[2] >> pp[3] >> pp[4]
            >> da0;
          if (ils != xls) cout << " Bad Ls" << endl;
          if (ilat != xlat) cout << " Bad Latitude " << ilat << " " << xlat << endl;
          if (ihgt != xhgt) cout << " Bad hgt " << ihgt << " " << xhgt << endl;
          // Adjust pressure and density 
          double dfa = DFA[0][od] + DFA[1][od] * xlat + DFA[2][od] * xlat * xlat
            + (DFA[3][od] + DFA[4][od] * xlat + DFA[5][od] * xlat * xlat) * sinLs
            + (DFA[6][od] + DFA[7][od] * xlat) * sinLs * sinLs;
          double dfb = DFB[0][od] + DFB[1][od] * xlat + DFB[2][od] * xlat * xlat
            + (DFB[3][od] + DFB[4][od] * xlat + DFB[5][od] * xlat * xlat) * sinLs
            + (DFB[6][od] + DFB[7][od] * xlat) * sinLs * sinLs;
          double dfm = dfa + dfb * da0;
          double ZZ = double(hgt) / 20.0;
          double dfc = 1.0 + DFC[0][od] * ZZ + DFC[1][od] * ZZ * ZZ + DFC[2][od] * ZZ * ZZ * ZZ;
          dfc = max(0.5, dfc);

          // Store tide amplitudes and phases 
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            mgcmTemperature[p][hgt][lat][ls][od] = tp[p];
            mgcmPressure[p][hgt][lat][ls][od] = pp[p];
          }
          mgcmDensity[hgt][lat][ls][od] = da0 * dfm / dfc;
          mgcmPressure[0][hgt][lat][ls][od] *= dfm / dfc;

          // Read ASCII version tide coefficient amplitudes and phases
          // for EW wind, and NS wind
          double up[PARAM_SIZE], vp[PARAM_SIZE];
          uvlo >> ils >> ihgt >> ilat
            >> up[0] >> up[1] >> up[2] >> up[3] >> up[4]
            >> vp[0] >> vp[1] >> vp[2] >> vp[3] >> vp[4];
          if (ils != xls) cout << " Bad Ls" << endl;
          if (ilat != xlat) cout << " Bad Latitude" << endl;
          if (ihgt != xhgt) cout << " Bad hgt" << endl;
          // Store tide amplitudes and phases
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            mgcmEWWind[p][hgt][lat][ls][od] = up[p];
            mgcmNSWind[p][hgt][lat][ls][od] = vp[p];
          }
        }
      }
    }
    // Close all input files
    tpdlo.close();
    uvlo.close();

    // Copy LS = 360 (at index LS_SIZE - 1) into LS = 0 (at index 0)
    for (size_t hgt = 0; hgt < MGCM_HEIGHT_SIZE; ++hgt) {
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        for (size_t od = 0; od < OD_SIZE; ++od) {
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            mgcmTemperature[p][hgt][lat][0][od] = mgcmTemperature[p][hgt][lat][LS_SIZE - 1][od];
            mgcmPressure[p][hgt][lat][0][od] = mgcmPressure[p][hgt][lat][LS_SIZE - 1][od];
            mgcmEWWind[p][hgt][lat][0][od] = mgcmEWWind[p][hgt][lat][LS_SIZE - 1][od];
            mgcmNSWind[p][hgt][lat][0][od] = mgcmNSWind[p][hgt][lat][LS_SIZE - 1][od];
          }
          mgcmDensity[hgt][lat][0][od] = mgcmDensity[hgt][lat][LS_SIZE - 1][od];
        }
      }
    }
  }

    cout << " Writing binary file MGCM_lower_data.bin." <<  endl;
    ofstream binaryFile;
    binaryFile.open("MGCM_lower_data.bin", ios::binary | ios::out);

    size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * OD_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmDensity);

    dataSize *= PARAM_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmTemperature);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmPressure);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmEWWind);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmNSWind);

    binaryFile.close();
  }

// Reads ASCII version MTGCM 80-240 km data and writes binary
void mtgcmbin(char version) {
  static const double Hf80[OD_SIZE][LS_SIZE][7] = {
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
  static const double AC[9] = {
      -2.8476225,  5.1614397,  11.808264,  1.5427184,  8.2490367,
      -6.7993216,  -7.6310943,  -7.1398243,  -11.384270
  };
  static const double BC[9] = {
      4.7110391,  -6.9787011,  -19.723706,  -2.5946215,  -13.403260,
      6.8682539,  15.810774,  9.5193664,  19.050826
  };
  static const string solact[2] = { "ls", "ms" };
  // Step through all solar activity levels
  for (size_t f107 = 0; f107 < MTGCM_F107_SIZE; ++f107) {
    ifstream baseHeight("zfht" + solact[f107] + version + ".txt");
    // Step through all dust optical depths
    for (size_t od = 0; od < OD_SIZE; ++od) {
      // Open ASCII version and binary version data files at 5m and 30m levels
      ifstream tpd("tpd" + solact[f107] + dust[od] + version + ".txt");
      ifstream uv("uv" + solact[f107] + dust[od] + version + ".txt");

      // Read (and ignore) the header line from the ASCII files
      string dummy;
      getline(tpd, dummy);
      getline(uv, dummy);

      // Step through all Ls values  (360 will be copied to 0 later)
      for (size_t ls = 1; ls < LS_SIZE; ls++) {
        int xls = ls * 30;
        double sinLs = sin(TO_RADIANS * xls);
        double cosLs = cos(TO_RADIANS * xls);
        // Step through all latitudes
        for (size_t lat = 0; lat < MTGCM_LAT_SIZE; ++lat) {
          double xlat = -87.5 + lat * 5.0;
          double zlat = min(max(xlat, -82.5), 82.5);
          double FitLat = Hf80[od][ls][0]
            + zlat * (Hf80[od][ls][1] + zlat * (Hf80[od][ls][2]
              + zlat * (Hf80[od][ls][3] + zlat * (Hf80[od][ls][4]
                + zlat * (Hf80[od][ls][5] + zlat * Hf80[od][ls][6])))));
          if (FitLat < 0.5) FitLat = 0.5;
          // Step through all heights
          for (size_t hgt = 0; hgt < MTGCM_HEIGHT_SIZE; ++hgt) {
            int xhgt = 80 + hgt * 5;
            // Read ASCII version tide coefficient amplitudes and phases
            // for temperature, pressure and density
            int ils, jhgt;
            double ylat, tp[PARAM_SIZE], pp[PARAM_SIZE], dp[PARAM_SIZE];
            tpd >> ils >> jhgt >> ylat
              >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4]
              >> pp[0] >> pp[1] >> pp[2] >> pp[3] >> pp[4]
              >> dp[0] >> dp[1] >> dp[2] >> dp[3] >> dp[4];
            if (ils != xls) cout << " Bad Ls" << endl;
            if (ylat != xlat) cout << " Bad Latitude" << endl;
            if (jhgt != xhgt) cout << " Bad hgt" << endl;
            // Store tide amplitudes and phases 
            for (size_t p = 0; p < PARAM_SIZE; ++p) {
              mtgcmTemperature[p][hgt][lat][ls][od][f107] = tp[p];
              mtgcmPressure[p][hgt][lat][ls][od][f107] = pp[p];
              mtgcmDensity[p][hgt][lat][ls][od][f107] = dp[p];
            }
            double zm80 = 5.0 * hgt / 100.0;
            if (zm80 > 0.55) zm80 = 0.55;
            double phi = zlat / 100.0;
            double phi2 = phi * phi;
            double Acalc = AC[0] + AC[1] * phi + AC[2] * phi2 + AC[3] * sinLs +
              AC[4] * cosLs + AC[5] * phi * sinLs + AC[6] * phi * cosLs +
              AC[7] * phi2 * sinLs + AC[8] * phi2 * cosLs;
            double Bcalc = BC[0] + BC[1] * phi + BC[2] * phi2 + BC[3] * sinLs +
              BC[4] * cosLs + BC[5] * phi * sinLs + BC[6] * phi * cosLs +
              BC[7] * phi2 * sinLs + BC[8] * phi2 * cosLs;
            double FitHgt = (1.0 + Acalc * zm80 + Bcalc * zm80 * zm80) * 0.9883;
            if (FitHgt < 0.5) FitHgt = 0.5;
            mtgcmDensity[0][hgt][lat][ls][od][f107] *= FitLat * FitHgt;
            mtgcmPressure[0][hgt][lat][ls][od][f107] *= FitLat * FitHgt;

            // Read ASCII version tide coefficient amplitudes and phases
            // for temperature, EW wind, and NS wind
            double up[PARAM_SIZE], vp[PARAM_SIZE];
            uv >> ils >> jhgt >> ylat
              >> up[0] >> up[1] >> up[2] >> up[3] >> up[4]
              >> vp[0] >> vp[1] >> vp[2] >> vp[3] >> vp[4];
            if (ils != xls) cout << " Bad Ls" << endl;
            if (ylat != xlat) cout << " Bad Latitude" << endl;
            if (jhgt != xhgt) cout << " Bad hgt" << endl;
            // Store tide amplitudes and phases 
            for (size_t p = 0; p < PARAM_SIZE; ++p) {
              mtgcmEWWind[p][hgt][lat][ls][od][f107] = up[p];
              mtgcmNSWind[p][hgt][lat][ls][od][f107] = vp[p];
            }
          }

          // Read thermosphere base heights
          double zfdust, ilat;
          int ils;
          baseHeight >> zfdust >> ils >> ilat
            >> mtgcmThermosphereBaseHeight[0][lat][ls][od][f107]
            >> mtgcmThermosphereBaseHeight[1][lat][ls][od][f107]
            >> mtgcmThermosphereBaseHeight[2][lat][ls][od][f107]
            >> mtgcmThermosphereBaseHeight[3][lat][ls][od][f107]
            >> mtgcmThermosphereBaseHeight[4][lat][ls][od][f107];

          int k1 = int((mtgcmThermosphereBaseHeight[0][lat][ls][od][f107] - 75.0) / 5.0) - 1;
          double ZFX = mtgcmThermosphereBaseHeight[0][lat][ls][od][f107];
          mtgcmThermosphereBaseHeight[0][lat][ls][od][f107] = ZFX + 5.0 * log(FitLat) /
            log(mtgcmDensity[0][k1][lat][ls][od][f107] / mtgcmDensity[0][k1 + 1][lat][ls][od][f107]);
        }
      }
      // Close all MGCM surface data input files
      tpd.close();
      uv.close();
    }
    baseHeight.close();
  }

  for (size_t f107 = 0; f107 < MTGCM_F107_SIZE; ++f107) {
    for (size_t od = 0; od < OD_SIZE; ++od) {
      for (size_t lat = 0; lat < MTGCM_LAT_SIZE; ++lat) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          for (size_t hgt = 0; hgt < MTGCM_HEIGHT_SIZE; ++hgt) {
            mtgcmTemperature[p][hgt][lat][0][od][f107] = mtgcmTemperature[p][hgt][lat][LS_SIZE - 1][od][f107];
            mtgcmPressure[p][hgt][lat][0][od][f107] = mtgcmPressure[p][hgt][lat][LS_SIZE - 1][od][f107];
            mtgcmDensity[p][hgt][lat][0][od][f107] = mtgcmDensity[p][hgt][lat][LS_SIZE - 1][od][f107];
            mtgcmEWWind[p][hgt][lat][0][od][f107] = mtgcmEWWind[p][hgt][lat][LS_SIZE - 1][od][f107];
            mtgcmNSWind[p][hgt][lat][0][od][f107] = mtgcmNSWind[p][hgt][lat][LS_SIZE - 1][od][f107];
          }
          mtgcmThermosphereBaseHeight[p][lat][0][od][f107] = mtgcmThermosphereBaseHeight[p][lat][LS_SIZE - 1][od][f107];
        }
      }
    }
  }

    cout << " Writing binary file MGCM_upper_data.bin." << endl;
    ofstream binaryFile;
    binaryFile.open("MGCM_upper_data.bin", ios::binary | ios::out);

    size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * OD_SIZE * MTGCM_F107_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmThermosphereBaseHeight);

    dataSize *= MTGCM_HEIGHT_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmTemperature);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmPressure);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmDensity);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmEWWind);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmNSWind);

    binaryFile.close();
  }

  void molabin(char version) {
    // Open ASCII version of the MOLA data
    ifstream mola("MOLATOPH.TXT");
    ifstream albedo1("albedo1.txt");

    // Read (and ignore) the header line from the ASCII files
    string dummy;
    getline(mola, dummy);
    getline(albedo1, dummy);


        // Step through all latitudes
        for (size_t lat = MOLA_LAT_SIZE - 2; lat > 0; --lat) {
          double xlat = -89.75 + (lat-1) * 0.5;
          // Step through all longitudes
          for (size_t lon = MOLA_LON_SIZE - 2; lon > 0; --lon) {
            double xlon = 360.0 - (0.25 + (lon-1) * 0.5);
            // Read ASCII version tide coefficient amplitudes and phases
            // for surface temperature
            int num;
            double ilat, ilon, rad, areoid, topo;
            mola >> ilon >> ilat >> rad >> areoid >> topo >> num;
            if (ilat != xlat) cout << " Bad Latitude" << endl;
            if (ilon != xlon) cout << " Bad Longitude" << endl;
            // Store the data
            areoRadius[lat][lon] = areoid * 1.0e-3;
            topoHeight[lat][lon] = topo * 1.0e-3;
            if (lon == 1) {
              areoRadius[lat][MOLA_LON_SIZE - 1] = areoid * 1.0e-3;
              topoHeight[lat][MOLA_LON_SIZE - 1] = topo * 1.0e-3;
            }
            else if (lon == MOLA_LON_SIZE - 2) {
              areoRadius[lat][0] = areoid * 1.0e-3;
              topoHeight[lat][0] = topo * 1.0e-3;
            }
          }
        }
        double aeroNorth = 0.0;
        double aeroSouth = 0.0;
        double topoNorth = 0.0;
        double topoSouth = 0.0;
        for (size_t lon = 1; lon < MOLA_LON_SIZE - 1; ++lon) {
          aeroNorth += areoRadius[MOLA_LAT_SIZE - 2][lon];
          aeroSouth += areoRadius[1][lon];
          topoNorth += topoHeight[MOLA_LAT_SIZE - 2][lon];
          topoSouth += topoHeight[1][lon];
        }
        aeroNorth /= MOLA_LON_SIZE - 2;
        aeroSouth /= MOLA_LON_SIZE - 2;
        topoNorth /= MOLA_LON_SIZE - 2;
        topoSouth /= MOLA_LON_SIZE - 2;
        for (size_t lon = 0; lon < MOLA_LON_SIZE; ++lon) {
          areoRadius[0][lon] = aeroSouth;
          topoHeight[0][lon] = topoSouth;
          areoRadius[MOLA_LAT_SIZE - 1][lon] = aeroNorth;
          topoHeight[MOLA_LAT_SIZE - 1][lon] = topoNorth;
        }

        // Step through all latitudes
        for (size_t lat = 1; lat < ALB_LAT_SIZE - 1; ++lat) {
          double xlat = -89.5 + (lat - 1);
          // Step through all longitudes
          for (size_t lon = 1; lon < ALB_LON_SIZE - 1; ++lon) {
            double xlon = 0.5 + (lon - 1);
            // Read ASCII version tide coefficient amplitudes and phases
            // for surface temperature
            double ilat, ilon, alb;
            albedo1 >> ilat >> ilon >> alb;
            if (ilat != xlat) cout << " Bad Latitude" << endl;
            if (ilon != xlon) cout << " Bad Longitude" << endl;
            // Store the data
            albedo[lat][lon] = alb;
            if (lon == ALB_LON_SIZE - 2) {
              albedo[lat][ALB_LON_SIZE - 1] = alb;
            }
            else if (lon == 1) {
              albedo[lat][0] = alb;
            }
          }
        }
        double albNorth = 0.0;
        double albSouth = 0.0;
        for (size_t lon = 1; lon < ALB_LON_SIZE - 1; ++lon) {
          albNorth += albedo[ALB_LAT_SIZE - 2][lon];
          albSouth += albedo[1][lon];
        }
        albNorth /= ALB_LON_SIZE - 2;
        albSouth /= ALB_LON_SIZE - 2;
        for (size_t lon = 0; lon < ALB_LON_SIZE; ++lon) {
          albedo[ALB_LAT_SIZE - 1][lon] = albNorth;
          albedo[0][lon] = albSouth;
        }

      // Close all MGCM surface data input files
      mola.close();
      albedo1.close();


    cout << " Writing binary file MOLA_data.bin." << endl;
    ofstream binaryFile;
    binaryFile.open("MOLA_data.bin", ios::binary | ios::out);

    size_t dataSize = MOLA_LAT_SIZE * MOLA_LON_SIZE * sizeof(double);
    binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
    binaryFile.write(reinterpret_cast<char*>(&areoRadius), dataSize);

    binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
    binaryFile.write(reinterpret_cast<char*>(&topoHeight), dataSize);

    dataSize = ALB_LAT_SIZE * ALB_LON_SIZE * sizeof(double);
    binaryFile.write(reinterpret_cast<char*>(&dataSize), sizeof(size_t));
    binaryFile.write(reinterpret_cast<char*>(&albedo), dataSize);

    binaryFile.close();
  }
