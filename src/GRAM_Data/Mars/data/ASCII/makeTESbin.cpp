#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

const double TO_RADIANS = (0.0174532925199432950);  // PI / 180

const size_t SURFACE_HEIGHT_SIZE = 3; // The number of height layers in the MGCM surface data (0, 0.005, 0.03 km).
const size_t SURFACE_LON_SIZE = 41;   // Tne number of longitude layers in the MGCM surface data (0 to 360 by 9 degrees).
const size_t MGCM_HEIGHT_SIZE = 30;   // The number of height layers in the MGCM data (-5 to 10 by 1km, 10 to 80 by 5 km).
const size_t MGCM_LAT_SIZE = 25;      // The number of latitude layers in the MGCM data (-90 to 90 by 7.5 degrees).
const size_t MTGCM_HEIGHT_SIZE = 33;  // The number of height layers in the MTGCM data (80 to 240 by 5 km).
const size_t MTGCM_LAT_SIZE = 36;     // The number of latitude layers in the MTGCM data (-87.7 to 87.5 by 5 degrees).
const size_t MTGCM_F107_SIZE = 3;     // The number of F10.7 layers at 1AU in the MTGCM data (70, 130 and 200).
const size_t PARAM_SIZE = 5;          // The number of tidal parameters.
const size_t LS_SIZE = 13;            // Longitude of the sun layers. From 0 to 360 by 30 degrees.
const size_t TES_YEAR_SIZE = 2;       // The number of TES years (1 and 2).
const size_t DUST_LON_SIZE = 40;
const size_t DUST_LAT_SIZE = 25;
const size_t DUST_LS_SIZE = 72;
double tes_tau[TES_YEAR_SIZE][DUST_LS_SIZE][DUST_LAT_SIZE][DUST_LON_SIZE];  //!< TES dust optical depth data.
double surfaceTemperatures[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double surfaceEWWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double surfaceNSWinds[PARAM_SIZE][SURFACE_HEIGHT_SIZE][MGCM_LAT_SIZE][SURFACE_LON_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double surfaceHeights[SURFACE_HEIGHT_SIZE] = { 0.0,  0.005,  0.03 };
double mgcmTemperature[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double mgcmPressure[MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double mgcmEWWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double mgcmNSWind[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double mgcmDensity[PARAM_SIZE][MGCM_HEIGHT_SIZE][MGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE] = {};
double mtgcmTemperature[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmPressure[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmDensity[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmEWWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmNSWind[PARAM_SIZE][MTGCM_HEIGHT_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
double mtgcmThermosphereBaseHeight[PARAM_SIZE][MTGCM_LAT_SIZE][LS_SIZE][TES_YEAR_SIZE][MTGCM_F107_SIZE] = {};
const string year[TES_YEAR_SIZE] = { "y1", "y2"};

void surfbin(char version);
void mgcmbin(char version);
void mtgcmbin(char version);
void dustbin(char version);


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

  cout << " Reading TES dust data" << endl;
  dustbin(version);

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
  for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
    // Open ASCII version and binary version data files at 5m and 30m levels
    ifstream level_0("sfc00" + year[tyr] + version + ".txt");
    ifstream level_5("sfc05" + year[tyr] + version + ".txt");
    ifstream level_30("sfc30" + year[tyr] + version + ".txt");

    // Read (and ignore) the header line from the ASCII files
    string dummy;
    getline(level_0, dummy);
    getline(level_5, dummy);
    getline(level_30, dummy);

    // Step through all Ls values  (360 will be copied to 0 later)
    for (size_t ls = 1; ls < LS_SIZE; ls++) {
      int xls = (int)ls * 30;
      // Step through all latitudes
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        double xlat = -90.0 + lat * 7.5;
        // Step through all longitudes
        for (size_t lon = SURFACE_LON_SIZE - 1; lon > 0; --lon) {
          int xlon = (int)lon * 9;
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
            surfaceTemperatures[p][0][lat][lon][ls][tyr] = tp[p];
            surfaceEWWinds[p][0][lat][lon][ls][tyr] = 0.0;
            surfaceNSWinds[p][0][lat][lon][ls][tyr] = 0.0;
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
            surfaceTemperatures[p][1][lat][lon][ls][tyr] = tp[p];
            surfaceEWWinds[p][1][lat][lon][ls][tyr] = up[p];
            surfaceNSWinds[p][1][lat][lon][ls][tyr] = vp[p];
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
            surfaceTemperatures[p][2][lat][lon][ls][tyr] = tp[p];
            surfaceEWWinds[p][2][lat][lon][ls][tyr] = up[p];
            surfaceNSWinds[p][2][lat][lon][ls][tyr] = vp[p];
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
    for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        for (size_t hgt = 0; hgt < SURFACE_HEIGHT_SIZE; ++hgt) {
          for (size_t lon = 1; lon < SURFACE_LON_SIZE; ++lon) {
            surfaceTemperatures[p][hgt][lat][lon][0][tyr] = surfaceTemperatures[p][hgt][lat][lon][LS_SIZE - 1][tyr];
            surfaceEWWinds[p][hgt][lat][lon][0][tyr] = surfaceEWWinds[p][hgt][lat][lon][LS_SIZE - 1][tyr];
            surfaceNSWinds[p][hgt][lat][lon][0][tyr] = surfaceNSWinds[p][hgt][lat][lon][LS_SIZE - 1][tyr];
            if (lon == SURFACE_LON_SIZE - 1) {
              for (size_t ls = 0; ls < LS_SIZE; ++ls) {
                surfaceTemperatures[p][hgt][lat][0][ls][tyr] = surfaceTemperatures[p][hgt][lat][lon][ls][tyr];
                surfaceEWWinds[p][hgt][lat][0][ls][tyr] = surfaceEWWinds[p][hgt][lat][lon][ls][tyr];
                surfaceNSWinds[p][hgt][lat][0][ls][tyr] = surfaceNSWinds[p][hgt][lat][lon][ls][tyr];
              }
            }
          }
        }
      }
    }
  }

  cout << " Writing binary file TES_surface_data.bin." << endl;
  ofstream binaryFile;
  binaryFile.open("TES_surface_data.bin", ios::binary | ios::out);

  size_t dataSize = PARAM_SIZE * SURFACE_HEIGHT_SIZE * MGCM_LAT_SIZE * SURFACE_LON_SIZE * LS_SIZE * TES_YEAR_SIZE;
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceTemperatures);
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceEWWinds);
  writeBinaryDataBlock(binaryFile, dataSize, (double*)surfaceNSWinds);

  binaryFile.close();
}

// Reads ASCII version MGCM 0 - 80 km data and writes binary version
void mgcmbin(char version) {
  // Step through all TES years
  for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
    // Open ASCII version and binary version data files at 5m and 30m levels
    ifstream tpdlo("tpdlo" + year[tyr] + version + ".txt");
    ifstream uvlo("uvlo" + year[tyr] + version + ".txt");

    // Read (and ignore) the header line from the ASCII files
    string dummy;
    getline(tpdlo, dummy);
    getline(uvlo, dummy);

    // Step through all Ls values  (360 will be copied to 0 later)
    for (size_t ls = 1; ls < LS_SIZE; ls++) {
      int xls = (int)ls * 30;
      // Step through all latitudes
      for (size_t lat = 0; lat < MGCM_LAT_SIZE; ++lat) {
        double xlat = -90.0 + lat * 7.5;
        // Step through all heights
        for (int hgt = MGCM_HEIGHT_SIZE - 1; hgt >= 0; --hgt) {
          int xhgt = hgt - 5;
          if (hgt > 15) xhgt = (hgt - 13) * 5;
          // Read ASCII version tide coefficient amplitudes and phases
          // for temperature, pressure, and density.
          int ils, ihgt;
          double ilat, tp[PARAM_SIZE], dp[PARAM_SIZE], pa0;
          tpdlo >> ils >> ihgt >> ilat
            >> tp[0] >> tp[1] >> tp[2] >> tp[3] >> tp[4]
            >> dp[0] >> dp[1] >> dp[2] >> dp[3] >> dp[4]
            >> pa0;
          if (ils != xls) cout << " Bad Ls" << endl;
          if (ilat != xlat) cout << " Bad Latitude" << endl;
          if (ihgt != xhgt) cout << " Bad hgt" << endl;
          // Store tide amplitudes and phases 
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            mgcmTemperature[p][hgt][lat][ls][tyr] = tp[p];
            mgcmDensity[p][hgt][lat][ls][tyr] = dp[p];
          }
          mgcmPressure[hgt][lat][ls][tyr] = pa0;

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
            mgcmEWWind[p][hgt][lat][ls][tyr] = up[p];
            mgcmNSWind[p][hgt][lat][ls][tyr] = vp[p];
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
        for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
          for (size_t p = 0; p < PARAM_SIZE; ++p) {
            mgcmTemperature[p][hgt][lat][0][tyr] = mgcmTemperature[p][hgt][lat][LS_SIZE - 1][tyr];
            mgcmDensity[p][hgt][lat][0][tyr] = mgcmDensity[p][hgt][lat][LS_SIZE - 1][tyr];
            mgcmEWWind[p][hgt][lat][0][tyr] = mgcmEWWind[p][hgt][lat][LS_SIZE - 1][tyr];
            mgcmNSWind[p][hgt][lat][0][tyr] = mgcmNSWind[p][hgt][lat][LS_SIZE - 1][tyr];
          }
          mgcmPressure[hgt][lat][0][tyr] = mgcmPressure[hgt][lat][LS_SIZE - 1][tyr];
        }
      }
    }
  }

    cout << " Writing binary file TES_lower_data.bin." <<  endl;
    ofstream binaryFile;
    binaryFile.open("TES_lower_data.bin", ios::binary | ios::out);

    size_t dataSize = MGCM_HEIGHT_SIZE * MGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmPressure);

    dataSize *= PARAM_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmTemperature);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmDensity);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmEWWind);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mgcmNSWind);

    binaryFile.close();
  }

// Reads ASCII version MTGCM 80-240 km data and writes binary
void mtgcmbin(char version) {
  static const string solact[MTGCM_F107_SIZE] = { "ls", "ms", "hs" };
  // Step through all solar activity levels
  for (size_t f107 = 0; f107 < MTGCM_F107_SIZE; ++f107) {
    ifstream baseHeight("zfTES" + solact[f107] + version + ".txt");
    // Step through all dust optical depths
    for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
      // Open ASCII version and binary version data files at 5m and 30m levels
      ifstream tpd("tpd" + solact[f107] + year[tyr] + version + ".txt");
      ifstream uv("uv" + solact[f107] + year[tyr] + version + ".txt");

      // Read (and ignore) the header line from the ASCII files
      string dummy;
      getline(tpd, dummy);
      getline(uv, dummy);

      // Step through all Ls values  (360 will be copied to 0 later)
      for (size_t ls = 1; ls < LS_SIZE; ls++) {
        int xls = (int)ls * 30;
        // Step through all latitudes
        for (size_t lat = 0; lat < MTGCM_LAT_SIZE; ++lat) {
          double xlat = -87.5 + lat * 5.0;
          // Step through all heights
          for (size_t hgt = 0; hgt < MTGCM_HEIGHT_SIZE; ++hgt) {
            int xhgt = 80 + (int)hgt * 5;
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
              mtgcmTemperature[p][hgt][lat][ls][tyr][f107] = tp[p];
              mtgcmPressure[p][hgt][lat][ls][tyr][f107] = pp[p];
              mtgcmDensity[p][hgt][lat][ls][tyr][f107] = dp[p];
            }

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
              mtgcmEWWind[p][hgt][lat][ls][tyr][f107] = up[p];
              mtgcmNSWind[p][hgt][lat][ls][tyr][f107] = vp[p];
            }
          }

          // Read thermosphere base heights
          double zfdust, ilat;
          int ils, nhgttop;
          baseHeight >> zfdust >> ils >> ilat
            >> mtgcmThermosphereBaseHeight[0][lat][ls][tyr][f107]
            >> mtgcmThermosphereBaseHeight[1][lat][ls][tyr][f107]
            >> mtgcmThermosphereBaseHeight[2][lat][ls][tyr][f107]
            >> mtgcmThermosphereBaseHeight[3][lat][ls][tyr][f107]
            >> mtgcmThermosphereBaseHeight[4][lat][ls][tyr][f107]
            >> nhgttop;
        }
      }
      // Close all MGCM surface data input files
      tpd.close();
      uv.close();
    }
    baseHeight.close();
  }

  for (size_t f107 = 0; f107 < MTGCM_F107_SIZE; ++f107) {
    for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
      for (size_t lat = 0; lat < MTGCM_LAT_SIZE; ++lat) {
        for (size_t p = 0; p < PARAM_SIZE; ++p) {
          for (size_t hgt = 0; hgt < MTGCM_HEIGHT_SIZE; ++hgt) {
            mtgcmTemperature[p][hgt][lat][0][tyr][f107] = mtgcmTemperature[p][hgt][lat][LS_SIZE - 1][tyr][f107];
            mtgcmPressure[p][hgt][lat][0][tyr][f107] = mtgcmPressure[p][hgt][lat][LS_SIZE - 1][tyr][f107];
            mtgcmDensity[p][hgt][lat][0][tyr][f107] = mtgcmDensity[p][hgt][lat][LS_SIZE - 1][tyr][f107];
            mtgcmEWWind[p][hgt][lat][0][tyr][f107] = mtgcmEWWind[p][hgt][lat][LS_SIZE - 1][tyr][f107];
            mtgcmNSWind[p][hgt][lat][0][tyr][f107] = mtgcmNSWind[p][hgt][lat][LS_SIZE - 1][tyr][f107];
          }
          mtgcmThermosphereBaseHeight[p][lat][0][tyr][f107] = mtgcmThermosphereBaseHeight[p][lat][LS_SIZE - 1][tyr][f107];
        }
      }
    }
  }

    cout << " Writing binary file TES_upper_data.bin." << endl;
    ofstream binaryFile;
    binaryFile.open("TES_upper_data.bin", ios::binary | ios::out);

    size_t dataSize = PARAM_SIZE * MTGCM_LAT_SIZE * LS_SIZE * TES_YEAR_SIZE * MTGCM_F107_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmThermosphereBaseHeight);

    dataSize *= MTGCM_HEIGHT_SIZE;
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmTemperature);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmPressure);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmDensity);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmEWWind);
    writeBinaryDataBlock(binaryFile, dataSize, (double*)mtgcmNSWind);

    binaryFile.close();
  }

// Reads ASCII version MGCM 0 - 80 km data and writes binary version
void dustbin(char version) {
  // Open ASCII version and binary version data files at 5m and 30m levels
  ifstream dust(string("TESdust") + version + ".txt");

  // Read (and ignore) the header line from the ASCII files
  string dummy;
  getline(dust, dummy);

  // Step through all TES years
  for (size_t tyr = 0; tyr < TES_YEAR_SIZE; ++tyr) {
    // Step through all Ls values  (360 will be copied to 0 later)
    for (size_t ls = 0; ls < DUST_LS_SIZE; ls++) {
      double xls = 2.5 + ls * 5;
      // Step through all latitudes
      for (size_t lat = 0; lat < DUST_LAT_SIZE; ++lat) {
        double xlat = -90.0 + lat * 7.5;
        // Step through all heights
        for (size_t lon = DUST_LON_SIZE; lon > 0; --lon) {
          int xlon = (int)lon * 9;
          // Read ASCII version tide coefficient amplitudes and phases
          // for temperature, pressure, and density.
          int iyr;
          double ils, ilat, ilon;
          if (lon == DUST_LON_SIZE) {
            dust >> iyr >> ils >> ilat >> ilon >> tes_tau[tyr][ls][lat][0];
          }
          else {
            dust >> iyr >> ils >> ilat >> ilon >> tes_tau[tyr][ls][lat][lon];
          }
          if (ils != xls) 
            cout << " Bad Ls" << endl;
          if (ilat != xlat) cout << " Bad Latitude" << endl;
          if (ilon != xlon) cout << " Bad Longitude" << endl;
        }
      }
    }
  }
  // Close all input files
  dust.close();

  cout << " Writing binary file TES_dust_data.bin." << endl;
  ofstream binaryFile;
  binaryFile.open("TES_dust_data.bin", ios::binary | ios::out);

  size_t dataSize = TES_YEAR_SIZE * DUST_LS_SIZE * DUST_LAT_SIZE * DUST_LON_SIZE * sizeof(double);
  binaryFile.write(reinterpret_cast<char*>(&tes_tau), dataSize);

  binaryFile.close();
}
