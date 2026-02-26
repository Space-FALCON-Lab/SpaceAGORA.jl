//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include "MERRA2.h"
#include "error_strings.h"

using namespace std;

// TODO comment new code

namespace GRAM {

const double MERRA2::pn[presSize]
    = {100000.0, 97500.0, 95000.0, 92500.0, 90000.0, 87500.0, 85000.0, 82500.0, 80000.0,
       77500.0, 75000.0, 72500.0, 70000.0, 65000.0, 60000.0, 55000.0, 50000.0, 45000.0, 40000.0,
       35000.0,  30000.0, 25000.0, 20000.0, 15000.0, 10000.0, 7000.0,  5000.0,  4000.0,
       3000.0,   2000.0,  1000.0,  700.0,   500.0,   400.0,   300.0,  200.0,  100.0,  
       70.0, 50.0, 40.0, 30.0, 10.0};

//! \brief Reads MERRA2 data from binary files.  
//!
//! This routine reads a MERRA-2 binary data file of any size.  The file size information is obtained
//! from the MERRA2info.txt file which must be present in data set. Since MERRA-2 data can be very large,
//! users can restrict reads by lat and lon.  Data arrays are allocated accordingly and populated.
//!
//! \b Inputs
//! \arg #M2Hour          
//! \arg #hour          
//! \arg #month          
//! \arg #M2Path          
//!
//! \returns The MERRA2 data arrays are populated.
void MERRA2::readM2File()
{
  if (initialized) {
    return;
  }
  
  // Get data sizes from the info file.
  readInfoFile();

  // Use bounds to compute sizes for dynamically allocated arrays.
  computeArraySizes();

  // Allocate memory for 1D arrays.
  allocateArrays();

  // Select the proper data file based on month and time of day code.
  buildDataFileName();

  // Notify the user
  if (CONSOLE_OUTPUT) {
    cout << "Reading MERRA2 data from " << (M2FileName) << endl;
  }

  ifstream file;
  file.open(M2FileName, ios::in | ios::binary | ios::ate);
  if (!file) {
    throw string(FILE_OPEN_ERROR_MESSAGE + M2FileName);
  }

  // Go to the beginning of the file.
  file.seekg(0, ios::beg);
  blockOfDoubles = new double[max(lonSize, lonSurfaceSize)];

  try {
    readBlock(file, dens, arraySize3D);
    readBlock(file, hgt, arraySize3D);
    readBlock(file, uwnd, arraySize3D);
    readBlock(file, temp, arraySize3D);
    readBlock(file, vwnd, arraySize3D);
    readBlock(file, uvcor, arraySize3D);
    readBlock(file, dewp, arraySize3D);
    readBlock(file, vprs, arraySize3D);
    readBlock(file, rhum, arraySize3D);
    readBlock(file, sden, arraySize3D);
    readBlock(file, spres, arraySize3D);
    readBlock(file, suwd, arraySize3D);
    readBlock(file, svwd, arraySize3D);
    readBlock(file, stmp, arraySize3D);
    readBlock(file, sdewp, arraySize3D);
    readBlock(file, svprs, arraySize3D);
    readBlock(file, srhum, arraySize3D);
    readBlock(file, spdav, arraySize3D);
    readBlock(file, spdsd, arraySize3D);
    readBlock(file, usfc, arraySize2DSurface);
    readBlock(file, susfc, arraySize2DSurface);
    readBlock(file, vsfc, arraySize2DSurface);
    readBlock(file, svsfc, arraySize2DSurface);
    readBlock(file, tsfc, arraySize2DSurface);
    readBlock(file, stsfc, arraySize2DSurface);
    readBlock(file, dsfc, arraySize2DSurface);
    readBlock(file, sdsfc, arraySize2DSurface);
    readBlock(file, psfc, arraySize2D);
    readBlock(file, spsfc, arraySize2D);
    readBlock(file, tdsfc, arraySize2D);
    readBlock(file, stdsfc, arraySize2D);
    readBlock(file, vpsfc, arraySize2D);
    readBlock(file, svpsfc, arraySize2D);
    readBlock(file, ghsfc, arraySize2D);
    readBlock(file, spdavsfc, arraySize2DSurface);
    readBlock(file, spdsdsfc, arraySize2DSurface);
    readBlock(file, uvcorsfc, arraySize2DSurface);
    readBlock(file, rhsfc, arraySize2D);
    readBlock(file, srhsfc, arraySize2D);
  }

  catch (const string& msg) {
    throw(FILE_READ_BINARY_ERROR_MESSAGE + M2FileName + "\n       " + msg);
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_BINARY_ERROR_MESSAGE + M2FileName + "\n       " + e.what());
  }

  // Close MERRA-2 data file
  file.close();

  delete[] blockOfDoubles;
  blockOfDoubles = nullptr;
  initialized = true;
}

//! \brief Read a block of binary data.
//!
//! Data can be restricted to a lat/lon rectangle by the user.  This routine reads only the data within
//! the specified rectangle.  For 3D arrays, all pressure levels are read in.
//!
//! \param file The input file stream.
//! \param[out] dataBlock A pre-allocated buffer.
//! \param blockSize The number of doubles to be read.
//!
//! \retval dataBlock Populated with the ingested data.
void MERRA2::readBlock(ifstream& file, greal* dataBlock, size_t blockSize)
{
  // Copy the doubles over to the data block (convert if necessary)
  if (blockSize == arraySize3D) {
    size_t j = 0;
    for (size_t presIndex = 0; presIndex < presSize; ++presIndex) {
      file.seekg(lonFullSize * latStartIndex * sizeof(double), ios::cur);
      for (size_t latIndex = 0; latIndex < latSize; ++latIndex) {
        for (int i = 0; i < 2; ++i) {
          if (lonData[i] > 0) {
            file.read((char*)(blockOfDoubles), lonData[i] * sizeof(double));
            for (size_t lonIndex = 0; lonIndex < lonData[i]; ++lonIndex) {
              size_t dataIndex = j;
              if (lonSkip[1] == 0) {
                if (i == 0) {
                  dataIndex += lonData[1];
                }
                else {
                  dataIndex -= lonData[0];
                }
              }
              dataBlock[dataIndex] = greal(blockOfDoubles[lonIndex]);
              j += 1;
            }
          }
          file.seekg(lonSkip[i] * sizeof(double), ios::cur);
        }
        if (lonSize > lonFullSize) {
          dataBlock[j] = dataBlock[j - lonFullSize];
          j += 1;
          dataBlock[j] = dataBlock[j - lonFullSize];
          j += 1;
        }
      }
      file.seekg((latFullSize - latStartIndex - latSize) * lonFullSize * sizeof(double), ios::cur);
    }
  }
  else if (blockSize == arraySize2D) {
    size_t j = 0;
    file.seekg(lonFullSize * latStartIndex * sizeof(double), ios::cur);
    for (size_t latIndex = 0; latIndex < latSize; ++latIndex) {
      for (int i = 0; i < 2; ++i) {
        if (lonData[i] > 0) {
          file.read((char*)(blockOfDoubles), lonData[i] * sizeof(double));
          for (size_t lonIndex = 0; lonIndex < lonData[i]; ++lonIndex) {
            size_t dataIndex = j;
            if (lonSkip[1] == 0) {
              if (i == 0) {
                dataIndex += lonData[1];
              }
              else {
                dataIndex -= lonData[0];
              }
            }
            dataBlock[dataIndex] = greal(blockOfDoubles[lonIndex]);
            j += 1;
          }
        }
        file.seekg(lonSkip[i] * sizeof(double), ios::cur);
      }
      if (lonSize > lonFullSize) {
        dataBlock[j] = dataBlock[j - lonFullSize];
        j += 1;
        dataBlock[j] = dataBlock[j - lonFullSize];
        j += 1;
      }
    }
    file.seekg((latFullSize - latStartIndex - latSize) * lonFullSize * sizeof(double), ios::cur);
  }
  else {
    size_t j = 0;
    file.seekg(lonSurfaceFullSize * latStartIndex * sizeof(double), ios::cur);
    for (size_t latIndex = 0; latIndex < latSize; ++latIndex) {
      for (int i = 0; i < 2; ++i) {
        if (lonSurfaceData[i] > 0) {
          file.read((char*)(blockOfDoubles), lonSurfaceData[i] * sizeof(double));
          for (size_t lonIndex = 0; lonIndex < lonSurfaceData[i]; ++lonIndex) {
            size_t dataIndex = j;
            if (lonSurfaceSkip[1] == 0) {
              if (i == 0) {
                dataIndex += lonSurfaceData[1];
              }
              else {
                dataIndex -= lonSurfaceData[0];
              }
            }
            dataBlock[dataIndex] = greal(blockOfDoubles[lonIndex]);
            j += 1;
          }
        }
        file.seekg(lonSurfaceSkip[i] * sizeof(double), ios::cur);
      }
      if (lonSurfaceSize > lonSurfaceFullSize) {
        dataBlock[j] = dataBlock[j - lonSurfaceFullSize];
        j += 1;
        dataBlock[j] = dataBlock[j - lonSurfaceFullSize];
        j += 1;
      }
    }
    file.seekg((latFullSize - latStartIndex - latSize) * lonSurfaceFullSize * sizeof(double), ios::cur);
  }
}

//! \brief Reads the MERRA2 info file for data sizes.
//!
//! Data can be restricted to a lat/lon rectangle by the user.  This routine reads only the data within
//! the specified rectangle.  For 3D arrays, all pressure levels are read in.
//!
//! \param #M2Path
//!
//! \retval #latFullSize
//! \retval #lonFullSize
//! \retval #lonSurfaceFullSize
void MERRA2::readInfoFile() {
  // Get the data info
  string infoFile = M2Path + "/MERRA2info.txt";

  // Read presSize, latFullSize, lonFullSize, lonSurfaceFullSize
  ifstream info;
  info.open(infoFile, ios::in);
  if (!info) {
    throw string(FILE_OPEN_ERROR_MESSAGE + infoFile);
  }
  string varName, equalSign;
  size_t pres;
  info >> varName;
  info >> varName >> equalSign >> pres;
  if (pres != presSize) {
    throw string("Error: Incorrect pressure levels found in\n"
                 "       MERRA-2 info file: " + infoFile + ".\n");
  }
  info >> varName >> equalSign >> latFullSize;
  info >> varName >> equalSign >> lonFullSize;
  info >> varName >> varName >> equalSign >> lonSurfaceFullSize;
  info.close();
}

//! \brief Determines array sizes based on lat/lon restrictions.
//!
//! Data can be restricted to a lat/lon rectangle by the user.  This routine determines the size of
//! each dimension of the 3D and 2Ddata arrays.
//!
//! \param #minimumLatitude
//! \param #maximumLatitude
//! \param #minimumLongitude
//! \param #maximumLongitude
//!
//! \retval #latSize
//! \retval #lonSize
//! \retval #lonSurfaceSize
//! \retval #arraySize3D
//! \retval #arraySize2D
//! \retval #arraySize2DSurface
void MERRA2::computeArraySizes() {
  if (maximumLatitude <= minimumLatitude) {
    throw string("Error: Incorrect latitude bounds.\n"
      "       MaximumLatitude must be larger than MinimumLatitude.\n");
  }
  if (maximumLatitude > 90.0_deg) {
    if (CONSOLE_OUTPUT) {
      maximumLatitude = 90.0_deg;
      cout << "Warning: MaximumLatitude too large.  Using 90 degrees.\n";
    }
  }
  if (minimumLatitude < -90.0_deg) {
    if (CONSOLE_OUTPUT) {
      minimumLatitude = -90.0_deg;
      cout << "Warning: MinimumLatitude too small.  Using -90 degrees.\n";
    }
  }

  // Compute array sizes based off data sizes and region of interest.
  double latStep = 180.0_deg / (latFullSize - 1);
  latStartIndex = (size_t)(floor((minimumLatitude + 90.0_deg) / latStep));
  size_t latEndIndex = (size_t)(ceil((maximumLatitude + 90.0_deg) / latStep));
  latSize = latEndIndex - latStartIndex + 1;

  if (maximumLongitude <= minimumLongitude) {
    throw string("Error: Incorrect longitude bounds.\n"
      "       MaximumLongitude must be larger than MinimumLongitude.\n");
  }
  if (minimumLongitude >= 0.0_deg) {
    if (maximumLongitude > 360.0_deg) {
      if (CONSOLE_OUTPUT) {
        maximumLongitude = 360.0_deg;
        cout << "Warning: MaximumLongitude too large.  Using 360 degrees.\n";
      }
    }
  }
  else {
    if (minimumLongitude < -180.0_deg) {
      if (CONSOLE_OUTPUT) {
        maximumLongitude = -180.0_deg;
        cout << "Warning: MinimumLongitude too small.  Using -180 degrees.\n";
      }
    }
    if (maximumLongitude > 360.0_deg + minimumLongitude) {
      if (CONSOLE_OUTPUT) {
        maximumLongitude = 360.0_deg + minimumLongitude;
        cout << "Warning: MaximumLongitude too large.  Using "
          << maximumLongitude << " degrees.\n";
      }
    }
  }

  double lonStep = 360.0_deg / lonFullSize;
  lonStartIndex = (size_t)(floor(wrapDegrees(minimumLongitude + 180.0_deg) / lonStep));
  size_t lonEndIndex = (size_t)(ceil(wrapDegrees(maximumLongitude + 180.0_deg) / lonStep));
  if (minimumLongitude < 180.0_deg) {
    lonShift = lonStartIndex * lonStep - 180.0_deg;
  }
  else {
    lonShift = lonStartIndex * lonStep + 180.0_deg;
  }
  if (lonStartIndex < lonEndIndex) {
    lonSize = lonEndIndex - lonStartIndex + 2;
    lonData[0] = 0;
    lonData[1] = lonSize;
    lonSkip[0] = lonStartIndex;
    lonSkip[1] = lonFullSize - lonSize - lonStartIndex;
  }
  else if (lonStartIndex > lonEndIndex) {
    lonSize = lonEndIndex + lonFullSize - lonStartIndex + 2;
    lonData[0] = lonStartIndex + lonSize - lonFullSize;
    lonData[1] = lonFullSize - lonStartIndex;
    lonSkip[0] = lonFullSize - lonSize;
    lonSkip[1] = 0;
  }
  else { // lonStartIndex == lonEndIndex
    lonSize = lonFullSize + 2;
    lonData[0] = lonStartIndex;
    lonData[1] = lonFullSize - lonStartIndex;
    lonSkip[0] = 0;
    lonSkip[1] = 0;
  }

  double lonSurfaceStep = 360.0_deg / lonSurfaceFullSize;
  lonSurfaceStartIndex = (size_t)(floor(wrapDegrees(minimumLongitude + 180.0_deg) / lonSurfaceStep));
  size_t lonSurfaceEndIndex = (size_t)(ceil(wrapDegrees(maximumLongitude + 180.0_deg) / lonSurfaceStep));
  if (minimumLongitude < 180.0_deg) {
    lonSurfaceShift = lonSurfaceStartIndex * lonSurfaceStep - 180.0_deg;
  }
  else {
    lonSurfaceShift = lonSurfaceStartIndex * lonSurfaceStep + 180.0_deg;
  }
  if (lonSurfaceStartIndex < lonSurfaceEndIndex) {
    lonSurfaceSize = lonSurfaceEndIndex - lonSurfaceStartIndex + 2;
    lonSurfaceData[0] = 0;
    lonSurfaceData[1] = lonSurfaceSize;
    lonSurfaceSkip[0] = lonSurfaceStartIndex;
    lonSurfaceSkip[1] = lonSurfaceFullSize - lonSurfaceSize - lonSurfaceStartIndex;
  }
  else if (lonSurfaceStartIndex > lonSurfaceEndIndex) {
    lonSurfaceSize = lonSurfaceEndIndex + lonSurfaceFullSize - lonSurfaceStartIndex + 2;
    lonSurfaceData[0] = lonSurfaceStartIndex + lonSurfaceSize - lonSurfaceFullSize;
    lonSurfaceData[1] = lonSurfaceFullSize - lonSurfaceStartIndex;
    lonSurfaceSkip[0] = lonSurfaceFullSize - lonSurfaceSize;
    lonSurfaceSkip[1] = 0;
  }
  else { // lonSurfaceStartIndex == lonSurfaceEndIndex
    lonSurfaceSize = lonSurfaceFullSize + 2;
    lonSurfaceData[0] = lonSurfaceStartIndex;
    lonSurfaceData[1] = lonSurfaceFullSize - lonSurfaceStartIndex;
    lonSurfaceSkip[0] = 0;
    lonSurfaceSkip[1] = 0;
  }

  // Size of the 3D arrays allocated as 1D.
  arraySize3D = presSize * latSize * lonSize;
  // Size of the 2D arrays allocated as 1D.
  arraySize2D = latSize * lonSize;
  // Size of the 2D surface data arrays allocated as 1D.
  arraySize2DSurface = latSize * lonSurfaceSize;
}

//! \brief Allocate data arrays.
//!
//! \param #arraySize3D
//! \param #arraySize2D
//! \param #arraySize2DSurface
//!
//! \returns Allocated data arrays.
void MERRA2::allocateArrays() {
  try {
    deleteArrays();
    temp = new greal[arraySize3D];
    dens = new greal[arraySize3D];
    dewp = new greal[arraySize3D];
    uwnd = new greal[arraySize3D];
    hgt = new greal[arraySize3D];
    vwnd = new greal[arraySize3D];
    stmp = new greal[arraySize3D];
    sden = new greal[arraySize3D];
    spres = new greal[arraySize3D];
    sdewp = new greal[arraySize3D];
    suwd = new greal[arraySize3D];
    svwd = new greal[arraySize3D];
    sprs = new greal[arraySize3D];
    rhum = new greal[arraySize3D];
    srhum = new greal[arraySize3D];
    vprs = new greal[arraySize3D];
    svprs = new greal[arraySize3D];
    spdav = new greal[arraySize3D];
    spdsd = new greal[arraySize3D];
    uvcor = new greal[arraySize3D];
    usfc = new greal[arraySize2DSurface];
    susfc = new greal[arraySize2DSurface];
    vsfc = new greal[arraySize2DSurface];
    svsfc = new greal[arraySize2DSurface];
    tsfc = new greal[arraySize2DSurface];
    stsfc = new greal[arraySize2DSurface];
    dsfc = new greal[arraySize2DSurface];
    sdsfc = new greal[arraySize2DSurface];
    tdsfc = new greal[arraySize2D];
    stdsfc = new greal[arraySize2D];
    vpsfc = new greal[arraySize2D];
    svpsfc = new greal[arraySize2D];
    psfc = new greal[arraySize2D];
    spsfc = new greal[arraySize2D];
    ghsfc = new greal[arraySize2D];
    spdavsfc = new greal[arraySize2DSurface];
    spdsdsfc = new greal[arraySize2DSurface];
    uvcorsfc = new greal[arraySize2DSurface];
    rhsfc = new greal[arraySize2D];
    srhsfc = new greal[arraySize2D];
  }
  catch (const std::bad_alloc& ex) {
    throw(MEMORY_ALLOCATION_ERROR_MESSAGE + "Attempting to allocate MERRA2 arrays. " + ex.what() + "\n");
  }
}

//! \brief Delete data arrays.
//!
//! \returns NULL data arrays.
void MERRA2::deleteArrays()
{
  if (temp != nullptr) {
    delete[] temp;
    delete[] dens;
    delete[] dewp;
    delete[] uwnd;
    delete[] hgt;
    delete[] vwnd;
    delete[] stmp;
    delete[] sden;
    delete[] sdewp;
    delete[] suwd;
    delete[] svwd;
    delete[] sprs;
    delete[] rhum;
    delete[] srhum;
    delete[] vprs;
    delete[] svprs;
    delete[] spdav;
    delete[] spdsd;
    delete[] uvcor;
    delete[] usfc;
    delete[] susfc;
    delete[] vsfc;
    delete[] svsfc;
    delete[] tsfc;
    delete[] stsfc;
    delete[] dsfc;
    delete[] sdsfc;
    delete[] tdsfc;
    delete[] stdsfc;
    delete[] vpsfc;
    delete[] svpsfc;
    delete[] ghsfc;
    delete[] spdavsfc;
    delete[] spdsdsfc;
    delete[] uvcorsfc;
    delete[] rhsfc;
    delete[] srhsfc;
    delete[] spres;
    temp = nullptr;
    dens = nullptr;
    dewp = nullptr;
    uwnd = nullptr;
    hgt = nullptr;
    vwnd = nullptr;
    stmp = nullptr;
    sden = nullptr;
    sdewp = nullptr;
    suwd = nullptr;
    svwd = nullptr;
    sprs = nullptr;
    rhum = nullptr;
    srhum = nullptr;
    vprs = nullptr;
    svprs = nullptr;
    spdav = nullptr;
    spdsd = nullptr;
    uvcor = nullptr;
    usfc = nullptr;
    susfc = nullptr;
    vsfc = nullptr;
    svsfc = nullptr;
    tsfc = nullptr;
    stsfc = nullptr;
    dsfc = nullptr;
    sdsfc = nullptr;
    tdsfc = nullptr;
    stdsfc = nullptr;
    vpsfc = nullptr;
    svpsfc = nullptr;
    ghsfc = nullptr;
    spdavsfc = nullptr;
    spdsdsfc = nullptr;
    uvcorsfc = nullptr;
    rhsfc = nullptr;
    srhsfc = nullptr;
    spres = nullptr;
  }
}

//! \brief Determine the MERRA-2 data file name.
//!
//! \param #month
//! \param #M2Hour
//! \param #M2Path
//!
//! \retval #M2FileName
void MERRA2::buildDataFileName()
{
  // Convert Month to a 2 digit string.
  string monthString = to_string(month);
  if (month < 10) {
    monthString = "0" + monthString;
  }

  // Day is divided into 3-hour intervals.
  int timeOfDayCode = M2Hour;
  if (timeOfDayCode == 0) {
    timeOfDayCode = 1 + int(round(hour / 3.0));
    if (timeOfDayCode > 8) {
      timeOfDayCode = 1;
    }
  }

  // Create the MERRA-2 file name.
  string M2CodeName;
  switch (timeOfDayCode) {
  case 1:
    M2CodeName = "00Z/MERRA2_3hr_00Z";
    break;
  case 2:
    M2CodeName = "03Z/MERRA2_3hr_03Z";
    break;
  case 3:
    M2CodeName = "06Z/MERRA2_3hr_06Z";
    break;
  case 4:
    M2CodeName = "09Z/MERRA2_3hr_09Z";
    break;
  case 5:
    M2CodeName = "12Z/MERRA2_3hr_12Z";
    break;
  case 6:
    M2CodeName = "15Z/MERRA2_3hr_15Z";
    break;
  case 7:
    M2CodeName = "18Z/MERRA2_3hr_18Z";
    break;
  case 8:
    M2CodeName = "21Z/MERRA2_3hr_21Z";
    break;
  case 9:
  default:
    M2CodeName = "All Mean/MERRA2All";
    break;
  }

  M2FileName = M2Path + "/" + M2CodeName + "_" + monthString + ".bin";
}

} // namespace GRAM