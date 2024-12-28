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
#include "NCEP.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

// Pressures (Pa).  Values for sea level (level=0) and surface(level = 1) are determined later.
const double NCEP::pn[presSize] = {
  999900.0, 999900.0, 100000.0, 92500.0, 85000.0, 70000.0, 60000.0, 50000.0, 40000.0, 30000.0,
  25000.0,  20000.0,  15000.0,  10000.0, 7000.0,  5000.0,  3000.0,  2000.0,  1000.0};


//! \brief Reads NCEP data from binary files.  
//!
//! Temperature (temp), density (dens),
//! dewpoint (dewp),  east-west wind (uwnd), north-south wind (vwnd),
//! geopotential height (geop), sea-level pressure (slpn), surface-level
//! pressure (sfcp), standard deviation of temperature (stmp), standard
//! deviation of density (sden), standard deviation of dewpoint (sdewp), 
//! standard deviation of east-west wind (suwd), standard deviation of 
//! north-south wind (svwd), standard deviation of pressure (sprs), standard
//! deviation of sea-level pressure (sslp), standard deviation of surface 
//! pressure, relative humidity (rhum), standard deviation of relative humidity
//! (srhum), vapor pressure (vprs), standard deviation of vapor pressure (svprs)
//! wind speed (spdav), standard deviation of wind speed (spdsd), u-v correlation
//! (uvcor) read into dynamic arrays using character pointer.  
//!
//! \b Inputs
//! \arg #NCEPHour          
//! \arg #NCEPYear          
//! \arg #hour          
//! \arg #month          
//! \arg #NCEPPath          
//!
//! \returns The NCEP data arrays are populated.
void NCEP::readNCEPFile()
{
  if (initialized) {
    return;
  }

  // Set the time of day code to use 
  // Code for UT hour of day if NCEP climatology is used : 
  // 1 = 00 UT, 2 = 06UT, 3 = 12UT, 4 = 18UT, 5 = all times of day combined, 
  // or 0 to use NCEP time-of-day based on input UTC hour (ihro)
  int timeOfDayCode = NCEPHour;
  if (timeOfDayCode == 0) {
    timeOfDayCode = 1 + int(round(hour / 6.0));
    if (timeOfDayCode > 4) {
      timeOfDayCode = 1;
    }
  }

  // Convert Month to a 2 digit string.
  string monthString = to_string(month);
  if (month < 10) {
    monthString = "0" + monthString;
  }

  // Create the NCEP file name.
  string NCEPFile = NCEPPath + "/" + "Nb" + to_string(NCEPYear) + monthString + ".bin";

  if (CONSOLE_OUTPUT) {
    cout << "Reading NCEP data from " << (NCEPFile) << endl;
  }

  // Try to open the file in folder NCEPpath.
  ifstream file;
  file.open(NCEPFile, ios::in | ios::binary | ios::ate);
  if (!file) {
    throw string(FILE_OPEN_ERROR_MESSAGE + NCEPFile);
  }

  // There are 19 3D arrays and 4 2D arrays.  Each is preceded and followed by an int.
  // So the total block size for an NCEP hour is:
  size_t totalBlockSize = (19 * arraySize3D + 4 * arraySize2D) * 8 + 46 * 4;

  // Go to the block corresponding to the NCEP hour.
  file.seekg((timeOfDayCode - 1) * totalBlockSize, ios::beg);

  blockOfDoubles = new double[arraySize3D];

  try {
    //Read temperature  
    readBlock(file, temp, arraySize3D);

    //Read density 
    readBlock(file, dens, arraySize3D);

    //Read dewpoint 
    readBlock(file, dewp, arraySize3D);

    //Read east-west wind 
    readBlock(file, uwnd, arraySize3D);

    //Read north-south wind 
    readBlock(file, vwnd, arraySize3D);

    //Read geopotential height 
    readBlock(file, geop, arraySize3D);

    //Read sea-level pressure
    readBlock(file, slpn, arraySize2D);

    //Read surface pressure
    readBlock(file, sfcp, arraySize2D);

    //Read standard deviation of temperature
    readBlock(file, stmp, arraySize3D);

    //Read standard deviation of density
    readBlock(file, sden, arraySize3D);

    //Read standard deviation of dewpoint
    readBlock(file, sdewp, arraySize3D);

    //Read standard deviation of east-west wind
    readBlock(file, suwd, arraySize3D);

    //Read standard deviation of north-south wind
    readBlock(file, svwd, arraySize3D);

    //Read standard deviation of pressure
    readBlock(file, sprs, arraySize3D);

    //Read standard deviation of sea-level pressure
    readBlock(file, sslp, arraySize2D);

    //Read standard deviation of surface pressure
    readBlock(file, ssfcp, arraySize2D);

    //Read relative humidity
    readBlock(file, rhum, arraySize3D);

    //Read standard deviation of relative humidity
    readBlock(file, srhum, arraySize3D);

    //Read vapor pressure
    readBlock(file, vprs, arraySize3D);

    //Read standard deviation of vapor pressure
    readBlock(file, svprs, arraySize3D);

    //Read wind speed
    readBlock(file, spdav, arraySize3D);

    //Read standard deviation of wind speed
    readBlock(file, spdsd, arraySize3D);

    //Read u-v correlation
    readBlock(file, uvcor, arraySize3D);
  }
  catch (const string& msg) {
    throw(FILE_READ_BINARY_ERROR_MESSAGE + NCEPFile
      + "\n       " + msg);
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_BINARY_ERROR_MESSAGE + NCEPFile
      + "\n       " + e.what());
  }

  //Close NCEP data file
  file.close();

  delete[] blockOfDoubles;
  blockOfDoubles = nullptr;
  initialized = true;
}

//! \brief Read a block of binary data.
//!
//! \param file The input file stream.
//! \param[out] dataBlock A pre-allocated buffer.
//! \param blockSize The number of doubles to be read.
//!
//! \retval dataBlock Populated with the ingested data.
void NCEP::readBlock(ifstream& file, greal* dataBlock, size_t blockSize)
{
  size_t totalBytes = sizeof(double) * blockSize;

  // Read the front block size tag
  size_t bytesFront = 0;
  file.read((char *)(&bytesFront), sizeof(int));

  // Read the data block of doubles.
  file.read((char *)(blockOfDoubles), totalBytes);

  // Read the back block size tag
  size_t bytesBack = 0;
  file.read((char *)(&bytesBack), sizeof(int));

  // Make sure the block sizes match.
  if (totalBytes != bytesFront || totalBytes != bytesBack) {
    throw ("Block sizes in file do not match.");
  }

  // Copy the doubles over to the data block (convert if necessary)
  if (sizeof(greal) == sizeof(double)) {
    memcpy(dataBlock, blockOfDoubles, totalBytes);
  }
  else {
    for (size_t i = 0; i < blockSize; ++i) {
      dataBlock[i] = greal(blockOfDoubles[i]);
    }
  }
}

} // namespace
