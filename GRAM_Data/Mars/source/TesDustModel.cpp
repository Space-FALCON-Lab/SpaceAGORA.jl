//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "TesDustModel.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool TesDustModel::isInitialized = false;
greal TesDustModel::tes_tau[MAP_YEAR_SIZE][LS_SIZE][LAT_SIZE][LON_SIZE];


//! \copydoc Atmosphere::Atmosphere()
TesDustModel::TesDustModel()
  : MarsDustModelBase()
{
  initializeData();
}

//! \fn  TesDustModel::TesDustModel(const TesDustModel& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TesDustModel::~TesDustModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void TesDustModel::setInputParameters(const MarsInputParameters& params)
{
  MarsDustModelBase::setInputParameters(params);

  mapYear = params.mapYear;
}

//! \brief Computes the dust optical depth (tau).
//!
//! Interpolates TES observed dust optical depth to current latitude, longitude, and Ls
//! for the specified map year.
//!
//! \b Inputs
//! \arg #position
//! \arg #longitudeSun
//! \arg #tes_tau
//!
//! \retval #dustOpticalDepth
void TesDustModel::updateDustOpticalDepth()
{
  // Step intervals
  greal stepLat = 180.0 / (LAT_SIZE - 1);
  greal stepLon = 360.0 / LON_SIZE;
  greal stepLs = 360.0 / LS_SIZE;     

  //=================================================================================================//
  //  define latitude parameters                                                                     //
  //-------------------------------------------------------------------------------------------------//
  // lower latitude data array index
  size_t lowerLatIndex = size_t(floor((position.latitude + 90.0_deg + stepLat) / stepLat)) - 1; 
  // ensure lower index is in bounds so index+1 exists
  if (lowerLatIndex >= LAT_SIZE - 1) {
    lowerLatIndex = LAT_SIZE - 2;
  }

  // upper latitude data array index
	size_t upperLatIndex = lowerLatIndex + 1;

  // latitude at lower index
  greal lowerLat = stepLat * upperLatIndex - (90.0_deg + stepLat); 
  // latitude relative displacement from value at lower index
  greal dLat = (position.latitude - lowerLat) / stepLat;      

  //=================================================================================================//
  //  define Ls parameters                                                                           //
  //-------------------------------------------------------------------------------------------------//
  // half step (for data offset, 1st value = 2.5)
  greal halfStepLs = stepLs / 2.0; 
  // lower Ls data array index (must be int to handle LS = 0)
	int lowerLsIndex = int(floor((longitudeSun + halfStepLs) / stepLs)) - 1;
  // upper Ls data array index
  int upperLsIndex = lowerLsIndex + 1;

  // Ls value at lower index
  greal lowerLs = stepLs * double(upperLsIndex) - halfStepLs;   

  // if lower index = 71 (max data array index) set upper index to 0
  // to handle cases near 360 degrees
  if (lowerLsIndex < 0 || lowerLsIndex == LS_SIZE - 1) {
    lowerLsIndex = LS_SIZE - 1;
    upperLsIndex = 0;   
  }

  // Ls relative displacement from value at lower index
  greal dLs = (longitudeSun - lowerLs) / stepLs;  

//=================================================================================================//
//  define longitude parameters                                                                    //
//-------------------------------------------------------------------------------------------------//
  // lower longitude data array index
  size_t lowerLonIndex = size_t(floor(position.getLongitude(WEST_POSITIVE) / stepLon));
  // if lower index=40 (for longitude=360) set index to 0 (for lon=0)
  if (lowerLonIndex == LON_SIZE) {
    lowerLonIndex = 0;
  }
  // upper longitude data array index
	size_t upperLonIndex = lowerLonIndex + 1;    

  // longitude value at lower index
  greal lowerLon = greal(lowerLonIndex) * stepLon;    
  // longitude relative displacement from value at lower index
  greal dLon = (position.getLongitude(WEST_POSITIVE) - lowerLon) / stepLon;

//=================================================================================================//
//  load TES tau values at corners of 3-d 'box' for Ls, latitude, and longitude dimensions         //
//-------------------------------------------------------------------------------------------------//
  // TES mapping year array index
  int yrx = mapYear - 1; 

  // declare arr3d structure to hold box values
  greal tauBox[2][2][2];     
	tauBox[0][0][0] = tes_tau[yrx][lowerLsIndex][lowerLatIndex][lowerLonIndex];
	tauBox[0][0][1] = tes_tau[yrx][lowerLsIndex][lowerLatIndex][upperLonIndex];
	tauBox[0][1][0] = tes_tau[yrx][lowerLsIndex][upperLatIndex][lowerLonIndex];
	tauBox[0][1][1] = tes_tau[yrx][lowerLsIndex][upperLatIndex][upperLonIndex];
	tauBox[1][0][0] = tes_tau[yrx][upperLsIndex][lowerLatIndex][lowerLonIndex];
  tauBox[1][0][1] = tes_tau[yrx][upperLsIndex][lowerLatIndex][upperLonIndex];
  tauBox[1][1][0] = tes_tau[yrx][upperLsIndex][upperLatIndex][lowerLonIndex];
  tauBox[1][1][1] = tes_tau[yrx][upperLsIndex][upperLatIndex][upperLonIndex];

//=================================================================================================//
//  compute and return interpolated dust optical depth value and exit routine                      //
//-------------------------------------------------------------------------------------------------//
  Interpolator interp3d(dLs, dLat, dLon);
	dustOpticalDepth = interp3d.linear(tauBox); 
}

//! \brief Reads the dust optical depth (tau) data from a binary file.
//!
//! The binary file contains a 4D-array of doubles containing the TES dust
//! optical depth (tau) values.
//!
//! \b Inputs
//! \arg #isInitialized
//!
//! \retval #tes_tau
void TesDustModel::initializeData()
{
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // The next two lines were used to read and convert the legacy file formats.
  //readLegacyData();
  //writeData();

  try {
    // Open a binary file stream.
    ifstream binaryFile;
    binaryFile.open(dataPath + tesDustFileName, ios::binary | ios::in);
    if (!binaryFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);  // catch and append file name
    }

    // Compute the size of the data block
    size_t dataSize = MAP_YEAR_SIZE * LS_SIZE * LAT_SIZE *LON_SIZE * sizeof(double);
    // Read the data block.
    binaryFile.read(reinterpret_cast<char*>(&tes_tau), dataSize);

    // Close the file stream.
    binaryFile.close();
  }
  catch (const std::string& msg) {
    throw msg + dataPath + tesDustFileName;
  }
}

#ifdef CONVERT_MARS_DATA
//! \brief Create a new format MOLA data file.
//!  
void TesDustModel::writeData()
{
  ofstream binaryFile;
   binaryFile.open(dataPath + "TES_dust_data.bin", ios::binary | ios::out);
  size_t dataSize = MAP_YEAR_SIZE * LS_SIZE * LAT_SIZE *LON_SIZE * sizeof(double);
  binaryFile.write(reinterpret_cast<char*>(&tes_tau), dataSize);
  binaryFile.close();
}

//! \brief Read the legacy data files for MOLA data.
//!  
void TesDustModel::readLegacyData()
{
  ifstream binaryFile;
   binaryFile.open(dataPath + "TESdust1.bin", ios::binary | ios::in);

  int bytes;
  binaryFile.read(reinterpret_cast<char*>(&bytes), sizeof(int));

  for (size_t lonx = 0; lonx < LON_SIZE; ++lonx) {
    for (size_t latx = 0; latx < LAT_SIZE; ++latx) {
      for (size_t lsx = 0; lsx < LS_SIZE; ++lsx) {
        for (size_t yrx = 0; yrx < MAP_YEAR_SIZE; ++yrx) {
          binaryFile.read(reinterpret_cast<char*>(&tes_tau[yrx][lsx][latx][lonx]), sizeof(double));
        }
      }
    }
  }

  binaryFile.read(reinterpret_cast<char*>(&bytes), sizeof(int));
  binaryFile.close();
}
#endif

} // namespace
