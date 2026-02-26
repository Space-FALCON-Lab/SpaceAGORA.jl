//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <array>
#include <string>
#include "gram.h"

namespace GRAM {

//! \brief This is a InputParameters.
//!
//! This is a InputParameters.
//! \ingroup CommonGRAM
class InputParameters
{
public:
  InputParameters();
  InputParameters(const InputParameters& orig)=default;
  virtual ~InputParameters() = default;
  InputParameters& operator=(const InputParameters& orig) = default;

  // Used by SpiceLoader
  std::string spicePath = DEFAULT_SPICE_PATH;  //!< For SpiceLoader via PerturbedAtmosphere (also SPICEDIR).
  std::string spiceLsk;
  std::string spicePck;
  std::string spiceVenus;
  std::string spiceEarth;
  std::string spiceMars;
  std::string spiceJupiter;
  std::string spiceSaturn;
  std::string spiceUranus;
  std::string spiceNeptune;
  std::string spiceTitan;

  // For models with data folders.
  std::string dataPath = "";                   //!< Data folder location.

  // ProfilePrinter parameters
  std::string listFileName = "LIST";              //!< For ProfilePrinter (also LSTFL).
  std::string columnFileName = "OUTPUT";          //!< For ProfilePrinter (also OUTFL).
  DensityPrintScale densityPrintScale = STANDARD; //!< For ProfilePrinter (also LOGSCALE).
  bool isEastLongitudePositiveOnOutput = true;    //!< For ProfilePrinter (also LONEAST).
  bool useLegacyOutputs = false;                  //!< For ProfilePrinter.
  size_t extraPrecision = 0;                      //!< For ProfilePrinter.

  // Used by AuxiliaryAtmosphere input files,
  // Trajectory input Files, and positional input parameters
  // (initialLongitude and deltaLongitude)
  bool isEastLongitudePositiveOnInput = true;     //!< For input classes (also LONEAST).

  // Perturbed Atmosphere parameters:
  // GramTime parameters used in PerturbedAtmosphere.
  GRAM_TIME_SCALE timeScale = UTC;     //!< For GramTime via PerturbedAtmosphere (also IUTC).
  GRAM_TIME_FRAME timeFrame = ERT;     //!< For GramTime via PerturbedAtmosphere (also IERT).
  int month = 1;                       //!< For GramTime via PerturbedAtmosphere (also MONTH).
  int day = 1;                         //!< For GramTime via PerturbedAtmosphere (also MDAY).
  int year = 2000;                     //!< For GramTime via PerturbedAtmosphere (also MYEAR).
  int hour = 0;                        //!< For GramTime via PerturbedAtmosphere (also IHOUR).
  int minute = 0;                      //!< For GramTime via PerturbedAtmosphere (also IMIN).
  greal seconds = 0.0;                 //!< For GramTime via PerturbedAtmosphere (also SEC).
  // Pertubation parameters
  int initialRandomSeed = 1001;              //!< For PerturbedAtmosphere (also NR1).
  greal densityPerturbationScale = 1.0;      //!< For PerturbedAtmosphere (also RPSCALE).
  greal ewWindPerturbationScale = 1.0;       //!< For PerturbedAtmosphere (also RPSCALE).
  greal nsWindPerturbationScale = 1.0;       //!< For PerturbedAtmosphere (also RPSCALE).
  greal verticalWindPerturbationScale = 1.0; //!< For PerturbedAtmosphere (also RPSCALE).
  greal minRelativeStepSize = 0.0;           //!< For PerturbedAtmosphere (also CORLMIN).

  // Ephemeris Parameters
  bool fastModeOn = false;

  // Trajectory/VerticalProfile parameters:
  std::string trajectoryFileName = ""; //!< For Trajectory (also TRAJFL).
  bool isPlanetoCentric = true;        //!< For Trajectory (also IPCLAT).
  bool initialPerturbationsUpdated = true;     //!< For Trajectory.
  int numberOfPositions = 21;          //!< For Trajectory (also NPOS).
  greal initialLatitude = 0.0;         //!< For Trajectory (also FLAT).
  greal initialLongitude = 0.0;        //!< For Trajectory (also FLON).
  greal initialHeight = 0.0;           //!< For Trajectory (also FHGT).
  greal deltaLatitude = 0.0;           //!< For Trajectory (also DELLAT).
  greal deltaLongitude = 0.0;          //!< For Trajectory (also DELLON).
  greal deltaHeight = 10.0;            //!< For Trajectory (also DELHGT).
  greal deltaTime = 0.0;               //!< For Trajectory (also DELTIME).

  // MonteCarlo parameters:
  int numberOfMonteCarloRuns = 1;      //!< For MonteCarlo (also NMONTE).

  // Auxiliary Atmosphere(profile) parameters:
  std::array<std::string, 10> auxiliaryAtmosphereFileName; //!< For AuxiliaryAtmosphere via AuxiliaryAdapter (also PROFILE).
  std::array<greal, 10> innerRadius = { 0.0 };             //!< For AuxiliaryAtmosphere via AuxiliaryAdapter (also PROFNEAR).
  std::array<greal, 10> outerRadius = { 0.0 };             //!< For AuxiliaryAtmosphere via AuxiliaryAdapter (also PROFFAR).

  // unused parameters
  int nvarx = 1;             //!< NVARX (not used)
  int nvary = 0;             //!< NVARY (not used)
  int iup = 0;               //!< IUP (not used)

  // FindDates Parameters
  bool findDates = false;            //!< Set to true to use FindDates routine
  greal targetLongitudeSun = 0.0;
  greal targetSolarTime = 0.0;
  // initialLongitude (see above)
  // all GramTime parameters (see above)
  // eastLongitudePositive (see above)
};

} // namespace
