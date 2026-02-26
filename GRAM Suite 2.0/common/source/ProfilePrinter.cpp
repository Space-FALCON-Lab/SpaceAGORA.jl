//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <thread>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "ProfilePrinter.h"
#include "ConstituentGas.h"

using namespace std;

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
//! In the default constructor the interpolation fraction defaults to 0.
ProfilePrinter::ProfilePrinter()
{
}

//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.
//! \param orig The source object to copy.
ProfilePrinter::ProfilePrinter(const ProfilePrinter& orig)
{
  densUnits = orig.densUnits;
  densityPrintScale = orig.densityPrintScale;
  eastLongitudePositive = orig.eastLongitudePositive;
  listFileName = orig.listFileName;
  columnFileName = orig.columnFileName;
  fileNamePrefix = orig.fileNamePrefix;
  outputStyle = orig.outputStyle;
  // DO NOT copy the output streams.
}

//! \fn  ProfilePrinter::~ProfilePrinter()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.
//! \param params The input parameters.
void ProfilePrinter::setInputParameters(const InputParameters& params)
{
  setEastLongitudePositive(params.isEastLongitudePositiveOnOutput);
  setColumnFileName(params.columnFileName);
  setListFileName(params.listFileName);
  setDensityPrintScale(params.densityPrintScale);
  setExtraPrecision(params.extraPrecision);
}

//! \brief Output the profile data for the specified atmosphere.
//!
//! For all of the selected print styles, this routine will print the provided Profile data for
//! the specified atmosphere.  This method will open files, print the data, and close the files.
//! If more than one set of Profile data is to be output, then the openOutput(), printData(), and 
//! closeOutput() methods should be called directly.
//! \param atmos The PerturbedAtmosphere which generated the data.
//! \param profile A list of ProfileData objects.
void ProfilePrinter::print(const PerturbedAtmosphere& atmos, const std::vector<ProfileData>& profile)
{
  // Open files, print data, close files.
  openOutput();
  printFileHeader(atmos);
  printSectionHeader(atmos);
  printData(profile);
  closeOutput();
}

void ProfilePrinter::openFile(std::ofstream& outputStream, const std::string& fileName) {
  outputStream.open(fileName);
  if (!outputFiles.empty()) {
    outputFiles += ", ";
  }
  outputFiles += fileName;
}

//! \brief Open files for output.
//!
//! Based on the outputStyle bit field, this method will open files for all selected print styles.
//! When adding a print style, create an ofstream member for the new style and open it here.
void ProfilePrinter::openOutput()
{
  outputFiles.clear();
  if ((outputStyle & GRAM_CSV_STYLE) && !columnFileName.empty()) {
    appendExtension(columnFileName, ".csv");
    openFile(gramCSVFile, fileNamePrefix + columnFileName);
  }
  if ((outputStyle & GRAM_MD_STYLE) && !listFileName.empty()) {
    appendExtension(listFileName, ".md");
    openFile(gramMDFile, fileNamePrefix + listFileName);
  }
  if ((outputStyle & GRAM_COLUMN_STYLE) && !columnFileName.empty()) {
    appendExtension(columnFileName, ".txt");
    openFile(gramColumnFile, fileNamePrefix + columnFileName);
  }
  if ((outputStyle & GRAM_LIST_STYLE) && !listFileName.empty()) {
    appendExtension(listFileName, ".txt");
    openFile(gramListFile, fileNamePrefix + listFileName);
  }
  if (outputStyle & GRAM_EPHEMERIS_STYLE) {
    openFile(gramEphemerisFile, fileNamePrefix + "Ephem.txt");
  }
  if (outputStyle & GRAM_DENSITY_STYLE) {
    openFile(gramDensityFile, fileNamePrefix + "Density.txt");
  }
  if (outputStyle & GRAM_PERTURB_STYLE) {
    openFile(gramPerturbFile, fileNamePrefix + "Perturb.txt");
  }
  if (outputStyle & GRAM_WINDS_STYLE) {
    openFile(gramWindsFile, fileNamePrefix + "Winds.txt");
  }
  if (outputStyle & GRAM_TPRESHGT_STYLE) {
    openFile(gramTPresHgtFile, fileNamePrefix + "TPresHgt.txt");
  }
  if (outputStyle & GRAM_SOUND_STYLE) {
    openFile(gramSoundFile, fileNamePrefix + "Sound.txt");
  }
}

//! \brief Print a file header for the selected print styles.
//!
//! A header is typically printed once at the top of each output file.  For columnar output, 
//! the header is a row of field names.
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printFileHeader(const PerturbedAtmosphere& atmos)
{
  if (outputStyle & GRAM_CSV_STYLE) {
    printGramCSVHeader(atmos);
  }
  if (outputStyle & GRAM_MD_STYLE) {
    printGramListMDHeader(atmos);
  }
  if (outputStyle & GRAM_COLUMN_STYLE) {
    printGramColumnHeader(atmos);
  }
  if (outputStyle & GRAM_LIST_STYLE) {
		printGramListHeader(atmos);
  }
  if (outputStyle & GRAM_EPHEMERIS_STYLE) {
		printGramEphemerisHeader(atmos);
  }
  if (outputStyle & GRAM_DENSITY_STYLE) {
		printGramDensityHeader(atmos);
  }
  if (outputStyle & GRAM_PERTURB_STYLE) {
		printGramPerturbHeader(atmos);
  }
  if (outputStyle & GRAM_WINDS_STYLE) {
		printGramWindsHeader(atmos);
  }
  if (outputStyle & GRAM_TPRESHGT_STYLE) {
		printGramTPresHgtHeader(atmos);
  }
  if (outputStyle & GRAM_SOUND_STYLE) {
		printGramSoundHeader(atmos);
  }
}

//! \brief Print a section header for the selected print styles.
//!
//! A section header is printed before each set of profile data.  A call to printSectionHeader()
//! should be paired with a call to printData().
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printSectionHeader(const PerturbedAtmosphere& atmos)
{
  if (outputStyle & GRAM_MD_STYLE) {
    mdCount = 1;
  }
  if (outputStyle & GRAM_LIST_STYLE) {
    gramListFile << "   Random seed = " << atmos.getSeed() << newline;
  }
}

//! \brief Print the profile data for the selected print styles.
//!
//! This method will loop through the list of Profile data and output that record
//! to each selected print style.  A call to printSectionHeader()
//! should be paired with a call to printData().
//! \param profile A vector of ProfileData.
void ProfilePrinter::printData(const std::vector<ProfileData>& profile)
{
  vector<thread> threads;
    if (outputStyle & GRAM_CSV_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramCSVStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_MD_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramListMDStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_COLUMN_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramColumnStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_LIST_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramListStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_EPHEMERIS_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramEphemerisStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_DENSITY_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramDensityStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_PERTURB_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramPerturbStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_WINDS_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramWindsStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_TPRESHGT_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramTPresHgtStyle(p);
        }
      }));
    }
    if (outputStyle & GRAM_SOUND_STYLE) {
      threads.push_back(thread([&]()
      {
        for (auto& p : profile) {
          printGramSoundStyle(p);
        }
      }));
    }
  for (auto& t : threads) {
    t.join();
  }
}

//! \brief Close files for output.
//!
//! Based on the outputStyle bit field, this method will close files for all selected print styles.
void ProfilePrinter::closeOutput()
{
  if (outputStyle & GRAM_CSV_STYLE) {
    gramCSVFile.close();
  }
  if (outputStyle & GRAM_MD_STYLE) {
    gramMDFile << "--------------------\n";
    gramMDFile << "## End of data" << newline;
    gramMDFile << "--------------------\n\n";
    gramMDFile.close();
  }
  if (outputStyle & GRAM_COLUMN_STYLE) {
    gramColumnFile.close();
  }
  if (outputStyle & GRAM_LIST_STYLE) {
    gramListFile.close();
  }
  if (outputStyle & GRAM_EPHEMERIS_STYLE) {
    gramEphemerisFile.close();
  }
  if (outputStyle & GRAM_DENSITY_STYLE) {
    gramDensityFile.close();
  }
  if (outputStyle & GRAM_PERTURB_STYLE) {
    gramPerturbFile.close();
  }
  if (outputStyle & GRAM_WINDS_STYLE) {
    gramWindsFile.close();
  }
  if (outputStyle & GRAM_TPRESHGT_STYLE) {
    gramTPresHgtFile.close();
  }
  if (outputStyle & GRAM_SOUND_STYLE) {
    gramSoundFile.close();
  }
}

//! \brief Scales density values based on densityPrintScale setting.
//!
//! The density values are scaled according to the densityPrintScale setting.  This method will alter
//! the values of the input AtmosphereState.
//! \param atmos An AtmosphereState.
void ProfilePrinter::scaleDensityForOutput(AtmosphereState& atmos)
{
	switch (densityPrintScale) {
	case LOG_10:
		atmos.density = log10(atmos.density);
		atmos.pressure = log10(atmos.pressure);
		atmos.lowDensity = log10(atmos.lowDensity);
		atmos.highDensity = log10(atmos.highDensity);
		atmos.perturbedDensity = log10(atmos.perturbedDensity);
    densUnits = "kg/m^3";
    break;
	case PERCENT_DEV:
	{
		atmos.density = atmos.lowDensityDeviation * 100.0;
		atmos.lowDensity = atmos.densityDeviation * 100.0;
		atmos.highDensity = atmos.highDensityDeviation * 100.0;
		atmos.perturbedDensity = atmos.perturbedDensityDeviation * 100.0;
		if (atmos.referencePressure <= 0.0) {
			atmos.pressure = -99.9;
		}
		else {
			atmos.pressure = 100.0 * (atmos.pressure - atmos.referencePressure) / atmos.referencePressure;
		}
    densUnits = "kg/m^3";
		break;
	}
	case KM:
		atmos.density = 1.0e9 * atmos.density;
		atmos.lowDensity = 1.0e9 * atmos.lowDensity;
		atmos.highDensity = 1.0e9 * atmos.highDensity;
		atmos.perturbedDensity = 1.0e9 * atmos.perturbedDensity;
		densUnits = "kg/km^3";
		break;
	default:
		break;
	}
}

//! \brief Header method for the GRAM_CSV_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramCSVHeader(const PerturbedAtmosphere& atmos)
{
  hasVerticalWinds = atmos.hasVerticalWindsInModel();
  gramCSVFile << "ElapsedTime_s,Height_km,Latitude_deg,Longitude" << eastWestChar()
    << "_deg,TotalRadius_km,LatitudeRadius_km,Gravity_ms2,"
    << "Temperature_K,Pressure_Pa,Density_kgm3,PressureScaleHeight_km,DensityScaleHeight_km,"
    << "SpeedOfSound_ms,PressureAtSurface_Pa,SigmaLevel,PressureAltitude_km,"
    << "ReferenceTemperature_K,ReferencePressure_Pa,ReferenceDensity_kgm3,ProfileWeight,"
    << "LowDensity_kgm3,HighDensity_kgm3,PerturbedDensity_kgm3,DensityPerturbation_pct,"
    << "DensityStandardDeviation_pct,PerturbedSpeedOfSound_ms,RelativeStepSize,"
    << "DensityDeviation_pct,LowDensityDeviation_pct,HighDensityDeviation_pct,PerturbedDensityDeviation_pct,"
    << "EWWind_ms,NSWind_ms,";
  if (hasVerticalWinds) {
    gramCSVFile << "VerticalWind_ms,";
  }
  gramCSVFile << "EWWindPerturbation_ms,NSWindPerturbation_ms,";
  if (hasVerticalWinds) {
    gramCSVFile << "VerticalWindPerturbation_ms,";
  }
  gramCSVFile << "PerturbedEWWind_ms,PerturbedNSWind_ms,";
  if (hasVerticalWinds) {
    gramCSVFile << "PerturbedVerticalWind_ms,";
  }
  gramCSVFile << "EWStandardDeviation_ms,NSStandardDeviation_ms,";
  if (hasVerticalWinds) {
    gramCSVFile << "VerticalStandardDeviation_ms,";
  }
  gramCSVFile << "LongitudeOfTheSun_deg,SubsolarLatitude_deg,SubsolarLongitude" << eastWestChar() << "_deg,LocalSolarTime_hr,"
    << "SolarZenithAngle_deg,OneWayLightTime_min,OrbitalRadius_AU,SecondsPerSol,"
    << "TotalNumberDensity_m3,SpecificGasConstant_JkgK,SpecificHeatRatio,AverageMolecularWeight,CompressibilityFactor,";
  printGasHeaders(gramCSVFile, atmos);

  // Print planet specific headers
  printPlanetCSVHeader(atmos);

  gramCSVFile << newline;
}

//! \brief Print method for the GRAM_CSV_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramCSVStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const AtmosphereState& atm = data.atmos;
  const EphemerisState& eph = data.ephem;

  char buffer[bufferSize];
  buffer[0] = 0;
  int iextra = (int)extra;
  int marker = 0;

  // Positional data
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, pos.elapsedTime,
    3 + iextra, pos.height,
    2 + iextra, pos.latitude,
    2 + iextra, pos.getLongitude(eastLongitudePositive),
    3 + iextra, pos.totalRadius,
    3 + iextra, pos.latitudeRadius,
    3 + iextra, pos.gravity);

  // Dynamics
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*e,%.*e,%.*f,%.*f,%.*f,%.*e,",
    1 + iextra, atm.temperature,
    3 + iextra, atm.pressure,
    3 + iextra, atm.density,
    3 + iextra, atm.pressureScaleHeight,
    3 + iextra, atm.densityScaleHeight,
    1 + iextra, atm.speedOfSound,
    3 + iextra, atm.pressureAtSurface);
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,%.*f,%.*e,%.*e,%.*f,",
    3 + iextra, atm.sigmaLevel,
    3 + iextra, atm.pressureAltitude,
    1 + iextra, atm.referenceTemperature,
    3 + iextra, atm.referencePressure,
    3 + iextra, atm.referenceDensity,
    1 + iextra, atm.profileWeight[0]);

  // Density Perturbations
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*e,%.*e,%.*e,%.*f,%.*f,%.*f,",
    3 + iextra, atm.lowDensity,
    3 + iextra, atm.highDensity,
    3 + iextra, atm.perturbedDensity,
    1 + iextra, atm.densityPerturbation * 100.0,
    2 + iextra, atm.densityStandardDeviation * 100.0,
    1 + iextra, atm.perturbedSpeedOfSound);
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    2 + iextra, atm.relativeStepSize,
    1 + iextra, atm.densityDeviation * 100.0,
    1 + iextra, atm.lowDensityDeviation * 100.0,
    1 + iextra, atm.highDensityDeviation * 100.0,
    1 + iextra, atm.perturbedDensityDeviation * 100.0,

    //Winds
    1 + iextra, atm.ewWind,
    1 + iextra, atm.nsWind);
  if (hasVerticalWinds) {
    marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,", 3 + iextra, atm.verticalWind);
  }
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,",
    1 + iextra, atm.ewWindPerturbation,
    1 + iextra, atm.nsWindPerturbation);
  if (hasVerticalWinds) {
    marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,", 3 + iextra, atm.verticalWindPerturbation);
  }
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,",
    1 + iextra, atm.perturbedEWWind,
    1 + iextra, atm.perturbedNSWind);
  if (hasVerticalWinds) {
    marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,", 3 + iextra, atm.perturbedVerticalWind);
  }
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,",
    1 + iextra, atm.ewStandardDeviation,
    1 + iextra, atm.nsStandardDeviation);
  if (hasVerticalWinds) {
    marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,", 3 + iextra, atm.verticalStandardDeviation);
  }

  // Ephemeris
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    2 + iextra, eph.longitudeSun,
    2 + iextra, eph.subsolarLatitude,
    2 + iextra, eph.getSubsolarLongitude(eastLongitudePositive),
    2 + iextra, eph.solarTime,
    2 + iextra, eph.solarZenithAngle,
    2 + iextra, eph.oneWayLightTime,
    2 + iextra, eph.orbitalRadius,
    1 + iextra, eph.secondsPerSol);

  // Gases
  marker += snprintf(buffer + marker, bufferSize - marker, "%.*e,%.*e,%.*e,%.*f,%.*f,",
    4 + iextra, atm.totalNumberDensity,
    3 + iextra, atm.specificGasConstant,
    3 + iextra, atm.specificHeatRatio,
    3 + iextra, atm.averageMolecularWeight,
    4 + iextra, atm.compressibilityFactor);
  gramCSVFile << buffer;

  printCSVGasData(atm.argon);
  printCSVGasData(atm.carbonDioxide);
  printCSVGasData(atm.carbonMonoxide);
  printCSVGasData(atm.dihydrogen);
  printCSVGasData(atm.dinitrogen);
  printCSVGasData(atm.dioxygen);
  printCSVGasData(atm.helium);
  printCSVGasData(atm.hydrogen);
  printCSVGasData(atm.methane);
  printCSVGasData(atm.nitrogen);
  printCSVGasData(atm.oxygen);
  printCSVGasData(atm.ozone);
  printCSVGasData(atm.nitrousOxide);
  printCSVGasData(atm.water);

  // Print planet specific values
  printPlanetCSVStyle(data);

  gramCSVFile << newline;
}

void ProfilePrinter::printCSVGasData( const ConstituentGas& gas)
{
  if (gas.isPresent) {
    char buffer[bufferSize];
    buffer[0] = 0;
    int iextra = (int)extra;
    snprintf(buffer, bufferSize, "%.*e,%.*f,%.*f,%.*f,",
     4 + iextra, gas.numberDensity,
     2 + iextra, gas.massFraction * 100.0,
     2 + iextra, gas.moleFraction * 100.0,
     2 + iextra, gas.averageMolecularWeight);
    gramCSVFile << buffer;
  }
}

//! \brief Header method for the GRAM_MD_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;

  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %-32s | %-26s | %-11s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.19s|%.34s|%.28s|%.13s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %-32s | %-26s | %-11d |\n",
    "Time Frame",
    ((timeFrame == ERT || timeFrame == EARTH_RECEIVE_TIME) ? "Earth Receive Time (ERT)" : "Planet Event Time (PET)"),
    "Initial Random Seed", atmos.getSeed());

  static const char *timeScaleString[3] = { "Coordinated Universal Time (UTC)",
                                            "Barycentric Dynamical Time (TDB)",
                                            "Terrestrial Dynamical Time (TDT)" };
  int timeScaleIndex = 0;
  switch (timeScale) {
  case COORDINATED_UNIVERSAL_TIME:
  case UTC:
    timeScaleIndex = 0;
    break;
  case BARYCENTRIC_DYNAMICAL_TIME:
  case TDB:
    timeScaleIndex = 1;
    break;
  case TERRESTRIAL_DYNAMICAL_TIME:
  case TDT:
    timeScaleIndex = 2;
    break;
  }

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %-32s | %-26s | %-11.3f |\n",
    "Time Scale", timeScaleString[timeScaleIndex],
    "Minimum Relative Step Size", atmos.getMinRelativeStepSize());

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %d/%d/%-*d | %-26s | %-11.2f |\n",
    "Start Date", month, day, 30 - (month < 10 ? 1 : 2) - (day < 10 ? 1 : 2), year,
    "Density Perturbation Scale", atmos.getDensityPerturbationScale());

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %02d:%02d:%05.2f                      | %-26s | %-11.2f |\n",
    "Start Time", hour, minute, seconds,
    "EW Wind Perturbation Scale", atmos.getEWWindPerturbationScale());
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-17s | %-32.6f | %-26s | %-11.2f |\n",
    "Julian Day", jdate,
    "NS Wind Perturbation Scale", atmos.getNSWindPerturbationScale());

  gramMDFile << buffer;

  printPlanetListMDHeader(atmos);

  gramMDFile << newline;
}

void ProfilePrinter::printHeightMDStyle(const ProfileData& data)
{
  char buffer[bufferSize];
  buffer[0] = 0;
  const Position& pos = data.position;
  snprintf(buffer, bufferSize, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
    "Height Above Ref. Ellipsoid (km)", pos.height,
    "Reference Radius (km)", pos.latitudeRadius);
  gramMDFile << buffer;
}

//! \brief Print method for the GRAM_MD_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramListMDStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "--------------------\n");
  marker += snprintf(buffer + marker, bufferSize - marker, "## Record #%zu\n", mdCount++);
  marker += snprintf(buffer + marker, bufferSize - marker, "--------------------\n\n");

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10s | %-33s | %-10s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.35s|%.12s|%.35s|%.12s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.2f |\n",
    "Elapsed Time (s)", pos.elapsedTime,
    "Elapsed Time (sols)", pos.elapsedTime / max(eph.secondsPerSol, (greal)1.0));

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  printHeightMDStyle(data);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Latitude (deg)", pos.latitude,
    "Local Solar Time (hrs)", eph.solarTime);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %10s%1c%-22s | %-10.2f | %-33s | %-10.2f |\n",
    "Longitude ", eastWestChar(), " (deg)", pos.getLongitude(eastLongitudePositive),
    "Longitude of the Sun (deg)", eph.longitudeSun);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Pressure Scale Height (km)", atm.pressureScaleHeight,
    "Orbital Radius (AU)", eph.orbitalRadius);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Density Scale Height (km)", atm.densityScaleHeight,
    "One Way Light Time (min)", eph.oneWayLightTime);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.1f | %-33s | %-10.2f |\n",
    "Temperature (K)", atm.temperature,
    "Subsolar Latitude (deg)", eph.subsolarLatitude);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3e | %19s%1c%-13s | %-10.2f |\n",
    "Pressure (Pa)", atm.pressure,
    "Subsolar Longitude ", eastWestChar(), " (deg)", eph.getSubsolarLongitude(eastLongitudePositive));
   
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Sigma Level", atm.sigmaLevel,
    "Solar Zenith Angle (km)", eph.solarZenithAngle);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Pressure Altitude (km)", atm.pressureAltitude,
    "Gravity (m/s^2)", pos.gravity);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.3f |\n",
    "Surface Pressure (Pa)", atm.pressureAtSurface,
    "Speed of Sound (m/s)", atm.speedOfSound);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.4f | %-33s | %-10.3f |\n",
    "Compressibility Factor (zeta)", atm.compressibilityFactor,
    "Specific Gas Constant (J/(kg K))", atm.specificGasConstant);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Specific Heat Ratio", atm.specificHeatRatio,
    "Profile Weight", atm.profileWeight[0]);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  // Override this function to append to the table above or to add tables.
  printPlanetListMDStyle(data);

  gramMDFile << newline;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10s | %-33s | %-10s |\n",
    "Density", "Low", "Average", "High");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.35s|%.12s|%.35s|%.12s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.4e | %-33.4e | %-10.4e |\n",
    ("Density (" + densUnits + ")").c_str(), 
    atm.lowDensity, atm.density, atm.highDensity);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.1f | %-33.1f | %-10.1f |\n",
    "Density Deviation (%)",
    atm.lowDensityDeviation * 100.0,
    atm.densityDeviation * 100.0,
    atm.highDensityDeviation * 100.0);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.4e | %-33s | %-10.1f |\n",
    ("Perturbed Density (" + densUnits + ")").c_str(), atm.perturbedDensity,
    "Perturbation (%)", atm.densityPerturbation * 100.0);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.2f |\n\n",
    "Perturbed Density Deviation (%)", atm.perturbedDensityDeviation * 100.0,
    "Perturbed Speed of Sound (m/s)", atm.perturbedSpeedOfSound);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-13s | %-13s | %-13s |\n",
    "Winds", "Mean", "Perturbation", "Perturbed");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.35s|%.15s|%.15s|%.15s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-13.1f | %-13.1f | %-13.1f |\n",
    "Eastward Wind  (m/s)", atm.ewWind, atm.ewWindPerturbation, atm.perturbedEWWind);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-13.1f | %-13.1f | %-13.1f |\n\n",
    "Northward Wind  (m/s)", atm.nsWind, atm.nsWindPerturbation, atm.perturbedNSWind);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-20s | %-22s | %-8s | %-8s | %-11s |\n",
    "Gases", "Number Density (#/m^3)", "Mass (%)", "Mole (%)", "Avg Mol Wgt");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.22s|%.24s|%.10s|%.10s|%.13s|\n",
    dashes, dashes, dashes, dashes, dashes);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  massTotal = 0.0;
  moleTotal = 0.0;

  printListGasData(atm.argon, "Argon (Ar)");
  printListGasData(atm.carbonDioxide, "Carbon Dioxide (CO2)");
  printListGasData(atm.carbonMonoxide, "Carbon Monoxide (CO)");
  printListGasData(atm.dihydrogen, "Dihydrogen (H2)");
  printListGasData(atm.dinitrogen, "Dinitrogen (N2)");
  printListGasData(atm.dioxygen, "Dioxygen (O2)");
  printListGasData(atm.helium, "Helium (He)");
  printListGasData(atm.hydrogen, "Hydrogen (H)");
  printListGasData(atm.methane, "Methane (CH4)");
  printListGasData(atm.nitrogen, "Nitrogen (N)");
  printListGasData(atm.oxygen, "Oxygen (O)");
  printListGasData(atm.ozone, "Ozone (O3)");
  printListGasData(atm.nitrousOxide, "Nitrous Oxide (N2O)");
  printListGasData(atm.water, "Water (H2O)");

  snprintf(buffer, bufferSize, "| %-20s | %-22.4e | %-8.1f | %-8.1f | %-11.2f |\n\n",
    "Total",
    atm.totalNumberDensity,
    massTotal * 100.0,
    moleTotal * 100.0,
    atm.averageMolecularWeight);

  gramMDFile << buffer;
}

void ProfilePrinter::printListGasData(const ConstituentGas& gas, const char* name)
{
  if (gas.isPresent) {
    char buffer[bufferSize];
    buffer[0] = 0;
    snprintf(buffer, bufferSize, "| %-20s | %-22.4e | %-8.1f | %-8.1f | %-11.2f |\n",
      name,
      gas.numberDensity,
      gas.massFraction * 100.0,
      gas.moleFraction * 100.0,
      gas.averageMolecularWeight);

    gramMDFile << buffer;

    massTotal += gas.massFraction;
    moleTotal += gas.moleFraction;
  }
}

//! \brief Header method for the GRAM_COLUMN_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
  gramColumnFile << "     Time   Height   Lat   Lon"
    << eastWestChar() << "  Denkgm3    Temp  EWind  NWind  sigD   Ls ";
  printGasHeaders(gramColumnFile, atmos, "%m");
  gramColumnFile << newline;
}

//! \brief Print method for the GRAM_COLUMN_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramColumnStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const AtmosphereState& atm = data.atmos;
  const EphemerisState& eph = data.ephem;
  gramColumnFile << fixed << setw(10) << setprecision(1) << pos.elapsedTime
    << setw(8) << setprecision(2) << pos.height
    << setw(7) << setprecision(2) << pos.latitude
    << setw(7) << setprecision(2) << pos.getLongitude(eastLongitudePositive)
    << scientific << setw(9) << setprecision(2) << atm.density
    << fixed << setw(7) << setprecision(1) << atm.temperature
    << setw(7) << setprecision(1) << atm.ewWind
    << setw(7) << setprecision(1) << atm.nsWind
    << setw(6) << setprecision(1) << atm.densityStandardDeviation * 100.0
    << setw(6) << setprecision(1) << eph.longitudeSun;
  if (atm.argon.isPresent)          gramColumnFile << setw(6) << setprecision(1) << atm.argon.massFraction * 100.0;
  if (atm.carbonDioxide.isPresent)  gramColumnFile << setw(6) << setprecision(1) << atm.carbonDioxide.massFraction * 100.0;
  if (atm.carbonMonoxide.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.carbonMonoxide.massFraction * 100.0;
  if (atm.dihydrogen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.dihydrogen.massFraction * 100.0;
  if (atm.dinitrogen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.dinitrogen.massFraction * 100.0;
  if (atm.dioxygen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.dioxygen.massFraction * 100.0;
  if (atm.helium.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.helium.massFraction * 100.0;
  if (atm.hydrogen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.hydrogen.massFraction * 100.0;
  if (atm.methane.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.methane.massFraction * 100.0;
  if (atm.nitrogen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.nitrogen.massFraction * 100.0;
  if (atm.oxygen.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.oxygen.massFraction * 100.0;
  if (atm.ozone.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.ozone.massFraction * 100.0;
  if (atm.nitrousOxide.isPresent) gramColumnFile << setw(6) << setprecision(1) << atm.nitrousOxide.massFraction * 100.0;
  gramColumnFile << newline;
}

//! \brief Header method for the GRAM_LIST_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
	const GramTime& time = atmos.getStartTime();
	GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
	GRAM_TIME_SCALE timeScale = time.getTimeScale();
	int year, month, day, hour, minute;
	double seconds;
	double jdate;
	time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
	time.getStartTime(timeScale, timeFrame, jdate);
	gramListFile << " " << atmos.getVersionString() << newline;
  gramListFile << " Input Time Frame is " << (timeFrame == ERT ? "Earth Receive Time (ERT)" : "Planet Event Time (PET)") << newline;
  gramListFile << " Input Time Scale is ";
	switch (timeScale) {
	case COORDINATED_UNIVERSAL_TIME:
	case UTC:
    gramListFile << "Coordinated Universal Time (UTC)";
		break;
	case BARYCENTRIC_DYNAMICAL_TIME:
	case TDB:
    gramListFile << "Barycentric Dynamical Time (TDB)";
		break;
	case TERRESTRIAL_DYNAMICAL_TIME:
	case TDT:
    gramListFile << "Terrestrial Dynamical Time (TDT)";
		break;
	}
  gramListFile << newline;
  gramListFile << " Date = " << month << "/" << day << "/" << year
		<< "  " << setw(2) << setfill('0') << hour << ":" << setw(2) << minute << ":" << setw(2) << rint(seconds) << setfill(' ')
		<< "   Julian Day " << fixed << jdate << newline;
  gramListFile << " Minimum Relative Step Size = " << atmos.getMinRelativeStepSize() << newline;
}

//! \brief Print method for the GRAM_LIST_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramListStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  gramListFile << " Time (rel. to T0) ="
    << fixed << setw(11) << setprecision(1) << pos.elapsedTime
    << " sec. (" << setw(7) << setprecision(3) << pos.elapsedTime / max(eph.secondsPerSol, (greal)1.0)
    << " sols)    Ls =" << setw(6) << setprecision(1) << eph.longitudeSun << newline;
  gramListFile << " Height Above Reference Ellipsoid =" << setw(9) << setprecision(3) << pos.height
    << " km          SZA =" << setw(7) << setprecision(2) << data.ephem.solarZenithAngle
    << " deg" << newline;
  gramListFile << " Total Radius (Ref Radius) = " << setw(9) << setprecision(3) << pos.totalRadius
    << " (" << setw(9) << setprecision(3) << pos.latitudeRadius
    << ") km     OWLT =" << setw(7) << setprecision(2) << eph.oneWayLightTime // RREF
    << " Min" << newline;
  gramListFile << " Scale Heights: H(p) =" << setw(6) << setprecision(2) << atm.pressureScaleHeight
    << "       H(rho) =" << setw(6) << setprecision(2) << atm.densityScaleHeight
    << " km    zeta =" << setw(7) << setprecision(4) << data.atmos.compressibilityFactor << newline;
  gramListFile << " Latitude = " << setw(7) << setprecision(2) << pos.latitude
    << "  degrees       Longitude = " << setw(7) << setprecision(2) << pos.longitude
    << " E (" << setw(7) << setprecision(2) << pos.getLongitudeDegrees(WEST_POSITIVE)
    << " W) deg." << newline;
  gramListFile << " Sun Latitude = " << setw(8) << setprecision(2) << eph.subsolarLatitude
    << " deg.      Orbital Radius  =" << setw(7) << setprecision(4) << eph.orbitalRadius
    << " AU" << newline;
  gramListFile << " Sun Longitude =" << setw(8) << setprecision(2) << eph.getSubsolarLongitude(eastLongitudePositive)
    << " " << eastWestChar() << " deg."  << "    Local True Solar Time =" << setw(6) << setprecision(2) << eph.solarTime
    << " Local hr" << newline;
  gramListFile << " Temperature = " << setw(7) << setprecision(1) << atm.temperature
    << " K      Pressure = " << scientific << setw(10) << setprecision(3) << atm.pressure
    << " N/m**2   profwgt =" << fixed << setw(6) << setprecision(3) << atm.profileWeight[0] << newline;
  gramListFile << " Density (Low, Avg., High) = " << scientific << setw(12) << setprecision(3) << atm.lowDensity
    << setw(12) << setprecision(3) << atm.density
    << setw(12) << setprecision(3) << atm.highDensity
    << " " << densUnits << newline;
  gramListFile << " Departure from Avg =     " << fixed << setw(10) << setprecision(1) << atm.lowDensityDeviation * 100.0
    << " %" << setw(10) << setprecision(1) << atm.densityDeviation * 100.0
    << " %" << setw(10) << setprecision(1) << atm.highDensityDeviation * 100.0
    << " %" << newline;
  gramListFile << " Tot.Dens. =" << scientific << setw(10) << setprecision(3) << atm.perturbedDensity
    << " " << densUnits <<
    "   Dens.Pert. =" << fixed << setw(7) << setprecision(2) << atm.densityPerturbation * 100.0
    << " % of mean  iupdate= " << iupdate << newline;
  gramListFile << " Eastward Wind  (Mean, Perturbed, Total) = " << setw(7) << setprecision(1) << atm.ewWind
    << setw(7) << setprecision(1) << atm.ewWindPerturbation
    << setw(7) << setprecision(1) << atm.perturbedEWWind
    << " m/s" << newline;
  gramListFile << " Northward Wind (Mean, Perturbed, Total) = " << setw(7) << setprecision(1) << atm.nsWind
    << setw(7) << setprecision(1) << atm.nsWindPerturbation
    << setw(7) << setprecision(1) << atm.perturbedNSWind
    << " m/s" << newline;
  gramListFile << "    AR          CO2           CO           H2           N2\n"
    << scientific
    << setw(11) << setprecision(3) << atm.argon.numberDensity
    << setw(13) << setprecision(3) << atm.carbonDioxide.numberDensity
    << setw(13) << setprecision(3) << atm.carbonMonoxide.numberDensity
    << setw(13) << setprecision(3) << atm.dihydrogen.numberDensity
    << setw(13) << setprecision(3) << atm.dinitrogen.numberDensity
    << setw(13) << setprecision(3) << atm.dioxygen.numberDensity
    << " #/m**3" << newline;
  gramListFile << fixed
    << setw( 7) << setprecision(3) << atm.argon.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.carbonDioxide.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.carbonMonoxide.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.dihydrogen.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.dinitrogen.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.dioxygen.massFraction * 100.0
    << "     % by mass" << newline;
  gramListFile << fixed
    << setw( 7) << setprecision(3) << atm.argon.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.carbonDioxide.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.carbonMonoxide.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.dihydrogen.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.dinitrogen.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.dioxygen.moleFraction * 100.0
    << "     % by volume" << newline;
  gramListFile << "    HE           H           CH4           N            O\n" << scientific
    << setw(11) << setprecision(3) << atm.helium.numberDensity
    << setw(13) << setprecision(3) << atm.hydrogen.numberDensity
    << setw(13) << setprecision(3) << atm.methane.numberDensity
    << setw(13) << setprecision(3) << atm.nitrogen.numberDensity
    << setw(13) << setprecision(3) << atm.oxygen.numberDensity
    << " #/m**3" << newline;
  gramListFile << fixed
    << setw( 7) << setprecision(3) << atm.helium.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.hydrogen.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.methane.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.nitrogen.massFraction * 100.0
    << setw(13) << setprecision(3) << atm.oxygen.massFraction * 100.0
    << "     % by mass" << newline;
  gramListFile << fixed
    << setw( 7) << setprecision(3) << atm.helium.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.hydrogen.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.methane.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.nitrogen.moleFraction * 100.0
    << setw(13) << setprecision(3) << atm.oxygen.moleFraction * 100.0
    << "     % by volume" << newline;
  gramListFile << " Total = " << scientific << setw(10) << setprecision(3) << atm.totalNumberDensity
    << " #/m**3     MolWgt = " << setw(6) << setprecision(3) << atm.averageMolecularWeight << newline;
  gramListFile << " ---------------------------------------------------------------------------\n";
}

//! \brief Header method for the GRAM_EPHEMERIS_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramEphemerisHeader(const PerturbedAtmosphere& atmos)
{
  gramEphemerisFile << "     Time   Height    Lat    Lon"
		<< eastWestChar() << "    Sslat  Sslon"
		<< eastWestChar() << "  LonSun   LST   OWLT      ORad\n";
}

//! \brief Print method for the GRAM_EPHEMERIS_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramEphemerisStyle(const ProfileData& data)
{
  const Position& pos = data.position;
//  const AtmosphereState& atm = data.atmos;
  const EphemerisState& eph = data.ephem;
  gramEphemerisFile << fixed << setw(10) << setprecision(1) << pos.elapsedTime
    << setw(8) << setprecision(2) << pos.height
    << setw(8) << setprecision(2) << pos.latitude
    << setw(8) << setprecision(2) << pos.getLongitude(eastLongitudePositive)
    << setw(8) << setprecision(2) << eph.subsolarLatitude
    << setw(8) << setprecision(2) << eph.getSubsolarLongitude(eastLongitudePositive)
    << setw(8) << setprecision(2) << eph.longitudeSun
    << setw(7) << setprecision(2) << eph.solarTime
    << setw(9) << setprecision(4) << eph.oneWayLightTime
    << setw(9) << setprecision(4) << eph.orbitalRadius
    << newline;
}

//! \brief Header method for the GRAM_DENSITY_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramDensityHeader(const PerturbedAtmosphere& atmos)
{
  gramDensityFile << "   Height        DENSLO     DENSAV     DENSHI    DENSTOT   Radius     Rref   Grav LOGSCALE profwgt\n";
}

//! \brief Print method for the GRAM_DENSITY_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramDensityStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;
  scaleDensityForOutput(atm);

  greal profwgt = data.atmos.profileWeight[0];
  if (profwgt < -999.0) profwgt = 0.0;

  gramDensityFile << fixed << setw(13) << setprecision(3) << data.position.height
    << scientific << setw(11) << setprecision(3) << atm.lowDensity
    << setw(11) << setprecision(3) << atm.density
    << setw(11) << setprecision(3) << atm.highDensity
    << setw(11) << setprecision(3) << atm.perturbedDensity
    << fixed << setw(9) << setprecision(2) << pos.totalRadius
    << setw(9) << setprecision(2) << pos.latitudeRadius
    << setw(7) << setprecision(3) << pos.gravity
    << setw(5) << densityPrintScale
    << setw(9) << setprecision(3) << profwgt
    << newline;
}

//! \brief Header method for the GRAM_PERTURB_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramPerturbHeader(const PerturbedAtmosphere& atmos)
{
  gramPerturbFile << "   Height      SigD  DensRand    corlim   SigU   SigV iupdate\n";
}

//! \brief Print method for the GRAM_PERTURB_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramPerturbStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  gramPerturbFile << fixed << setw(13) << setprecision(3) << data.position.height
    << setw(6) << setprecision(2) << atm.densityStandardDeviation * 100.0
    << setw(10) << setprecision(3) << atm.densityPerturbation * 100.0
    << scientific << setw(10) << setprecision(3) << data.atmos.relativeStepSize
    << fixed << setw(7) << setprecision(2) << atm.ewStandardDeviation
    << setw(7) << setprecision(2) << atm.nsStandardDeviation
    << setw(5) << iupdate
    << newline;
}

//! \brief Header method for the GRAM_WINDS_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramWindsHeader(const PerturbedAtmosphere& atmos)
{
  gramWindsFile << "   Height      EWmean  EWpert   EWtot  NSmean  NSpert   NStot iupdate\n";
}

//! \brief Print method for the GRAM_WINDS_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramWindsStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  gramWindsFile << fixed << setw(13) << setprecision(3) << data.position.height
    << setw(8) << setprecision(2) << atm.ewWind
    << setw(8) << setprecision(2) << atm.ewWindPerturbation
    << setw(8) << setprecision(2) << atm.perturbedEWWind
    << setw(8) << setprecision(2) << atm.nsWind
    << setw(8) << setprecision(2) << atm.nsWindPerturbation
    << setw(8) << setprecision(2) << atm.perturbedNSWind
    << setw(5) << iupdate
    << newline;
}

//! \brief Header method for the GRAM_TPRESHGT_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramTPresHgtHeader(const PerturbedAtmosphere& atmos)
{
  gramTPresHgtFile << "   Height       Temp      Pres   TdegC   Pres_mb    Hrho  Hpres MolWt";
  printGasHeaders(gramTPresHgtFile, atmos, "%v");
  gramTPresHgtFile << " LOGSCALE\n";
}

//! \brief Print method for the GRAM_TPRESHGT_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramTPresHgtStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);

  gramTPresHgtFile << fixed << setw(13) << setprecision(3) << data.position.height
                  << setw(7)  << setprecision(1) << atm.temperature
    << scientific << setw(11) << setprecision(3) << atm.pressure
    << fixed      << setw(7)  << setprecision(1) << atm.temperature - 273.15
    << scientific << setw(11) << setprecision(3) << atm.pressure * 0.01
    << fixed      << setw(7)  << setprecision(2) << atm.densityScaleHeight
                  << setw(7)  << setprecision(2) << atm.pressureScaleHeight
                  << setw(6)  << setprecision(2) << atm.averageMolecularWeight;
  if (atm.argon.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.argon.moleFraction * 100.0;
  if (atm.carbonDioxide.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.carbonDioxide.moleFraction * 100.0;
  if (atm.carbonMonoxide.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.carbonMonoxide.moleFraction * 100.0;
  if (atm.dihydrogen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.dihydrogen.moleFraction * 100.0;
  if (atm.dinitrogen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.dinitrogen.moleFraction * 100.0;
  if (atm.dioxygen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.dioxygen.moleFraction * 100.0;
  if (atm.helium.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.helium.moleFraction * 100.0;
  if (atm.hydrogen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.hydrogen.moleFraction * 100.0;
  if (atm.methane.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.methane.moleFraction * 100.0;
  if (atm.nitrogen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.nitrogen.moleFraction * 100.0;
  if (atm.oxygen.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.oxygen.moleFraction * 100.0;
  if (atm.ozone.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.ozone.moleFraction * 100.0;
  if (atm.nitrousOxide.isPresent) gramTPresHgtFile << setw(6) << setprecision(1) << atm.nitrousOxide.moleFraction * 100.0;
  gramTPresHgtFile <<  setw(5)  << densityPrintScale << newline;
}

//! \brief Header method for the GRAM_SOUND_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGramSoundHeader(const PerturbedAtmosphere& atmos)
{
  gramSoundFile << "     Time   Height   Lat   Lon"
    << eastWestChar() << "     Temp     Pres     Denkgm3    SOS     PSOS ";
  printGasHeaders(gramSoundFile, atmos, "%m");
  gramSoundFile << newline;
}

//! \brief Print method for the GRAM_SOUND_STYLE.
//!
//! \param data The ProfileData to be printed.
void ProfilePrinter::printGramSoundStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const AtmosphereState& atm = data.atmos;
  gramSoundFile << fixed << setw(10) << setprecision(1) << pos.elapsedTime
    << setw(8) << setprecision(2) << pos.height
    << setw(7) << setprecision(2) << pos.latitude
    << setw(7) << setprecision(2) << pos.getLongitude(eastLongitudePositive)
    << fixed << setw(10) << setprecision(3) << atm.temperature
    << scientific << setw(10) << setprecision(3) << atm.pressure
    << scientific << setw(10) << setprecision(3) << atm.density
    << fixed << setw(8) << setprecision(1) << atm.speedOfSound
  << fixed << setw(8) << setprecision(1) << atm.perturbedSpeedOfSound;
  if (atm.argon.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.argon.massFraction * 100.0;
  if (atm.carbonDioxide.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.carbonDioxide.massFraction * 100.0;
  if (atm.carbonMonoxide.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.carbonMonoxide.massFraction * 100.0;
  if (atm.dihydrogen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.dihydrogen.massFraction * 100.0;
  if (atm.dinitrogen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.dinitrogen.massFraction * 100.0;
  if (atm.dioxygen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.dioxygen.massFraction * 100.0;
  if (atm.helium.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.helium.massFraction * 100.0;
  if (atm.hydrogen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.hydrogen.massFraction * 100.0;
  if (atm.methane.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.methane.massFraction * 100.0;
  if (atm.nitrogen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.nitrogen.massFraction * 100.0;
  if (atm.oxygen.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.oxygen.massFraction * 100.0;
  if (atm.ozone.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.ozone.massFraction * 100.0;
  if (atm.nitrousOxide.isPresent) gramSoundFile << setw(6) << setprecision(1) << atm.nitrousOxide.massFraction * 100.0;
  gramSoundFile << newline;
}

//! \brief A utility method for printing column style Gas headers.
//!
//! This method prints heading information for gases that are present in the specified atmosphere.
//! \param stream The output stream.
//! \param atmos The PerturbedAtmosphere which generated the data.
//! \param units An two character string for the units of the value (%m or %v).
void ProfilePrinter::printGasHeaders(std::ostream& stream, const PerturbedAtmosphere& atmos, const char units[4])
{
  if (atmos.isGasPresent(ARGON))           stream << "  AR" << units;
  if (atmos.isGasPresent(CARBON_DIOXIDE))  stream << " CO2" << units;
  if (atmos.isGasPresent(CARBON_MONOXIDE)) stream << "  CO" << units;
  if (atmos.isGasPresent(DIHYDROGEN))      stream << "  H2" << units;
  if (atmos.isGasPresent(DINITROGEN))      stream << "  N2" << units;
  if (atmos.isGasPresent(DIOXYGEN))        stream << "  O2" << units;
  if (atmos.isGasPresent(HELIUM))          stream << "  He" << units;
  if (atmos.isGasPresent(HYDROGEN))        stream << "   H" << units;
  if (atmos.isGasPresent(METHANE))         stream << " CH4" << units;
  if (atmos.isGasPresent(NITROGEN))        stream << "   N" << units;
  if (atmos.isGasPresent(OXYGEN))          stream << "   O" << units;
  if (atmos.isGasPresent(OZONE))           stream << "  O3" << units;
  if (atmos.isGasPresent(NITROUS_OXIDE))   stream << " N2O" << units;
  if (atmos.isGasPresent(WATER))           stream << " H2O" << units;
}

//! \brief A utility method for printing column style Gas headers.
//!
//! This method prints heading information for gases that are present in the specified atmosphere.
//! \param stream The output stream.
//! \param atmos The PerturbedAtmosphere which generated the data.
void ProfilePrinter::printGasHeaders(std::ostream& stream, const PerturbedAtmosphere& atmos)
{
  if (atmos.isGasPresent(ARGON))           stream << "Arnd_m3,Armass_pct,Armole_pct,Aramw,";
  if (atmos.isGasPresent(CARBON_DIOXIDE))  stream << "CO2nd_m3,CO2mass_pct,CO2mole_pct,CO2amw,";
  if (atmos.isGasPresent(CARBON_MONOXIDE)) stream << "COnd_m3,COmass_pct,COmole_pct,COamw,";
  if (atmos.isGasPresent(DIHYDROGEN))      stream << "H2nd_m3,H2mass_pct,H2mole_pct,H2amw,";
  if (atmos.isGasPresent(DINITROGEN))      stream << "N2nd_m3,N2mass_pct,N2mole_pct,N2amw,";
  if (atmos.isGasPresent(DIOXYGEN))        stream << "O2nd_m3,O2mass_pct,O2mole_pct,O2amw,";
  if (atmos.isGasPresent(HELIUM))          stream << "Hend_m3,Hemass_pct,Hemole_pct,Heamw,";
  if (atmos.isGasPresent(HYDROGEN))        stream << "Hnd_m3,Hmass_pct,Hmole_pct,Hamw,";
  if (atmos.isGasPresent(METHANE))         stream << "CH4nd_m3,CH4mass_pct,CH4mole_pct,CH4amw,";
  if (atmos.isGasPresent(NITROGEN))        stream << "Nnd_m3,Nmass_pct,Nmole_pct,Namw,";
  if (atmos.isGasPresent(OXYGEN))          stream << "Ond_m3,Omass_pct,Omole_pct,Oamw,";
  if (atmos.isGasPresent(OZONE))           stream << "O3nd_m3,O3mass_pct,O3mole_pct,O3amw,";
  if (atmos.isGasPresent(NITROUS_OXIDE))   stream << "N2Ond_m3,N2Omass_pct,N2Omole_pct,N2Oamw,";
  if (atmos.isGasPresent(WATER))           stream << "H2Ond_m3,H2Omass_pct,H2Omole_pct,H2Oamw,";
}

void ProfilePrinter::printDates(GramTime gramTime[3], greal lonSun[3], greal tlst[3])
{
  int year, mon, day, hr, min;
  greal secs;
  string labels[3] = { "@desired Solar Time", "@desired LS" , "@desired Solar Time" };
  gramTime[1].getTime(gramTime[0].getTimeScale(), gramTime[0].getTimeFrame(), year, mon, day, hr, min, secs);

  cout << "                          One Way    Solar   Longitude\n";
  cout << "Date       Time " << gramTime[0].getScaleString() << "_" << gramTime[0].getFrameString()
    << "  Light Time  Time    of the Sun\n";
  double seconds;
  for (int i = 0; i < 3; ++i) {
    gramTime[i].getTime(gramTime[i].getTimeScale(), gramTime[i].getTimeFrame(), year, mon, day, hr, min, seconds);
    cout << setfill('0') << year << "/"
      << setw(2) << mon << "/"
      << setw(2) << day << " "
      << setw(2) << hr << ":"
      << setw(2) << min << ":"
      << fixed << setw(5) << setprecision(2) << seconds << "   "
      << setfill(' ')
      << setw(8) << gramTime[i].getOneWayLightTime() << "  "
      << setw(8) << setprecision(4) << tlst[i] << "   "
      << setw(8) << setprecision(4) << lonSun[i] << "  "
      << labels[i] << newline;
  }
  cout << endl;

}

void ProfilePrinter::appendExtension(std::string& fileName, const std::string& ext)
{
  // Remove the current extension, if any (must match txt, csv, or md).  
  size_t dotPosition = fileName.find_last_of('.');
  if (dotPosition != fileName.npos) {
    string oldExtension = fileName.substr(dotPosition + 1);
    if (oldExtension == "txt" || oldExtension == "csv" || oldExtension == "md") {
      fileName = fileName.substr(0, dotPosition);
    }
  }
  // Append the new extension.
  fileName += ext;
}


//! \fn void ProfilePrinter::setEastLongitudePositive(bool flag = true)
//! \brief Designate the output convention for longitudes.
//!
//! \param flag True for east positive, false for west positive.

//! \fn void ProfilePrinter::setWestLongitudePositive()
//! \brief Designate the output convention for longitudes as west positive.

//! \fn void ProfilePrinter::setStyle(size_t styleFlag)
//! \brief Designate the desired print styles.
//!
//! The selection of print styles is accomplised by performing a bitwise-or, "|",
//! of the style masks.  For example: 
//! \code 
//! size_t styles = GRAM_DENSITY_STYLE | GRAM_PERTURB_STYLE | GRAM_WINDS_STYLE; 
//! printer.setStyle(styles); 
//! \endcode
//! \param styleFlag A bitfield composed of style masks.

//! \fn void ProfilePrinter::setDensityPrintScale(DensityPrintScale scale)
//! \brief Designate the output scale for density values.
//!
//! \param scale A DensityPrintScale option.

//! \fn void ProfilePrinter::setListFileName(const std::string& fileName)
//! \brief Designate the file name to use for List style output.
//!
//! \param fileName A file name.

//! \fn void ProfilePrinter::setColumnFileName(const std::string& fileName)
//! \brief Designate the file name to use for Column style output.
//!
//! \param fileName A file name.

//! \fn void ProfilePrinter::setFileNamePrefix(const std::string& prefix)
//! \brief Designate a prefix that will be prepended to all output file names.
//!
//! \param prefix A prefix string (using characters that form legitimate file names).

//! \fn void ProfilePrinter::eastWestChar(bool flip = false)
//! \brief Prints E/W longitude designation.
//!
//! \param flip False to print E/W, true to print W/E.

} // namespace
