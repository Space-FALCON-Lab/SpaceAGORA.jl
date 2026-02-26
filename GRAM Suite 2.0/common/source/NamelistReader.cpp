//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include <iostream>
#include "NamelistReader.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
NamelistReader::NamelistReader() {}

//! \fn NamelistReader::NamelistReader(const NamelistReader& orig)
//! \brief Copying this object is discouraged.

//! \fn  NamelistReader::~NamelistReader()
//! \copydoc Atmosphere::~Atmosphere()

void NamelistReader::tryGetSpicePath(std::string &spicePath)
{
  string dataPath;
  tryGetSpicePath(spicePath, dataPath);
}

void NamelistReader::tryGetSpicePath(std::string &spicePath, std::string& dataPath)
{
  InputParameters params;
  tryGetSpicePath(params);
  spicePath = params.spicePath;
  dataPath = params.dataPath;
}



void NamelistReader::tryGetSpicePath(InputParameters& params)
{
  params.spicePath = DEFAULT_SPICE_PATH;
  string fileName = "spice.txt";
  ifstream spiceFile(fileName);
  if (spiceFile.good()) {
    spiceFile.close();

    // Load up the itemMap from the namelist file.
    readFile(fileName);

    // Get the common input parameters
    getItem("SPICEDIR", params.spicePath);
    getItem("SPICEPATH", params.spicePath);
    processPath(params.spicePath);

    getItem("DATADIR", params.dataPath);
    getItem("DATAPATH", params.dataPath);
    processPath(params.dataPath);

    getItem("SPICELSK", params.spiceLsk);
    processPath(params.spiceLsk);
    getItem("SPICEPCK", params.spicePck);
    processPath(params.spicePck);
    getItem("SPICEVENUS", params.spiceVenus);
    processPath(params.spiceVenus);
    getItem("SPICEEARTH", params.spiceEarth);
    processPath(params.spiceEarth);
    getItem("SPICEMARS", params.spiceMars);
    processPath(params.spiceMars);
    getItem("SPICEJUPITER", params.spiceJupiter);
    processPath(params.spiceJupiter);
    getItem("SPICESATURN", params.spiceSaturn);
    processPath(params.spiceSaturn);
    getItem("SPICEURANUS", params.spiceUranus);
    processPath(params.spiceUranus);
    getItem("SPICENEPTUNE", params.spiceNeptune);
    processPath(params.spiceNeptune);
    getItem("SPICETITAN", params.spiceTitan);
    processPath(params.spiceTitan);
  }
}

//! \brief Reads the namelist file into an item map.
//!
//! This method parses the namelist file into name-value pairs where both the name and
//! the value are strings.  The pairs are stored in a map for further parsing by getParameters().
//! \param fileName The full path/name of the namelist file.
void NamelistReader::readFile(const std::string &fileName)
{
  ifstream namelistFile(fileName);
  if (!namelistFile) {
    throw string(FILE_OPEN_ERROR_MESSAGE + fileName);
  }

  try {
    std::string lineBuffer = "";

    // Search for the start of the name list denoted by a $
    while (!namelistFile.eof() && lineBuffer.find("$") == string::npos) {
      getline(namelistFile, lineBuffer);
      if (namelistFile.bad()) {
        throw string("Unable to read the file. Current buffer: ") + string(lineBuffer);
      }
    }

    // Process the name list until the $ (end) is reached
    while (!namelistFile.eof()) {
      getline(namelistFile, lineBuffer);
      if (namelistFile.bad()) {
        throw string("Unable to read the file. Current buffer: ") + string(lineBuffer);
      }
      if (lineBuffer.find("$") != string::npos) {
        break;
      }
      parseLine(lineBuffer);
    }
  }
  catch (string msg) {
    warnings << "Error: Namelist file error." << endl;
    warnings << "       " << msg << endl;
    warnings << "       " << fileName << endl;
    warnings << "       The program will attempt to continue.\n" << endl;
  }
}

void NamelistReader::parseLine(const std::string& lineBuffer)
{
  std::string upperBuffer = lineBuffer;
  std::transform(upperBuffer.begin(), upperBuffer.end(), upperBuffer.begin(), ::toupper);
  // Find the start of the item name
  size_t itemStart = upperBuffer.find_first_not_of(" \t");
  if (itemStart != string::npos) {
    // Find the end of the item name
    size_t itemEnd = upperBuffer.find_first_of("= \t", itemStart + 1);
    if (itemEnd != string::npos) {
      // Save the item name
      string itemName = upperBuffer.substr(itemStart, itemEnd - itemStart);
      // Find the equal sign
      size_t equals = lineBuffer.find("=", itemEnd);
      if (equals != string::npos) {
        // Find the start of the right hand side
        itemStart = lineBuffer.find_first_not_of("' \t", equals + 1);
        if (itemStart != string::npos) {
          // Find the end of the right hand side
          if (upperBuffer[itemStart - 1] == '\'') {
            // End for a quoted string
            itemEnd = upperBuffer.find_first_of("'", itemStart);
          }
          else {
            // End for all other types
            itemEnd = upperBuffer.find_first_of(" \t", itemStart + 1);
          }
          // Extract the right hand side
          string itemValue = lineBuffer.substr(itemStart, itemEnd - itemStart);
          // Add this into the item map
          itemMap[itemName] = itemValue;
        }
        // Trap the case of PARAMETER = '' and replace with null string.
        else if (lineBuffer.substr(equals + 1).find("''") != string::npos) {
          itemMap[itemName] = "null";
        }
        else {
          warnings << "Error: Namelist file error." << endl;
          warnings << "       Cannot find right had side of assignment: " << lineBuffer << endl;
          warnings << "       The program will attempt to continue.\n" << endl;
        }
      }
    }
  }
}

//! \brief Gets input parameters from the parsed item map.
//!
//! This method searches the item map for key strings.  If the strings are found, then the
//! item value is stored in the InputParameters object.
//! \param inputParameters An InputParameters object.
void NamelistReader::getParameters(InputParameters &inputParameters)
{
  // File and path parameters
  getItem("SPICEDIR", inputParameters.spicePath);
  getItem("SPICEPATH", inputParameters.spicePath);
  processPath(inputParameters.spicePath);

  getItem("SPICELSK", inputParameters.spiceLsk);
  processPath(inputParameters.spiceLsk);
  getItem("SPICEPCK", inputParameters.spicePck);
  processPath(inputParameters.spicePck);
  getItem("SPICEVENUS", inputParameters.spiceVenus);
  processPath(inputParameters.spiceVenus);
  getItem("SPICEEARTH", inputParameters.spiceEarth);
  processPath(inputParameters.spiceEarth);
  getItem("SPICEMARS", inputParameters.spiceMars);
  processPath(inputParameters.spiceMars);
  getItem("SPICEJUPITER", inputParameters.spiceJupiter);
  processPath(inputParameters.spiceJupiter);
  getItem("SPICESATURN", inputParameters.spiceSaturn);
  processPath(inputParameters.spiceSaturn);
  getItem("SPICEURANUS", inputParameters.spiceUranus);
  processPath(inputParameters.spiceUranus);
  getItem("SPICENEPTUNE", inputParameters.spiceNeptune);
  processPath(inputParameters.spiceNeptune);
  getItem("SPICETITAN", inputParameters.spiceTitan);
  processPath(inputParameters.spiceTitan);

  getItem("LSTFL", inputParameters.listFileName);
  getItem("LISTFILENAME", inputParameters.listFileName);
  processPath(inputParameters.listFileName);

  getItem("OUTFL", inputParameters.columnFileName);
  getItem("COLUMNFILENAME", inputParameters.columnFileName);
  processPath(inputParameters.columnFileName);

  getItem("TRAJFL", inputParameters.trajectoryFileName);
  getItem("TRAJECTORYFILENAME", inputParameters.trajectoryFileName);
  processPath(inputParameters.trajectoryFileName);

  getItem("DATADIR", inputParameters.dataPath);
  getItem("DATAPATH", inputParameters.dataPath);
  processPath(inputParameters.dataPath);

  // I/O parameters
  getItem("LONEAST", inputParameters.isEastLongitudePositiveOnInput);
  getItem("LONEW", inputParameters.isEastLongitudePositiveOnInput);
  getItem("EASTLONGITUDEPOSITIVE", inputParameters.isEastLongitudePositiveOnInput);
  inputParameters.isEastLongitudePositiveOnOutput = inputParameters.isEastLongitudePositiveOnInput;
  getItem("USELEGACYOUTPUT", inputParameters.useLegacyOutputs);
  getItem("USELEGACYOUTPUTS", inputParameters.useLegacyOutputs);
  int precision = 0;
  getItem("EXTRAPRECISION", precision);
  inputParameters.extraPrecision = static_cast<size_t>(clampSize(precision, 15));

  // Time parameters
  int flag = 0;
  getItem("IERT", flag);
  getItem("TIMEFRAME", flag);
  inputParameters.timeFrame = (flag == 0 ? PET : ERT);

  getItem("IUTC", flag);
  getItem("TIMESCALE", flag);
  if (flag == 0)
    inputParameters.timeScale = TDT;
  else if (flag == 2)
    inputParameters.timeScale = TDB;
  else
    inputParameters.timeScale = UTC;

  getItem("MONTH", inputParameters.month);
  getItem("MDAY", inputParameters.day);
  getItem("DAY", inputParameters.day);
  getItem("MYEAR", inputParameters.year);
  getItem("YEAR", inputParameters.year);
  if (inputParameters.year < 70) {
    inputParameters.year += 2000;
  }
  else if (inputParameters.year < 100) {
    inputParameters.year += 1900;
  }
  getItem("IHOUR", inputParameters.hour);
  getItem("IHR", inputParameters.hour);
  getItem("HOUR", inputParameters.hour);
  getItem("IMIN", inputParameters.minute);
  getItem("MINUTE", inputParameters.minute);
  getItem("SEC", inputParameters.seconds);
  getItem("SECONDS", inputParameters.seconds);

  string msg;
  if (dateHasError(inputParameters, msg)) {
    throw(msg + "\n       This is an unrecoverable error.");
  }

  // Pertubation parameters
  getItem("NR1", inputParameters.initialRandomSeed);
  getItem("INITIALRANDOMSEED", inputParameters.initialRandomSeed);

  greal factor = -1.0;
  getItem("RPSCALE", factor);
  getItem("PERTURBATIONSCALES", factor);
  if (factor >= 0.0) {
    inputParameters.densityPerturbationScale = factor;
    inputParameters.ewWindPerturbationScale = factor;
    inputParameters.nsWindPerturbationScale = factor;
    inputParameters.verticalWindPerturbationScale = factor;
  }
  getItem("DENSITYPERTURBATIONSCALE", inputParameters.densityPerturbationScale);
  getItem("EWWINDPERTURBATIONSCALE", inputParameters.ewWindPerturbationScale);
  getItem("NSWINDPERTURBATIONSCALE", inputParameters.nsWindPerturbationScale);
  getItem("VERTICALWINDPERTURBATIONSCALE", inputParameters.verticalWindPerturbationScale);

  getItem("CORLMIN", inputParameters.minRelativeStepSize);
  getItem("MINRELATIVESTEPSIZE", inputParameters.minRelativeStepSize);
  getItem("MINIMUMRELATIVESTEPSIZE", inputParameters.minRelativeStepSize);

  // Ouput parameters
  getItem("NVARX", inputParameters.nvarx); // not used
  getItem("NVARY", inputParameters.nvary); // not used

  int logScale = 0;
  getItem("LOGSCALE", logScale);
  getItem("DENSITYPRINTSCALE", logScale);
  switch (logScale) {
  case 1:
    inputParameters.densityPrintScale = LOG_10;
    break;
  case 2:
    inputParameters.densityPrintScale = PERCENT_DEV;
    break;
  case 3:
    inputParameters.densityPrintScale = KM;
    break;
  default:
    inputParameters.densityPrintScale = STANDARD;
    break;
  }

  // Input Trajectory parameters
  getItem("NPOS", inputParameters.numberOfPositions);
  getItem("NUMBEROFPOSITIONS", inputParameters.numberOfPositions);

  getItem("FLAT", inputParameters.initialLatitude);
  getItem("INITIALLATITUDE", inputParameters.initialLatitude);

  getItem("FLON", inputParameters.initialLongitude);
  getItem("INITIALLONGITUDE", inputParameters.initialLongitude);

  getItem("FHGT", inputParameters.initialHeight);
  getItem("INITIALHEIGHT", inputParameters.initialHeight);

  getItem("DELLAT", inputParameters.deltaLatitude);
  getItem("DELTALATITUDE", inputParameters.deltaLatitude);

  getItem("DELLON", inputParameters.deltaLongitude);
  getItem("DELTALONGITUDE", inputParameters.deltaLongitude);

  getItem("DELHGT", inputParameters.deltaHeight);
  getItem("DELTAHEIGHT", inputParameters.deltaHeight);

  getItem("DELTIME", inputParameters.deltaTime);
  getItem("DELTATIME", inputParameters.deltaTime);

  getItem("IPCLAT", inputParameters.isPlanetoCentric);
  getItem("ISPLANETOCENTRIC", inputParameters.isPlanetoCentric);
  getItem("ISGEOCENTRIC", inputParameters.isPlanetoCentric);
  getItem("INITIALPERTURBATIONSUPDATED", inputParameters.initialPerturbationsUpdated);

  // Monte Carlo parameters
  getItem("NMONTE", inputParameters.numberOfMonteCarloRuns);
  getItem("NUMBEROFMONTECARLORUNS", inputParameters.numberOfMonteCarloRuns);

  // Auxiliary Profile parameters
  getItem("PROFILE", inputParameters.auxiliaryAtmosphereFileName[0]);
  getItem("AUXILIARYATMOSPHEREFILENAME", inputParameters.auxiliaryAtmosphereFileName[0]);

  getItem("PROFNEAR", inputParameters.innerRadius[0]);
  getItem("INNERRADIUS", inputParameters.innerRadius[0]);

  getItem("PROFFAR", inputParameters.outerRadius[0]);
  getItem("OUTERRADIUS", inputParameters.outerRadius[0]);

  getItem("FASTMODEON", inputParameters.fastModeOn);

  getItem("IUP", inputParameters.iup); // not used

  // FindDates parameters
  getItem("FINDDATES", inputParameters.findDates);
  getItem("TARGETLONGITUDESUN", inputParameters.targetLongitudeSun);
  getItem("TARGETSOLARTIME", inputParameters.targetSolarTime);
  // initialLongitude (see above)
  // all GramTime parameters (see above)
}

//! \brief Checks the parsed item map for unknown entries.
//!
//! This method checks each entry in the item map against the list of valid key names.
//! An exception is thrown if any entries are invalid.
void NamelistReader::checkForErrors()
{
  string invalid;
  sort(validKeyNames.begin(), validKeyNames.end());
  for (auto &item : itemMap) {
    if (!binary_search(validKeyNames.begin(), validKeyNames.end(), item.first)) {
      if (!invalid.empty()) {
        invalid += ", ";
      }
      invalid += item.first;
    }
  }
  if (!invalid.empty()) {
    warnings << "Error: Namelist file error." << endl;
    warnings << "       Invalid items found in namelist file." << endl;
    warnings << "       Invalid items: " << invalid << endl;
    warnings << "       The program will attempt to continue.\n" << endl;
  }

  // If console output is on, emit warnings.
  // The call to tellp() says the warnings are not empty if not at the beginning.
  if (CONSOLE_OUTPUT && warnings.tellp() != 0) {
    cerr << warnings.str() << endl;
  }
}

//! \brief Looks for a key in the nameline item map.
//!
//! This method will search for the name/key of an item in the item map.  If the
//! item is found, then the value string is returned.  Otherwise, an empty string is returned.
//! \returns The value string or an empty string.
const std::string &NamelistReader::findItem(const std::string &itemName)
{
  // Keep an empty string around for "not found" items
  static const string emptyString = "";

  // Search for the item
  auto iter = itemMap.find(itemName);

  // if not found, return an empty string
  if (iter == itemMap.end()) {
    return emptyString;
  }

  // return the found item;
  return iter->second;
}

//! \brief Parses value string for an integer value.
//!
//! If the item name is found in the item map, then this method attempts to parse
//! it into the appropriate format.
//! \param itemName The key name in the item map.
//! \param value The parsed value.
//! \returns True if value was parsed.
bool NamelistReader::getItem(const std::string &itemName, int &value)
{
  validKeyNames.push_back(itemName);
  const std::string &stringValue = findItem(itemName);
  if (!stringValue.empty()) {
    try {
      value = stoi(stringValue);
      return true;
    }
    catch (const invalid_argument &) {
      warnings << "Error: Namelist file error." << endl;
      warnings << "       Invalid integer assignment: " << itemName << " = " << stringValue << endl;
      warnings << "       The program will attempt to continue.\n" << endl;
    }
  }
  return false;
}

//! \brief Parses value string for an double value.
//!
//! If the item name is found in the item map, then this method attempts to parse
//! it into the appropriate format.
//! \param itemName The key name in the item map.
//! \param value The parsed value.
//! \returns True if value was parsed.
bool NamelistReader::getItem(const std::string &itemName, double &value)
{
  validKeyNames.push_back(itemName);
  const std::string &stringValue = findItem(itemName);
  if (!stringValue.empty()) {
    try {
      value = stod(stringValue);
      return true;
    }
    catch (const invalid_argument &) {
      warnings << "Error: Namelist file error." << endl;
      warnings << "       Invalid double assignment: " << itemName << " = " << stringValue << endl;
      warnings << "       The program will attempt to continue.\n" << endl;
    }
  }
  return false;
}

//! \brief Parses value string for an float value.
//!
//! If the item name is found in the item map, then this method attempts to parse
//! it into the appropriate format.
//! \param itemName The key name in the item map.
//! \param value The parsed value.
//! \returns True if value was parsed.
bool NamelistReader::getItem(const std::string &itemName, float &value)
{
  validKeyNames.push_back(itemName);
  const std::string &stringValue = findItem(itemName);
  if (!stringValue.empty()) {
    try {
      value = stof(stringValue);
      return true;
    }
    catch (const invalid_argument &) {
      warnings << "Error: Namelist file error." << endl;
      warnings << "       Invalid float assignment: " << itemName << " = " << stringValue << endl;
      warnings << "       The program will attempt to continue.\n" << endl;
    }
  }
  return false;
}

//! \brief Parses value string for an boolean value.
//!
//! If the item name is found in the item map, then this method attempts to parse
//! it into the appropriate format.
//! \param itemName The key name in the item map.
//! \param value The parsed value.
//! \returns True if value was parsed.
bool NamelistReader::getItem(const std::string &itemName, bool &value)
{
  validKeyNames.push_back(itemName);
  const std::string &stringValue = findItem(itemName);
  if (!stringValue.empty()) {
    try {
      int intValue = stoi(stringValue);
      value = (intValue != 0);
      return true;
    }
    catch (const invalid_argument &) {
      warnings << "Error: Namelist file error." << endl;
      warnings << "       Invalid integer assignment: " << itemName << " = " << stringValue << endl;
      warnings << "       The program will attempt to continue.\n" << endl;
    }
  }
  return false;
}

//! \brief Parses value string for an string value.
//!
//! If the item name is found in the item map, then this method attempts to parse
//! it into the appropriate format.
//! \param itemName The key name in the item map.
//! \param value The parsed value.
//! \returns True if value was parsed.
bool NamelistReader::getItem(const std::string &itemName, std::string &value)
{
  validKeyNames.push_back(itemName);
  const std::string &stringValue = findItem(itemName);
  if (!stringValue.empty()) {
    if (stringValue != "null") {
      value = stringValue;
      return true;
    }
    else {
      value.clear();
    }
  }
  return false;
}

bool NamelistReader::dateHasError(const InputParameters &inputParameters, string msg)
{
  msg = "Error: Namelist file error.\n";
  bool dateError = false;
  // Validate the start time
  int month = inputParameters.month;
  int day = inputParameters.day;
  if (month < 1 || month > 12) {
    msg += "       Invalid month (1-12): " + to_string(month) + "\n";
    dateError = true;
  }
  // Thirty days hath September...
  if (month == 9 || month == 4 || month == 6 || month == 11) {
    if (day < 1 || day > 30) {
      msg += "       Invalid day (1-30): " + to_string(day) + "\n";
      dateError = true;
    }
  }
  // All the rest have 31...
  else if (month != 2) {
    if (day < 1 || day > 31) {
      msg += "       Invalid day (1-31): " + to_string(day) + "\n";
      dateError = true;
    }
  }
  // Excepting February alone....
  else {
    bool leapyear = false;
    int year = inputParameters.year;
    if (year / 4.0 - (int)(year / 4.0) == 0.0) {
      leapyear = true;
      if (year / 400.0 - (int)(year / 400.0) == 0.0) {
        leapyear = false;
      }
    }
    if (day < 1 || (leapyear && day > 29) || (!leapyear && day > 28)) {
      msg += "       Invalid day (1-2";
      msg += (leapyear ? "9" : "8");
      msg += "): " + to_string(day) + "\n";
      dateError = true;
    }
  }
  if (inputParameters.hour < 0 || inputParameters.hour > 23) {
    msg += "       Invalid hour (0-23): " + to_string(inputParameters.hour) + "\n";
    dateError = true;
  }
  if (inputParameters.minute < 0 || inputParameters.minute > 59) {
    msg += "       Invalid minute (0-59): " + to_string(inputParameters.minute) + "\n";
    dateError = true;
  }
  if (inputParameters.seconds < 0.0 || inputParameters.seconds >= 60.0) {
    msg += "       Invalid seconds (0.0 to less than 60.0): " + to_string(inputParameters.seconds) + "\n";
    dateError = true;
  }
  return dateError;
}

void NamelistReader::processPath(std::string &path)
{
  replace(path.begin(), path.end(), '\\', '/');
}

} // namespace GRAM
