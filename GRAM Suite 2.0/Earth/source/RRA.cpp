//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include "RRA.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
RRA::RRA()
  : EarthCommon(this)
{
  atmos.setPlanetSpecificMetrics(earthAtmos);
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
RRA::RRA(const RRA& orig)
  : EarthCommon(this)
{
  // Copy input parameters.
  useRRA = orig.useRRA;
  rraYear = orig.rraYear;
  outerRadius = orig.outerRadius;
  innerRadius = orig.innerRadius;
  ewWindPerturbationScale = orig.ewWindPerturbationScale;
  nsWindPerturbationScale = orig.nsWindPerturbationScale;
  densityPerturbationScale = orig.densityPerturbationScale;
  pressurePerturbationScale = orig.pressurePerturbationScale;
  temperaturePerturbationScale = orig.temperaturePerturbationScale;
  month = orig.month;
  rraPath = orig.rraPath;
  rraSiteList = orig.rraSiteList;

  // If using RRA, copy site list data.
  // RRA site data is not copied.  It will be read if needed.
  if (useRRA) {
    siteNameList = orig.siteNameList;
    siteGeodeticLatitudeList = orig.siteGeodeticLatitudeList;
    siteGeocentricLatitudeList = orig.siteGeocentricLatitudeList;
    siteLongitudeList = orig.siteLongitudeList;
    siteSurfaceHeightList = orig.siteSurfaceHeightList;
    siteMaxHeightList = orig.siteMaxHeightList;
  }
}

//! \fn  RRA::~RRA()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void RRA::setInputParameters(const EarthInputParameters& params)
{
  useRRA = params.useRRA;
  outerRadius = params.rraOuterRadius;
  innerRadius = params.rraInnerRadius;
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
  densityPerturbationScale = params.randomPerturbationScale;
  pressurePerturbationScale = params.randomPerturbationScale;
  temperaturePerturbationScale = params.randomPerturbationScale;
  rraYear = params.rraYear;
  month = params.month;
  rraSiteList = params.rraSiteList;
  rraPath = params.rraPath;
  readSiteList();
}

//! \fn  RRA::setRRASiteList(const std::string& fileName)
//! \brief Set the name of the RRA site list file.
//!
//! By default, this file name is "rrasites.txt".  It is assumed to be located in
//! the RRAPath folder.
//!
//! \param fileName   The name of the RRA site list file.

//! \fn  RRA::setMonth(int mon)
//! \brief Set the month for reading RRA data files.
//!
//! This parameter must be set prior to reading the RRA data.  Once set, the
//! month parameter should not be altered.
//!
//! \param mon   The month (1-12) for RRA data.

//! \brief Set RRA lookup parameters.
//!
//! These parameters control which RRA dataset is used and how close a position
//! must be to location of an RRA site.
//!
//! \param year   The RRAYearType defines which RRA data set to use.
//! \param inner  The inner fairing radius \units{\text{degrees}}.
//! \param outer  The outer fairing radius \units{\text{degrees}}.
void RRA::setRRAParameters(RRAYearType year, greal inner, greal outer)
{
  rraYear = year;
  innerRadius = inner;
  outerRadius = outer;
}

//! \brief Activate the RRA capability.
//!
//! The RRA capability is off by default.  In order to use RRA data, 
//! this method must be called after setting the RRA parameters.
//!
//! \param useFlag True to enable RRA, false to disable.
void RRA::setUseRRA(bool useFlag)
{
  useRRA = useFlag;
  readSiteList();
}

//! \brief Set the input atmosphere state.
//!
//! The RRA routines will fair an existing atmosphere state with that of a nearby
//! RRA site when appropriate.  This method supplies the existing atmosphere
//! state reference to the RRA routines.  The RRA will modify this atmosphere 
//! state with the faired values upon update.
//!
//! \param state   The input atmosphere state.
void RRA::setInputState(AtmosphereState& state)
{
  inState = &state;
}

//! \fn  RRA::setEpsilon(greal eps)
//! \brief Set the epsilon used for computing vapor density.
//!
//! \param eps   The epsilon used for computing vapor density.

//! \fn  RRA::setPpmToND(greal ppm2nd)
//! \brief Set the ppm to number density conversion factor.
//!
//! \param ppm2nd   The ppm to number density conversion factor.

//! \fn  RRA::isWeighted()
//! \brief Logical for when an RRA site weight exists.
//!
//! \returns True if RRA is active, a site is located, and the current weight is non-zero.

//! \fn  RRA::isActive()
//! \brief Logical for when the RRA object is activated.
//!
//! \returns True if RRA is active.

//! \fn  RRA::getRRAYear()
//! \brief Get the selected RRAYearType.
//!
//! \returns The selected RRAYearType.

//! \brief Returns the fairing heights for the current RRA site.
//!
//! The top of an RRA column is smoothed with existing atmospheres using
//! over the region returned by this method.
//! This method assumes an RRA site has been found, typically by calling
//! the update() method. If not, lower and upper are set to -DBL_MAX.
//!
//! \param[out] lower A greal.
//! \param[out] upper A greal.
//!
//! \retval lower The lower height of the fairing region \units{km}.
//! \retval upper The upper height of the fairing region \units{km}.
void RRA::getFairingHeights(greal &lower, greal &upper) const
{
  if (siteIndex >= 0) {
    lower = siteMaxHeightList[siteIndex] - 5;
    upper = siteMaxHeightList[siteIndex];
  }
  else {
    lower = upper = -DBL_MAX;
  }
}

//! \brief Read the RRA site list.
//!
//! The RRA site list contains the name, location, and height details about all of
//! the available RRA sites.  This method reads and stores that information.
//!
//! \b Inputs:
//! \arg #rraPath
//! \arg #rraSiteList
//! \arg #rraYear
//!
//! \returns The following vectors are populated: #siteNameList, #siteLongitudeList,
//!          #siteGeodeticLatitudeList, #siteGeocentricLatitudeList, #siteSurfaceHeightList, #siteMaxHeightList.
void RRA::readSiteList()
{
  if (!isActive() || !siteNameList.empty()) {
    return;
  }

  ifstream rrasites;
  rrasites.open(rraPath  + "/" + rraSiteList);
	
	// File open error check
  if (!rrasites) {
    throw string(FILE_OPEN_ERROR_MESSAGE + rraPath + "/" + rraSiteList);
  }

  // Convert the year type to a string
  std::string rraYearString;
  switch (rraYear) {
  case RRA_1983:
		rraYearString = "83";
    break;
  case RRA_2006:
    rraYearString = "06";
    break;
  case RRA_2013:
    rraYearString = "13";
    break;
  case RRA_2019:
		rraYearString = "19";
    break;
	}

  // Clear out all lists.
  siteNameList.clear();
  siteLongitudeList.clear();
  siteSurfaceHeightList.clear();
  siteMaxHeightList.clear();
  siteGeocentricLatitudeList.clear();

  // Ignore the header line
  string notUsed;
  getline(rrasites, notUsed);

	// Read file to end of file
  while (!rrasites.eof()) {

    greal longitudeEast, surfaceHeightMeters, maxHeight;
    string rraCode, rraSiteYear, wmo;
    double geodeticLatitude, geocentricLatitude;

    // Get the first string.  It should be the RRA code.
    rrasites >> rraCode;

    // If the line is blank.  Then we must be at the end of the list.  Stop the loop.
    if (rraCode.empty()) {
      break;
    }

    // Read more data from the line.
    rrasites >> rraSiteYear
             >> geodeticLatitude 
             >> geocentricLatitude // not used
             >> longitudeEast 
             >> surfaceHeightMeters 
             >> maxHeight;

    // The rest of the line is not used
    getline(rrasites, notUsed);

    // If the year matches, then store the data
    if (rraSiteYear.substr(2, 3) == rraYearString) {
      siteNameList.push_back(rraCode + rraSiteYear.substr(2, 3));
      siteLongitudeList.push_back(longitudeEast);
      siteSurfaceHeightList.push_back(surfaceHeightMeters);
      siteMaxHeightList.push_back(maxHeight);
      siteGeodeticLatitudeList.push_back(geodeticLatitude);

      //Compute and store geocentric latitude
      Position sitePosition;
      sitePosition.isPlanetoCentric = false;
      sitePosition.latitude = geodeticLatitude;
      sitePosition.convertToPlanetocentric(polarRadius, equatorialRadius);
      siteGeocentricLatitudeList.push_back(sitePosition.latitude);
    }
  }

	//Close RRA site data file to re-use unit for RRA data files
	rrasites.close();
}

//! \brief Parses the first line of an RRA file.
//!
//! Parses the first line of an RRA file and verified the sites longitude and latitude.
//!
//! \param rraDataFile The RRA file stream (already open).
//! \param fileName    The RRA file name (for error messages).
void RRA::parseTopLine(std::ifstream& rraDataFile, const std::string& fileName) {
  string topLine;
  getline(rraDataFile, topLine);
  // Parse out the latitude
  size_t start = 15;
  size_t end = topLine.find_first_of(' ', start);
  start = topLine.find_first_not_of(' ', end);
  end = topLine.find_first_of("NS", start);
  greal zlat = stod(topLine.substr(start, end - start));
  string ns = topLine.substr(end, 1);

  // Parse out the longitude
  start = topLine.find_first_not_of(' ', end + 1);
  end = topLine.find_first_of("EW", start);
  greal zlon = stod(topLine.substr(start, end - start));
  string ew = topLine.substr(end, 1);

  // Make north latitudes positive and south latitudes negative or
  // write error message and stop otherwise
  if (ns == "S") {
    zlat *= -1.0;
  }
  else if (ns != "N") {
    throw(FILE_PARSE_ERROR_MESSAGE + fileName
      + "\n       Bad latitude code on line: " + topLine);
  }

  // Make east longitudes positive and west longitudes negative or
  // write error message and stop otherwise
  if (ew == "W") {
    zlon *= -1.0;
  }
  else if (ew != "E") {
    throw(FILE_PARSE_ERROR_MESSAGE + fileName
      + "\n       Bad longitude code on line: " + topLine);
  }

  // Write error message and stop if RRA latitude from file does not agree with
  // expected latitude.
  if (abs(siteLatitude - zlat) > 0.005) {
    throw(FILE_PARSE_ERROR_MESSAGE + fileName
      + "\n       Bad latitude on line: " + topLine);
  }

  // Write error message and stop if RRA longitude from file does not agree with
  // expected longitude
  if (abs(siteLongitude - zlon) > 0.005) {
    throw(FILE_PARSE_ERROR_MESSAGE + fileName
      + "\n       Bad longitude on line: " + topLine);
  }
}

//! \brief Seeks the desired month's data within the RRA file.
//!
//! After the first line of the file has been parsed, this function is used
//! to position the input stream at the start of the desired month's data.
//!
//! \b Inputs:
//! \arg #month

//! \param rraDataFile The RRA file stream (already open, first line parsed).
void RRA::parseToMonth(std::ifstream& rraDataFile)
{
  // Initialize month read from file and previous height read 
  // (first data in file are annual mean values, month = 0)
  int m1;
  if (siteName.substr(3, 4) == "83") {
    m1 = 0;
  }
  else {
    m1 = 1;
  }

  string lineInput;
  getline(rraDataFile, lineInput);

  // Skip all months before selected month
  while (m1 < month) {
    // Read a line
    getline(rraDataFile, lineInput);

    // Find the index of the first non-blank and non-one character
    size_t firstIndex = lineInput.find_first_not_of(" 1");

    // If it is a T, then were are at the next T A B L E header line
    if (lineInput[firstIndex] == 'T') {
      // Increment the month counter.
      ++m1;
    }
  }

  // Get and ignore the next four header lines
  getline(rraDataFile, lineInput);
  getline(rraDataFile, lineInput);
  getline(rraDataFile, lineInput);
  getline(rraDataFile, lineInput);
}

//! \brief Read the RRA winds (T1) file.
//!
//! Reads Range Reference Atmosphere (RRA) wind components (u and v) and
//! standard deviations (su and sv) for site at geocentric latitude = xlat
//! longitude = xlon from RRA data file T1yyrra.txt.  Month is m.
//!
//! \b Inputs:
//! \arg #rraPath
//! \arg #siteName
//! \arg #month
//!
//! \returns The following vectors are populated: z1, u, su, v, sv, ruvt,
//!           avspd, sdspd.
void RRA::readrra1()
{
  if (CONSOLE_OUTPUT) {
    cout << "Reading RRA data: " << siteName << '\n';
  }

  // Clear existing data (if any)
  z1.clear();
  u.clear();
  su.clear();
  v.clear();
  sv.clear();
  ruvt.clear();
  avspd.clear();
  sdspd.clear();

  // Open file for reading RRA data
  string fileName = rraPath + "T1" + siteName + ".txt";
  ifstream rraDataFile;
  rraDataFile.open(fileName);

  // File open error check
  if (!rraDataFile) {
    throw string(FILE_OPEN_ERROR_MESSAGE + fileName);
  }

  try {
    // Read the first line for lat-lon
    parseTopLine(rraDataFile, fileName);

    // Read (and ignore) up to the desired month
    parseToMonth(rraDataFile);

    // Highest RRA data allowed
    greal ztop = siteMaxHeightList[siteIndex];

    surfaceIndex1 = 0;

    // Loop through the data lines until we hit the maximum height  
    greal zcurr = 0.0;
    while (zcurr < ztop) {
      // Read a line of data
      greal uwn, usd, ruv, vwn, vsd, ws, sws, skew, nobs;
      rraDataFile >> zcurr >> uwn >> usd >> ruv >> vwn >> vsd >> ws >> sws >> skew >> nobs;

      // The rational behind removing height zero: All surface heights are greater than zero
      // at all RRA sites.  Not all T1 files have data for height zero.  When a T1 file does
      // have a height zero line, all values are zero on that line.  It is possible that a position
      // will fall within RRA limits and have an above ground height that is less than the RRA
      // site's surface height.  Since no wind data exists below the RRA surface level, the surface
      // level wind metrics will be used.  The cases where the position is actually below the surface
      // will be handled elsewhere in the program.

      // Ignore lines where nobs is too small, or data is invalid
      if (nobs > 10.0 && (abs(uwn) > 0.0 || abs(vwn) > 0.0)) {

        // Set minimum allowed standard deviations
        if (usd < 0.01) {
          usd = 0.01;
        }
        if (vsd < 0.01) {
          vsd = 0.01;
        }

        // Load RRA array values
        z1.push_back(zcurr);
        u.push_back(uwn);
        su.push_back(usd);
        v.push_back(vwn);
        sv.push_back(vsd);
        ruvt.push_back(ruv);
        avspd.push_back(ws);
        sdspd.push_back(sws);
      }

      // If we are within one meter of the site list surface height, save the index.
      if (abs(1000.0 * zcurr - siteSurfaceHeightList[siteIndex]) < 1.0) {
        surfaceIndex1 = int(z1.size() - 1);
      }
    }
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_TEXT_ERROR_MESSAGE + fileName
      + "\n       " + e.what());
  }

  // Close the RRA input file (to re-use the RRA unit number)
  rraDataFile.close();
}

//! \brief Read the RRA pressre, density, and temperature (T2) file.
//!
//! Reads Range Reference Atmosphere (RRA) pressure (p), temperature (t), density (d)
//! and standard deviations (sp, st, sd) for site at latitude = xlat, longitude = xlon
//! from RRA data file T2yyrra.txt.  Month is m.
//!
//! \b Inputs:
//! \arg #rraPath
//! \arg #siteName
//! \arg #month
//!
//! \returns The following vectors are populated: z2, p, sp, t, st, d, sd.
void RRA::readrra2()
{
  // Clear existing data (if any)
  z2.clear();
  p.clear();
  sp.clear();
  t.clear();
  st.clear();
  d.clear();
  sd.clear();

  // Open file for reading RRA data
  string fileName = rraPath + "T2" + siteName + ".txt";
  ifstream rraDataFile;
  rraDataFile.open(fileName);

  // File open error check
  if (!rraDataFile) {
    throw(FILE_OPEN_ERROR_MESSAGE + fileName);
  }

  try {
    // Read the first line for lat-lon
    parseTopLine(rraDataFile, fileName);

    // Read (and ignore) up to the desired month
    parseToMonth(rraDataFile);

    surfaceIndex2 = 0;

    // Highest RRA data allowed is 70 km; max number = 300
    greal ztop = siteMaxHeightList[siteIndex];

    // Loop through the data lines until we hit the maximum height
    greal zcurr = 0.0;
    while (zcurr < ztop) {
      // Read a line of data
      greal pr, spr, skp, tmp, stmp, skt, dns, sdns, skd, xnp, xnt, xnd;
      rraDataFile >> zcurr >> pr >> spr >> skp >> tmp >> stmp >> skt
        >> dns >> sdns >> skd >> xnp >> xnt >> xnd;

      // The rational behind removing height zero: All surface heights are greater than zero
      // at all RRA sites.  The sub-surface T2 data is trusted for 1983 and 2006 data files. 
      // The only sub-surface heights appearing in most files is a height zero line.  A few of
      // the 1983 files have multiple sub-surface heights.  So for 2013 and 2019 data, the
      // zero heights will be removed.
      
      // It is possible that a valid position will fall within RRA limits and have an above ground 
      // height that is less than the RRA site's surface height.  The cases where the position 
      // is actually below the surface will be handled elsewhere in the program.

      // Accept sub-surface heights for 1983 and 2006 data only.
      // Ignore any data lines with too few observed data points or invalid data points
      if (
       // (rraYear == RRA_1983 || rraYear == RRA_2006 || zcurr > 0.0) &&
        xnp > 10.0 && xnt > 10.0 && xnd > 10.0 && tmp <= 999.0 && stmp <= 99.0 && pr > 0.0
        && dns > 0.0 && tmp > 0.0) {

        // Set minimum allowed standard deviations
        if (spr < 0.001) {
          spr = 0.001;
        }
        if (stmp < 0.01) {
          stmp = 0.01;
        }
        if (sdns < 0.001) {
          sdns = 0.001;
        }

        // Ensure pressure standard deviations are not too large
        if (spr > 0.3 * pr) {
          spr = 0.3 * pr;
          if (CONSOLE_OUTPUT) {
            cout << " Large RRA pressure sigma (reset) " << zcurr << "  "
              << pr << "  " << spr << endl;
          }
        }

        // Ensure temperature standard deviations are not too large
        if (stmp > 0.3 * tmp) {
          stmp = 0.3 * tmp;
          if (CONSOLE_OUTPUT) {
            cout << " Large RRA temp. sigma (reset) " << zcurr << "  "
              << tmp << "  " << stmp << endl;
          }
        }

        // Ensure density standard deviations are not too large
        if (sdns > 0.3 * dns) {
          sdns = 0.3 * dns;
          if (CONSOLE_OUTPUT) {
            cout << " Large RRA density sigma (reset) " << zcurr << "  "
              << dns << "  " << sdns << endl;
          }
        }

        // Load RRA array values
        z2.push_back(zcurr);
        p.push_back(pr * 100.0);
        sp.push_back(spr / pr);
        t.push_back(tmp);
        st.push_back(stmp / tmp);
        d.push_back(dns / 1000.0);
        sd.push_back(sdns / dns);

        // If we are within one meter of the site list surface height, save the index.
        if (abs(1000.0 * zcurr - siteSurfaceHeightList[siteIndex]) < 1.0) {
          surfaceIndex2 = int(z2.size() - 1);
        }
      }
    }
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_TEXT_ERROR_MESSAGE + fileName
      + "\n       " + e.what());
  }

  // Close the RRA input file (to re-use the RRA unit number)
  rraDataFile.close();
}


//! \brief Read the RRA vapor pressre, and dewpoint (T3) file.
//!
//! Read Range Reference Atmosphere (RRA) vapor pressure (vp), dewpoint
//! temperature (td) and standard deviations (svp, std) for a site at latitude
//! xlat, longitude = xlon from RRA data file T3yyrra.txt.  Month is m.
//!
//! \b Inputs:
//! \arg #rraPath
//! \arg #siteName
//! \arg #month
//!
//! \returns The following vectors are populated: z3, vp, svp, td, std.
void RRA::readrra3()
{
  // Clear existing data (if any)
  z3.clear();
  vp.clear();
  svp.clear();
  td.clear();
  std.clear();

  // Open file for reading RRA data
  string fileName = rraPath + "T3" + siteName + ".txt";
  ifstream rraDataFile;
  rraDataFile.open(fileName);

  // File open error check
  if (!rraDataFile.is_open()) {
    throw(FILE_OPEN_ERROR_MESSAGE + fileName);
  }

  try {
    // Read the first line for lat-lon
    parseTopLine(rraDataFile, fileName);

    // Read (and ignore) up to the desired month
    parseToMonth(rraDataFile);

    // Read and ignore one more header line;
    string lineInput;
    getline(rraDataFile, lineInput);

    surfaceIndex3 = 0;

    // Highest RRA data
    greal ztop = 30.0;

    // Loop through the data lines until we hit the maximum height
    greal zcurr = 0.0;
    while (zcurr < ztop) {
      // Read a line of data
      greal vpr, svpr, skvp, tvir, stvir, sktv, tdp, stdp, sktd, dum, xntp, xntt;
      rraDataFile >> zcurr >> vpr >> svpr >> skvp >> tvir >> stvir >> sktv
        >> tdp >> stdp >> sktd >> xntp >> xntt;
      // These files have an extra field (not used)
      if ((siteName == "fad83") || (siteName == "nel83") || (siteName == "shm83")
        || (siteName == "thu83") || (siteName == "wak83")) {
        rraDataFile >> dum;
      }

      // The rational behind removing height zero: All surface heights are greater than zero
      // at all RRA sites.  The sub-surface T3 data is trusted for 1983 and 2006 data files. 
      // The only sub-surface heights appearing in most files is a height zero line.  A few of
      // the 1983 files have multiple sub-surface heights.  So for 2013 and 2019 data, the
      // zero heights will be removed.

      // It is possible that a valid position will fall within RRA limits and have an above ground 
      // height that is less than the RRA site's surface height.  The cases where the position 
      // is actually below the surface will be handled elsewhere in the program.

      // Accept sub-surface heights for 1983 and 2006 data only.
      // Ignore any data lines with invalid data points
      if (
        //(rraYear == RRA_1983 || rraYear == RRA_2006 || zcurr > 0.0) &&
        tdp <= 999.9 && stdp <= 99.9 && vpr <= 99.9 && svpr <= 99.9
        && tdp > 0.0 && vpr > 0.0 && xntp > 10.0 && xntt > 10.0) {

        // Set minimum allowed standard deviations
        if (svpr < 0.001) {
          svpr = 0.001;
        }
        if (stdp < 0.01) {
          stdp = 0.01;
        }

        // Ensure vapor pressure standard deviations are not too large
        if (svpr > 1.08 * vpr) {
          svpr = 1.08 * vpr;
          if (CONSOLE_OUTPUT) {
            cout << "Large RRA vapor pressure sigma (reset)"
              << "  " << zcurr << "  " << vpr << "  " << svpr << endl;
          }
        }

        // Ensure dewpoint temperature standard deviations are not too large
        if (stdp > 0.3 * tdp) {
          stdp = 0.3 * tdp;
          if (CONSOLE_OUTPUT) {
            cout << "Large RRA dewpoint sigma (reset)"
              << "  " << zcurr << "  " << tdp << "  " << stdp << endl;
          }
        }

        // Load RRA array values
        z3.push_back(zcurr);
        vp.push_back(vpr * 100.0);
        svp.push_back(svpr * 100.0);
        td.push_back(tdp);
        std.push_back(stdp);

        // If we are within one meter of the site list surface height, save the index.
        if (abs(1000.0 * zcurr - siteSurfaceHeightList[siteIndex]) < 1.0) {
          surfaceIndex3 = int(z3.size() - 1);
        }
      }
    }
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_TEXT_ERROR_MESSAGE + fileName
      + "\n       " + e.what());
  }

  // Close the RRA input file (to re-use the RRA unit number)
  rraDataFile.close();
}

//! \brief Get a height based weight for RRA to GRAM data transitions.
//!
//! Gets the weighting factor for smooth altitude transition from Range Reference
//! Atmosphere (RRA) data to GRAM data.
//! The computed weight factor is 1 up to the next to last (top) entry in the 
//! supplied height list. The weight is 0 if the height is above the top
//! RRA height. Linear variation is used between the top and next to top height.
//!
//! \param h The height of the GRAM data \units{km}.
//! \param z A list of increasing heights \units{km}.
//!
//! \returns A weighting factor.
greal RRA::getHeightWeight(greal h, const std::vector<greal> &z)
{
  // Assumptions (checking occurs in debug mode only)
  // There must be at least 2 data points.
  assert(z.size() > 1);

  greal heightwt_out;
  size_t last = z.size() - 1;

  // Height weight factor = 1 up to RRA height next to top, 0 if above top
  // RRA height, linear variation between top and next to top height
  if (h <= z[last - 1]) {
    heightwt_out = 1.0;
  }
  else if (h >= z[last]) {
    heightwt_out = 0.0;
  }
  else {
    heightwt_out = (z[last] - h) / (z[last] - z[last - 1]);
  }

  return heightwt_out;
}

//! \brief Get index from height array.
//!
//! Find array index number for Range Reference Atmosphere (RRA) height array (z) for input altitude (h).
//! The index returned is the lower bound index of the segment bracketing the height.
//!
//! \param h The height of the GRAM data \units{km}.
//! \param z A list of increasing heights \units{km}.
//!
//! \returns The lower bound index of the segment bracketing the height.
int RRA::findHeightIndex(greal h, const std::vector<greal> &z)
{
  size_t heightIndex = 0;

  // Finds array index for which z[i] < h < z[i+1]
  while (heightIndex < z.size() - 1 && h > z[heightIndex]) {
    if (h <= z[heightIndex + 1]) {
      break;
    }
    ++heightIndex;
  }

  return int(heightIndex);
}

//! \brief Horizontal weight factor for Range Reference Atmosphere (RRA) data.
//!
//! Inputs are delta, the latitude-longitude distance from location to RRA
//! site, sitelim, the maximum lat-lon distance allowed for adjustment to
//! take place, and sitenear.  The lat-lon limit within which RRA data is
//! used with full weight of 1.
//!
//! \b Inputs:
//! \arg #innerRadius
//! \arg #outerRadius
//!
//! \param delta The latitude-longitude distance from location to RRA site \units{\text{degrees}}.
//!
//! \returns A weighting factor.
greal RRA::getWeightFactor(greal delta)
{
  greal weightfact_out;

  if (delta < innerRadius) {
    // RRA weight factor = 1 if location is near RRA site
    weightfact_out = 1.0;
  }
  else if (delta > outerRadius) {
    // RRA weight factor = 0 if location outside sitelim
    weightfact_out = 0.0;
  }
  else {
    // RRA weight factor decrease to 0 (as cosine squared) between 'sitenear'
    // limit and distance = sitelim
    weightfact_out = square(cos(PI * (delta - innerRadius) / (2.0 * (outerRadius - innerRadius))));
  }

  return weightfact_out;
}

//! \brief Locate the closest RRA site.
//!
//! This method looks for the closest RRA site to the current position.  The site
//! must be within the outerRadius lat-lon limit.  If a site is found, the 
//! horizontal weighting factor is computed.
//!
//! \b Inputs:
//! \arg #latitude
//! \arg #longitude
//! \arg #siteGeocentricLatitudeList
//! \arg #siteLongitudeList
//! \arg #siteNameList
//!
//! \retval siteIndex   Index of the site in the RRA site list.
//! \retval siteWeight  The horizontal weighting factor.
//! \retval siteName    Name of the RRA site (or "GRM" if no site).
void RRA::locateNearestSite()
{
  greal minAngle = 999.0;
  siteName = "GRM";
  siteIndex = NO_SITE;
  siteWeight = 0.0;
  siteLatitude = 0.0;
  siteLongitude = 0.0;

  //Step through the RRA sites
  for (size_t i = 0; i < siteNameList.size(); i++) {

    //Compute arc angle of given lat-lon from RRA site
    greal angleToSiteCenter = getArcAngle(latitude, longitude, siteGeocentricLatitudeList[i], siteLongitudeList[i]);

    //Find minimum radius that is within lat-lon limit (if any)
    if ((angleToSiteCenter <= outerRadius) && (angleToSiteCenter < minAngle)) {
      minAngle = angleToSiteCenter;
      siteIndex = int(i);
    }
  }

  // If a site was found, set the weight, name and location.
  if (siteIndex != NO_SITE) {
    //Horizontal weighting factor for RRA site
    siteWeight = getWeightFactor(minAngle);

    //RRA site name year & code
    siteName = siteNameList[siteIndex];

    // These are used to verify RRA data files.
    siteLatitude = siteGeodeticLatitudeList[siteIndex];
    siteLongitude = siteLongitudeList[siteIndex];
  }
}

//! \brief Interface for the primary atmosphere computations.
//!
//! Based on a given height, geocentric latitude, and longitude
//! determine if GRAM values need to be adjusted by Range Reference
//! Atmosphere (RRA) values (given location must be within latitude-longitude
//! "radius" for adjustment to be used).  If more than one RRA
//! site is within adjustment range, the RRA site with minimum distance from
//! the given location is used.  Computes modified atmospheric values, based
//! on weighted average of original GRAM value and RRA value.
void RRA::update()
{
  // Only update if RRA is activated.
  if (!isActive()) {
    return;
  }

  // Find the RRA site nearest to the current position.
  locateNearestSite();

  // Set references to the incoming atmosphere state data.
  AtmosphereState& inAtmos = *inState;
  EarthAtmosphereState& inEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();
  // Yes, these are the same object, but it reads better this way.
  AtmosphereState& outAtmos = *inState;
  EarthAtmosphereState& outEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();

  // Set initial values of windCorrelation and rra weights.
  windCorrelation = inEarth.windCorrelation;
  rraWeight1 = 0.0;
  rraWeight2 = 0.0;
  rraWeight3 = 0.0;

  // Bypass RRA modification if all RRA sites are outside the lat-lon limit range
  if (siteIndex == NO_SITE) {
    siteName = "GRM";
    inEarth.rraWeight = 0.0;
    inEarth.rraSiteName[siteName.copy(inEarth.rraSiteName, 6)] = 0;
    return;
  }

  // Read RRA data unless this is the RRA site that was last read in
  if (siteIndex != siteReadIndex) {

    // Read Table 1 RRA data (wind components and standard deviations)
    readrra1();

    // Read Table 2 RRA data (pressure, density, temperature, and their standard deviations)
    readrra2();

    // Read Table 3 RRA data (vapor pressure and dewpoint temperature and their standard deviations)
    readrra3();

    // Store RRA site number of site just read in
    siteReadIndex = siteIndex;
  }

  // Modify incoming GRAM wind values if in height range
  if (height <= z1.back()) {
    // Find height index values for vertical interpolation
    int n1 = findHeightIndex(height, z1);
    int n2 = n1 + 1;

    // Create a height-based interpolater which does not extrapolate below
    // the first data point.  In most cases, findHeightIndex() will return
    // an index so that z1[n1] <= height <= z1[n2].  However, if height < z1[0],
    // then we will use z1[0] as the interpolation height.
    Interpolator interpHeight;
    interpHeight.makeFraction(z1[n1], z1[n2], max(height, z1[n1]));

    // Do vertical interpolation to get wind components
    ewWind = interpHeight.linear(u[n1], u[n2]);
    nsWind = interpHeight.linear(v[n1], v[n2]);
    windCorrelation = interpHeight.linear(ruvt[n1], ruvt[n2]);
    if (height <= z1[surfaceIndex1]) {
      windCorrelation = ruvt[surfaceIndex1];
    }

    // Do vertical interpolation to get wind standard deviations
    ewStandardDeviation = interpHeight.linear(su[n1], su[n2]);
    nsStandardDeviation = interpHeight.linear(sv[n1], sv[n2]);

    // Do vertical interpolation to get speed statistics
    windSpeed = interpHeight.linear(avspd[n1], avspd[n2]);
    windSpeedStandardDeviation = interpHeight.linear(sdspd[n1], sdspd[n2]);

    if (height <= z1[surfaceIndex1]) {
      ewWind = u[surfaceIndex1];
      nsWind = v[surfaceIndex1];
      windSpeed = avspd[surfaceIndex1];
      windSpeedStandardDeviation = sdspd[surfaceIndex1];
      ewStandardDeviation = su[surfaceIndex1] * ewWindPerturbationScale;
      nsStandardDeviation = sv[surfaceIndex1] * nsWindPerturbationScale;
    }
    else {
      ewStandardDeviation = ewStandardDeviation * ewWindPerturbationScale;
      nsStandardDeviation = nsStandardDeviation * nsWindPerturbationScale;
    }

    // Total RRA weighting factor = products of height weight and horizontal
    // weight factors
    rraWeight1 = siteWeight * getHeightWeight(height, z1);

    // Create an interolator for incoming and RRA data.
    Interpolator interp1(rraWeight1);

    // Revise GRAM winds and standard deviations by weighted average with RRA data
    outAtmos.ewWind = interp1.linear(inAtmos.ewWind, ewWind);
    outAtmos.nsWind = interp1.linear(inAtmos.nsWind, nsWind);
    outAtmos.ewStandardDeviation = interp1.linear(inAtmos.ewStandardDeviation, ewStandardDeviation);
    outAtmos.nsStandardDeviation = interp1.linear(inAtmos.nsStandardDeviation, nsStandardDeviation);
    outEarth.windCorrelation = interp1.linear(inEarth.windCorrelation, windCorrelation);
  }
  else {
    rraWeight1 = 0.0;
  }

  // Below 27 km, blend in surface values.
  if (height <= SURFACE_LIMIT) {
    // For winds, note that index 0 is the index of the surface data.
    ewWindAtSurface = u[surfaceIndex1];
    nsWindAtSurface = v[surfaceIndex1];
    ewWindSDAtSurface = su[surfaceIndex1];
    nsWindSDAtSurface = sv[surfaceIndex1];
    windCorrelationAtSurface = ruvt[surfaceIndex1];
    windSpeedAtSurface = avspd[surfaceIndex1];
    windSpeedSDAtSurface = sdspd[surfaceIndex1];

    // Create an interolator for incoming and RRA data based on distance to the RRA location.
    Interpolator interpS(siteWeight);
    outEarth.ewWindAtSurface = interpS.linear(inEarth.ewWindAtSurface, ewWindAtSurface);
    outEarth.nsWindAtSurface = interpS.linear(inEarth.nsWindAtSurface, nsWindAtSurface);
    outEarth.ewWindSDAtSurface = interpS.linear(inEarth.ewWindSDAtSurface, ewWindSDAtSurface);
    outEarth.nsWindSDAtSurface = interpS.linear(inEarth.nsWindSDAtSurface, nsWindSDAtSurface);
    outEarth.windCorrelationAtSurface = interpS.linear(inEarth.windCorrelationAtSurface, windCorrelationAtSurface);
    outEarth.windSpeedAtSurface = interpS.linear(inEarth.windSpeedAtSurface, windSpeedAtSurface);
    outEarth.windSpeedSDAtSurface = interpS.linear(inEarth.windSpeedSDAtSurface, windSpeedSDAtSurface);
    // Surface height is in the Position class. 
    surfaceHeight = interpS.linear(surfaceHeight, z1[surfaceIndex1]);
  }

  // Modify incoming GRAM pressure, density and temperature values with RRA values if in height range
  if (height <= z2.back()) {

    // Find height index values for vertical interpolation
    int n1 = findHeightIndex(height, z2);
    int n2 = n1 + 1;

    if (height <= z2[surfaceIndex2]) {
      pressure = p[surfaceIndex2];
      density = d[surfaceIndex2];
      temperature = t[surfaceIndex2];
      pressureStandardDeviation = sp[surfaceIndex2] * pressurePerturbationScale;
      densityStandardDeviation = sd[surfaceIndex2] * densityPerturbationScale;
      temperatureStandardDeviation = st[surfaceIndex2] * temperaturePerturbationScale;
    }
    else {
      // Create a height-based interpolater.
      // Since height > z2[0], no extrapolation will occur below the first data point.
      Interpolator interpHeight;
      interpHeight.makeFraction(z2[n1], z2[n2], height); 

      // Temperature: linear interpolation
      temperature = interpHeight.linear(t[n1], t[n2]);

      // Gas constant: linear interpolation
      greal r1 = p[n1] / (d[n1] * t[n1]);
      greal r2 = p[n2] / (d[n2] * t[n2]);
      greal r = interpHeight.linear(r1, r2);

      // Pressure: logarithmic interpolation if temperature gradient = 0
      if (abs(t[n2] - t[n1]) <= 0.001) {
        pressure = interpHeight.log(p[n1], p[n2]);
      }
      // Power law interpolation on pressure
      else {
        greal expon = log(p[n2] / p[n1]) / log(t[n1] / t[n2]);
        greal base = (t[n1] / temperature);
        pressure = p[n1] * pow(base, expon);
      }

      // Density from perfect gas law
      density = pressure / (r * temperature);

      // Standard deviations: linear interpolation
      pressureStandardDeviation = interpHeight.linear(sp[n1], sp[n2]) * pressurePerturbationScale;
      densityStandardDeviation = interpHeight.linear(sd[n1], sd[n2]) * densityPerturbationScale;
      temperatureStandardDeviation = interpHeight.linear(st[n1], st[n2]) * temperaturePerturbationScale;
    }

    // Compute height weighting factor for RRA values
    greal heightWeight = getHeightWeight(height, z2);

    // Total RRA weighting factor = product of height weight and horizontal weight factors
    rraWeight2 = siteWeight * heightWeight;

    // Create an interolator for incoming and RRA data.
    Interpolator interp2(rraWeight2);

    // Revise GRAM pressure, density, temperature and standard deviations
    // by weighted average with RRA data
    presBeforeRRA = inAtmos.pressure;
    tempBeforeRRA = inAtmos.temperature;
    outAtmos.pressure = interp2.linear(inAtmos.pressure, pressure);
    outAtmos.density = interp2.linear(inAtmos.density, density);
    outAtmos.temperature = interp2.linear(inAtmos.temperature, temperature);
    outAtmos.pressureStandardDeviation = interp2.linear(inAtmos.pressureStandardDeviation, pressureStandardDeviation);
    outAtmos.densityStandardDeviation = interp2.linear(inAtmos.densityStandardDeviation, densityStandardDeviation);
    outAtmos.temperatureStandardDeviation = interp2.linear(inAtmos.temperatureStandardDeviation, temperatureStandardDeviation);

    // Compute height weighting factor for surface RRA values
    if (height <= SURFACE_LIMIT) {
      // Create an interolator for incoming and RRA data based on distance to the RRA location.
      Interpolator interpS(siteWeight);

      // Set the RRA surface values.
      pressureSDAtSurface = sp[surfaceIndex2];
      densitySDAtSurface = sd[surfaceIndex2];
      pressureAtSurface = p[surfaceIndex2];
      densityAtSurface = d[surfaceIndex2];
      temperatureAtSurface = t[surfaceIndex2];
      temperatureSDAtSurface = st[surfaceIndex2];

      // If incoming surface pressure/density data exists, then smooth with RRA data.
      if ((inAtmos.pressureAtSurface > 0.0) && (inEarth.densityAtSurface > 0.0)) {

        outAtmos.pressureAtSurface = interpS.log(inAtmos.pressureAtSurface, pressureAtSurface);
        outEarth.densityAtSurface = interpS.log(inEarth.densityAtSurface, densityAtSurface);
        outEarth.pressureSDAtSurface = interpS.linear(inEarth.pressureSDAtSurface, pressureSDAtSurface);
        outEarth.densitySDAtSurface = interpS.linear(inEarth.densitySDAtSurface, densitySDAtSurface);
      }
      else {
        // No incoming surface pressure/density data.
        outAtmos.pressureAtSurface = 0.0;
        outEarth.densityAtSurface = 0.0;
        outEarth.pressureSDAtSurface = 0.0;
        outEarth.densitySDAtSurface = 0.0;
      }

      // Always smooth incoming temperature with RRA temperatue at surface
      outEarth.temperatureAtSurface = interpS.linear(inEarth.temperatureAtSurface, temperatureAtSurface);
      outEarth.temperatureSDAtSurface = interpS.linear(inEarth.temperatureSDAtSurface, temperatureSDAtSurface * temperatureAtSurface);
    }
  }
  else {
    rraWeight2 = 0.0;
  }

  // Modify GRAM vapor pressure and dewpoint temperature values with RRA values
  // if in height range
  if (height <= z3.back()) {

    //Find height index values for vertical interpolation
    int n1 = findHeightIndex(height, z3);
    int n2 = n1 + 1;

    // Create a height-based interpolater which does not extrapolate below
    // the first data point.  In most cases, findHeightIndex() will return
    // an index so that z3[n1] <= height <= z3[n2].  However, if height < z3[0],
    // then we will use z3[0] as the interpolation height.
    Interpolator interpHeight;
    interpHeight.makeFraction(z3[n1], z3[n2], max(height, z3[n1]));

    // Linear interpolation
    vaporPressure = interpHeight.linear(vp[n1], vp[n2]);
    dewPoint = interpHeight.linear(td[n1], td[n2]);

    // For sub-surface heights, use the surface values for standard deviations.
    if (height <= z3[surfaceIndex3]) {
      vaporPressureSD = svp[surfaceIndex3];
      dewPointSD = std[surfaceIndex3];
    }
    else {
      //Do vertical interpolations to get standard deviations
      vaporPressureSD = interpHeight.linear(svp[n1], svp[n2]);
      dewPointSD = interpHeight.linear(std[n1], std[n2]);
    }

    //Compute height weighting factor for RRA values
    greal heightWeight = getHeightWeight(height, z3);
    rraWeight3 = siteWeight * heightWeight;
  }
  else {
    rraWeight3 = 0.0;
  } 

  // Not used and not output
  //if (height <= SURFACE_LIMIT) { //***
  //  dewPointSurface = td[surfaceIndex3]; //***
  //  dewPointSDSurface = std[surfaceIndex3]; //***
  //} //***

  // For reporting, save the PDT weight.
  rraWeight = rraWeight2;

  //Set rra_id to "GRM", weight to 0 if no RRA adjustment done
  if ((height >= z1.back()) || (height >= z2.back())) {
    siteName = "GRM";
    rraWeight = 0.0;
  }

  // Copy values to the incoming state.
  outEarth.rraWeight = rraWeight;
  outEarth.rraSiteName[siteName.copy(inEarth.rraSiteName, 6)] = 0;
}

//! \brief Update water related metrics.
//!
//! Fair the GRAM with the RRA values for vapor pressure, vapor density,
//! dewpoint, relative humidity, and their standard deviations.  The mole
//! fraction and number density for water is also updated.
//!
//! \b Inputs:
//! \arg #rraWeight3
//! \arg #height
//! \arg #inState
//! \arg #atmos (RRA metrics)
//! \arg #presBeforeRRA
//! \arg #ppmToND
//!
//! \returns #inState
void RRA::updatePsychrometrics()
{
  if (isActive() && rraWeight3 > 0.0) {
    // Set references to the incoming atmosphere state data.
    AtmosphereState& inAtmos = *inState;
    EarthAtmosphereState& inEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();
    // Yes, these are the same object, but it reads better this way.
    AtmosphereState& outAtmos = *inState;
    EarthAtmosphereState& outEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();

    // Create an interpolator for smoothing incoming psychrometric data with RRA data.
    Interpolator interp3(rraWeight3);

    // Don't modify values above the data range.
    if (height <= z3.back()) {
      outEarth.vaporPressure = interp3.linear(inEarth.vaporPressure, vaporPressure);
      outEarth.vaporPressureSD = interp3.linear(inEarth.vaporPressureSD, vaporPressureSD);
      outEarth.dewPoint = interp3.linear(inEarth.dewPoint, dewPoint);
      outEarth.dewPointSD = interp3.linear(inEarth.dewPointSD, dewPointSD);

      double dest = 0.5 * d2edt2(temperature, pressure) * temperatureStandardDeviation;
      relativeHumidity = 100.0 * vaporPressure / (wexler(temperature, pressure) + dest);
      relativeHumidity = clamp(relativeHumidity, 3.0, 100.0);
      outEarth.relativeHumidity = interp3.linear(inEarth.relativeHumidity, relativeHumidity);
      outEarth.vaporDensity = epsilon * inEarth.vaporPressure / (dryAirGasConstant * tempBeforeRRA);

      if (vaporPressure > 0.0) {
        double sigest = dedt(tempBeforeRRA, presBeforeRRA) * inAtmos.temperatureStandardDeviation
                        * tempBeforeRRA / inEarth.vaporPressure;
        double srhf = 1.43 + 0.315 * log10(pressure / 1.0e+05);
        double corr = 0.85 - 0.7 * pow(cos(toRadians(latitude)), 6.0);
        relativeHumiditySD = relativeHumidity
              * sqrt(square(inEarth.vaporPressureSD / inEarth.vaporPressure) + square(sigest)
                     - 2.0 * corr * sigest * inEarth.vaporPressureSD / inEarth.vaporPressure)
              / srhf;
        relativeHumiditySD = min(50.0, relativeHumiditySD);
        relativeHumiditySD = min(1.05 * relativeHumidity, relativeHumiditySD);
        relativeHumiditySD = max(0.05 * relativeHumidity, relativeHumiditySD);
        outEarth.relativeHumiditySD = interp3.linear(inEarth.relativeHumiditySD, relativeHumiditySD);
        outEarth.vaporDensitySD = inEarth.vaporDensity
                                 * sqrt(square(inEarth.vaporPressureSD / inEarth.vaporPressure)
                                        + square(inEarth.dewPointSD / inEarth.dewPoint));
      }

      outAtmos.water.moleFraction = inEarth.vaporPressure / (presBeforeRRA - inEarth.vaporPressure);
      outAtmos.water.numberDensity = inAtmos.water.moleFraction * 1.0e6 * ppmToND;
    }

    // Not used and not output
    //if (height <= SURFACE_LIMIT) {
    //  outEarth.dewPointSurface = grmwt * inEarth.dewPoint + rrawt * dewPointSurface;
    //  outEarth.dewPointSDSurface = grmwt * inEarth.dewPointSD + rrawt * dewPointSDSurface;
    //}
  }
}

//! \brief Update wind speed and its standard deviation.
//!
//! Fair the GRAM with the RRA values for wind speed and its standard deviation.
//!
//! \b Inputs:
//! \arg #rraWeight1
//! \arg #siteWeight
//! \arg #height
//! \arg #inState
//! \arg #atmos (RRA metrics)
//!
//! \returns #inState
void RRA::updateWindSpeed()
{
  // Set references to the incoming atmosphere state data.
  AtmosphereState& inAtmos = *inState;
  EarthAtmosphereState& inEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();
  // Yes, these are the same object, but it reads better this way.
  EarthAtmosphereState& outEarth = inAtmos.getPlanetSpecificMetrics<EarthAtmosphereState>();

  Interpolator interp1(rraWeight1);

  if (siteWeight > 0.0 && height <= z1.back()) {
    outEarth.windSpeed = interp1.linear(inEarth.windSpeed, windSpeed);
    outEarth.windSpeedStandardDeviation = interp1.linear(inEarth.windSpeedStandardDeviation, windSpeedStandardDeviation);
  }
}

} // namespace
