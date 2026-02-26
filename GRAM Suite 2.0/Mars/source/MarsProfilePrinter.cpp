//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "MarsProfilePrinter.h"
#include "MarsInputParameters.h"
#include "MarsAtmosphere.h"
#include "Ephemeris.h"

using namespace std;

namespace GRAM {

//! \copydoc ProfilePrinter::ProfilePrinter()
MarsProfilePrinter::MarsProfilePrinter()
{
}

//! \fn MarsProfilePrinter::MarsProfilePrinter(const MarsProfilePrinter& orig)
//! \copydoc ProfilePrinter::ProfilePrinter(const ProfilePrinter& orig)

//! \fn MarsProfilePrinter::~MarsProfilePrinter()
//! \copydoc ProfilePrinter::~ProfilePrinter()

//! \brief Header method for planet specific data in the GRAM_CSV_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printPlanetCSVHeader(const PerturbedAtmosphere& atmos)
{
  gramCSVFile << "TemperatureDaily_K,PressureDaily_Pa,DensityDaily_kgm3,EWWindDaily_ms,NSWindDaily_ms,"
    << "TemperatureMin_K,TemperatureMax_K,DensityMin_kgm3,DensityMax_kgm3,"
    << "PlanetoGraphicHeight_km,PlanetoGraphicLatitude_deg,ReferenceHeight_km,ReferenceRadius_km,"
    << "GroundTemperature_K,ThermosphereBaseHeight_km,ThermosphereBaseTemperature_K,ExosphericTemperature_K,"
    << "F1PeakHeight_km,Albedo,HeightOffset_km,LocalHeightOffset_km,"
    << "DustOpticalDepth,DustColumnArealDensity_kgm2,DustMixingRatio,DustMassDensity_ugm3,DustNumberDensity_m3,"
    << "IceIsPresent,WavePerturbation_pct";
}

//! \brief Print method for planet specific data the GRAM_CSV_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printPlanetCSVStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  gramCSVFile << fixed
    << setprecision(1 + extra) << atmMars.temperatureDaily << comma
    << scientific
    << setprecision(3 + extra) << atmMars.pressureDaily << comma
    << setprecision(3 + extra) << atmMars.densityDaily << comma
    << fixed
    << setprecision(1 + extra) << atmMars.ewWindDaily << comma
    << setprecision(1 + extra) << atmMars.nsWindDaily << comma
    << setprecision(1 + extra) << atmMars.temperatureMin << comma
    << setprecision(1 + extra) << atmMars.temperatureMax << comma
    << scientific
    << setprecision(3 + extra) << atmMars.densityMin << comma
    << setprecision(3 + extra) << atmMars.densityMax << comma
    << fixed
    << setprecision(3 + extra) << atmMars.planetoGraphicHeight << comma
    << setprecision(2 + extra) << atmMars.planetoGraphicLatitude << comma
    << setprecision(3 + extra) << atmMars.referenceHeight << comma
    << setprecision(3 + extra) << atmMars.referenceRadius << comma
    << setprecision(1 + extra) << atmMars.groundTemperature << comma
    << setprecision(3 + extra) << atmMars.thermosphereBaseHeight << comma
    << setprecision(1 + extra) << atmMars.thermosphereBaseTemperature << comma
    << setprecision(1 + extra) << atmMars.exosphericTemperature << comma
    << setprecision(3 + extra) << atmMars.f1PeakHeight << comma
    << setprecision(3 + extra) << atmMars.albedo << comma
    << setprecision(3 + extra) << atmMars.heightOffset << comma
    << setprecision(3 + extra) << atmMars.localHeightOffset << comma
    << setprecision(3 + extra) << atmMars.dustOpticalDepth << comma
    << setprecision(3 + extra) << atmMars.dustColumnArealDensity << comma
//    << scientific
    << setprecision(3 + extra) << atmMars.dustMixingRatio << comma
//    << fixed
    << setprecision(3 + extra) << atmMars.dustMassDensity << comma
    << setprecision(3 + extra) << atmMars.dustNumberDensity << comma
    << atmMars.iceIsPresent << comma
    << setprecision(1 + extra) << atmMars.wavePerturbation * 100.0;
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void MarsProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const MarsAtmosphere& mars = static_cast<const MarsAtmosphere&>(atmos);

  gramMDFile << "# " << mars.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

//! \brief Print method for the GRAM_MD_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printPlanetListMDHeader(const PerturbedAtmosphere& atmos)
{
  const MarsInputParameters& mars = static_cast<const MarsInputParameters&>(atmos.getInputParameters());
  isMolaHeights = mars.isMolaHeights;

  Ephemeris ephem;
  ephem.setBody(MARS);
  ephem.setTime(atmos.getStartTime());
  ephem.update();

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int mark = 0;

  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32s | %-26s | %-11.2f |\n",
    "Input Heights", (mars.isPlanetoCentric ? "Planetocentric" : "Planetographic"),
    "Ref Ellips Equat Rad (km)", atmos.getEquatorialRadius());
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32s | %-26s | %-11.2f |\n",
    "Relative to", (isMolaHeights ? "MOLA areoid" : "Reference ellipsoid"),
    "Ref Ellips Polar Rad (km)", atmos.getPolarRadius());
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32.4f | %-26s | %-11.2f |\n",
    "Dust Nu", mars.dustNu,
    "Pert Wave Length Scale", mars.perturbationWaveLengthScale);
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32.2e | %-26s | %-11.2f |\n",
    "Dust Diameter (m)", mars.dustDiameter,
    "Mean Winds Scale", mars.meanWindsScale);
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32.1f | %-26s | %-11.2f |\n",
    "Dust Den (kg/m^3)", mars.dustDensity,
    "Boundary Layer Winds Scale", mars.boundaryLayerWindsScale);
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-17s | %-32.1f | %-26s | %-11.1f |\n\n",
    "F10.7 Flux 1 AU", mars.F107,
    "F10.7 Flux at Mars ", mars.F107 / square(ephem.getEphemerisState().orbitalRadius));

  gramMDFile << buffer;
  buffer[0] = 0;
  mark = 0;

  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27s | %-26s | %-11s |\n",
    "Field", "Value", "Field", "Value");

  mark += snprintf(buffer + mark, bufferSize - mark, "|%.24s|%.29s|%.28s|%.13s|\n",
    dashes, dashes, dashes, dashes);

  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27.2f | %-26s | %-11.2f |\n",
    "Dust Storm LS (deg)", mars.stormLongitudeSun,
    "Storm Duration (deg)", mars.stormDuration);
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27.2f | %-26s | %-11.2f |\n",
    "Storm Intensity", mars.stormIntensity,
    "Storm MaxRadius (km)", mars.stormMaxRadius);
   
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27.2f | %-26s | %-11.2f |\n",
    "Storm Latitude (deg)", mars.stormLatitude,
    "Storm Longitude (deg)", mars.stormLongitude);
   
  ostringstream trio;
  trio << fixed << setprecision(3) << mars.waveAmplitude1 << ", " << mars.wavePhase1 << ", " << mars.wavePhase1Rate;
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27s | %-26s | %-11.2f |\n",
    "Wave 1 (A, P, R)", trio.str().c_str(),
    "Wave Date", mars.waveDate);
   
  trio.clear();
  trio.str("");
  trio << fixed << setprecision(3) << mars.waveAmplitude2 << ", " << mars.wavePhase2 << ", " << mars.wavePhase2Rate;
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27s | %-26s | %-11.2f |\n",
    "Wave 2 (A, P, R)", trio.str().c_str(),
    "Wave Mean Offset", mars.waveMeanOffset);
   
  trio.clear();
  trio.str("");
  trio << fixed << setprecision(3) << mars.waveAmplitude3 << ", " << mars.wavePhase3 << ", " << mars.wavePhase3Rate;
  mark += snprintf(buffer + mark, bufferSize - mark, "| %-22s | %-27s | %-26s | %-11.2f |\n",
    "Wave 3 (A, P, R)", trio.str().c_str(),
    "Wave Scale", mars.waveScale);

  gramMDFile << buffer;
}

//! \brief Print height for the GRAM_MD_STYLE.
//!
//! This routine prints height relative to the MOLA areoid or the reference ellipsoid.
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printHeightMDStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  char buffer[bufferSize];
  buffer[0] = 0;
  if (isMolaHeights) {
    snprintf(buffer, bufferSize, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
      "Height Above MOLA Areoid (km)", pos.height,
      "MOLA Areoid Radius (km)", pos.latitudeRadius);
  }
  else {
    const MarsAtmosphereState& atmMars = data.atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();
    snprintf(buffer, bufferSize, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
      "Height Above Ref. Ellipsoid (km)", atmMars.referenceHeight,
      "Reference Ellipsoid Radius (km)", atmMars.referenceRadius);
  }
  gramMDFile << buffer;
}

//! \brief Print method for the GRAM_MD_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printPlanetListMDStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  Position pos = data.position;

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
   "Height Above MOLA Surface (km)", pos.height - pos.surfaceHeight,
   "Topographic Height (km)", pos.surfaceHeight);

  if (isMolaHeights) {
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
      "Height Above Ref. Ellipsoid (km)", atmMars.referenceHeight,
      "Reference Ellipsoid Radius (km)", atmMars.referenceRadius);
  }
  else {
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
      "Height Above MOLA Areoid (km)", pos.height,
      "MOLA Areoid Radius (km)", pos.latitudeRadius);
  }

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Planetographic Height (km)", atmMars.planetoGraphicHeight,
    "Planetographic Latitude (deg)", atmMars.planetoGraphicLatitude);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.1f | %-33s | %-10.3f |\n",
    "Ground Temperature (K)", atmMars.groundTemperature,
    "Height Offset (km)", atmMars.heightOffset);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Surface Albedo", atmMars.albedo,
    "Local Height Offset (km)", atmMars.localHeightOffset);

  if (pos.height > 80.0) {
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.1f | %-33s | %-10.3f |\n",
      "F1 Peak Height (km)", atmMars.f1PeakHeight,
      "Thermosphere Base Height (km)", atmMars.thermosphereBaseHeight);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.1f | %-33s | %-10.1f |\n",
      "Exospheric Temperature (K)", atmMars.exosphericTemperature,
      "Thermosphere Base Temperature (K)", atmMars.thermosphereBaseTemperature);
  }

  snprintf(buffer + marker, bufferSize - marker, "\n");
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10s | %-33s | %-10s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.35s|%.12s|%.35s|%.12s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.4e |\n",
    "Daily Mean Temperature (K)", atmMars.temperatureDaily,
    "Daily Mean Density (kg/m^3)", atmMars.densityDaily);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.4e |\n",
    "Min Daily Temperature (K)", atmMars.temperatureMin,
    "Min Daily Density (kg/m^3)", atmMars.densityMin);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.4e |\n",
    "Max Daily Temperature (K)", atmMars.temperatureMax,
    "Max Daily Density (kg/m^3)", atmMars.densityMax);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10d |\n\n",
    "Daily Mean Pressure (Pa)", atmMars.pressureDaily,
    "Ice Is Present", atmMars.iceIsPresent);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10s | %-33s | %-10s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.35s|%.12s|%.35s|%.12s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3e |\n",
    "Dust Optical Depth", atmMars.dustOpticalDepth,
    "Dust Column Areal Density(kg/m^2)", atmMars.dustColumnArealDensity);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.3e |\n",
    "Dust Mass Density (ug/m^2)", atmMars.dustMassDensity,
    "Dust Mixing Ratio (dust/air)", atmMars.dustMixingRatio);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.4f |\n",
    "Dust Number Density (#/m^3)", atmMars.dustNumberDensity,
    "cos(Solar Zenith Angle)", cos(toRadians(data.ephem.solarZenithAngle)));

  gramMDFile << buffer;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void MarsProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
  gramColumnFile << "     Time   Height   Lat   Lon"
    << eastWestChar() << "    Denkgm3   Temp  EWind  NWind  sigD  Ls   Dust  LTST CO2%m  N2%m  Ar%m  O2%m  CO%m   O%m  He%m  H2%m   H%m H2O%m  DensP\n";
}

//! \copydoc ProfilePrinter::printGramColumnStyle()
void MarsProfilePrinter::printGramColumnStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  //  scaleDensityForOutput(atm);

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%9.0f.%7.2f%7.2f%8.2f%10.3e%7.1f%7.1f%7.1f%6.1f%6.1f%5.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%7.3f\n",
    pos.elapsedTime, pos.height, pos.latitude, pos.getLongitude(eastLongitudePositive),
    atm.density, atm.temperature, atm.ewWind, atm.nsWind, atm.densityStandardDeviation * 100.0,
    eph.longitudeSun, atmMars.dustOpticalDepth, eph.solarTime,
    atm.carbonDioxide.massFraction * 100.0,
    atm.dinitrogen.massFraction * 100.0,
    atm.argon.massFraction * 100.0,
    atm.dioxygen.massFraction * 100.0,
    atm.carbonMonoxide.massFraction * 100.0,
    atm.oxygen.massFraction * 100.0,
    atm.helium.massFraction * 100.0,
    atm.dihydrogen.massFraction * 100.0,
    atm.hydrogen.massFraction * 100.0,
    atm.water.massFraction * 100.0,
    atm.perturbedDensity / atm.density);

  gramColumnFile << buffer;
}

//! \copydoc ProfilePrinter::printGramListHeader()
void MarsProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
  const MarsInputParameters& params = static_cast<const MarsInputParameters&>(atmos.getInputParameters());
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);
  gramListFile << " " << MarsAtmosphere::getVersionString() << newline;
  gramListFile << " Input time is " << (timeFrame == ERT ? "Earth-Receive Time (ERT)" : "Planet Event Time (PET)") << newline;
  gramListFile << " Input time is ";
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
  gramListFile << "\n Date = " << month << "/" << day << "/" << year
    << "   Julian Day = " << fixed << jdate
    << " Time = " << setw(2) << setfill('0') << hour << ":" << setw(2) << minute << ":" << setw(2) << rint(seconds) << setfill(' ') << newline;
  gramListFile << " Input heights are planeto-" << (params.isPlanetoCentric ? "centric" : "graphic") 
               << ", relative to " << (params.isMolaHeights ? "MOLA areoid" : "reference ellipsoid") << newline;
  gramListFile << " Reference ellipsoid radii (km): Equator = " << setw(7) << setprecision(2) << atmos.getEquatorialRadius() 
               << " Pole = " << setw(7) << setprecision(2) << atmos.getPolarRadius() << newline;
  gramListFile << " Output heights are planeto-centric, except as noted." << newline;
  gramListFile << " Longitude & ephemeris use IAU 2000 rotational system." << newline;
  if (params.stormLongitudeSun > 0.0) {
    if (abs(params.stormMaxRadius) > 0.0) {
      gramListFile << " Local scale dust storm, starting at Ls = " << setw(5) << setprecision(1) << params.stormLongitudeSun
        << " deg.,  Intensity =" << setw(4) << params.stormIntensity << newline;
      gramListFile << "  with duration = " << setw(5) << params.stormDuration << " degrees of Ls angle." << newline;
      gramListFile << " Max. radius = " << setw(7) << params.stormMaxRadius
        << " km,  At Lat-Lon = " << setw(6) << setprecision(2) << params.stormLatitude
        << " N,  " << setw(6) << params.stormLongitude 
        << " W (" << setw(6) << 360.0 - params.stormLongitude << " E)" << newline;
    }
    else {
      gramListFile << " Global scale dust storm, starting at Ls = " << setw(5) << setprecision(1) << params.stormLongitudeSun
        << " deg.,  Intensity =" << setw(4) << params.stormIntensity << newline;
      gramListFile << "  with duration = " << setw(5) << params.stormDuration << " degrees of Ls angle." << newline;
    }
  }
  gramListFile << " F10.7 flux = " << setw(5) << setprecision(1) << params.F107 << " (1 AU) " << newline;
  if (params.mapYear == 0) {
    gramListFile << " Dust optical depth from NAMELIST input" << newline;
  }
  else {
    gramListFile << " Dust optical depth vs lat and Ls from TES Mapping Year " << params.mapYear << " data" << newline;
  }
  gramListFile << " Dustnu =" << setw(7) << setprecision(4) << params.dustNu 
    << "   Dustdiam =" << setw(6) << setprecision(2) << params.dustDiameter 
    << " E-6 meters   Dustdens =" << setw(8) << setprecision(1) << params.dustDensity << " kg/m**3" << newline;
  gramListFile << "   Random seed =  " << atmos.getSeed() 
    << "  Dens.Pert.Scale Factor = " << setprecision(2) << atmos.getDensityPerturbationScale() 
    << "   corlmin = " << setprecision(3) << atmos.getMinRelativeStepSize() << newline;
  gramListFile << "   Wind.Pert.Scale Factor =" << setw(6) << setprecision(2) << params.ewWindPerturbationScale
    << "    Wavelength Scale Factor =" << setw(6) << params.perturbationWaveLengthScale << newline;
  gramListFile << "   Mean Wind Scale Factor =" << setw(6) << params.meanWindsScale
    << "    Slope Wind Scale Factor =" << setw(6) << params.boundaryLayerWindsScale << newline;
  gramListFile << " A0,A1,phi1,A2,phi2,A3,phi3=" << setw(6) << setprecision(3) << params.waveMeanOffset
    << setw(6) << setprecision(3) << params.waveAmplitude1 << setw(7) << setprecision(1) << params.wavePhase1
    << setw(6) << setprecision(3) << params.waveAmplitude2 << setw(7) << setprecision(1) << params.wavePhase2
    << setw(6) << setprecision(3) << params.waveAmplitude3 << setw(7) << setprecision(1) << params.wavePhase3
    << newline;
  if (params.waveDate > 0.0) {
    gramListFile << " Traveling wave phases initialized at Julian DAY =" << setw(13) << setprecision(3) << params.waveDate << newline;
    gramListFile << " at phase rates (phi1dot,phi2dot,phi3dot, deg/day)=" << setw(8) << params.wavePhase1Rate
      << setw(8) << params.wavePhase2Rate << setw(8) << params.wavePhase3Rate << newline;
  }
  gramListFile << "   Wave Scale =" << setw(8) << setprecision(1) << params.waveScale 
    << " km.    Wave phases are in degrees of " << eastWestName() << " Longitude" << newline;
}

//! \copydoc ProfilePrinter::printGramListStyle()
void MarsProfilePrinter::printGramListStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  scaleDensityForOutput(atm);
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  string densUnits = "kg/m**3";

  greal fmm = data.atmos.minMaxFactor;
  if (fmm < -900.0) fmm = 2.0;

  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker,
    " Time (rel. to T0) =%10.1f sec. (%8.3f sols)  Ls =%6.1f  Dust =%5.2f\n",
    pos.elapsedTime, pos.elapsedTime / eph.secondsPerSol, eph.longitudeSun, atmMars.dustOpticalDepth);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Height Above MOLA (or Surface) =%8.3f km (%8.3f km)  OWLT =%6.2f Min\n",
    pos.height, pos.height - pos.surfaceHeight, eph.oneWayLightTime);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Topographic Height =%9.3f km   Radius (Areoid) =%9.3f (%8.3f) km\n",
    pos.surfaceHeight, pos.totalRadius, pos.latitudeRadius);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Hgt Above Ellipsoid = %7.3f km   Scale Hgt H(p)=%7.2f H(rho)=%7.2f km\n",
    pos.totalRadius - atmMars.referenceRadius, atm.pressureScaleHeight, atm.densityScaleHeight);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Height Offset Parameters:   ibougher =%2d    Local Height Offset =%7.3f km\n",
    ibougher, atmMars.localHeightOffset);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Planeto-Centric Lat = %7.2f deg  Longitude = %7.2f W (%7.2f E) deg.\n",
    pos.latitude, pos.getLongitude(WEST_POSITIVE), pos.getLongitude(EAST_POSITIVE));
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Planeto-Graphic Lat = %7.2f deg  Planeto-Graphic Hgt (Ellps)= %8.3f km\n",
    atmMars.planetoGraphicLatitude, atmMars.planetoGraphicHeight);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Planeto-Cent Sun Lat = %6.2f deg  Mars Orbital Radius =%6.3f AU\n",
    eph.subsolarLatitude, eph.orbitalRadius);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Sun Longitude = %7.2f deg.W      Local True Solar Time = %6.2f Mars hrs\n",
    eph.getSubsolarLongitude(WEST_POSITIVE), eph.solarTime);

  gramListFile << buffer;
  buffer[0] = 0;
  marker = 0;

  if (pos.height > 80.0) {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " Exospheric Temp. = %6.1f K        Tbase = %6.1f K    Zbase = %6.1f km\n",
      atmMars.exosphericTemperature, atmMars.thermosphereBaseTemperature, atmMars.thermosphereBaseHeight);
    marker += snprintf(buffer + marker, bufferSize - marker,
      " Solar Zenith Angle =%6.1f deg     F1 Peak = %6.1f km\n",
      eph.solarZenithAngle, atmMars.f1PeakHeight);
  }
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Temperature = %7.1f K        Pressure = %9.3e N/m**2   profwgt =%6.3f\n",
    atm.temperature, atm.pressure, atm.profileWeight[0]);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Density (Low, Avg., High) =%12.3e%12.3e%12.3e %s\n",
    atm.lowDensity, atm.density, atm.highDensity, densUnits.c_str());
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Departure, COSPAR NH Mean =%9.1f %%%10.1f %%%10.1f %%  iupdate = %1d\n",
    atm.lowDensityDeviation * 100.0, atm.densityDeviation * 100.0, atm.highDensityDeviation * 100.0, iupdate);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Tot.Dens. =%10.3e %7s   Dens.Pert. =%7.2f%% Wave =%5.2f%% of mean \n",
    atm.perturbedDensity, densUnits.c_str(), atm.densityPerturbation * 100.0, atmMars.wavePerturbation * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Eastward Wind (Mean,Perturbed,Total) = %6.1f%7.1f%7.1f m/s    VertWind\n",
    atm.ewWind, atm.ewWindPerturbation, atm.perturbedEWWind);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Northward Wind(Mean,Perturbed,Total) = %6.1f%7.1f%7.1f m/s%8.1f m/s\n",
    atm.nsWind, atm.nsWindPerturbation, atm.perturbedNSWind, atm.perturbedVerticalWind);

  gramListFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker,
    " CO2=%9.3e N2=%9.3e Ar=%9.3e O2=%9.3e CO=%9.3e #/m**3\n",
    atm.carbonDioxide.numberDensity, atm.dinitrogen.numberDensity, atm.argon.numberDensity,
    atm.dioxygen.numberDensity, atm.carbonMonoxide.numberDensity);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%14.3f%13.3f%13.3f%13.3f%13.3f %% by mass\n",
    atm.carbonDioxide.massFraction * 100.0, atm.dinitrogen.massFraction * 100.0,
    atm.argon.massFraction * 100.0, atm.dioxygen.massFraction * 100.0, atm.carbonMonoxide.massFraction * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%14.3f%13.3f%13.3f%13.3f%13.3f %% by volume\n",
    atm.carbonDioxide.moleFraction * 100.0, atm.dinitrogen.moleFraction * 100.0,
    atm.argon.moleFraction * 100.0, atm.dioxygen.moleFraction * 100.0, atm.carbonMonoxide.moleFraction * 100.0);

  gramListFile << buffer;
  buffer[0] = 0;
  marker = 0;

  if (pos.height > 80.0) {
    marker += snprintf(buffer + marker, bufferSize - marker,
      "   O=%9.3e He=%9.3e H2=%9.3e  H=%9.3e   Total=%9.3e #/m**3\n",
      atm.oxygen.numberDensity,  atm.helium.numberDensity,
      atm.dihydrogen.numberDensity, atm.hydrogen.numberDensity, atm.totalNumberDensity);
    marker += snprintf(buffer + marker, bufferSize - marker,
      "%14.3f%13.3f%13.3f%13.3f %% by mass  MolWgt=%6.3f\n",
      atm.oxygen.massFraction * 100.0, atm.helium.massFraction * 100.0, atm.dihydrogen.massFraction * 100.0,
      atm.hydrogen.massFraction * 100.0, atm.averageMolecularWeight);
    marker += snprintf(buffer + marker, bufferSize - marker,
      "%14.3f%13.3f%13.3f%13.3f %% volume (mole) fraction\n",
      atm.oxygen.moleFraction * 100.0, atm.helium.moleFraction * 100.0,
      atm.dihydrogen.moleFraction * 100.0, atm.hydrogen.moleFraction * 100.0);
  }
  else {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " H2O=%9.3e  O=%9.3e He=%9.3e H2=%9.3e   Total=%9.3e #/m**3\n",
      atm.water.numberDensity, atm.oxygen.numberDensity, atm.helium.numberDensity,
      atm.dihydrogen.numberDensity, atm.totalNumberDensity);
    marker += snprintf(buffer + marker, bufferSize - marker,
      "%14.3f%13.3f%13.3f%13.3f %% by mass  MolWgt=%6.3f\n",
      atm.water.massFraction * 100.0, atm.oxygen.massFraction * 100.0,
      atm.helium.massFraction * 100.0, atm.dihydrogen.massFraction * 100.0, atm.averageMolecularWeight);
    marker += snprintf(buffer + marker, bufferSize - marker,
      "%14.3f%13.3f%13.3f%13.3f %% volume (mole) fraction\n",
      atm.water.moleFraction * 100.0, atm.oxygen.moleFraction * 100.0,
      atm.helium.moleFraction * 100.0, atm.dihydrogen.moleFraction * 100.0);
  }
  marker += snprintf(buffer + marker, bufferSize - marker,
    " -----------------------------------------------------------------------------\n");

  gramListFile << buffer;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void MarsProfilePrinter::printGramTPresHgtHeader(const PerturbedAtmosphere& atmos)
{
  gramTPresHgtFile << "   Height       Temp      Pres   TdegC   Pres_mb    Hrho  Hpres MolWt TerHgt Tgrnd  Areoid  dAreoid CO2%v  N2%v  Ar%v  O2%v  CO%v   O%v  He%v  H2%v   H%v  H2O%v LOGSCALE\n";
}

//! \copydoc ProfilePrinter::printGramTPresHgtStyle()
void MarsProfilePrinter::printGramTPresHgtStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  const Position& pos = data.position;

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%13.3f%7.1f%11.3e%7.1f%11.3e%8.2f%8.2f%6.2f%7.3f%6.1f%9.3f%8.3f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%5d\n",
    data.position.height, atm.temperature, atm.pressure,
    atm.temperature - 273.15, atm.pressure * 0.01,
    atm.densityScaleHeight, atm.pressureScaleHeight,
    atm.averageMolecularWeight,
    pos.surfaceHeight, atmMars.groundTemperature,
    pos.latitudeRadius,  pos.latitudeRadius - atmMars.referenceRadius,
    atm.carbonDioxide.moleFraction * 100.0,
    atm.dinitrogen.moleFraction * 100.0,
    atm.argon.moleFraction * 100.0,
    atm.dioxygen.moleFraction * 100.0,
    atm.carbonMonoxide.moleFraction * 100.0,
    atm.oxygen.moleFraction * 100.0,
    atm.helium.moleFraction * 100.0,
    atm.dihydrogen.moleFraction * 100.0,
    atm.hydrogen.moleFraction * 100.0,
    atm.water.moleFraction * 100.0,
    densityPrintScale);

  gramTPresHgtFile << buffer;
}

//! \brief Header method for the GRAM_DENSITY_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printGramDensityHeader(const PerturbedAtmosphere& atmos)
{
  gramDensityFile << "   Height        DENSLO     DENSAV     DENSHI    DENSTOT  DustOD  Radius   Grav RadAU LOGSCALE hgtoffset ibougher MapYear profwgt\n";
  const MarsInputParameters& params = static_cast<const MarsInputParameters&>(atmos.getInputParameters());
  ibougher = params.offsetModel;
  mapyear = params.mapYear;
}

//! \brief Print method for the GRAM_DENSITY_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printGramDensityStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  scaleDensityForOutput(atm);

  greal profwgt = data.atmos.profileWeight[0];
  if (profwgt < -999.0) profwgt = 0.0;

  if (pos.height < 10.0) {
    gramDensityFile << fixed << setw(9) << setprecision(4) << pos.height;
  }
  else if (pos.height < 100.0) {
    gramDensityFile << fixed << setw(9) << setprecision(3) << pos.height;
  }
  else {
    gramDensityFile << fixed << setw(9) << setprecision(2) << pos.height;
  }
  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%15.3e%11.3e%11.3e%11.3e%7.4f%9.3f%7.3f%6.3f%4d%13.3f%6d%8d%9.3f\n",
    atm.lowDensity,
    atm.density,
    atm.highDensity,
    atm.perturbedDensity,
    atmMars.dustOpticalDepth,
    pos.totalRadius,
    pos.gravity,
    eph.orbitalRadius,
    densityPrintScale,
    atmMars.localHeightOffset,
    ibougher,
    mapyear,
    profwgt);

  gramDensityFile << buffer;
}

//! \brief Header method for the GRAM_WINDS_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printGramWindsHeader(const PerturbedAtmosphere& atmos)
{
  gramWindsFile << "   Height      EWmean  EWpert   EWtot  NSmean  NSpert   NStot  VWpert iupdate\n";
}

//! \brief Print method for the GRAM_WINDS_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printGramWindsStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  greal height = data.position.height;
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;
  char buffer[bufferSize];
  buffer[0] = 0;

  int prec = 2;
  if (height < 10.0) prec = 4;
  else if (height < 100.0) prec = 3;

  snprintf(buffer, bufferSize,
    "%9.*f%12.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%3d\n",
    prec, data.position.height,
    atm.ewWind, atm.ewWindPerturbation, atm.perturbedEWWind,
    atm.nsWind, atm.nsWindPerturbation, atm.perturbedNSWind,
    atm.perturbedVerticalWind, iupdate);

  gramWindsFile << buffer;
}

//! \brief Header method for the GRAM_PERTURB_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printGramPerturbHeader(const PerturbedAtmosphere& atmos)
{
  gramPerturbFile << "   Height      SigD  DensRand  DensWave     Densp    corlim   SigU   SigW iupdate\n";
}

//! \brief Print method for the GRAM_PERTURB_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printGramPerturbStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;
  char buffer[bufferSize];
  buffer[0] = 0;
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  int prec = 0;
  if (pos.height < 10.0) {
    prec = 4;
  }
  else if (pos.height < 100.0) {
    prec = 3;
  }
  else {
    prec = 2;
  }

  snprintf(buffer, bufferSize,
    "%9.*f%10.2f%10.3f%10.3f%10.3f%10.3e%7.2f%7.2f%5d\n",
    prec, pos.height,
    atm.densityStandardDeviation * 100.0, atm.densityPerturbation * 100.0,
    0.0, atm.densityPerturbation * 100.0,
    data.atmos.relativeStepSize, atm.ewStandardDeviation,
    atm.verticalStandardDeviation, iupdate);

  gramPerturbFile << buffer;
}

//! \brief Header method for the GRAM_PERTURB_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printDayDataHeader(const PerturbedAtmosphere& atmos)
{
  dayDataFile << "    Height    TempDay  PresDay    DensDay   EWwnDay NSwnDay Tempmin Tempmax   Densmin    Densmax  LOGSCALE  DensAV\n";
}

//! \brief Print method for the GRAM_PERTURB_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printDayDataStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  const Position& pos = data.position;

  // If the temp is zero, then no daily data is available at that height.
  if (atmMars.temperatureDaily == 0.0) {
    return;
  }

  int prec = 0;
  if (pos.height < 10.0) {
    prec = 4;
  }
  else if (pos.height < 100.0) {
    prec = 3;
  }
  else {
    prec = 2;
  }

  char buffer[bufferSize];
  buffer[0] = 0;
  snprintf(buffer, bufferSize,
    "%9.*f%11.1f%11.3e%11.3e%8.1f%8.1f%8.1f%8.1f%11.3e%11.3e%5d%14.3e\n",
    prec, pos.height,
    atmMars.temperatureDaily, atmMars.pressureDaily, atmMars.densityDaily,
    atmMars.ewWindDaily, atmMars.nsWindDaily,
    atmMars.temperatureMin, atmMars.temperatureMax,
    atmMars.densityMin, atmMars.densityMax,
    densityPrintScale, atm.density);

  dayDataFile << buffer;
}

//! \brief Header method for the GRAM_PERTURB_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printMarsRadHeader(const PerturbedAtmosphere& atmos)
{
  marsRadFile << "   Height       alb      mu0 Dareaden  Dmixrat  Dmasden  Dnumden Ice\n";
}

//! \brief Print method for the GRAM_PERTURB_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printMarsRadStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  const Position& pos = data.position;
  const EphemerisState& ephem = data.ephem;

  int prec = 0;
  if (pos.height < 10.0) {
    prec = 4;
  }
  else if (pos.height < 100.0) {
    prec = 3;
  }
  else {
    prec = 2;
  }

  char buffer[bufferSize];
  buffer[0] = 0;
  snprintf(buffer, bufferSize,
    "%9.*f%10.3f%9.5f%9.2e%9.2e%9.2e%9.2e%3d\n",
    prec, pos.height, atmMars.albedo,
    cos(toRadians(ephem.solarZenithAngle)),
    atmMars.dustColumnArealDensity,
    atmMars.dustMixingRatio,
    atmMars.dustMassDensity,
    atmMars.dustNumberDensity,
    atmMars.iceIsPresent);

  marsRadFile << buffer;
}

//! \brief Header method for the GRAM_PERTURB_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void MarsProfilePrinter::printThrmDataHeader(const PerturbedAtmosphere& atmos)
{
  thrmDataFile << "   Height       Tbase   Zbase  F1peak  MolWgt   Texos  hgtoffset ibougher\n";
}

//! \brief Print method for the GRAM_PERTURB_STYLE.
//!
//! \param data The ProfileData to be printed.
void MarsProfilePrinter::printThrmDataStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const MarsAtmosphereState& atmMars = atm.getPlanetSpecificMetrics<MarsAtmosphereState>();
  const Position& pos = data.position;

  int prec = 0;
  if (pos.height < 10.0) {
    prec = 4;
  }
  else if (pos.height < 100.0) {
    prec = 3;
  }
  else {
    prec = 2;
  }

  char buffer[bufferSize];
  buffer[0] = 0;
  snprintf(buffer, bufferSize,
    "%9.*f%12.1f%8.1f%8.1f%8.2f%8.1f%10.3f%6d\n",
    prec, pos.height,
    atmMars.thermosphereBaseTemperature, atmMars.thermosphereBaseHeight,
    atmMars.f1PeakHeight, atm.averageMolecularWeight,
    atmMars.exosphericTemperature, atmMars.heightOffset, ibougher);

  thrmDataFile << buffer;
}

//! \copydoc ProfilePrinter::openOutput()
void MarsProfilePrinter::openOutput()
{
  ProfilePrinter::openOutput();
  if (outputStyle & MARS_DAY_STYLE) {
    openFile(dayDataFile, fileNamePrefix + "DayData.txt");
  }
  if (outputStyle & MARS_RAD_STYLE) {
    openFile(marsRadFile, fileNamePrefix + "MarsRad.txt");
  }
  if (outputStyle & MARS_THRM_STYLE) {
    openFile(thrmDataFile, fileNamePrefix + "ThrmData.txt");
  }
}

//! \copydoc ProfilePrinter::printFileHeader()
void MarsProfilePrinter::printFileHeader(const PerturbedAtmosphere& atmos)
{
  ProfilePrinter::printFileHeader(atmos);
  if (outputStyle & MARS_DAY_STYLE) {
    printDayDataHeader(atmos);
  }
  if (outputStyle & MARS_RAD_STYLE) {
    printMarsRadHeader(atmos);
  }
  if (outputStyle & MARS_THRM_STYLE) {
    printThrmDataHeader(atmos);
  }
}

//! \copydoc ProfilePrinter::printData()
void MarsProfilePrinter::printData(const std::vector<ProfileData>& profile)
{
  ProfilePrinter::printData(profile);
  if (outputStyle & MARS_DAY_STYLE) {
    for (auto& p : profile) {
      printDayDataStyle(p);
    }
  }
  if (outputStyle & MARS_RAD_STYLE) {
    for (auto& p : profile) {
      printMarsRadStyle(p);
    }
  }
  if (outputStyle & MARS_THRM_STYLE) {
    for (auto& p : profile) {
      printThrmDataStyle(p);
    }
  }
}

//! \copydoc ProfilePrinter::closeOutput()
void MarsProfilePrinter::closeOutput()
{
  ProfilePrinter::closeOutput();
  if (outputStyle & MARS_DAY_STYLE) {
    dayDataFile.close();
  }
  if (outputStyle & MARS_RAD_STYLE) {
    marsRadFile.close();
  }
  if (outputStyle & MARS_THRM_STYLE) {
    thrmDataFile.close();
  }
}

} // namespace
