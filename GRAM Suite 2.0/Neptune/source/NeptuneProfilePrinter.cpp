//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//
// Adapted from Neptune-GRAM 2004 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include "NeptuneProfilePrinter.h"
#include "NeptuneAtmosphere.h"

using namespace std;

namespace GRAM {

//! \copydoc ProfilePrinter::ProfilePrinter()
NeptuneProfilePrinter::NeptuneProfilePrinter()
{
}

//! \fn NeptuneProfilePrinter::NeptuneProfilePrinter(const NeptuneProfilePrinter& orig)
//! \copydoc ProfilePrinter::ProfilePrinter(const ProfilePrinter& orig)

//! \fn NeptuneProfilePrinter::~NeptuneProfilePrinter()
//! \copydoc ProfilePrinter::~ProfilePrinter()

//! \brief Header method for planet specific data in the GRAM_CSV_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void NeptuneProfilePrinter::printPlanetCSVHeader(const PerturbedAtmosphere& atmos)
{
  gramCSVFile << "MinMaxFactor";
}

//! \brief Print method for planet specific data the GRAM_CSV_STYLE.
//!
//! \param data The ProfileData to be printed.
void NeptuneProfilePrinter::printPlanetCSVStyle(const ProfileData& data)
{
  gramCSVFile << fixed << setprecision(3 + extra) << data.atmos.minMaxFactor;
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void NeptuneProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const NeptuneAtmosphere& neptune = static_cast<const NeptuneAtmosphere&>(atmos);

  gramMDFile << "# " << neptune.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

//! \brief Print method for the GRAM_MD_STYLE.
//!
//! \param data The ProfileData to be printed.
void NeptuneProfilePrinter::printPlanetListMDStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  char buffer[bufferSize];
  buffer[0] = 0;
  snprintf(buffer, bufferSize, "| %-33s | %-10.3f | %-33s | %-10s |\n",
    "Min/Max Factor", atm.minMaxFactor,
    " ",  " ");
  gramMDFile << buffer;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void NeptuneProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
  gramColumnFile << "     Time   Height   Lat   Lon"
    << eastWestChar() << "  Denkgm3    Temp  EWind  NWind  sigD  Ls     Fmm  H2%m He%m CH4%m N2%m";
  gramColumnFile << newline;
}

//! \copydoc ProfilePrinter::printGramColumnStyle()
void NeptuneProfilePrinter::printGramColumnStyle(const ProfileData& data)
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
    << setw(6) << setprecision(1) << eph.longitudeSun
    << setw(6) << setprecision(2) << atm.minMaxFactor;  // fmm
  gramColumnFile << setw(7) << setprecision(3) << atm.dihydrogen.massFraction * 100.0;
  gramColumnFile << setw(7) << setprecision(3) << atm.helium.massFraction * 100.0;
  gramColumnFile << setw(7) << setprecision(3) << atm.methane.massFraction * 100.0;
  gramColumnFile << setw(7) << setprecision(3) << atm.dinitrogen.massFraction * 100.0;
  gramColumnFile << newline;
}

//! \copydoc ProfilePrinter::printGramListHeader()
void NeptuneProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);
  gramListFile << " " << NeptuneAtmosphere::getVersionString() << newline;
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
    gramListFile << "   Random seed =  " << atmos.getSeed() << "   Scale factor = " << setprecision(1) << atmos.getDensityPerturbationScale() << "   corlmin = " << setprecision(3) << atmos.getMinRelativeStepSize() << newline;
}

//! \copydoc ProfilePrinter::printGramListStyle()
void NeptuneProfilePrinter::printGramListStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, 
    " Time (rel. to T0) =%11.1f sec. (%7.3f Neptune Days)    Ls =%6.1f\n", 
    pos.elapsedTime, pos.elapsedTime / max(eph.secondsPerSol, (greal)1.0), eph.longitudeSun);
  marker += snprintf(buffer + marker, bufferSize - marker, 
    " Height Above Reference Ellipsoid =%9.3f km          Fminmax =%7.3f\n", 
    pos.height, atm.minMaxFactor);
  marker += snprintf(buffer + marker, bufferSize - marker, 
    " Total Radius (Ref Radius) = %9.3f (%9.3f) km   OWLT =%7.2f Min\n", 
    pos.totalRadius, pos.latitudeRadius, eph.oneWayLightTime);
  marker += snprintf(buffer + marker, bufferSize - marker, 
    " Scale Heights: H(p) =%6.2f       H(rho) =%6.2f km    zeta =%7.4f\n",
    atm.pressureScaleHeight, atm.densityScaleHeight, data.atmos.compressibilityFactor);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Latitude = %7.2f  degrees       Longitude = %7.2f %1c (%7.2f %1c) deg.\n",
    pos.latitude, pos.getLongitude(eastLongitudePositive), eastWestChar(), 
    pos.getLongitude(!eastLongitudePositive), eastWestChar(true));
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Sun Latitude = %8.2f deg.      Neptune Orbital Radius  =%7.3f AU\n",
    eph.subsolarLatitude, eph.orbitalRadius);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Sun Longitude =%8.2f deg.%1c     Local True Solar Time =%6.2f Neptune hr\n",
    eph.getSubsolarLongitude(eastLongitudePositive), eastWestChar(), eph.solarTime);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Temperature = %7.1f K           Pressure = %12.3e N/m**2   profwgt =%6.3f\n",
    atm.temperature, atm.pressure, atm.profileWeight[0]);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Density (Low, Avg., High) = %12.3e%12.3e%12.3e %s\n",
    atm.lowDensity, atm.density, atm.highDensity, densUnits.c_str());
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Departure from Neptune Avg =     %6.1f %%%10.1f %%%10.1f %%\n",
    atm.lowDensityDeviation * 100.0, atm.densityDeviation * 100.0, atm.highDensityDeviation * 100.0);

  gramListFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker,
    " Tot.Dens. =%10.3e %6s      Dens.Pert. =%7.2f %% of mean  iupdate= %1d\n",
    atm.perturbedDensity, densUnits.c_str(), atm.densityPerturbation * 100.0, iupdate);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Eastward Wind  (Mean, Perturbed, Total) = %7.1f%7.1f%7.1f m/s\n",
    atm.ewWind, atm.ewWindPerturbation, atm.perturbedEWWind);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Northward Wind (Mean, Perturbed, Total) = %7.1f%7.1f%7.1f m/s\n",
    atm.nsWind, atm.nsWindPerturbation, atm.perturbedNSWind);
  if (atm.dinitrogen.moleFraction != 0) {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " N2 =%10.3e #/m**3 %9.3f%% by mass %7.3f%% mole fraction\n",
      atm.dinitrogen.numberDensity, atm.dinitrogen.massFraction * 100.0, atm.dinitrogen.moleFraction * 100.0);
  }
  marker += snprintf(buffer + marker, bufferSize - marker,
    " H2 =%10.3e   He =%10.3e  CH4=%10.3e  Total=%10.3e #/m**3\n",
    atm.dihydrogen.numberDensity, atm.helium.numberDensity, atm.methane.numberDensity, atm.totalNumberDensity);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%11.3f%% %14.3f%% %13.3f%%  by mass      MolWgt=%6.3f\n",
    atm.dihydrogen.massFraction * 100.0, atm.helium.massFraction * 100.0, atm.methane.massFraction * 100.0, atm.averageMolecularWeight);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%11.3f%% %14.3f%% %13.3f%%  mole (or volume) fraction\n",
    atm.dihydrogen.moleFraction * 100.0, atm.helium.moleFraction * 100.0, atm.methane.moleFraction * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " ---------------------------------------------------------------------------\n");
  gramListFile << buffer;
}

} // namespace
