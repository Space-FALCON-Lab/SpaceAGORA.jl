//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include "VenusProfilePrinter.h"
#include "VenusAtmosphere.h"

using namespace std;

namespace GRAM {

//! \copydoc ProfilePrinter::ProfilePrinter()
VenusProfilePrinter::VenusProfilePrinter()
{
}

//! \fn VenusProfilePrinter::VenusProfilePrinter(const VenusProfilePrinter& orig)
//! \copydoc ProfilePrinter::ProfilePrinter(const ProfilePrinter& orig)

//! \fn VenusProfilePrinter::~VenusProfilePrinter()
//! \copydoc ProfilePrinter::~ProfilePrinter()

//! \brief Header method for planet specific data in the GRAM_CSV_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void VenusProfilePrinter::printPlanetCSVHeader(const PerturbedAtmosphere& atmos)
{
  gramCSVFile << "HeightAboveSurface_km,SurfaceHeight_km";
}

//! \brief Print method for planet specific data the GRAM_CSV_STYLE.
//!
//! \param data The ProfileData to be printed.
void VenusProfilePrinter::printPlanetCSVStyle(const ProfileData& data)
{
  gramCSVFile << fixed
    << setprecision(3 + extra) << data.position.height - data.position.surfaceHeight << comma
    << setprecision(3 + extra) << data.position.surfaceHeight;
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void VenusProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const  VenusAtmosphere&  venus = static_cast<const  VenusAtmosphere&>(atmos);

  gramMDFile << "# " << venus.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

//! \brief Print method for the GRAM_MD_STYLE.
//!
//! \param data The ProfileData to be printed.
void VenusProfilePrinter::printPlanetListMDStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  char buffer[bufferSize];

  snprintf(buffer, bufferSize, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Height Above Surface (km)", pos.height - pos.surfaceHeight,
    "Surface Height Above Raduis (km)", pos.surfaceHeight);
  gramMDFile << buffer;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void VenusProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
  gramColumnFile << "     Time   Height   Lat    Lon"
    << eastWestChar() << "  Denkgm3    Temp  EWind  NWind  sigD   Ls    SZA CO2%m  N2%m  O%m CO%m He%m  N%m   H%m\n";
}

//! \copydoc ProfilePrinter::printGramColumnStyle()
void VenusProfilePrinter::printGramColumnStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%10.1f%8.2f%7.2f%8.2f%9.2e%7.1f%7.1f%7.1f%6.1f%6.1f%6.1f%6.1f%6.1f%5.1f%5.1f%5.1f%5.1f%6.1f\n",
    pos.elapsedTime, pos.height, pos.latitude, pos.getLongitude(eastLongitudePositive),
    atm.density, atm.temperature, atm.ewWind, atm.nsWind, atm.densityStandardDeviation * 100.0, 
    eph.longitudeSun, data.ephem.solarZenithAngle,
    atm.carbonDioxide.massFraction * 100.0, 
    atm.dinitrogen.massFraction * 100.0,
    atm.oxygen.massFraction * 100.0, 
    atm.carbonMonoxide.massFraction * 100.0,
    atm.helium.massFraction * 100.0, 
    atm.nitrogen.massFraction * 100.0,
    atm.hydrogen.massFraction * 100.0);

  gramColumnFile << buffer;
}


//! \copydoc ProfilePrinter::printGramListHeader()
void VenusProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);
  gramListFile << " " << VenusAtmosphere::getVersionString() << newline;
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
    << "  Julian Day = " << fixed << jdate
    << " Time = " << setw(2) << setfill('0') << hour << ":" << setw(2) << minute << ":" << setw(2) << rint(seconds) << setfill(' ') << newline;
  gramListFile << "   Random seed =  " << atmos.getSeed() << "   Scale factor = " << setprecision(1) << atmos.getDensityPerturbationScale() << "   corlmin = " << setprecision(3) << atmos.getMinRelativeStepSize() << newline;
}

//! \copydoc ProfilePrinter::printGramListStyle()
void VenusProfilePrinter::printGramListStyle(const ProfileData& data)
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
    " Time (rel. to T0) =%11.1f sec. (%7.3f Venus Days)    Ls =%6.1f\n",
    pos.elapsedTime, pos.elapsedTime/eph.secondsPerSol, eph.longitudeSun);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Height Above Reference Ellipsoid =%9.3f km          SZA =%7.2f deg\n",
    pos.height, data.ephem.solarZenithAngle);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Total Radius (Ref Radius) = %9.3f (%9.3f) km   OWLT =%7.2f Min\n",
    pos.totalRadius, pos.latitudeRadius, eph.oneWayLightTime);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Scale Heights: H(p) =%6.2f       H(rho) =%6.2f km    zeta =%7.4f\n",
    atm.pressureScaleHeight, atm.densityScaleHeight, data.atmos.compressibilityFactor);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Latitude = %7.2f  degrees       Longitude = %7.2f E (%7.2f W) deg.\n",
    pos.latitude, pos.getLongitude(EAST_POSITIVE), pos.getLongitude(WEST_POSITIVE));
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Sun Latitude = %8.2f deg.      Venus Orbital Radius  =%7.4f AU\n",
    eph.subsolarLatitude, eph.orbitalRadius);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Sun Longitude =%8.2f deg.E     Local True Solar Time =%6.2f Venus hr\n",
    eph.getSubsolarLongitude(true), eph.solarTime);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Temperature = %7.1f K      Pressure =%10.3e N/m**2   profwgt =%6.3f\n",
    atm.temperature, atm.pressure, atm.profileWeight[0]);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Density (Low, Avg., High) = %12.3e%12.3e%12.3e %s\n",
    atm.lowDensity, atm.density, atm.highDensity, densUnits.c_str());
  marker += snprintf(buffer + marker, bufferSize - marker,
    " Departure from Venus Avg = %10.1f %%%10.1f %%%10.1f %%\n",
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
  marker += snprintf(buffer + marker, bufferSize - marker,
    " CO2 =%10.3e N2 =%10.3e  O =%10.3e CO =%10.3e #/m**3\n",
    atm.carbonDioxide.numberDensity, atm.dinitrogen.numberDensity,
    atm.oxygen.numberDensity, atm.carbonMonoxide.numberDensity);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%16.3f%15.3f%15.3f%15.3f %% by mass\n",
    atm.carbonDioxide.massFraction * 100.0, atm.dinitrogen.massFraction * 100.0, 
    atm.oxygen.massFraction * 100.0, atm.carbonMonoxide.massFraction * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%16.3f%15.3f%15.3f%15.3f %% by volume\n",
    atm.carbonDioxide.moleFraction * 100.0, atm.dinitrogen.moleFraction * 100.0,
    atm.oxygen.moleFraction * 100.0, atm.carbonMonoxide.moleFraction * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "  He =%10.3e  N =%10.3e  H =%10.3e   Total=%10.3e #/m**3\n",
    atm.helium.numberDensity, atm.nitrogen.numberDensity, atm.hydrogen.numberDensity, atm.totalNumberDensity);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%16.3f%15.3f%15.3f %% by mass       MolWgt=%6.3f\n",
    atm.helium.massFraction * 100.0, atm.nitrogen.massFraction * 100.0,
    atm.hydrogen.massFraction * 100.0, atm.averageMolecularWeight);
  marker += snprintf(buffer + marker, bufferSize - marker,
    "%16.3f%15.3f%15.3f %% volume (or mole) fraction\n",
    atm.helium.moleFraction * 100.0, atm.nitrogen.moleFraction * 100.0, atm.hydrogen.moleFraction * 100.0);
  marker += snprintf(buffer + marker, bufferSize - marker,
    " ---------------------------------------------------------------------------\n");

  gramListFile << buffer;
}

//! \copydoc ProfilePrinter::printGramTPresHgtHeader()
void VenusProfilePrinter::printGramTPresHgtHeader(const PerturbedAtmosphere& atmos)
{
  gramTPresHgtFile << "   Height       Temp      Pres   TdegC   Pres_mb    Hrho  Hpres MolWt CO2%v N2%v  O%v CO%v He%v  N%v   H%v LOGSCALE\n";
}

//! \copydoc ProfilePrinter::printGramTPresHgtStyle()
void VenusProfilePrinter::printGramTPresHgtStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);

  if (data.position.height < 10.0) {
    gramTPresHgtFile << fixed << setw(9) << setprecision(4) << data.position.height;
  }
  else if (data.position.height < 100.0) {
    gramTPresHgtFile << fixed << setw(9) << setprecision(3) << data.position.height;
  }
  else {
    gramTPresHgtFile << fixed << setw(9) << setprecision(2) << data.position.height;
  }
  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%11.1f%11.3e%7.1f%11.3e%7.2f%7.2f%6.2f%6.1f%5.1f%5.1f%5.1f%5.1f%5.1f%6.1f%5d\n",
    atm.temperature,
    atm.pressure,
    atm.temperature - 273.15,
    atm.pressure * 0.01,
    atm.densityScaleHeight,
    atm.pressureScaleHeight,
    atm.averageMolecularWeight,
    atm.carbonDioxide.moleFraction * 100.0,
    atm.dinitrogen.moleFraction * 100.0,
    atm.oxygen.moleFraction * 100.0,
    atm.carbonMonoxide.moleFraction * 100.0,
    atm.helium.moleFraction * 100.0,
    atm.nitrogen.moleFraction * 100.0,
    atm.hydrogen.moleFraction * 100.0,
    densityPrintScale);

    gramTPresHgtFile << buffer;
}

} // namespace
