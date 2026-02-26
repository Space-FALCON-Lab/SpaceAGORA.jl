//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include "TitanProfilePrinter.h"
#include "TitanAtmosphere.h"

using namespace std;

namespace GRAM {

TitanProfilePrinter::TitanProfilePrinter()
{
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void TitanProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const TitanAtmosphere& titan = static_cast<const TitanAtmosphere&>(atmos);

  gramMDFile << "# " << titan.getVersionString() << '\n';
  gramMDFile << '\n';

  ProfilePrinter::printGramListMDHeader(atmos);
}

//! \brief Header method for planet specific data in the GRAM_CSV_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void TitanProfilePrinter::printPlanetCSVHeader(const PerturbedAtmosphere& atmos)
{
  gramCSVFile << "MinMaxFactor";
}

//! \brief Print method for planet specific data the GRAM_CSV_STYLE.
//!
//! \param data The ProfileData to be printed.
void TitanProfilePrinter::printPlanetCSVStyle(const ProfileData& data)
{
  gramCSVFile << fixed << setprecision(3 + extra) << data.atmos.minMaxFactor;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void TitanProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
    gramColumnFile << "     Time   Height   Lat   Lon"
      << eastWestChar() << "  Denkgm3    Temp  EWind  NWind  sigD   Ls    FMM N2%m CH4%m Ar%m\n";
}


void TitanProfilePrinter::printGramColumnStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  const AtmosphereState& atm = data.atmos;
//  scaleDensityForOutput(atm);

  greal fmm = data.atmos.minMaxFactor;
  if (fmm < -900.0) fmm = 2.0;

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%10.1f%8.2f%7.2f%7.2f%9.2e%7.1f%7.1f%7.1f%6.1f%6.1f%6.2f%5.1f%5.1f%5.1f\n",
    pos.elapsedTime, pos.height, pos.latitude, pos.getLongitude(eastLongitudePositive),
    atm.density, atm.temperature, atm.ewWind, atm.nsWind,
    atm.densityStandardDeviation * 100.0,
    eph.longitudeSun, fmm,
    atm.dinitrogen.massFraction * 100.0,
    atm.methane.massFraction * 100.0,
    atm.argon.massFraction * 100.0);

  gramColumnFile << buffer;
}

//! \copydoc ProfilePrinter::printGramListHeader()
void TitanProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);
  gramListFile << " " << TitanAtmosphere::getVersionString() << newline;
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

void TitanProfilePrinter::printGramListStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const EphemerisState& eph = data.ephem;
  AtmosphereState atm = data.atmos;
  scaleDensityForOutput(atm);
  int iupdate = 0;
  if (atm.updateStatus == NO_UPDATES) iupdate = -1;
  if (atm.updateStatus == PERTS_UPDATED) iupdate = 1;

  densUnits = "kg/m**3";

  greal fmm = data.atmos.minMaxFactor;
  if (fmm < -900.0) fmm = 2.0;


  gramListFile  << " Time (rel. to T0) ="
          << fixed << setw(11) << setprecision(1) << pos.elapsedTime
          << " sec. (" << setw(7)  << setprecision(3) << pos.elapsedTime/eph.secondsPerSol
          << " Titan Days)      Ls =" << setw(6)  << setprecision(1) << eph.longitudeSun
          << newline;
  gramListFile  << " Height Above Reference Ellipsoid =" << setw(9)  << setprecision(3) << pos.height
          << " km          Fminmax =" << setw(7)  << setprecision(3) << fmm
          << newline;
  gramListFile  << " Total Radius (Ref Radius) = " << setw(8)  << setprecision(3) << pos.totalRadius
          << " (" << setw(8)  << setprecision(3) << pos.latitudeRadius
          << ") km     OWLT =" << setw(7)  << setprecision(2) << eph.oneWayLightTime // RREF
          << " Min\n";
  gramListFile  << " Scale Heights: H(p) =" << setw(6)  << setprecision(2) << atm.pressureScaleHeight
          << "       H(rho) =" << setw(6)  << setprecision(2) << atm.densityScaleHeight
          << " km\n";
  gramListFile  << " Latitude = " << setw(7)  << setprecision(2) << pos.latitude
          << "  degrees       Longitude = " << setw(7)  << setprecision(2) << pos.getLongitude(eastLongitudePositive)
          << " " << eastWestChar() << " (" << setw(7)  << setprecision(2) << pos.getLongitude(!eastLongitudePositive)
          << " " << eastWestChar(true) << ") deg.\n";
  gramListFile  << " Sun Latitude = " << setw(8)  << setprecision(2) << eph.subsolarLatitude
          << " deg.      Titan Orbital Radius  =" << setw(7)  << setprecision(3) << eph.orbitalRadius
          << " AU\n";
  gramListFile  << " Sun Longitude =" << setw(8)  << setprecision(2) << eph.getSubsolarLongitude(eastLongitudePositive)
          << " deg." << eastWestChar() <<"     Local True Solar Time = " << setw(6)  << setprecision(2) << eph.solarTime
          << " Titan hr\n";
  gramListFile  << " Temperature = " << setw(7)  << setprecision(1) << atm.temperature
          << " K        Pressure = " << scientific << setw(10)  << setprecision(3) << atm.pressure
          << " N/m**2   profwgt =" << fixed << setw(6)  << setprecision(3) << atm.profileWeight[0]
          << newline;
  gramListFile  << " Density (Low, Avg., High) = " << scientific << setw(12)  << setprecision(3) << atm.lowDensity
          << setw(12)  << setprecision(3) << atm.density
          << setw(12)  << setprecision(3) << atm.highDensity
          << " " << densUnits << newline;
  gramListFile  << " Departure from Yelle Avg  =" << fixed << setw(10)  << setprecision(1) << atm.lowDensityDeviation * 100.0
          << " %" << setw(10)  << setprecision(1) << atm.densityDeviation * 100.0
          << " %" << setw(10)  << setprecision(1) << atm.highDensityDeviation * 100.0 << " %\n";
  gramListFile  << " Tot.Dens. =" << scientific << setw(10)  << setprecision(3) << atm.perturbedDensity
          << " " << densUnits <<
          "    Dens.Pert. =" << fixed << setw(7)  << setprecision(2) << atm.densityPerturbation * 100.0
          << " % of mean  iupdate= "<< iupdate << newline;
  gramListFile  << " Eastward Wind  (Mean, Perturbed, Total) = " << setw(7)  << setprecision(1) << atm.ewWind
          << setw(7)  << setprecision(1) << atm.ewWindPerturbation
          << setw(7)  << setprecision(1) << atm.perturbedEWWind
          << " m/s\n";
  gramListFile  << " Northward Wind (Mean, Perturbed, Total) = " << setw(7)  << setprecision(1) << atm.nsWind
          << setw(7)  << setprecision(1) << atm.nsWindPerturbation
          << setw(7)  << setprecision(1) << atm.perturbedNSWind
          << " m/s\n";
  gramListFile  << " N2 =" << scientific << setw(10)  << setprecision(3) << atm.dinitrogen.numberDensity
          << "  CH4 =" << setw(10)  << setprecision(3) << atm.methane.numberDensity
          << "  Ar=" << setw(10)  << setprecision(3) << atm.argon.numberDensity
          << "  Total=" << setw(10)  << setprecision(3) << atm.totalNumberDensity
          << " #/m**3\n";
  gramListFile  << fixed << setw(11)  << setprecision(3) << atm.dinitrogen.massFraction * 100.0
          << "%" << setw(15)  << setprecision(3) << atm.methane.massFraction * 100.0
          << "%" << setw(14)  << setprecision(3) << atm.argon.massFraction * 100.0
          << "%  by mass       MolWgt=" << setw(7)  << setprecision(3) << atm.averageMolecularWeight
          << newline;
  gramListFile  << fixed << setw(11)  << setprecision(3) << atm.dinitrogen.moleFraction * 100.0
          << "%" << setw(15)  << setprecision(3) << atm.methane.moleFraction * 100.0
          << "%" << setw(14)  << setprecision(3) << atm.argon.moleFraction * 100.0
          << "%  mole (or volume) fraction\n";
  gramListFile  << " ---------------------------------------------------------------------------\n";
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void TitanProfilePrinter::printGramTPresHgtHeader(const PerturbedAtmosphere& atmos)
{
  gramTPresHgtFile << "   Height       Temp      Pres   TdegC   Pres_mb    Hrho  Hpres MolWt N2%v CH4%v Ar%v LOGSCALE\n";
}

void TitanProfilePrinter::printGramTPresHgtStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%13.3f%7.1f%11.3e%7.1f%11.3e%7.2f%7.2f%6.2f%5.1f%5.1f%5.1f%5d\n",
    data.position.height,
    atm.temperature,
    atm.pressure,
    atm.temperature - 273.15,
    atm.pressure * 0.01,
    atm.densityScaleHeight,
    atm.pressureScaleHeight,
    atm.averageMolecularWeight,
    atm.dinitrogen.moleFraction * 100.0,
    atm.methane.moleFraction * 100.0,
    atm.argon.moleFraction * 100.0,
    densityPrintScale);

  gramTPresHgtFile << buffer;
}

} // namespace
