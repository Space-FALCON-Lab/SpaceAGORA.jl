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
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "EarthProfilePrinter.h"
#include "EarthInputParameters.h"
#include "EarthAtmosphere.h"
#include "EarthAtmosphereState.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EarthProfilePrinter::EarthProfilePrinter()
{
}

//! \fn EarthProfilePrinter::EarthProfilePrinter(const EarthProfilePrinter& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn EarthProfilePrinter::~EarthProfilePrinter()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc ProfilePrinter::setInputParameters(const InputParameters& params)
void EarthProfilePrinter::setInputParameters(const InputParameters& params)
{
  ProfilePrinter::setInputParameters(params);

  const EarthInputParameters& eParams = static_cast<const EarthInputParameters&>(params);
  if (!eParams.speciesPath.empty()) {
    speciesFileName = eParams.speciesPath;
  }
  if (!eParams.boundaryLayerPath.empty()) {
    bltestFileName = eParams.boundaryLayerPath;
  }
}

//! \copydoc ProfilePrinter::printSectionHeader()
void EarthProfilePrinter::printSectionHeader(const PerturbedAtmosphere& atmos)
{
  if (outputStyle & GRAM_MD_STYLE) {
    mdCount = 1;
  }
  //if (outputStyle & GRAM_LIST_STYLE) {
  //  gramListFile << " --------- ------ -------  ----------  ---------- ------  ------ ------ ----- ----" << '\n';
  //  gramListFile << "   Random seed = " << atmos.getSeed() << '\n';
  //}
}

//! \copydoc ProfilePrinter::printGramCSVHeader()
void EarthProfilePrinter::printPlanetCSVHeader(const PerturbedAtmosphere& atmos)
{
  gramCSVFile << "TemperaturePerturbation_pct,TemperatureStandardDeviation_pct,PerturbedTemperature_K,"
    << "PressurePerturbation_pct,PressureStandardDeviation_pct,PerturbedPressure_Pa,"
    << "PresPertSmall_pct,DensPertSmall_pct,TempPertSmall_pct,EWWindPertSmall_ms,NSWindPertSmall_ms,"
    << "PresSDSmall_pct,DensSDSmall_pct,TempSDSmall_pct,EWWindSDSmall_ms,NSWindSDSmall_ms,"
    << "PresPertLarge_pct,DensPertLarge_pct,TempPertLarge_pct,EWWindPertLarge_ms,NSWindPertLarge_ms,"
    << "PresSDLarge_pct,DensSDLarge_pct,TempSDLarge_pct,EWWindSDLarge_ms,NSWindSDLarge_ms,"
    << "DewPoint_K,DewPointSD_pct,VaporPressure_Pa,VaporPressureSD_pct,VaporDensity_kgm3,VaporDensitySD_pct,"
    << "RelativeHumidity_pct,RelativeHumiditySD_pct,WindSpeed_ms,WindSpeedStandardDeviation_pct,"
    << "WindCorrelation,WindSpeedAtSurface_ms,WindSpeedSDAtSurface_pct,"
    << "TemperatureAtSurface_K,TemperatureSDAtSurface_pct,PressureSDAtSurface_Pa,DensitySDAtSurface_pct,"
    << "DensityAtSurface_kgm3,EWWindAtSurface_ms,NSWindAtSurface_ms,EWWindSDAtSurface_pct,NSWindSDAtSurface_pct,"
    << "WindCorrelationAtSurface_pct,SurfaceHeight_km,GeodeticLatitude_deg,RRAWeight,SurfaceRoughness_m,"
    << "NetRadiationIndex,Stability,InverseLength_m,FrictionVelocity_ms,BVFrequencySquare_s2,MetersAboveSurface_m,"
    << "SigmaRatio,SigmaW_ms,BoundaryLayerDepth_m,NeutralBoundaryLayerDepth_m,PerturbedWindSpeedAtSurface_ms,UnstableBLFactor,"
    << "SolarDays,SolarHourAngle_deg,SolarElevation_deg,ElevationAtMidnight_deg,ElevationAtNoon_deg,LandCode,SeverityLevel";
}

//! \copydoc ProfilePrinter::printGramCSVStyle()
void EarthProfilePrinter::printPlanetCSVStyle(const ProfileData& data)
{
  const AtmosphereState& atm = data.atmos;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();
  const EphemerisState& ephem = data.ephem;

  char buffer[bufferSize];
  buffer[0] = 0;
  int iextra = (int)extra;
  int marker = 0;
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*e,%.*f,%.*f,%.*e,",
    1 + iextra, earth.temperaturePerturbation * 100.0,
    2 + iextra, atm.temperatureStandardDeviation * 100.0,
    3 + iextra, earth.perturbedTemperature,
    1 + iextra, earth.pressurePerturbation * 100.0,
    2 + iextra, atm.pressureStandardDeviation * 100.0,
    3 + iextra, earth.perturbedPressure);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.presPertSmall   * 100.0,
    1 + iextra, earth.densPertSmall   * 100.0,
    1 + iextra, earth.tempPertSmall   * 100.0,
    1 + iextra, earth.ewWindPertSmall,
    1 + iextra, earth.nsWindPertSmall);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.presSDSmall   * 100.0,
    1 + iextra, earth.densSDSmall   * 100.0,
    1 + iextra, earth.tempSDSmall   * 100.0,
    1 + iextra, earth.ewWindSDSmall,
    1 + iextra, earth.nsWindSDSmall);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.presPertLarge   * 100.0,
    1 + iextra, earth.densPertLarge   * 100.0,
    1 + iextra, earth.tempPertLarge   * 100.0,
    1 + iextra, earth.ewWindPertLarge,
    1 + iextra, earth.nsWindPertLarge);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.presSDLarge   * 100.0,
    1 + iextra, earth.densSDLarge   * 100.0,
    1 + iextra, earth.tempSDLarge   * 100.0,
    1 + iextra, earth.ewWindSDLarge,
    1 + iextra, earth.nsWindSDLarge);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*e,%.*e,%.*e,",
    2 + iextra, earth.dewPoint,
    1 + iextra, earth.dewPointSD,
    3 + iextra, earth.vaporPressure,
    3 + iextra, earth.vaporPressureSD,
    3 + iextra, earth.vaporDensity);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.vaporDensitySD,
    1 + iextra, earth.relativeHumidity,
    1 + iextra, earth.relativeHumiditySD,
    1 + iextra, earth.windSpeed,
    1 + iextra, earth.windSpeedStandardDeviation,
    1 + iextra, earth.windCorrelation);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    2 + iextra, earth.windSpeedAtSurface,
    2 + iextra, earth.windSpeedSDAtSurface,
    2 + iextra, earth.temperatureAtSurface,
    2 + iextra, earth.temperatureSDAtSurface * 100.0,
    2 + iextra, earth.pressureSDAtSurface * 100.0,
    2 + iextra, earth.densitySDAtSurface * 100.0);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*e,%.*f,%.*f,%.*f,%.*f,",
    3 + iextra, earth.densityAtSurface,
    2 + iextra, earth.ewWindAtSurface,
    2 + iextra, earth.nsWindAtSurface,
    2 + iextra, earth.ewWindSDAtSurface,
    2 + iextra, earth.nsWindSDAtSurface);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.windCorrelationAtSurface,
    3 + iextra, data.position.surfaceHeight,
    2 + iextra, earth.geodeticLatitude,
    1 + iextra, earth.rraWeight,
    5 + iextra, earth.surfaceRoughness);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    2 + iextra, earth.netRadiationIndex,
    2 + iextra, earth.stability,
    5 + iextra, earth.inverseLength,
    3 + iextra, earth.frictionVelocity,
    5 + iextra, earth.BVFrequencySquare);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,",
    1 + iextra, earth.metersAboveSurface,
    3 + iextra, earth.sigmaRatio,
    2 + iextra, earth.sigmaW,
    1 + iextra, earth.boundaryLayerDepth,
    1 + iextra, earth.neutralBoundaryLayerDepth);
  marker += snprintf(buffer + marker,bufferSize - marker, "%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,%.*f,",
    2 + iextra, earth.perturbedWindSpeedAtSurface,
    3 + iextra, earth.unstableBLFactor,
    3 + iextra, earth.solarDays,
    2 + iextra, ephem.solarHourAngle,
    1 + iextra, earth.solarElevation,
    1 + iextra, earth.elevationAtMidnight,
    1 + iextra, earth.elevationAtNoon);
  gramCSVFile << buffer << earth.landCode << comma << earth.severityLevel;
}

//! \copydoc ProfilePrinter::printGramListMDHeader()
void EarthProfilePrinter::printGramListMDHeader(const PerturbedAtmosphere& atmos)
{
  const EarthAtmosphere& earth = static_cast<const EarthAtmosphere&>(atmos);
  const GramTime& time = atmos.getStartTime();
  const EarthInputParameters& params = static_cast<const EarthInputParameters&>(atmos.getInputParameters());

  gramMDFile << "# " << earth.getVersionString() << '\n';
  gramMDFile << '\n';

  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds, jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15s | %-28s | %-15s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker, bufferSize - marker, "|%.30s|%.17s|%.30s|%.17s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %d/%d/%-*d | %-28s | %-15d |\n",
    "Start Date", month, day, 13 - (month < 10 ? 1 : 2) - (day < 10 ? 1 : 2), year,
    "Initial Random Seed", earth.getSeed());
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %02d:%02d:%05.2f     | %-28s | %-15.2f |\n",
    "Start Time", hour, minute, seconds,
    "Random Perturbation Scale", params.randomPerturbationScale);
   
  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15.6f | %-28s | %-15.2f |\n",
    "Julian Day", jdate,
    "Hor Wind Perturbation Scale", params.horizontalWindPerturbationScale);
   
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  string thermo;
  switch (params.thermosphereModel) {
  case EG_MSIS:
    thermo = "MSIS";
    break;
  case EG_MET:
    thermo = "MET";
    break;
  case EG_JB2008:
    thermo = "JB2008";
    break;
  }

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15s | %-28s | %-15.2f |\n",
    "Thermosphere Model", thermo.c_str(),
    "Ver Wind Perturbation Scale", params.verticalWindPerturbationScale);
   
  if (params.useNCEP) {
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15.2f | %-28s | %-15d |\n",
                  "Daily F10.7", params.dailyF10, "NCEP Global Climatology Data", params.NCEPYear);
  }
  else {
	marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15.2f | %-28s | %-15s |\n",
                "Daily F10.7", params.dailyF10, "MERRA-2 Global Climatology", " ");
  }

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15.2f | %-28s | %-15d |\n",
    "Mean F10.7", params.meanF10,
    params.useNCEP ? "NCEP Hour" : "MERRA-2 Hour", params.useNCEP ? params.NCEPHour : params.M2Hour);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15.2f | %-28s | %-15s |\n",
    "AP Index", params.ap,
    "Patchy Turbulence Option", (params.patchy ? "On" : "Off"));

  string rrayr;
  if (params.useRRA) {
    switch (params.rraYear) {
    case RRA_1983:
      rrayr = "1983";
      break;
    case RRA_2006:
      rrayr = "2006";
      break;
    case RRA_2013:
      rrayr = "2013";
      break;
    case RRA_2019:
      rrayr = "2019";
      break;
    }
  }
  else {
    rrayr = "not used";
  }

  char rraradii[16];
  snprintf(rraradii, 16, "%.1f, %.1f", params.rraInnerRadius, params.rraOuterRadius);

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-15s | %-28s | %-15s |\n\n",
    "Range Reference Atmosphere", rrayr.c_str(),
    "RRA Inner, Outer Radii", rraradii);

  gramMDFile << buffer;

  if (params.useNCEP) {
    gramMDFile << "- NCEP Path: " << params.NCEPPath << '\n';
  }
  else {
    gramMDFile << "- MERRA-2 Path: " << params.M2Path << '\n';
  }

  if (params.useRRA) {
    gramMDFile << "- RRA Path: " << params.rraPath << '\n';
  }
  gramMDFile << '\n';

}

//! \copydoc ProfilePrinter::printGramListMDStyle()
void EarthProfilePrinter::printGramListMDStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();
  const EphemerisState& eph = data.ephem;
  const Position& pos = data.position;

  static const char dashes[51] = "--------------------------------------------------";
  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "--------------------\n");
  marker += snprintf(buffer + marker,bufferSize - marker, "## Record #%zu\n", mdCount++);
  marker += snprintf(buffer + marker,bufferSize - marker, "--------------------\n\n");

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10s | %-33s | %-10s |\n",
    "Field", "Value", "Field", "Value");

  marker += snprintf(buffer + marker,bufferSize - marker, "|%.35s|%.12s|%.35s|%.12s|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.2f |\n",
    "Elapsed Time (s)", pos.elapsedTime,
    "Local Solar Time (hrs)", eph.solarTime);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.1f |\n",
    "Height Above Ref. Ellipsoid (km)",  pos.height,
    "Reference Radius (km)",  pos.latitudeRadius);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Latitude (deg)", pos.latitude,
    "Geodetic Latitude (deg)", earth.geodeticLatitude);
    
  marker += snprintf(buffer + marker,bufferSize - marker, "| %10s%1c%-22s | %-10.2f | %-33s | %-10.2f |\n",
    "Longitude ", eastWestChar(), " (deg)", pos.getLongitude(eastLongitudePositive),
    "Longitude of the Sun (deg)", eph.longitudeSun);
    
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.2f | %19s%1c%-13s | %-10.2f |\n",
    "Subsolar Latitude (deg)", eph.subsolarLatitude,
    "Subsolar Longitude ", eastWestChar(), " (deg)", eph.getSubsolarLongitude(eastLongitudePositive));
    
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Pressure Scale Height (km)", atm.pressureScaleHeight,
    "Orbital Radius (AU)", eph.orbitalRadius);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.2f |\n",
    "Density Scale Height (km)", atm.densityScaleHeight,
    "Solar Zenith Angle (deg)", eph.solarZenithAngle);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Sigma Level", atm.sigmaLevel,
    "Gravity (m/s^2)", pos.gravity);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Pressure Altitude (km)", atm.pressureAltitude,
    "Speed of Sound (m/s)", atm.speedOfSound);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.4f | %-33s | %-10.3f |\n",
    "Compressibility Factor (zeta)", atm.compressibilityFactor,
    "Perturbed Speed of Sound (m/s)", atm.perturbedSpeedOfSound);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10.3f |\n",
    "Specific Heat Ratio", atm.specificHeatRatio,
    "Profile Weight", atm.profileWeight[0]);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3f | %-33s | %-10s |\n",
    "Specific Gas Constant (J/(kg K))", atm.specificGasConstant,
    "RRA Site Name", earth.rraSiteName);
   
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10d | %-33s | %-10.3f |\n",
    "Severity Level", earth.severityLevel,
    "RRA Weight", earth.rraWeight);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.2f |\n",
    "Vapor Pressure (Pa)", earth.vaporPressure,
    "Vapor Pressure SD (%)", earth.vaporPressureSD);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.2f |\n",
    "Vapor Density (kg/m^3)", earth.vaporDensity,
    "Vapor Density SD (%)", earth.vaporDensitySD);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.2f |\n",
    "Dew Point (K)", earth.dewPoint,
    "Dew Point SD (%)", earth.dewPointSD);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.2f | %-33s | %-10.2f |\n",
    "Relative Humidity (%)", earth.relativeHumidity,
    "Relative Humidity SD (%)", earth.relativeHumiditySD);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.2f |\n",
    "Wind Speed (m/s)", earth.windSpeed,
    "Wind Speed SD (%)", earth.windSpeedStandardDeviation);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.2f |\n",
    "Wind Speed at Surface (m/s)", earth.windSpeedAtSurface,
    "Wind Speed SD at Surface (%)", earth.windSpeedSDAtSurface);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-33s | %-10.3e | %-33s | %-10.2f |\n\n",
    "Wind Correlation", earth.windCorrelation,
    "Wind Correlation at Surface", earth.windCorrelationAtSurface);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19s | %-19s | %-19s |\n",
    "Field", "Pressure (Pa)", "Density (kg/m^3)", "Temperature (K)");

  marker += snprintf(buffer + marker,bufferSize - marker, "|%.31s|%.20s:|%.20s:|%.20s:|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.3e | %-19.3e | %-19.2f |\n",
    "Mean", atm.pressure, atm.density, atm.temperature);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.3e | %-19.3e | %-19.2f |\n",
    "Total Perturbed", earth.perturbedPressure, atm.perturbedDensity, earth.perturbedTemperature);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.3e | %-19.3e | %-19.2f |\n",
    "Reference Mean (US-76)", atm.referencePressure, atm.referenceDensity, atm.referenceTemperature);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.3e | %-19.3e | %-19.2f |\n",
    "At the Surface", atm.pressureAtSurface, earth.densityAtSurface, earth.temperatureAtSurface);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Perturbation (%)", 
    earth.pressurePerturbation * 100.0, 
    atm.densityPerturbation * 100.0, 
    earth.temperaturePerturbation * 100.0);
   
  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Deviations (%)", 
    (atm.pressure - atm.referencePressure) / atm.referencePressure * 100.0,
    atm.densityDeviation * 100.0, 
    (atm.temperature - atm.referenceTemperature) / atm.referenceTemperature * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Perturbed Deviations (%)",
    (earth.perturbedPressure - atm.referencePressure) / atm.referencePressure  * 100.0,
    atm.perturbedDensityDeviation * 100.0,
    (earth.perturbedTemperature - atm.referenceTemperature) / atm.referenceTemperature * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Small-Scale Perturbations (%)",
    earth.presPertSmall * 100.0,
    earth.densPertSmall * 100.0,
    earth.tempPertSmall * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Small-Scale Std. Devs. (%)",
    earth.presSDSmall * 100.0,
    earth.densSDSmall * 100.0,
    earth.tempSDSmall * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Large-Scale Perturbations (%)",
    earth.presPertLarge * 100.0,
    earth.densPertLarge * 100.0,
    earth.tempPertLarge * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Large-Scale Std. Devs. (%)",
    earth.presSDLarge * 100.0,
    earth.densSDLarge * 100.0,
    earth.tempSDLarge * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Standard Deviations (%)",
    atm.pressureStandardDeviation * 100.0,
    atm.densityStandardDeviation * 100.0,
    atm.temperatureStandardDeviation * 100.0);
   
  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n\n",
    "Std. Devs. at the Surface (%)",
    earth.pressureSDAtSurface * 100.0,
    earth.densitySDAtSurface * 100.0,
    earth.temperatureSDAtSurface * 100.0);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19s | %-19s | %-19s |\n",
    "Field", "E/W Wind (m/s)", "N/S Wind (m/s)", "Vertical Wind (m/s)");

  marker += snprintf(buffer + marker,bufferSize - marker, "|%.31s|%.20s:|%.20s:|%.20s:|\n",
    dashes, dashes, dashes, dashes);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Mean", atm.ewWind, atm.nsWind, atm.verticalWind);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Total Perturbed", atm.perturbedEWWind, atm.perturbedNSWind, atm.perturbedVerticalWind);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Perturbation", atm.ewWindPerturbation, atm.nsWindPerturbation, atm.verticalWindPerturbation);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n",
    "At the Surface", earth.ewWindAtSurface, earth.nsWindAtSurface);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n",
    "Small-Scale Perturbations", earth.ewWindPertSmall, earth.nsWindPertSmall);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n",
    "Small-Scale Std. Devs.", earth.ewWindSDSmall, earth.nsWindSDSmall);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n",
    "Large-Scale Perturbations", earth.ewWindPertLarge, earth.nsWindPertLarge);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n",
    "Large-Scale Std. Devs.", earth.ewWindSDLarge, earth.nsWindSDLarge);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f | %-19.2f |\n",
    "Standard Deviations", atm.ewStandardDeviation, atm.nsStandardDeviation, atm.verticalStandardDeviation);

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-29s | %-19.2f | %-19.2f |                     |\n\n",
    "Std. Devs. at the Surface", earth.ewWindSDAtSurface, earth.nsWindSDAtSurface);

  gramMDFile << buffer;
  buffer[0] = 0;
  marker = 0;

  marker += snprintf(buffer + marker,bufferSize - marker, "| %-20s | %-22s | %-8s | %-8s | %-11s |\n",
    "Gases", "Number Density (#/m^3)", "Mass (%)", "Mole (%)", "Avg Mol Wgt");

  marker += snprintf(buffer + marker,bufferSize - marker, "|%.22s|%.24s|%.10s|%.10s|%.13s|\n",
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

  marker += snprintf(buffer + marker, bufferSize - marker, "| %-20s | %-22.4e | %-8.1f | %-8.1f | %-11.2f |\n\n",
     "Total", 
     atm.totalNumberDensity,
     massTotal * 100.0,
     moleTotal * 100.0,
     atm.averageMolecularWeight);

  int indexAboveBoundaryLayer = 1;
  if (pos.surfaceHeight > 2.0) {
    indexAboveBoundaryLayer = 2;
  }
  if (pos.height < 5.0 * indexAboveBoundaryLayer) {
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10s | %-38s | %-10s |\n",
      "Field", "Value", "Field", "Value");

    marker += snprintf(buffer + marker, bufferSize - marker, "|%.30s|%.12s|%.40s|%.12s|\n",
      dashes, dashes, dashes, dashes);

    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10d | %-38s | %-10.2f |\n",
      "Land Code", earth.landCode,
      "Solar Days", earth.solarDays);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.5f | %-38s | %-10.2f |\n",
      "Surface Roughness (m)", earth.surfaceRoughness,
      "Solar Hour Angle (deg)", eph.solarHourAngle);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.2f | %-38s | %-10.2f |\n",
      "Net Radiation Index", earth.netRadiationIndex,
      "Solar Elevation (deg)", earth.solarElevation);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.2f | %-38s | %-10.2f |\n",
      "Stability", earth.stability,
      "Elevation At Midnight (deg)", earth.elevationAtMidnight);
     
    gramMDFile << buffer;
    buffer[0] = 0;
    marker = 0;

    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.5f | %-38s | %-10.2f |\n",
      "Inverse Length (1/m)", earth.inverseLength,
      "Elevation At Noon (deg)", earth.elevationAtNoon);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.5f | %-38s | %-10.2f |\n",
      "BV Frequency Square (1/s^2)", earth.BVFrequencySquare,
      "Local Solar Time (hours)", eph.solarTime);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.3f | %-38s | %-10.3f |\n",
      "Friction Velocity (m/s)", earth.frictionVelocity,
      "Perturbed Wind Speed At Surface (m/s)", earth.perturbedWindSpeedAtSurface);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.3f | %-38s | %-10.1f |\n",
      "Sigma Ratio", earth.sigmaRatio,
      "Meters Above Surface (m)", earth.metersAboveSurface);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.3f | %-38s | %-10.1f |\n",
      "Sigma W (m/s)", earth.sigmaW,
      "Boundary Layer Depth (m)", earth.boundaryLayerDepth);
     
    marker += snprintf(buffer + marker, bufferSize - marker, "| %-28s | %-10.3f | %-38s | %-10.1f |\n\n",
      "Unstable BL Factor", earth.unstableBLFactor,
      "Neutral Boundary Layer Depth (m)", earth.neutralBoundaryLayerDepth);
  }

  gramMDFile << buffer;
}

//! \copydoc ProfilePrinter::printGramColumnHeader()
void EarthProfilePrinter::printGramColumnHeader(const PerturbedAtmosphere& atmos)
{
  gramColumnFile << "      Time    Hgtkm GeocenLat  Lon(" << eastWestName()
    << ")   DensMean   PresMean    Tmean   EWmean   NSmean"
    << "   DensPert   PresPert    Tpert   EWpert  NSpert SDden% SDprs% SDtemK SDuwnd SDvwnd SDwwnd Wpert SpdAvg"
    << "  SpdStd SOSmean SOSpert Sev" << '\n';
}


//! \copydoc ProfilePrinter::printGramColumnStyle()
void EarthProfilePrinter::printGramColumnStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  const AtmosphereState& atm = data.atmos;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%10.2f%9.3f%10.5f%11.5f%12.4e%12.4e%8.2f%8.2f%8.2f%12.4e%12.4e%8.2f%8.2f%8.2f"
    "%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%8.2f%8.2f%4d\n",
    pos.elapsedTime, pos.height, pos.latitude, wrapDegrees180(pos.getLongitude(eastLongitudePositive)),
    atm.density, atm.pressure, atm.temperature, atm.ewWind, atm.nsWind,
    atm.perturbedDensity, earth.perturbedPressure, earth.perturbedTemperature,
    atm.perturbedEWWind, atm.perturbedNSWind,
    atm.densityStandardDeviation * 100.0, atm.pressureStandardDeviation * 100.0,
    atm.temperatureStandardDeviation * atm.temperature,
    atm.ewStandardDeviation, atm.nsStandardDeviation, atm.verticalStandardDeviation,
    atm.perturbedVerticalWind, earth.windSpeed, earth.windSpeedStandardDeviation,
    atm.speedOfSound, atm.perturbedSpeedOfSound, earth.severityLevel);

  gramColumnFile << buffer;
}

//! \copydoc ProfilePrinter::printGramListHeader()
void EarthProfilePrinter::printGramListHeader(const PerturbedAtmosphere& atmos)
{
  const EarthInputParameters& params = static_cast<const EarthInputParameters&>(atmos.getInputParameters());
  const GramTime& time = atmos.getStartTime();
  GRAM_TIME_FRAME timeFrame = time.getTimeFrame();
  GRAM_TIME_SCALE timeScale = time.getTimeScale();
  int year, month, day, hour, minute;
  double seconds;
  double jdate;
  time.getStartTime(timeScale, timeFrame, year, month, day, hour, minute, seconds);
  time.getStartTime(timeScale, timeFrame, jdate);
  string rraYear;
  switch (params.rraYear) {
  case RRA_1983:
    rraYear = "1983";
    break;
  case RRA_2006:
    rraYear = "2006";
    break;
  case RRA_2013:
    rraYear = "2013";
    break;
  case RRA_2019:
    rraYear = "2019";
    break;
  }
  string thermo;
  switch (params.thermosphereModel) {
  case EG_MSIS:
    thermo = "MSIS";
    break;
  case EG_MET:
    thermo = "MET";
    break;
  case EG_JB2008:
    thermo = "JB2008";
    break;
  }
  gramListFile << " **** Earth Global Reference Atmospheric Model - " << atmos.getVersionString() << " ****\n"
    << " MM/DD/YYYY = " << month << "/" << day << "/" << year
    << "   HH:MM:SS(UTC) = " << setw(2) << setfill('0') << hour << ":" << setw(2) << minute << ":" << setw(2) << rint(seconds) << setfill(' ')
    << "   Julian Day " << fixed << setprecision(4) << jdate 
    << "\nF10.7 = " << setprecision(1) << params.dailyF10 << " Mean F10.7 = " << params.meanF10 << " ap Index = " << params.ap
    << "\nRange Reference Atmosphere (RRA) data available for " << rraYear
    << "\nNCEP Global Climatology Data:  POR = " << params.NCEPYear << " NCEPhr = " << params.NCEPHour
    << "\nPath = " << params.NCEPPath
    << "\nThermospheric conditions from " << thermo << " model\n"
    << "1st Random No. = " << params.initialRandomSeed 
    << " Random Scale Factors = " << setprecision(2) << params.randomPerturbationScale
    << " " << params.horizontalWindPerturbationScale << " " << params.verticalWindPerturbationScale
    << "\nPatchy Turbulence Option = " << (params.patchy ? "On" : "Off")
    << "\n\n Mean-76 and Total-76 are percent deviations from 1976 US Standard Atmosphere.\n"
    << " Other deviations in percent are with respect to mean values.  RH is relative\n"
    << " humidity in percent.  Zeroes for H2O indicate no estimate available.\n"
    << " E-W wind positive toward East; N-S wind positive toward North.\n\n"
    << "  Height/  GcLat    Long.                        Tempera-               Vert.\n"
    << "  Radius   GdLat   [E+W-]  Pressure   Density/     ture/   E-W    N-S   Wind\n"
    << "   (km)    (deg)    (deg)  /Vap.Pr.   Vap.Dens.    Dewpt.  Wind   Wind  (m/s)\n"
    << "  Time_sec  /RRA   /Weight (Nt/m**2)  (kg/m**3)    (K)    (m/s)  (m/s)  RH(%)\n";
}

//! \copydoc ProfilePrinter::printGramListStyle()
void EarthProfilePrinter::printGramListStyle(const ProfileData& data)
{
  const Position& pos = data.position;
  AtmosphereState atm = data.atmos;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();

  densUnits = "kg/m**3";

  char buffer[bufferSize];
  buffer[0] = 0;
  
  int marker = snprintf(buffer, bufferSize,
    " --------- ------ -------  ----------  ---------- ------  ------ ------ ----- ----\n");

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %8.3f %7.3f %8.3f %10.3e %11.3e %6.1f %6.1f %6.1f %6.3f Mean\n",
    pos.height,
    pos.latitude,
    wrapDegrees180(pos.getLongitude(eastLongitudePositive)),
    atm.pressure,
    atm.density,
    atm.temperature,
    atm.ewWind,
    atm.nsWind,
    atm.verticalWind);

  if (atm.referencePressure*atm.referenceTemperature > 0.0) {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " %8.3f %7.3f          %7.2f%%   %8.2f%%  %6.2f%%                      M-76\n",
      pos.totalRadius,
      earth.geodeticLatitude,
      (atm.pressure - atm.referencePressure) / atm.referencePressure * 100.0,
      atm.densityDeviation * 100.0,
      (atm.temperature - atm.referenceTemperature) / atm.referenceTemperature * 100.0);
  }
  else {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " %8.3f %7.3f          %7.2f%%   %8.2f%%  %6.2f%%                      M-76\n",
      pos.totalRadius,
      earth.geodeticLatitude, 0.0, 0.0, 0.0);
  }

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %9.1f %6s %8.4f %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f        ranS\n",
    pos.elapsedTime,
    earth.rraSiteName,
    earth.rraWeight,
    earth.presPertSmall * 100.0,
    earth.densPertSmall * 100.0,
    earth.tempPertSmall * 100.0,
    earth.ewWindPertSmall,
    earth.nsWindPertSmall);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %25s %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f        sigS\n",
    " ",
    earth.presSDSmall * 100.0,
    earth.densSDSmall * 100.0,
    earth.tempSDSmall * 100.0,
    earth.ewWindSDSmall,
    earth.nsWindSDSmall);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %-25s %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f        ranL\n",
    "Wind and SoS",
    earth.presPertLarge * 100.0,
    earth.densPertLarge * 100.0,
    earth.tempPertLarge * 100.0,
    earth.ewWindPertLarge,
    earth.nsWindPertLarge);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %-25s %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f        sigL\n",
    "------------",
    earth.presSDLarge * 100.0,
    earth.densSDLarge * 100.0,
    earth.tempSDLarge * 100.0,
    earth.ewWindSDLarge,
    earth.nsWindSDLarge);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " Ruv = %6.3f              %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f %6.2f ranT\n",
    earth.windCorrelation,
    earth.pressurePerturbation * 100.0,
    atm.densityPerturbation * 100.0,
    earth.temperaturePerturbation * 100.0,
    atm.ewWindPerturbation,
    atm.nsWindPerturbation,
    atm.verticalWindPerturbation);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " SpdAv=%6.1f              %7.2f%%   %8.2f%%  %6.2f%% %6.1f %6.1f %6.2f sigT\n",
    earth.windSpeed,
    atm.pressureStandardDeviation * 100.0,
    atm.densityStandardDeviation * 100.0,
    atm.temperatureStandardDeviation * 100.0,
    atm.ewStandardDeviation,
    atm.nsStandardDeviation,
    atm.verticalStandardDeviation);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " SpdSd=%6.1f             %11.3e %11.3e %6.1f %6.1f %6.1f %6.2f Tot.\n",
    earth.windSpeedStandardDeviation,
    earth.perturbedPressure,
    atm.perturbedDensity,
    earth.perturbedTemperature,
    atm.perturbedEWWind,
    atm.perturbedNSWind,
    atm.perturbedVerticalWind);

  if (atm.referencePressure*atm.referenceTemperature > 0.0) {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " SoSav=%6.1f              %7.2f%%   %8.2f%%  %6.2f%%                      T-76\n",
      atm.speedOfSound,
      100.0 * (earth.perturbedPressure - atm.referencePressure) / atm.referencePressure,
      atm.perturbedDensityDeviation * 100.0,
      100.0 * (earth.perturbedTemperature - atm.referenceTemperature) / atm.referenceTemperature);
  }
  else {
    marker += snprintf(buffer + marker, bufferSize - marker,
      " SoSav=%6.1f              %7.2f%%   %8.2f%%  %6.2f%%                      T-76\n",
      atm.speedOfSound, 0.0, 0.0, 0.0);
  }

  marker += snprintf(buffer + marker, bufferSize - marker,
    " SoSpt=%6.1f             %11.3e %11.3e %6.1f              %6.1f%% H2O\n",
    atm.perturbedSpeedOfSound,
    earth.vaporPressure,
    earth.vaporDensity,
    earth.dewPoint,
    earth.relativeHumidity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    " %24s %11.3e %11.3e %6.1f              %6.1f%% sigH\n",
    " ",
    earth.vaporPressureSD,
    earth.vaporDensitySD,
    earth.dewPointSD,
    earth.relativeHumiditySD);

  gramListFile << buffer;
}


//! \brief Header method for the EARTH_SPECIES_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void EarthProfilePrinter::printSpeciesHeader(const PerturbedAtmosphere& atmos)
{
  speciesFile << "   Height   GcLat    Long.     Conc.     #Dens.   |    Conc.     #Dens.     Species\n"
                 "    (km)    (deg)    (deg)    (ppmv)    (#/m^3)   |    (ppmv)   (#/m^3)\n";
}

//! \brief Print method for the EARTH_SPECIES_STYLE.
//!
//! \param data The ProfileData to be printed.
void EarthProfilePrinter::printSpeciesStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;

  char buffer[bufferSize];
  buffer[0] = 0;
  int marker = snprintf(buffer, bufferSize,
    " %8.3f %7.3f %8.3f %11.3e %10.3e | %10.3e %10.3e  H2O |  O3\n",
    pos.height,
    pos.latitude,
    wrapDegrees180(pos.getLongitude(eastLongitudePositive)),
    atm.water.moleFraction * 1.0e6,
    atm.water.numberDensity,
    atm.ozone.moleFraction * 1.0e6,
    atm.ozone.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | %10.3e %10.3e  N2O |  CO\n",
    " ",
    atm.nitrousOxide.moleFraction * 1.0e6,
    atm.nitrousOxide.numberDensity,
    atm.carbonMonoxide.moleFraction * 1.0e6,
    atm.carbonMonoxide.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | %10.3e %10.3e  CH4 | CO2\n",
    " ",
    atm.methane.moleFraction * 1.0e6,
    atm.methane.numberDensity,
    atm.carbonDioxide.moleFraction * 1.0e6,
    atm.carbonDioxide.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | %10.3e %10.3e   N2 |  O2\n",
    " ",
    atm.dinitrogen.moleFraction * 1.0e6,
    atm.dinitrogen.numberDensity,
    atm.dioxygen.moleFraction * 1.0e6,
    atm.dioxygen.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | %10.3e %10.3e    O |  Ar\n",
    " ",
    atm.oxygen.moleFraction * 1.0e6,
    atm.oxygen.numberDensity,
    atm.argon.moleFraction * 1.0e6,
    atm.argon.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | %10.3e %10.3e   He |   H\n",
    " ",
    atm.helium.moleFraction * 1.0e6,
    atm.helium.numberDensity,
    atm.hydrogen.moleFraction * 1.0e6,
    atm.hydrogen.numberDensity);

  marker += snprintf(buffer + marker, bufferSize - marker,
    "%26s %11.3e %10.3e | MW=%6.3f  %10.3e    N | Tot\n",
    " ",
    atm.nitrogen.moleFraction * 1.0e6,
    atm.nitrogen.numberDensity,
    atm.averageMolecularWeight,
    atm.totalNumberDensity);
    
  speciesFile << buffer << " -------- ------- --------  ---------- ----------   ---------- ----------  ---- ---- \n";
}

//! \brief Header method for the EARTH_BLTEST_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void EarthProfilePrinter::printBoundaryLayerHeader(const PerturbedAtmosphere& atmos)
{
  bltestFile << "   Day    LST    Hgt    Lat      Lon   hsrf spdsrf"
    " LC     z0   Elmn   El   Elmd     sha blfct  nri   S"
    "       ool  ustar    BVfsq    hN   hbl   chb  sigrat"
    "  swb   swh spdavsrf sdpsdsrf   tsrf  stsrf\n";
}

//! \brief Print method for the EARTH_BLTEST_STYLE.
//!
//! \param data The ProfileData to be printed.
void EarthProfilePrinter::printBoundaryLayerStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();
  const EphemerisState& ephem = data.ephem;

  int indexAboveBoundaryLayer = 1;
  if (pos.surfaceHeight > 2.0)
    indexAboveBoundaryLayer = 2;
  if (pos.height < 5.0 * indexAboveBoundaryLayer) {

    char buffer[bufferSize];
    buffer[0] = 0;

    snprintf(buffer, bufferSize,
      "%7.3f%6.2f%7.3f%8.3f%9.3f%6.3f%7.2f%3d%8.5f%6.1f%6.1f%6.1f%8.2f%6.3f%6.2f%5.2f"
      "%9.5f%6.3f%9.5f%6.0f%6.0f%6.0f%7.3f%6.2f%6.2f%9.2f%9.2f%7.2f%7.2f\n",
      earth.solarDays, ephem.solarTime, pos.height, pos.latitude, wrapDegrees180(pos.getLongitude(eastLongitudePositive)),
      pos.surfaceHeight, earth.perturbedWindSpeedAtSurface, earth.landCode, earth.surfaceRoughness,
      earth.elevationAtMidnight, earth.solarElevation, earth.elevationAtNoon, ephem.solarHourAngle,
      earth.unstableBLFactor, earth.netRadiationIndex, earth.stability, earth.inverseLength,
      earth.frictionVelocity, earth.BVFrequencySquare, earth.neutralBoundaryLayerDepth,
      earth.boundaryLayerDepth, earth.metersAboveSurface, earth.sigmaRatio, earth.sigmaW,
      atm.verticalStandardDeviation, earth.windSpeedAtSurface, earth.windSpeedSDAtSurface,
      earth.temperatureAtSurface, earth.temperatureSDAtSurface);

    bltestFile << buffer;
  }
}

//! \brief Header method for the EARTH_CORR_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void EarthProfilePrinter::printCorrHeader(const PerturbedAtmosphere& atmos)
{
  corrFile << "    Hgtkm GeocenLat  Lon(East)   DensPert   PresPert    Tpert"
    "   EWpert   NSpert   Wpert   DensMean   PresMean    Tmean   EWmean   NSmean  Wmean"
    "  SDden% SDprs% SDtemK SDuwnd SDvwnd SDwwnd\n";
}

//! \brief Print method for the EARTH_CORR_STYLE.
//!
//! \param data The ProfileData to be printed.
void EarthProfilePrinter::printCorrStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%9.3f%10.5f%11.5f%12.4e%12.4e%8.2f%8.2f%8.2f%8.2f%12.4e%12.4e%8.2f%8.2f%8.2f%8.3f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f\n",
    pos.height, pos.latitude, pos.getLongitude(eastLongitudePositive),
    atm.perturbedDensity, earth.perturbedPressure, earth.perturbedTemperature,
    atm.perturbedEWWind, atm.perturbedNSWind, atm.perturbedVerticalWind,
    atm.density, atm.pressure, atm.temperature, atm.ewWind, atm.nsWind, atm.verticalWind,
    atm.densityStandardDeviation * 100.0, atm.pressureStandardDeviation * 100.0,
    atm.temperatureStandardDeviation * atm.temperature,
    atm.ewStandardDeviation, atm.nsStandardDeviation, atm.verticalStandardDeviation);

    corrFile << buffer;
}

//! \brief Header method for the EARTH_TRAJ_STYLE.
//!
//! \param atmos The PerturbedAtmosphere which generated the data.
void EarthProfilePrinter::printTrajHeader(const PerturbedAtmosphere& atmos)
{
  trajFile << "    Hgtkm GeocenLat  Lon(East)  PresPert     DensPert   Tpert    EWpert  NSpert   Wpert\n";
}

//! \brief Print method for the EARTH_TRAJ_STYLE.
//!
//! \param data The ProfileData to be printed.
void EarthProfilePrinter::printTrajStyle(const ProfileData& data)
{
  AtmosphereState atm = data.atmos;
  const Position& pos = data.position;
  const EarthAtmosphereState& earth = atm.getPlanetSpecificMetrics<EarthAtmosphereState>();

  char buffer[bufferSize];
  buffer[0] = 0;

  snprintf(buffer, bufferSize,
    "%9.3f%10.5f%11.5f%12.4e%12.4e%8.2f%8.2f%8.2f%8.2f\n",
    pos.height, pos.latitude, wrapDegrees180(pos.getLongitude(eastLongitudePositive)),
    earth.perturbedPressure, atm.perturbedDensity, earth.perturbedTemperature,
    atm.perturbedEWWind, atm.perturbedNSWind, atm.perturbedVerticalWind);

    trajFile << buffer;
}

//! \copydoc ProfilePrinter::openOutput()
void EarthProfilePrinter::openOutput()
{
  ProfilePrinter::openOutput();

  if ((outputStyle & EARTH_SPECIES_STYLE) && !speciesFileName.empty()) {
    appendExtension(speciesFileName, ".txt");
    openFile(speciesFile, fileNamePrefix + speciesFileName);
  }
  if ((outputStyle & EARTH_BLTEST_STYLE) && !bltestFileName.empty()) {
    appendExtension(bltestFileName, ".txt");
    openFile(bltestFile, fileNamePrefix + bltestFileName);
  }
  if ((outputStyle & EARTH_CORR_STYLE) && !corrFileName.empty()) {
    appendExtension(corrFileName, ".txt");
    openFile(corrFile, fileNamePrefix + corrFileName);
  }
  if ((outputStyle & EARTH_CORR1_STYLE) && !corrFileName.empty()) {
    //appendExtension(corrFileName, "1.txt");
    openFile(corrFile, fileNamePrefix + "Profiles1.txt");
  }
  if ((outputStyle & EARTH_TRAJ_STYLE) && !trajFileName.empty()) {
    appendExtension(trajFileName, ".txt");
    openFile(trajFile, fileNamePrefix + trajFileName);
  }
}

//! \copydoc ProfilePrinter::printFileHeader()
void EarthProfilePrinter::printFileHeader(const PerturbedAtmosphere& atmos)
{
  ProfilePrinter::printFileHeader(atmos);
  if (outputStyle & EARTH_SPECIES_STYLE) {
    printSpeciesHeader(atmos);
  }
  if (outputStyle & EARTH_BLTEST_STYLE) {
    printBoundaryLayerHeader(atmos);
  }
  if (outputStyle & (EARTH_CORR_STYLE | EARTH_CORR1_STYLE)) {
    printCorrHeader(atmos);
  }
  if (outputStyle & EARTH_TRAJ_STYLE) {
    printTrajHeader(atmos);
  }
}

//! \copydoc ProfilePrinter::printData()
void EarthProfilePrinter::printData(const std::vector<ProfileData>& profile)
{
  ProfilePrinter::printData(profile);
  if (outputStyle & EARTH_SPECIES_STYLE) {
    for (auto& p : profile) {
      printSpeciesStyle(p);
    }
  }
  if (outputStyle & EARTH_BLTEST_STYLE) {
    for (auto& p : profile) {
      printBoundaryLayerStyle(p);
    }
  }
  if (outputStyle & (EARTH_CORR_STYLE | EARTH_CORR1_STYLE)) {
    for (auto& p : profile) {
      printCorrStyle(p);
    }
  }
  if (outputStyle & EARTH_TRAJ_STYLE) {
    for (auto& p : profile) {
      printTrajStyle(p);
    }
  }
}

//! \copydoc ProfilePrinter::closeOutput()
void EarthProfilePrinter::closeOutput()
{
  ProfilePrinter::closeOutput();
  if (outputStyle & EARTH_SPECIES_STYLE) {
    speciesFile.close();
  }
  if (outputStyle & EARTH_BLTEST_STYLE) {
    bltestFile.close();
  }
  if (outputStyle & (EARTH_CORR_STYLE | EARTH_CORR1_STYLE)) {
    corrFile.close();
  }
  if (outputStyle & EARTH_TRAJ_STYLE) {
    trajFile.close();
  }
}

} // namespace
