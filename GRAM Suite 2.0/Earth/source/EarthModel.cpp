//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "EarthModel.h"
#include "EarthInputParameters.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EarthModel::EarthModel()
  : Atmosphere(), EarthCommon(this)
{
  // Set the Earth specific metrics.
  atmos.setPlanetSpecificMetrics(earthAtmos);
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
EarthModel::EarthModel(const EarthModel& orig)
  : Atmosphere(orig), EarthCommon(this), ncep(orig.ncep), merra2(orig.merra2), map(orig.map), msis(orig.msis), met(orig.met), jb2008(orig.jb2008)
{
  // Note that all classes, NCEP, MAP, etc., are in the (copy) constructor list above.
  // Data not requiring a copy constructor (plain old data), is copied below.
  // Not all items below need to be copied, but it is simpler to maintain by copying everything.
  earthAtmos = orig.earthAtmos;
  year = orig.year;
  month = orig.month;
  thermosphereModel = orig.thermosphereModel;
  ewWindPerturbationScale = orig.ewWindPerturbationScale;
  nsWindPerturbationScale = orig.nsWindPerturbationScale;
  lowerLowAtmosFairingHeight = orig.lowerLowAtmosFairingHeight;
  upperLowAtmosFairingHeight = orig.upperLowAtmosFairingHeight;
  ppm = orig.ppm;
  ppmToND = orig.ppmToND;
  epsilon = orig.epsilon;
  waterSD = orig.waterSD;
  ppmWaterLow = orig.ppmWaterLow;
  vaporPressureLow = orig.vaporPressureLow;
  vaporPressureSDLow = orig.vaporPressureSDLow;
  relativeHumidityLow = orig.relativeHumidityLow;
  relativeHumiditySDLow = orig.relativeHumiditySDLow;
  temperatureSDLow = orig.temperatureSDLow;
  windCorrelationLow = orig.windCorrelationLow;
  dtdz = orig.dtdz;
  dmdz = orig.dmdz;
  useMERRA2 = orig.useMERRA2;

  // Set the Earth specific metrics.
  atmos.setPlanetSpecificMetrics(earthAtmos);
}

//! \fn  EarthModel::~EarthModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void EarthModel::setInputParameters(const EarthInputParameters& params)
{
  // Pass parameters on the the model classes.
  useMERRA2 = !params.useNCEP;
  if (useMERRA2) {
    merra2.setInputParameters(params);
  }
  else {
    ncep.setInputParameters(params);
  }
  map.setInputParameters(params);
  msis.setInputParameters(params);
  met.setInputParameters(params);
  jb2008.setInputParameters(params);

  thermosphereModel = params.thermosphereModel;
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
  year = params.year;
  month = params.month;
}

//! \brief Set solar flux and index parameters.
//!
//! Passes dailyF10, meanF10, and ap to thermosphere models.
//! \param params Input parameters.
void EarthModel::setSolarParameters(const EarthInputParameters& params)
{
  msis.setInputParameters(params);
  met.setInputParameters(params);
  jb2008.setInputParameters(params);
}

//! \brief Set JB2008 parameters.
//!
//! Passes dailyS10, meanS10, dailyXM10, meanXM10, dailyY10, meanY10, and dstdtc to the JB2008 thermosphere model.
//! \param params Input parameters.
void EarthModel::setJB2008Parameters(const EarthInputParameters& params)
{
  jb2008.setInputParameters(params);
}

//! \brief Set perturbations scales for PDT and winds.
//!
//! \param params Input parameters.
void EarthModel::setPerturbationScales(const EarthInputParameters& params)
{
  map.setInputParameters(params);
  ncep.setInputParameters(params);
  merra2.setInputParameters(params);
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
}

//! \brief Sets the day of the year.
//!
//! \param doy Day of the year.
//! \param jDay Julian day.
void EarthModel::setDayOfYear(greal doy, greal jDay)
{
  met.setDayOfYear(doy);
  msis.setDayOfYear(doy);
  jb2008.setDayOfYear(doy);
  jb2008.setJulianDay(jDay);
}

//! \fn EarthModel::setThermosphereModel(ThermosphereModelType model)
//! \copydoc EarthAtmosphere::setThermosphereModel()

//! \fn EarthModel::setUseNCEP(bool flag)
//! \brief Control use of NCEP model over MERRA-2 model.
//!
//! \param flag Set to true to use NCEP model instead of MERRA-2.

//! \fn EarthModel::setNCEPParameters(int NCEPYear, int NCEPHour)
//! \copydoc EarthAtmosphere::setNCEPParameters()

//! \fn EarthModel::getEpsilon()
//! \brief Gets the ratio of molecular weights of water to air.
//!
//! \returns Ratio of molecular weights of water to air (epsilon).

//! \fn EarthModel::getPpmToND()
//! \brief Gets the ppm to number density conversion factor.
//!
//! \returns The ppm to number density conversion factor.

//! \copydoc Atmosphere::update()
void EarthModel::update()
{
  windSpeed = 0.0;

  updateSurface(atmos, false);

  // If height is low, call the getHeights member function to get the height at
  // the lowest two pressure levels to use for fairing low and middle atmospheres.
  if (useMERRA2) {
    if (height < 80.0_km) {
      merra2.getHeights(latitude, longitude, lowerLowAtmosFairingHeight, upperLowAtmosFairingHeight);
    }
    else {
      // Above the MERRA2 data set, so set this below the current height
      upperLowAtmosFairingHeight = 0.0_km;
    }
  }
  else {
    if (height < 40.0_km) {
      ncep.getHeights(latitude, longitude, lowerLowAtmosFairingHeight, upperLowAtmosFairingHeight);
    }
    else {
      // Above the NCEP data set, so set this below the current height
      upperLowAtmosFairingHeight = 0.0_km;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // If height <= height of lowest pressure level, then use lower atmosphere data
  /////////////////////////////////////////////////////////////////////////////////////////
  temperatureSDLow = 0.0;
  if (height <= upperLowAtmosFairingHeight) {
    // Compute lower atmosphere data values
    AtmosphereState atmosLow;
    EarthAtmosphereState earthLow;
    atmosLow.setPlanetSpecificMetrics(earthLow);
    if (useMERRA2) {
      updateMERRA2(atmosLow, dtdz);
    }
    else {
      updateNCEP(atmosLow, dtdz);
    }

    // Set these now as they do not fair with MAP
    dewPoint = earthLow.dewPoint;
    dewPointSD = earthLow.dewPointSD;
    windSpeed = earthLow.windSpeed;
    windSpeedStandardDeviation = earthLow.windSpeedStandardDeviation;

    // Save these for computing species concentrations.
    temperatureSDLow = atmosLow.temperatureStandardDeviation;
    relativeHumidityLow = earthLow.relativeHumidity;
    relativeHumiditySDLow = earthLow.relativeHumiditySD;
    vaporPressureLow = earthLow.vaporPressure;
    vaporPressureSDLow = earthLow.vaporPressureSD;

    // Get mean and sigma water vapor volume mixing ratio.
    ppmWaterLow = 1.0e+6 * vaporPressureLow / (atmosLow.pressure - vaporPressureLow);
    waterSD = ppmWaterLow * vaporPressureSDLow / vaporPressureLow;

    // No fairing occurs below the lower fairing height.
    if (height < lowerLowAtmosFairingHeight) {
      atmos.pressure = atmosLow.pressure;
      atmos.density = atmosLow.density;
      atmos.temperature = atmosLow.temperature;
      atmos.ewWind = atmosLow.ewWind;
      atmos.nsWind = atmosLow.nsWind;
      atmos.verticalWind = atmosLow.verticalWind;
    }
    else {
      // Fair between lower and middle-atmosphere values if height 
      // is between the two lower atmosphere fairing heights.

      // Call map member function for calculating values from MAP data.
      AtmosphereState atmosMAP;
      greal dtdzMAP;
      updateMAP(atmosMAP, dtdzMAP);

      // TODO why is PW using 0 here?
      // fair between lower atmosphere and MAP
      fair(1, lowerLowAtmosFairingHeight, atmosLow, upperLowAtmosFairingHeight, atmosMAP, height, atmos);

      // If the MAP and lower atmosphere temperatures differ, 
      // then fair the temperature gradient.
      if (atmosMAP.temperature != atmosLow.temperature) {
        greal ffac = (temperature - atmosLow.temperature)
                     / (atmosMAP.temperature - atmosLow.temperature);
        dtdz += ffac * (dtdzMAP - dtdz);
      }

      // Calculate standard deviation of winds from MAP data.
      // Used to calculate windSpeed and windSpeedStandardDeviation.
      greal ewWindSDMap = 0.0, nsWindSDMap = 0.0;
      map.setPosition(position);
      map.getWindStandardDeviations(ewWindSDMap, nsWindSDMap);

      Interpolator cosInterp;
      cosInterp.makeCosineSquaredFraction(lowerLowAtmosFairingHeight, upperLowAtmosFairingHeight, height);
      greal ewWindSD
          = cosInterp.linear(atmosLow.ewStandardDeviation, ewWindSDMap) * ewWindPerturbationScale;
      greal nsWindSD
          = cosInterp.linear(atmosLow.nsStandardDeviation, nsWindSDMap) * nsWindPerturbationScale;

      // Calculate a faired value for wind speed and wind speed standard deviation.
      greal FS = pow(ewWindSD, 2.0) + pow(nsWindSD, 2.0);
      windSpeed
          = cosInterp.linear(windSpeed, sqrt(pow(ewWind, 2.0) + pow(nsWind, 2.0) + 0.605 * FS));
      windSpeedStandardDeviation = cosInterp.linear(windSpeedStandardDeviation, sqrt(0.395 * FS));
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////
  // Use MET, MSIS, or JB2008 thermosphere model if height above lowerThermoFairingHeight, 90 km.
  /////////////////////////////////////////////////////////////////////////////////////////
  else if (height >= lowerThermoFairingHeight) {
    AtmosphereState atmosTherm;
    // Call thermosphere model based on user-selected input.
    switch (thermosphereModel) {
    case EG_MET:
      updateMET(atmosTherm, dtdz, dmdz);
      break;
    case EG_MSIS:
      updateMSIS(atmosTherm, dtdz, dmdz);
      break;
    case EG_JB2008:
      updateJB2008(atmosTherm, dtdz, dmdz);
      break;
    }

    // Use only thermosphere model data if height above upperThermoFairingHeight, 120 km.
    if (height >= upperThermoFairingHeight) {
      atmos.pressure = atmosTherm.pressure;
      atmos.density = atmosTherm.density;
      atmos.temperature = atmosTherm.temperature;
      atmos.ewWind = atmosTherm.ewWind;
      atmos.nsWind = atmosTherm.nsWind;
      atmos.verticalWind = atmosTherm.verticalWind;
    }
    else {
      // Fair between thermosphere and middle-atmosphere values
      // if height between lowerThermoFairingHeight and upperThermoFairingHeight.
      AtmosphereState atmosMAP;
      greal dtdzMAP;
      updateMAP(atmosMAP, dtdzMAP);

      // Fair by height between atmosMAP and atmosTherm. Result goes to atmos.
      fair(1, lowerThermoFairingHeight, atmosMAP, upperThermoFairingHeight, atmosTherm, height, atmos);

      // Interpolate dtdz using temperature.
      if ((atmosMAP.temperature - atmosTherm.temperature) != 0.0) {
        greal ffac = (temperature - atmosTherm.temperature)
                     / (atmosMAP.temperature - atmosTherm.temperature);
        dtdz = dtdz + ffac * (dtdzMAP - dtdz);
      }
    }
  }
  else {
    // Use only middle-atmosphere values if height between 
    // upperLowerAtmosFairingHeight and lowerThermoFairingHeight.
    updateMAP(atmos, dtdz);
  }

  getStandardDeviations(atmos, false);

  // Calculate species concentration data.
  updateSpeciesConcentrations();

  // Calculate wind speed if data is not available.
  if (windSpeed <= 0) {
    // Calculate wind speed if data is not available.
    greal ewSD, nsSD;
    map.setPosition(position);
    map.getWindStandardDeviations(ewSD, nsSD);
    ewSD *= ewWindPerturbationScale;
    nsSD *= nsWindPerturbationScale;
    greal FS = square(ewSD) + square(nsSD);
    windSpeed = sqrt(square(ewWind) + square(nsWind) + 0.605 * FS);
    windSpeedStandardDeviation = sqrt(0.3950 * FS);
  }

}

//! \brief Compute surface values.
//!
//! The lower atmosphere model is used to compute surface metrics.  Initial values are computed by
//! setting the init argument to \c true.
//!
//! \param[out] inState The AtmosphereState to populate.
//! \param init  Set to \c true if values are being initialized.
//!
//! \retval inState Has been populated with surface values.
void EarthModel::updateSurface(AtmosphereState& inState, bool init)
{
  // Get lower atmosphere data at surface height
  if (useMERRA2) {
    merra2.setPosition(position);
    merra2.updateSurface();
  }
  else {
    ncep.setPosition(position);
    ncep.updateSurface();
  }
  // Get a reference to the lower atmosphere state and the earth specific values.
  const AtmosphereState& surfaceState = useMERRA2 ? merra2.getAtmosphereState() : ncep.getAtmosphereState();
  const EarthAtmosphereState& surfaceEarth = surfaceState.getPlanetSpecificMetrics<EarthAtmosphereState>();

  // Copy the lower atmosphere values into the inState's surface metrics.
  EarthAtmosphereState& earth = inState.getPlanetSpecificMetrics<EarthAtmosphereState>();
  inState.pressureAtSurface = surfaceState.pressure;
  earth.densityAtSurface = surfaceState.density;
  earth.temperatureAtSurface = surfaceState.temperature;
  earth.temperatureSDAtSurface = surfaceState.temperatureStandardDeviation;
  earth.ewWindAtSurface = surfaceState.ewWind;
  earth.nsWindAtSurface = surfaceState.nsWind;
  earth.ewWindSDAtSurface = surfaceState.ewStandardDeviation;
  earth.nsWindSDAtSurface = surfaceState.nsStandardDeviation;
  earth.windCorrelationAtSurface = surfaceEarth.windCorrelation;
  earth.windSpeedAtSurface = surfaceEarth.windSpeed;
  earth.windSpeedSDAtSurface = surfaceEarth.windSpeedStandardDeviation;
  if (init) {
    earth.temperatureSDAtSurface = earth.temperatureAtSurface * earth.temperatureSDAtSurface;
  }
  else {
    earth.temperatureSDAtSurface = earth.temperatureAtSurface * sqrt(earth.temperatureSDAtSurface);
  }
}

//! \brief Compute standard deviations and surface values.
//!
//! The lower atmosphere and MAP models are used to compute standard deviations based on height.
//! The lower atmosphere model is used to compute surface metrics.  Initial values are computed by
//! setting the init argument to \c true.
//!
//! \param[out] inState The AtmosphereState to populate.
//! \param init  Set to \c true if values are being initialized.
//!
//! \retval inState Has been populated with standard deviations and surface values.
void EarthModel::getStandardDeviations(AtmosphereState& inState, bool init)
{
  // Get a reference to the earth specific metrics in the inState.
  EarthAtmosphereState& inEarth = inState.getPlanetSpecificMetrics<EarthAtmosphereState>();

  // Get the fairing heights at the lowest two pressure levels.
  if (useMERRA2) {
    merra2.getHeights(latitude, longitude, lowerLowAtmosFairingHeight, upperLowAtmosFairingHeight);
  }
  else {
    ncep.getHeights(latitude, longitude, lowerLowAtmosFairingHeight, upperLowAtmosFairingHeight);
  }
  // Get a reference to the lower atmosphere state and the earth specific values.
  const AtmosphereState& lowState = useMERRA2 ? merra2.getAtmosphereState() : ncep.getAtmosphereState();
  const EarthAtmosphereState& lowEarth = lowState.getPlanetSpecificMetrics<EarthAtmosphereState>();

  if (init) {
    updateSurface(inState, init);
  }

  // Get lower atmosphere data at current height
  if (init || height <= upperLowAtmosFairingHeight) {
    // If initializing, then the lower atmosphere needs to be computed.
    if (init) {
      if (useMERRA2) {
        merra2.setPosition(position);
        merra2.update();
      }
      else {
        ncep.setPosition(position);
        ncep.update();
      }
    }

    // Copy the lower atmosphere values into the inState's metrics.
    inState.pressureStandardDeviation = lowState.pressureStandardDeviation;
    inState.densityStandardDeviation = lowState.densityStandardDeviation;
    inState.temperatureStandardDeviation = lowState.temperatureStandardDeviation;
    inState.ewStandardDeviation = lowState.ewStandardDeviation;
    inState.nsStandardDeviation = lowState.nsStandardDeviation;
    inEarth.windCorrelation = lowEarth.windCorrelation;
  }
  else {
    inEarth.windCorrelation = 0.0;
  }

  if (height > lowerLowAtmosFairingHeight) {
    // Use MAP model
    map.setPosition(position);
    // Get wind standard deviation.
    map.getWindStandardDeviations(inState.ewStandardDeviation, inState.nsStandardDeviation);
    // Get interpolated standard deviations of p, d, t
    map.getStandardDeviations(inState.pressureStandardDeviation, 
                              inState.densityStandardDeviation,
                              inState.temperatureStandardDeviation);

    // Fair between lower atmosphere and MAP data
    if (height < upperLowAtmosFairingHeight) {
      greal spj, sdj, stj, suj, svj, swj;
      fair(2,
        lowerLowAtmosFairingHeight, lowState.pressureStandardDeviation, lowState.densityStandardDeviation, 
           lowState.temperatureStandardDeviation, lowState.ewStandardDeviation, 
           lowState.nsStandardDeviation, inState.verticalStandardDeviation,
        upperLowAtmosFairingHeight, inState.pressureStandardDeviation, inState.densityStandardDeviation,
           inState.temperatureStandardDeviation, inState.ewStandardDeviation,
           inState.nsStandardDeviation, inState.verticalStandardDeviation, 
        height, spj, sdj, stj, suj, svj, swj);
      inState.pressureStandardDeviation = spj;
      inState.densityStandardDeviation = sdj;
      inState.temperatureStandardDeviation = stj;
      inState.ewStandardDeviation = suj;
      inState.nsStandardDeviation = svj;
    }
  }
   if (init) {
    return;
  }

  inState.pressureStandardDeviation = sqrt(abs(inState.pressureStandardDeviation));
  inState.densityStandardDeviation = sqrt(abs(inState.densityStandardDeviation));
  inState.temperatureStandardDeviation = sqrt(abs(inState.temperatureStandardDeviation));
}

//! \brief Calculate pressure scale height and density scale height.
//!
//! Computes the pressure scale height (km) and density scale height (km) for the input
//! pressure, density, and temperature.
//!
//! \param pres The current pressure (Pa).
//! \param dens The current density (kg/m^3).
//! \param temp the current temperature (K).
//! \param[out] pressureScaleHeight A greal reference.
//! \param[out] densityScaleHeight A greal reference.
//!
//! \retval pressureScaleHeight The pressure scale height (km).
//! \retval densityScaleHeight The density scale height (km).
void EarthModel::getScaleHeights(greal pres, greal dens, greal temp, 
                                 greal& pressureScaleHeight, greal& densityScaleHeight)
{
  pressureScaleHeight = pres / (dens * gravity);
  if (height <= upperLowAtmosFairingHeight) {
    dtdz = dtdz / 1000.0;
  }
  if (height <= lowerThermoFairingHeight) {
    densityScaleHeight = pressureScaleHeight / (1.0 + pressureScaleHeight * dtdz / temp);
  }
  else {
    densityScaleHeight = pressureScaleHeight
                         / (1.0 + pressureScaleHeight * dtdz / temp
                            - pressureScaleHeight * dmdz / averageMolecularWeight);
  }
  pressureScaleHeight *= 1.0e-3;
  densityScaleHeight *= 1.0e-3;
}

//! \brief Uses the MAP model to compute the atmosphere state.
//!
//! The MAP model is updated with the current position and updated.  The provided AtmosphereState is 
//! populated with pressure, density, temperature, and winds.  The temperature gradient is also returned.
//!
//! \param[out] atmosMAP An AtmosphereState reference.
//! \param[out] dtdzMAP A greal reference.
//!
//! \retval atmosMAP The AtmosphereState with metrics populated.
//! \retval dtdzMAP The temperature gradient (K/km).
void EarthModel::updateMAP(AtmosphereState& atmosMAP, greal& dtdzMAP)
{
  // Update the MAP model for the current position to compute mean pressure, density, temperture, 
  // east-west wind, north-south wind, and vertical wind from the middle atmosphere program.
  map.setPosition(position);
  map.update();
  const AtmosphereState& mapState = map.getAtmosphereState();

  dtdzMAP = map.getTemperatureGradient();
  atmosMAP.pressure = mapState.pressure;
  atmosMAP.temperature = mapState.temperature;
  atmosMAP.density = mapState.density;
  atmosMAP.ewWind = mapState.ewWind;
  atmosMAP.nsWind = mapState.nsWind;
  atmosMAP.verticalWind = mapState.verticalWind;
}

//! \brief Uses the MET model to compute the atmosphere state.
//!
//! The MET model is updated with the current position and updated.  The provided AtmosphereState is 
//! populated with pressure, density, temperature, and winds.  The gas densities are populated.
//! The temperature gradient and the the molecular weight gradient are also returned.
//!
//! \param[out] atmosMET An AtmosphereState reference.
//! \param[out] dtdzMET A greal reference.
//! \param[out] dmdzMET A greal reference.
//!
//! \retval atmosMET The AtmosphereState with metrics populated.
//! \retval dtdzMET The temperature gradient (K/km).
//! \retval dmdzMET The molecular weight gradient (1/km).
void EarthModel::updateMET(AtmosphereState& atmosMET, greal& dtdzMET, greal& dmdzMET)
{
  // Update the MET model for the current position to calculate mean pressure,
  // density, temperature, east-west wind, north-south wind, and vertical wind.
  // Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
  met.setEphemerisState(ephem);
  met.setPosition(position);
  met.update();
  const AtmosphereState& metState = met.getAtmosphereState();

  dtdzMET = met.getTemperatureGradient();
  dmdzMET = met.getMolecularWeightGradient();
  atmosMET.pressure = metState.pressure;
  atmosMET.density = metState.density;
  atmosMET.temperature = metState.temperature;
  atmosMET.ewWind = metState.ewWind;
  atmosMET.nsWind = metState.nsWind;
  atmosMET.verticalWind = metState.verticalWind;

  // Put number densities into the earthModel atmos
  argon.numberDensity = metState.argon.numberDensity;
  helium.numberDensity = metState.helium.numberDensity;
  hydrogen.numberDensity = metState.hydrogen.numberDensity;
  dioxygen.numberDensity = metState.dioxygen.numberDensity;
  dinitrogen.numberDensity = metState.dinitrogen.numberDensity;
  oxygen.numberDensity = metState.oxygen.numberDensity;
  averageMolecularWeight = metState.averageMolecularWeight;
}

//! \brief Uses the MSIS model to compute the atmosphere state.
//!
//! The MSIS model is updated with the current position and updated.  The provided AtmosphereState is 
//! populated with pressure, density, temperature, and winds.  The gas densities are populated.
//! The temperature gradient and the the molecular weight gradient are also returned.
//!
//! \param[out] atmosMSIS An AtmosphereState reference.
//! \param[out] dtdzMSIS A greal reference.
//! \param[out] dmdzMSIS A greal reference.
//!
//! \retval atmosMSIS The AtmosphereState with metrics populated.
//! \retval dtdzMSIS The temperature gradient (K/km).
//! \retval dmdzMSIS The molecular weight gradient (1/km).
void EarthModel::updateMSIS(AtmosphereState& atmosMSIS, greal& dtdzMSIS, greal& dmdzMSIS)
{
  // Update the MSIS model for the current position to calculate mean pressure,
  // density, temperature, east-west wind, north-south wind, and vertical wind.  
  // Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
  msis.setEphemerisState(ephem);
  msis.setPosition(position);
  msis.update();
  const AtmosphereState& msisState = msis.getAtmosphereState();

  dtdzMSIS = msis.getTemperatureGradient();
  dmdzMSIS = msis.getMolecularWeightGradient();
  atmosMSIS.pressure = msisState.pressure;
  atmosMSIS.density = msisState.density;
  atmosMSIS.temperature = msisState.temperature;
  atmosMSIS.ewWind = msisState.ewWind;
  atmosMSIS.nsWind = msisState.nsWind;
  atmosMSIS.verticalWind = msisState.verticalWind;

  // Put number densities into the earthModel atmos
  argon.numberDensity = msisState.argon.numberDensity;
  helium.numberDensity = msisState.helium.numberDensity;
  hydrogen.numberDensity = msisState.hydrogen.numberDensity;
  dioxygen.numberDensity = msisState.dioxygen.numberDensity;
  dinitrogen.numberDensity = msisState.dinitrogen.numberDensity;
  oxygen.numberDensity = msisState.oxygen.numberDensity;
  nitrogen.numberDensity = msisState.nitrogen.numberDensity;
  averageMolecularWeight = msisState.averageMolecularWeight;
}

//! \brief Uses the JB2008 model to compute the atmosphere state.
//!
//! The JB2008 model is updated with the current position and updated.  The provided AtmosphereState is 
//! populated with pressure, density, temperature, and winds.  The gas densities are populated.
//! The temperature gradient and the the molecular weight gradient are also returned.
//!
//! \param[out] atmosJB An AtmosphereState reference.
//! \param[out] dtdzJB A greal reference.
//! \param[out] dmdzJB A greal reference.
//!
//! \retval atmosJB The AtmosphereState with metrics populated.
//! \retval dtdzJB The temperature gradient (K/km).
//! \retval dmdzJB The molecular weight gradient (1/km).
void EarthModel::updateJB2008(AtmosphereState& atmosJB, greal& dtdzJB, greal& dmdzJB)
{
  // Update the JB2008 model for the current position to calculate mean pressure,
  // density, temperature, east-west wind, north-south wind, and vertical wind.  
  // Atmospheric constituents ar, he, h, o2, n2, and o are also calculated.
  jb2008.setEphemerisState(ephem);
  jb2008.setPosition(position);
  jb2008.update();
  const AtmosphereState& jbState = jb2008.getAtmosphereState();

  dtdzJB = jb2008.getTemperatureGradient();
  dmdzJB = jb2008.getMolecularWeightGradient();
  atmosJB.pressure = jbState.pressure;
  atmosJB.density = jbState.density;
  atmosJB.temperature = jbState.temperature;
  atmosJB.ewWind = jbState.ewWind;
  atmosJB.nsWind = jbState.nsWind;
  atmosJB.verticalWind = jbState.verticalWind;

  // Put number densities into the earthModel atmos
  argon.numberDensity = jbState.argon.numberDensity;
  helium.numberDensity = jbState.helium.numberDensity;
  hydrogen.numberDensity = jbState.hydrogen.numberDensity;
  dioxygen.numberDensity = jbState.dioxygen.numberDensity;
  dinitrogen.numberDensity = jbState.dinitrogen.numberDensity;
  oxygen.numberDensity = jbState.oxygen.numberDensity;
  nitrogen.numberDensity = jbState.nitrogen.numberDensity;
  averageMolecularWeight = jbState.averageMolecularWeight;
}

//! \brief Uses the NCEP model to compute the atmosphere state.
//!
//! The NCEP model is updated with the current position and updated.  The provided AtmosphereState is 
//! populated with pressure, density, temperature, and winds.  The gas densities are populated.
//! The temperature gradient and the the molecular weight gradient are also returned.
//!
//! \param[out] atmosNCEP An AtmosphereState reference.
//! \param[out] dtdzNCEP A greal reference.
//!
//! \retval atmosNCEP The AtmosphereState with metrics populated.
//! \retval dtdzNCEP The temperature gradient (K/km).
void EarthModel::updateNCEP(AtmosphereState& atmosNCEP, greal &dtdzNCEP)
{
  //The member function ncedmd from the NCEP class used to calculate mean and standard
  //deviation for pressure, density, temperature, east-west wind, north-south wind, 
  //vertical wind, dewpoint, vapor pressure, relative humidity and wind speed.  
  //The u-v correlation is also provided from NCEP data.
  ncep.setPosition(position);
  ncep.update();
  atmosNCEP = ncep.getAtmosphereState();
  dtdzNCEP = ncep.getTemperatureGradient();
}

void EarthModel::updateMERRA2(AtmosphereState& atmosM2, greal& dtdzM2)
{

  merra2.setPosition(position);
  merra2.update();
  atmosM2 = merra2.getAtmosphereState();
  dtdzM2 = merra2.getTemperatureGradient();
}

//! \brief Calculates species concentrations.
//!
//! This method calculates the species concentrations.  Below 120 km, the MAP nad NCEP models
//! are used.  Above 120 km, thermosphere quantities are assumed.
//!
//! \b Inputs
//! \arg #pressure          
//! \arg #temperature          
//! \arg #density          
//! \arg #averageMolecularWeight          
//! \arg #position          
//! \arg #vaporPressureLow
//! \arg #relativeHumidityLow
//! \arg #temperatureStandardDeviation          
//! \arg #gases          
//!
//! returns The mole fraction for all gases and psychrometric quantities.
void EarthModel::updateSpeciesConcentrations()
{
  constexpr greal arppm = 9.34e+03;
  constexpr greal heppm = 5.2;
  greal mb300 = 30000.0_pa;    // NCEP
  if (useMERRA2) {
    mb300 = 30.0_pa;
  }

  greal airMW = averageMolecularWeight;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Heights < 120 then NCEP/LaRC/MAP/AFGL values or mixed values (Thermosphere-NCEP/LaRC/MAP/AFGL)
  //////////////////////////////////////////////////////////////////////////////////////////////////
  if (height < upperThermoFairingHeight) {
    // Evaluate concentrations
    map.setPosition(position);
    map.getSpeciesConcentrations(year, pressure, ppmWaterLow, waterSD, oxygen.numberDensity, ppm);

    // Mean H2O parameters from NCEP if pressure > 300 mb
    if (pressure >= mb300) {
      vaporPressure = vaporPressureLow;
      relativeHumidity = relativeHumidityLow;
      airMW = dryAirMolWt;
      epsilon = water.averageMolecularWeight / airMW;
      vaporDensity = epsilon * vaporPressure / (dryAirGasConstant * temperature);
    }
    // Mean H2O parameters from LaRC/MAP/AFGL concentrations (ppm) if pressure < 300 mb
    else {
      // Molecular weight of air
      if (height <= lowerThermoFairingHeight) {
        airMW = dryAirMolWt;
      }
      else if (height >= upperThermoFairingHeight) {
        airMW = averageMolecularWeight;
      }
      else {
        // Intermediate value if height = 90-120 km
        airMW = dryAirGasConstant * density * temperature * dryAirMolWt / pressure;
      }
      // Mixing ratio, RH, dewpoint, vapor pressure, vapor density
      epsilon = water.averageMolecularWeight / airMW;
      greal mixingRatio = getMixingRatio(temperature, pressure);
      greal waterMole = ppm.water * 1.0e-6;
      relativeHumidity = 0.0;
      if (mixingRatio > 0.0) {
        relativeHumidity = (waterMole / mixingRatio) * (mixingRatio + epsilon) 
                           / (waterMole + epsilon);
      }
      dewPoint = getDewPoint(temperature, relativeHumidity, pressure);
      vaporPressure = waterMole * pressure / (waterMole + epsilon);
      relativeHumidity = relativeHumidity * 100.0;
      vaporDensity = epsilon * vaporPressure / (dryAirGasConstant * temperature);

      if (pressure > 10000.0) {
        vaporPressureSD = dewPointSD * dedt(dewPoint, pressure);
        greal corr = 0.8;
        greal seRatio = vaporPressureSD / vaporPressure;
        greal stRatio = temperatureStandardDeviation * dedt(temperature, pressure)
                        / (100.0 * vaporPressure);
        relativeHumiditySD = relativeHumidity
              * sqrt(square(seRatio) + square(seRatio) - 2 * corr * seRatio * stRatio);
        vaporDensitySD = vaporDensity * sqrt(square(seRatio) + temperatureStandardDeviation);
      }
    }

    // Conversion factor from ppm to number density
    ppmToND = AVOGADRO * density / (1000.0 * airMW);

    // Water vapor concentrations
    water.numberDensity = ppm.water * ppmToND;

    // For standard deviation of water vapor use either:
    // (a) NCEP data (or RH extrapolation), surface to 100 mb
    // (b) MAP vol 31 sigma/mean for 100-0.01 mb
    // (c) A sigma/mean ratio of 0.36 (average of values in Table 1 of
    // Harries, Rev. Geophys. Space Phys., vol 14(4), p 565, 1976) for
    // heights above 0.01 mb level

    if (pressure <= 10000.0) {
      if (pressure < 1.0) {
        // case with sigma/mean = 0.36 case
        vaporPressureSD = 0.36 * vaporPressure;
      }
      else {
        vaporPressureSD = (waterSD / ppm.water) * vaporPressure;
      }

      // MAP sigma case
      if (vaporPressure <= 1.0e-9) {
        vaporPressureSD = 0.0;
        dewPointSD = 0.0;
        relativeHumiditySD = 0.0;
        vaporDensitySD = 0.0;
      }
      else {
        dewPointSD = vaporPressureSD / dedt(dewPoint, pressure);
        greal corr = 0.8;
        greal seRatio = vaporPressureSD / vaporPressure;
        greal stRatio = temperatureSDLow * dedt(temperature, pressure)
                        / (100.0 * vaporPressure);
        relativeHumiditySD = relativeHumidity
              * sqrt(square(seRatio) + square(stRatio) - 2.0 * corr * seRatio * stRatio);
        vaporDensitySD = vaporDensity * sqrt(square(seRatio) + square((dewPointSD / dewPoint)));
      }
    }
    else {
      // NCEP of RH-extrapolated case
      greal seRatio;
      if (pressure >= mb300) {
        vaporPressureSD = vaporPressureSDLow;
        seRatio = vaporPressureSD / vaporPressure;
        relativeHumiditySD = relativeHumiditySDLow;
      }
      else {
        vaporPressureSD = dewPointSD * dedt(dewPoint, pressure);
        seRatio = vaporPressureSD / vaporPressure;
        greal stRatio = temperatureSDLow * dedt(temperature, pressure)
                        / (100.0 * vaporPressure);
        greal corr = 0.8;
        relativeHumiditySD = vaporDensity
              * sqrt(square(seRatio) + square(stRatio) - 2.0 * corr * seRatio * stRatio);
      }
      vaporDensitySD = vaporDensity
            * sqrt(square(seRatio) + square((dewPointSD / dewPoint)));
    }

    // Convert oxygen number density to ppm
    ppm.oxygen = oxygen.numberDensity / ppmToND;

    // Zero if ch<90 km; fixed He and A; AFGL O2 and N2, ch and N=0
    if (height <= lowerThermoFairingHeight) {
      hydrogen.numberDensity = 0.0;
      ppm.hydrogen = 0.0;
      ppm.helium = heppm;
      ppm.argon = arppm;
      helium.numberDensity = heppm * ppmToND;
      argon.numberDensity = arppm * ppmToND;
      nitrogen.numberDensity = 0.0;
      ppm.nitrogen = 0.0;
    }
    else {
      // Thermosphere Ar, He if h = 90-120 km
      ppm.argon = argon.numberDensity / ppmToND;
      ppm.helium = helium.numberDensity / ppmToND;
      ppm.hydrogen = hydrogen.numberDensity / ppmToND;
      ppm.nitrogen = nitrogen.numberDensity / ppmToND;
      greal ppmo2 = dioxygen.numberDensity / ppmToND;
      greal ppmn2 = dinitrogen.numberDensity / ppmToND;
      // Faired values for O2 and N2 if h = 90-120 km
      Interpolator interpHeight;
      interpHeight.makeFraction(lowerThermoFairingHeight, upperThermoFairingHeight, height);
      ppm.dioxygen = interpHeight.linear(ppm.dioxygen, ppmo2);
      ppm.dinitrogen = interpHeight.linear(ppm.dinitrogen, ppmn2);
      dioxygen.numberDensity = ppm.dioxygen * ppmToND;
      dinitrogen.numberDensity = ppm.dinitrogen * ppmToND;
    }

    // AFGL concentrations for other species at all heights < 120
    ozone.numberDensity = ppm.ozone * ppmToND;
    nitrousOxide.numberDensity = ppm.nitrousOxide * ppmToND;
    carbonMonoxide.numberDensity = ppm.carbonMonoxide * ppmToND;
    methane.numberDensity = ppm.methane * ppmToND;
    carbonDioxide.numberDensity = ppm.carbonDioxide * ppmToND;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Heights >= 120 km, concentrations (=0.0) in pure thermosphere range
  //////////////////////////////////////////////////////////////////////////////////////////////////
  else {
    vaporPressure = 0.0;
    vaporDensity = 0.0;
    dewPoint = 0.0;
    relativeHumidity = 0.0;
    relativeHumiditySD = 0.0;
    vaporPressureSD = 0.0;
    vaporDensitySD = 0.0;
    dewPointSD = 0.0;
    ppm.water = 0.0;
    ppm.ozone = 0.0;
    ppm.nitrousOxide = 0.0;
    ppm.carbonMonoxide = 0.0;
    ppm.methane = 0.0;
    ppm.carbonDioxide = 0.0;

    water.numberDensity = 0.0;
    ozone.numberDensity = 0.0;
    nitrousOxide.numberDensity = 0.0;
    carbonMonoxide.numberDensity = 0.0;
    methane.numberDensity = 0.0;
    carbonDioxide.numberDensity = 0.0;

    // Conversion factor from ppm to number density
    ppmToND = AVOGADRO * density / (1000.0 * averageMolecularWeight);
    ppm.dinitrogen = dinitrogen.numberDensity / ppmToND;
    ppm.dioxygen = dioxygen.numberDensity / ppmToND;
    ppm.oxygen = oxygen.numberDensity / ppmToND;
    ppm.argon = argon.numberDensity / ppmToND;
    ppm.helium = helium.numberDensity / ppmToND;
    ppm.hydrogen = hydrogen.numberDensity / ppmToND;
    ppm.nitrogen = nitrogen.numberDensity / ppmToND;
  }
  nitrogen.moleFraction = ppm.nitrogen * 1.0e-6;
  dinitrogen.moleFraction = ppm.dinitrogen * 1.0e-6;
  dioxygen.moleFraction = ppm.dioxygen * 1.0e-6;
  oxygen.moleFraction = ppm.oxygen * 1.0e-6;
  argon.moleFraction = ppm.argon * 1.0e-6;
  helium.moleFraction = ppm.helium * 1.0e-6;
  hydrogen.moleFraction = ppm.hydrogen * 1.0e-6;
  ozone.moleFraction = ppm.ozone * 1.0e-6;
  nitrousOxide.moleFraction = ppm.nitrousOxide * 1.0e-6;
  carbonMonoxide.moleFraction = ppm.carbonMonoxide * 1.0e-6;
  methane.moleFraction = ppm.methane * 1.0e-6;
  carbonDioxide.moleFraction = ppm.carbonDioxide * 1.0e-6;
  water.moleFraction = ppm.water * 1.0e-6;
}

//! \brief Calculates dew point temperature.
//!
//! The method returns dew point temperature (K) as a function of temperature T (K), 
//! relative humidity rh (0-1), and pressure p (Pa), from the inverse of the Buck 4 formula;
//! See Table 2 of Elliott and Gaffen, Bull. Amer. Meterol. Soc., 72(10), 1507, 1991.
//!
//! \param t Temperature (K).
//! \param rh Relative humidity (% between 0 and 1).
//! \param p Pressure (Pa).
//!
//! \returns Dew point temperature (K).
greal EarthModel::getDewPoint(greal t, greal rh, greal p)
{
  constexpr greal e0 = 6.1121;
  constexpr greal a = 1.0 / 227.3;
  constexpr greal b = 18.729;
  constexpr greal c = 257.87;
  constexpr greal aa = 1.0007;
  constexpr greal bb = 3.46e-08;
  constexpr greal ckf = 273.15;

  greal w = wexler(t, p);

  // Assumptions (checking occurs in debug mode only)
  assert((rh > 0.0) && (w > 0.0));

  greal alph = log(rh * w / (100.0 * e0 * (aa + bb * p)));
  greal tdbuck = 2.0 * alph * c / ((b - alph) + sqrt((pow((b - alph), 2) - 4.0 * alph * a * c)));

  // Formula gives dewpoint in C.  Convert to Kelvin.
  return tdbuck + ckf;
}

//! \brief Calculates the water vapor mixing ratio.
//!
//! This method returns the water vapor volume mixing ratio (m**3 vapor/m**3 dry air).
//! 
//! \param tx Dew point temperature (K).
//! \param p Pressure (Pa).
//!
//! \returns Water vapor volume mixing ratio (K).
greal EarthModel::getMixingRatio(greal tx, greal p)
{
  greal mixingRatio = 0.0;

  // Get vapor pressure in (N/m**2)
  greal vaporPres = wexler(tx, p);

  // Avoid cases that would give mixing ratio > 1
  if (vaporPres < p / 2.0) {
    mixingRatio = vaporPres / (p - vaporPres);
  }
  else {
    mixingRatio = 1.0;
  }

  return mixingRatio;
}

} // namespace