//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "MarsAtmosphere.h"
#include "MOLATopography.h"
#include "SlopeWindsModel.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool MarsAtmosphere::dataLoaded = false;
std::vector<MarsAtmosphere::WaveData> MarsAtmosphere::waveDataList;


//! \copydoc Atmosphere::Atmosphere()
MarsAtmosphere::MarsAtmosphere()
  : MarsCommon(this)
{
  atmos.setPlanetSpecificMetrics(marsAtmos);
  hasVerticalWinds = true;
  gramBody = MARS;
  marsgcmPtr = new MarsGCM();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
MarsAtmosphere::MarsAtmosphere(const MarsAtmosphere& orig)
  : PerturbedAtmosphere(orig), MarsCommon(this)
{
  hasVerticalWinds = true;
  gramBody = MARS;
  inputParameters = orig.inputParameters;
  marsgcmPtr = new MarsGCM(*orig.marsgcmPtr);
  refPosition = orig.refPosition;
  isPlanetoCentric = orig.isPlanetoCentric;
  isMolaHeights = orig.isMolaHeights;
  computeMinMax = orig.computeMinMax;
  waveScale = orig.waveScale;
  dustNu = orig.dustNu;
  dustDensity = orig.dustDensity;
  dustDiameter = orig.dustDiameter;
  perturbationWaveLengthScale = orig.perturbationWaveLengthScale;
  waveData = orig.waveData;
  marsAtmos = orig.marsAtmos;
  atmos.setPlanetSpecificMetrics(marsAtmos);
}

//! \copydoc Atmosphere::~Atmosphere()
MarsAtmosphere::~MarsAtmosphere()
{
  if (marsgcmPtr != NULL) {
    delete marsgcmPtr;
  }
}

//! \copydoc Atmosphere::getVersionString()
const std::string& MarsAtmosphere::getVersionString()
{
  static const string version = "MarsGRAM 2021b :: " + Atmosphere::getVersionString();
  return version;
}

//! \fn void MarsAtmosphere::cp(greal t)
//! \brief Computes heat capacity given temperature.
//!
//! This is a polynomial fit of the heat capacity curve given temperature.
//!
//! \param t Temperature \units{K}.
//! \returns The heat capacity \units{J/K}.

//! \copydoc MarsGCM::setMapYear(int year)
void MarsAtmosphere::setMapYear(int year)
{
  inputParameters.mapYear = clampSize(year, 3);
  marsgcmPtr->setMapYear(inputParameters.mapYear);
}

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.  It also causes the (re)initializaion
//! of the random number generator by setting the initial random seed.
//! \param params The input parameters.
void MarsAtmosphere::setInputParameters(const MarsInputParameters& params)
{
  // Save the input parameters
  inputParameters = params;

  inputParameters.heightAboveSurface = clamp(inputParameters.heightAboveSurface, 0.0_km, 4.5_km);
  inputParameters.boundaryLayerWindsScale = max(0.0, inputParameters.boundaryLayerWindsScale);
  inputParameters.mgcmMaxDustLevel = clamp(inputParameters.mgcmMaxDustLevel, 0.1, 1.0);
  inputParameters.mgcmMinDustLevel = clamp(inputParameters.mgcmMinDustLevel, 0.1, inputParameters.mgcmMaxDustLevel);
  inputParameters.stormDuration = clamp(inputParameters.stormDuration, 12.0, 48.0);
  inputParameters.stormIntensity = clamp(inputParameters.stormIntensity, 0.0, 3.0);
  if (inputParameters.stormMaxRadius < 0.0 || inputParameters.stormMaxRadius > 10000.0) {
    inputParameters.stormMaxRadius = 0.0;
  }
  inputParameters.waveScale = clamp(inputParameters.waveScale, 10.0, 10000.0);
  if (inputParameters.waveDate <= 0.0) {
    inputParameters.waveDate = 0.0;
    inputParameters.wavePhase1Rate = 0.0;
    inputParameters.wavePhase2Rate = 0.0;
    inputParameters.wavePhase3Rate = 0.0;
  }

  if (inputParameters.equatorialRadius != 0.0) {
    equatorialRadius = inputParameters.equatorialRadius;
  }
  else {
    inputParameters.equatorialRadius = equatorialRadius;
  }

  if (inputParameters.polarRadius != 0.0) {
    polarRadius = inputParameters.polarRadius;
  }
  else {
    inputParameters.polarRadius = polarRadius;
  }

  // Set the common perturbation parameters
  PerturbedAtmosphere::setInputParameters(inputParameters);
  setSeed(inputParameters.initialRandomSeed);

  // Add auxiliary atmospheres (if present)
  AuxiliaryAdapter::setInputParameters(inputParameters);

  // Set the data path before creating/initializing any Mars object
  setDataPath(inputParameters.dataPath);

  // Pass parameters to models.
  marsgcmPtr->setInputParameters(inputParameters);

  // Load wave data (if necessary) or set wave parameters
  initializeWaveData(inputParameters);

  // Set parameters.
  waveScale = inputParameters.waveScale;
  isPlanetoCentric = inputParameters.isPlanetoCentric;
  isMolaHeights = inputParameters.isMolaHeights;
  computeMinMax = inputParameters.computeMinMax;
  dustNu = inputParameters.dustNu;
  dustDensity = inputParameters.dustDensity;
  dustDiameter = inputParameters.dustDiameter;
  setPerturbationWaveLengthScale(inputParameters.perturbationWaveLengthScale);
}

//! \brief Set the path to Mars data folder.
//!
//! Set the path to the folder that contains the Mars binary data files for
//! MGCM models, TES models, and topographic heights.
//! \param path A valid Windows or Linux file path.
void MarsAtmosphere::setDataPath(const std::string& path)
{
  char slash = '/';
  dataPath = path;
  replace(dataPath.begin(), dataPath.end(), '\\', slash);
  if (dataPath[dataPath.size() - 1] != slash) {
    dataPath += slash;
  }
}

//! \brief Set the radii for the reference ellipsoid.
//!
//! Set the equatorial and polar radii for the reference ellipsiod.  The values will override
//! the internal constants, #equatorialRadiusDefault and #polarRadiusDefault, set in MarsCommon.h.
//! \param eqrad The equatorial radius \units{km}.
//! \param prad The polar radius \units{km}.
void MarsAtmosphere::setPlanetaryRadii(greal eqrad, greal prad)
{
  if (eqrad > 0.0 && prad > 0.0) {
    equatorialRadius = inputParameters.equatorialRadius = eqrad;
    polarRadius = inputParameters.polarRadius = prad;
  }
}

//! \brief Set the height offset model.
//!
//! Height offsets can be used for two purposes: (1) to control the smoothness of 
//! the transition at 80 km altitude between MGCM data and MTGCM data, or (2) as
//! another means (besides wave multipliers) to adjust MTGCM data for better
//! agreement with observations, such as obtained during aerobraking operations.
//! \param model The #MarsOffsetModel option.
//! \param hgtOffset The (optional) height offset required by the #MARS_CONSTANT and #MARS_SEASONAL options.
void MarsAtmosphere::setHeightOffsetModel(MarsOffsetModel model, greal hgtOffset)
{
  inputParameters.offsetModel = model;
  inputParameters.constantHeightOffset = hgtOffset;
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the height above surface offset.
//!
//! This height offset is added to the height whenever height becomes too negative (less than -8.7 km).
//! The offset must be between 0 and 4.5 km.
//! \param heightAboveSurface The height above surface \units{ km }.
void MarsAtmosphere::setHeightAboveSurface(greal heightAboveSurface)
{
  inputParameters.heightAboveSurface = clamp(heightAboveSurface, 0.0_km, 4.5_km);
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the MGCM dust levels.
//!
//! For the MGCM model, dust optical depths may be controlled using this method.  For a constant 
//! dust optical depth. set \p constantDustLevel between 0.1 and 3.0.  For a seasonal effect, set the 
//! \p constantDustLevel to 0 and set \p minDustLevel (at Ls = 90) and \p maxDustLevel (at Ls = 270).
//! \param constantDustLevel The optical depth of background dust level, 0.1 to 3.0, or use 0 for
//!        assumed seasonal variation of background dust.
//! \param minDustLevel The minimum seasonal dust tau (>= 0.1) if input constantDustLevel = 0.
//! \param maxDustLevel The maximum seasonal dust tau (<= 1.0) if input constantDustLevel = 0.
void MarsAtmosphere::setMGCMDustLevels(greal constantDustLevel, greal minDustLevel, greal maxDustLevel)
{
  inputParameters.mgcmConstantDustLevel = constantDustLevel;
  inputParameters.mgcmMaxDustLevel = clamp(maxDustLevel, 0.1, 1.0);
  inputParameters.mgcmMinDustLevel = clamp(minDustLevel, 0.1, inputParameters.mgcmMaxDustLevel);
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the perturbation wave length scale.
//!
//! This scaling factor is applied to the vertical and horizontal perturbation scales.
//! \param scale The scale factor for perturbation wave lengths (0.1-10).
void MarsAtmosphere::setPerturbationWaveLengthScale(greal scale)
{
  scale = clamp(scale, 0.1, 10.0);
  perturbationWaveLengthScale = inputParameters.perturbationWaveLengthScale = scale;
}

//! \brief Set the 10.7 cm solar flux.
//!
//! Sets the 10.7 cm solar flux (10**-22 W/cm**2 at 1 AU).
//! \param f107 The 10.7 cm solar flux (10**-22 W/cm**2 at 1 AU).
void MarsAtmosphere::setF107(greal f107)
{
  inputParameters.F107 = f107;
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the use of MOLA topographic heights.
//!
//! Set to true if input heights are relative to the MOLA areoid. 
//! Otherwise, input heights are relative to reference ellipsoid.
//! \param isMola True if input heights are relative to the MOLA areoid.
void MarsAtmosphere::setMOLAHeights(bool isMola)
{
  isMolaHeights = inputParameters.isMolaHeights = isMola;
}

//! \brief Controls the computation of daily min/max values.
//!
//! Set to true if daily min/max values are to be computed.
//! \param minMax True if daily min/max values are to be computed.
void MarsAtmosphere::setMinMax(bool minMax)
{
  computeMinMax = inputParameters.computeMinMax = minMax;
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the parameters of a dust storm.
//!
//! A dust storm can be simulated by specifying the start time and duration via longitude of the sun.
//! One must also provide the center of the storm as a lat/lon pair, the radius of the storm, and
//! the storm intensity.
//! \param lonSun The starting Ls value \units{\text{degrees}} for a dust storm (0 = none).
//! \param duration The duration \units{\text{in LS degrees}} for the dust storm.
//! \param intensity The dust storm intensity (0.0 - 3.0).
//! \param maxRadius The max. radius \units{ km } of dust storm (0 or >10000 = global).
//! \param lat The latitude \units{\text{degrees}} of the center of the dust storm.
//! \param lon The longitude \units{\text{degrees, East positive}} of the center of the dust storm.
void MarsAtmosphere::setDustStorm(greal lonSun, greal duration, greal intensity, greal maxRadius, greal lat, greal lon)
{
  inputParameters.stormLongitudeSun = lonSun;
  inputParameters.stormLatitude = lat;
  inputParameters.stormLongitude = lon;
  inputParameters.stormDuration = clamp(duration, 12.0, 48.0);
  inputParameters.stormIntensity = clamp(intensity, 0.0, 3.0);
  inputParameters.stormMaxRadius = maxRadius;
  if (inputParameters.stormMaxRadius < 0.0 || inputParameters.stormMaxRadius > 10000.0) {
    inputParameters.stormMaxRadius = 0.0;
  }
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the dust density parameters.
//!
//! These dust density parameters are used to compute dust density variables by methods  
//! of Haberle et al., J. Geophys. Res., 104, 8957 (1999).  Outputs affected
//! include dustMixingRatio, dustMassDensity, and dustNumberDensity.
//! \param nu The parameter for vertical distribution of dust density.
//! \param diameter The dust particle diameter (micrometers, assumed monodisperse).
//! \param dens The dust particle density \units{kg/m**3}.
void MarsAtmosphere::setDustDensity(greal nu, greal diameter, greal dens)
{
  dustNu = inputParameters.dustNu = nu;
  dustDiameter = inputParameters.dustDiameter = diameter;
  dustDensity = inputParameters.dustDensity = dens;
}

//! \brief Set the exospheric temperature parameters for the Stewart model.
//!
//! The exospheric temperature offset and the factor will be applied to the initial (uncorrected) exospheric temperature.
//! The formula used is \f$ \mathrm{exosphericTemperature} = \mathrm{initialExoTemperature}\: e^{(0.16\ \mathrm{exosphericTemperatureFactor})} + \mathrm{exosphericTemperatureOffset}\f$.
//! \param offset The adjustment for exospheric temperature \units{K}.
//! \param factor The standard deviation for thermosphere variation (-3.0 to +3.0).
void MarsAtmosphere::setExosphericTemperature(greal offset, greal factor)
{
  inputParameters.exosphericTemperatureOffset = offset;
  inputParameters.exosphericTemperatureFactor = factor;
  marsgcmPtr->setInputParameters(inputParameters);
}

//! \brief Set the default values for the wave perturbation parameters.
//!
//! The wave perturbation factor typically depends on parameters supplied in a wave data file.
//! If no wave data file is supplied, then the default parameters are used.
//! \param date The Julian date for (primary) peak(s) of wave (0 for no traveling component).
//! \param scale The vertical scale \units{km} of longitude-dependent wave damping
//!              at altitudes below 100 km(10 <= scale <= 10, 000 km).
//! \param mean The mean term of longitude-dependent wave multiplier for density.
//! \param a1 The amplitude of the wave-1 component of longitude-dependent wave multiplier for density.
//! \param p1 The phase of the wave-1 component of longitude-dependent wave multiplier 
//!           (East positive longitude in degrees).
//! \param r1 The rate of longitude movement (East positive degrees per day) for the wave-1 component.
//! \param a2 The amplitude of the wave-2 component.
//! \param p2 The phase of the wave-2 component.
//! \param r2 The rate of longitude movement for the wave-2 component.
//! \param a3 The amplitude of the wave-3 component.
//! \param p3 The phase of the wave-3 component.
//! \param r3 The rate of longitude movement for the wave-3 component.
void MarsAtmosphere::setWaveDefaults(greal date, greal scale, greal mean, greal a1, greal p1, greal r1, greal a2, greal p2, greal r2, greal a3, greal p3, greal r3)
{
  inputParameters.waveDate = date;
  inputParameters.waveScale = clamp(scale, 10.0, 10000.0);
  inputParameters.waveMeanOffset = mean;
  inputParameters.waveAmplitude1 = a1;
  inputParameters.wavePhase1 = p1;
  inputParameters.wavePhase1Rate = r1;
  inputParameters.waveAmplitude2 = a2;
  inputParameters.wavePhase2 = p2;
  inputParameters.wavePhase2Rate = r2;
  inputParameters.waveAmplitude3 = a3;
  inputParameters.wavePhase3 = p3;
  inputParameters.wavePhase3Rate = r3;
  if (inputParameters.waveDate <= 0.0) {
    inputParameters.waveDate = 0.0;
    inputParameters.wavePhase1Rate = 0.0;
    inputParameters.wavePhase2Rate = 0.0;
    inputParameters.wavePhase3Rate = 0.0;
  }
  initializeWaveData(inputParameters);
}

//! \brief Set the wave perturbation parameter file.
//!
//! Specifies the path to the optional file containing time-dependent wave coefficient data.
//! See setWaveDefaults() for a list of parameters preceded by an elasped time \units{\text{seconds}}.
//! \param waveFile The full or relative path to the wave file.
void MarsAtmosphere::setWaveFile(const std::string& waveFile)
{
  inputParameters.waveFile = waveFile;
  initializeWaveData(inputParameters);
}

//! \brief Set the wind scale factors.
//!
//! These scales are applied to the nominal mean and boundary layer winds.
//! \param meanWinds The scale factor for mean winds.
//! \param boundaryLayerWinds The scale factor for boundary layer slope winds (0 = none).
void MarsAtmosphere::setWindScales(greal meanWinds, greal boundaryLayerWinds)
{
  inputParameters.meanWindsScale = meanWinds;
  inputParameters.boundaryLayerWindsScale = max(0.0, boundaryLayerWinds);
}

//! \brief Interface for the primary atmosphere computations.
//!
//! This routine controls the computation of the atmospheric state for the current position.
//! The ephemeris state and the atmosphere state are updated. 
//! The state is updated by the auxiliary atmospheres, if present.
//! Wave perturbations factors are applied.
//! Then perturbations are computed prior to computing a few final metrics.
//!
//! \b Inputs
//! \arg #position      
//!
//! \returns The AtmosphereState is populated.
void MarsAtmosphere::update()
{
  MarsGCM& marsgcm = *marsgcmPtr;

  // Update the time object with the current elapsed time.
  time.setElapsedTime(elapsedTime);

  // Update ephemeris values.
  updateEphemeris();

  // Update the wave perturbation model.
  updateWavePertubation();

  // Set the current position and ephemeris in the Mars model.
  marsgcm.setWavePerturbation(wavePerturbation);
  marsgcm.setEphemerisState(ephem);
  marsgcm.setPosition(position);
  // Compute the atmosphere state.
  marsgcm.update();
  greal wp = wavePerturbation;
  atmos = marsgcm.getAtmosphereState();
  wavePerturbation = wp;

  // Update the atmosphere state with any auxilary atmosphere data.
  updateAuxiliaryAtmospheres(position, atmos);

  // Save the mean density (without wave perturbations) for the perturbation model. 
  greal meanDensity = density;
  greal meanEWWind = ewWind;
  greal meanNSWind = nsWind;

  // Should this be in the MGCM model?
  updateSpecificHeatRatio();
  updateSpeedOfSound();
  updateSlopeWinds();

  // Apply wave perturbation to pressure and density values 
  applyWavePerturbation();

  // Apply perturbations.
  updatePerturbations(meanDensity, meanEWWind, meanNSWind);

  // Limit perturbed winds to approximately the speed of sound over sqrt(2).
  greal sosLimit = 0.7 * speedOfSound;
  if (abs(perturbedEWWind) > sosLimit) {
    perturbedEWWind = copysign(sosLimit, perturbedEWWind);
    ewWindPerturbation = perturbedEWWind - ewWind;
  }
  if (abs(perturbedNSWind) > sosLimit) {
    perturbedNSWind = copysign(sosLimit, perturbedNSWind);
    nsWindPerturbation = perturbedNSWind - nsWind;
  }

  // If height is too negative, use user input height.
  if (height < -8.7_km) {
    height = surfaceHeight + inputParameters.heightAboveSurface;
    totalRadius = latitudeRadius + height;
  }

  // Update metrics that depend on atmosphere state.
  updateMetrics();
  
  // Compute areodetic position values (for output only)
  planetodeticPosition = position;
  planetodeticPosition.convertToPlanetodetic(polarRadius, equatorialRadius);
  marsAtmos.planetoGraphicHeight = planetodeticPosition.height;
  marsAtmos.planetoGraphicLatitude = planetodeticPosition.latitude;
  marsAtmos.referenceRadius = refPosition.latitudeRadius;
  marsAtmos.referenceHeight = refPosition.height;
}

//! \brief Updates the position to be relative to the MOLA areoid.
//!
//! This routine updates the current position to be relative to the MOLA areoid. The routine also
//! updates #refPosition to be the current position relative to the reference ellipsoid.
//! The flag #isMolaHeights is used to declare the state of the incoming #position data.  When true,
//! the input position is relative to the MOLA areoid.  When false, the input position is relative
//! to the reference ellipsoid.
//!
//! \b Inputs
//! \arg #isMolaHeights      
//! \arg #position      
//!
//! \retval #position
//! \retval #refPosition
//! \retval #albedo
void MarsAtmosphere::updatePosition()
{
  if (isMolaHeights && !isPlanetoCentric) {
    throw string("Error: Positions relative to the MOLA areoid must be planetocentric.\n"
                 "       Please verify the IsMolaHeights and IsPlanetoCentric parameters.");
  }

  // Position will be wrt the mola areoid.
  // RefPosition will be wrt the reference ellipsoid.
  refPosition = position;

  // Get radius at latitude to ref. ellipsoid.  (oldref)
  greal er2 = equatorialRadius * equatorialRadius;
  greal pr2 = polarRadius * polarRadius;
  greal tlat2 = pow(tan(toRadians(latitude)), 2);
  greal xx = er2 * pr2 / (pr2 + er2 * tlat2);
  greal yy = xx * tlat2;
  refPosition.latitudeRadius = sqrt(xx + yy);

  // Get radius at latitude to mola areoid.  (rref)
  MOLATopography mola;
  mola.setAreoidRadiusCallback(getAreoidRadiusCallback);
  mola.setTopographicHeightCallback(getTopographicHeightCallback);
  mola.setCallbackData(callbackDataPointer);
  latitudeRadius = mola.getAreoidRadius(latitude, longitude);

  // Large heights are assumed to include the radius.
  //if (height > 3000.0_km) {
  //  // Fix totalRadius (same for mola and ref)
  //  refPosition.totalRadius = totalRadius = height;
  //  // Fix heights
  //  height = totalRadius - latitudeRadius;
  //  refPosition.height = refPosition.totalRadius - refPosition.latitudeRadius;
  //}
  //else { // normal heights
  if (isMolaHeights) {
    // Get total radius from mola coordinates.
    refPosition.totalRadius = totalRadius = height + latitudeRadius;
    // Find height above ref. ellipsoid.
    refPosition.height = refPosition.totalRadius - refPosition.latitudeRadius;
  }
  else {
    // Get total radius from ref. coordinates.
    totalRadius =  refPosition.totalRadius = refPosition.latitudeRadius + height;
    // Find height above mola areoid.
    height = totalRadius - latitudeRadius;
  }
  //}

  // Get MOLA quantities.
  surfaceHeight = mola.getTopographicHeight(latitude, longitude);

  // Update gravity based on updated positions.
  gravity = getGravity(latitude, latitudeRadius, height);
  refPosition.gravity = getGravity(latitude, refPosition.latitudeRadius, refPosition.height);
}

//! \copydoc Atmosphere::updatePressureAtSurface
void MarsAtmosphere::updatePressureAtSurface()
{
  // pressureAtSurface has already been calculated in MarsGCM update
}

//! \copydoc Atmosphere::updateMetrics
void MarsAtmosphere::updateMetrics()
{
  Atmosphere::updateMetrics();

  MOLATopography mola;
  albedo = mola.getAlbedo(latitude, longitude);

  // Compute dust density variables by methods of Haberle et al., 
  // Icarus, 50, 322 (1982) and Haberle et al., J. Geophys. Res., 104, 8957 (1999)                                            
  // Dust column areal density (kg/m**2)                          
  dustColumnArealDensity = 5.0e-3 * marsAtmos.dustOpticalDepth;

  // Dust mixing ratio (kg dust/kg air) at surface                
  greal qsurf = dustColumnArealDensity * gravity / (0.994 * exp(-dustNu) * pressureAtSurface);

  // Dust mixing ratio at current position and pressure           
  greal expfact = dustNu * (1.0 - 1.0 / sigmaLevel);
  if (expfact > -85.0) {
    dustMixingRatio = qsurf * exp(expfact);
  }
  else {
    dustMixingRatio = 0.0;
  }

  // Dust mass density (micrograms dust / m**3)                   
  dustMassDensity = 1.0e9 * dustMixingRatio * density;

  // Dust number density (number dust particles / m**3)           
  dustNumberDensity = dustMassDensity / (5.23599e-10 * dustDensity * pow(dustDiameter, 3));

}

//! \brief Applies the wave perturbation factor to pressure and density metrics.
//!
//! This routine applies the wave perturbation factor to pressure and density metrics.
//!
//! \b Inputs
//! \arg #atmos      
//! \arg #wavePerturbation      
//!
//! \retval #atmos
void MarsAtmosphere::applyWavePerturbation()
{
  pressure *= wavePerturbation;
  pressureAtSurface *= wavePerturbation;
  density *= wavePerturbation;
  marsAtmos.pressureDaily *= wavePerturbation;
  marsAtmos.densityDaily *= wavePerturbation;
  marsAtmos.densityMin *= wavePerturbation;
  marsAtmos.densityMax *= wavePerturbation;

  wavePerturbation -= 1.0;  // for reporting
}

//! \brief Updates the wave perturbation factor.
//!
//!  Computes relative density and pressure perturbation due to wave model for longitude dependent
//!  (terrain - fixed) traveling or standing waves.  The coefficients for the wave model are input
//!  in the namelist file or can be read in from a wavefile. Use of the wavefile allows for the   
//!  coefficientws to change as a function of time during a given simulation.                     
//!
//! \b Inputs
//! \arg #waveDataList      
//! \arg #elapsedTime      
//! \arg #position      
//! \arg #time      
//!
//! \retval #wavePerturbation
void MarsAtmosphere::updateWavePertubation()
{
  // If a wave file was loaded, then find the wave coefficients.                
  if (dataLoaded) {
    // Loop through wave data entries.
    size_t i;
    for (i = 0; i < waveDataList.size() - 1; ++i) {
      // Find the entries with times that bound the current time.
      if (elapsedTime < waveDataList[i+1].waveTime && elapsedTime >= waveDataList[i].waveTime) {
        break;
      }
    }
    // This is either the lower bounding entry, or the last entry in the list.
    waveData = waveDataList[i];
  }

  double abslat = abs(latitude);  

  // Compute the height factor.
	double heightFactor;                            
  if (height < 100.0) {
    heightFactor = exp((height - 100.0) / waveScale);
  }
  else {
    heightFactor = 1.0;
  }

  // Days between current julian date and waveDate
  greal julianDate;
  time.getTime(inputParameters.timeScale, inputParameters.timeFrame, julianDate);
  double deltaDays;
  if (waveData.waveDate <= 0.0) {
    deltaDays = 0.0;
  }
  else {
    deltaDays = julianDate - waveData.waveDate;
  }

  //  compute wave perturbation as sum of mean and wave 1, 2, and 3 terms using west positive longitudes
  wavePerturbation = heightFactor * (waveData.waveMeanOffset
    + waveData.waveAmplitude1 * cos(TO_RADIANS * (position.getLongitude(WEST_POSITIVE) - waveData.wavePhase1 - waveData.wavePhase1Rate * deltaDays))
    + waveData.waveAmplitude2 * cos(2.0 * TO_RADIANS * (position.getLongitude(WEST_POSITIVE) - waveData.wavePhase2 - waveData.wavePhase2Rate * deltaDays))
    + waveData.waveAmplitude3 * cos(3.0 * TO_RADIANS * (position.getLongitude(WEST_POSITIVE) - waveData.wavePhase3 - waveData.wavePhase3Rate * deltaDays)) - 1.0);

  //  make adjustments near the poles
  if (abslat >= 85.0_deg) { 
    double poleFactor = square(cos(18.0 * toRadians(abslat - 85.0_deg))); 
    wavePerturbation = poleFactor * (wavePerturbation - (waveData.waveMeanOffset - 1.0) * heightFactor)
      + heightFactor * (waveData.waveMeanOffset - 1.0);
  }

  //  ensure wavepert is in proper range (-0.9 to 9.0)
  wavePerturbation = 1.0 + clamp(wavePerturbation, -0.9, 9.0);
} 

//! \copydoc Atmosphere::updateReferenceValues()
void MarsAtmosphere::updateReferenceValues()
{
  // Get the reference values for the current height from the Neptune minmax model.
  marsgcmPtr->getReferenceValues(height, referenceTemperature, referencePressure, referenceDensity);
}

//********************************* Perturbation Factor Methods ********************************//

//! \brief Get high and low perturbation factors.
//!
//! \param[out] pertLow A greal.
//! \param[out] pertHigh A greal.
//! \retval pertLow The low perturbation factor.
//! \retval pertHigh The high perturbation factor.
void MarsAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
{
  static const greal maxHeight = marsgcmPtr->getMaxHeight();
  const greal surfacePerturbationPercent = 0.02;
  // Thermospheric perturbation sigmas assumed to be 45%
  const greal thermospherePerturbationPercent = 0.45;

  if (height <= maxHeight) {
    greal perturbationPercent = 0.0;
    if (height <= marsgcmPtr->getSurfaceFairingHeight()) {
      perturbationPercent = surfacePerturbationPercent;
    }
    else if (height < 100.0) {
      perturbationPercent = max(surfacePerturbationPercent, 0.30 * exp((height - 100.0) / 35.0));
    }
    else {
      perturbationPercent = min(thermospherePerturbationPercent, 0.30 + 0.005 * (height - 100.0));
    }
    pertHigh = 1.0 + perturbationPercent;
  }
  else {
    pertHigh = 1.0 + thermospherePerturbationPercent;
  }
  pertLow = 1.0 / pertHigh;
}

//! \brief Get vertical and horizontal scale parameters.
//!
//! \param[out] verticalScale A greal.
//! \param[out] horizontalScale A greal.
//! \retval verticalScale The vertical scale parameter.
//! \retval horizontalScale The horizontal scale parameter.
void MarsAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
{
  verticalScale = 8.0 * perturbationWaveLengthScale;
  horizontalScale = min(600.0, 30.0 + 0.01875 * position.height * position.height);
  horizontalScale *= perturbationWaveLengthScale;
}

//! \brief Get wind standard deviations.
//!
//! \param[out] ewStdDev A greal.
//! \param[out] nsStdDev A greal.
//! \param[out] vertStdDev A greal.
//! \retval ewStdDev The east/west wind standard deviation.
//! \retval nsStdDev The north/south wind standard deviation.
//! \retval vertStdDev The vertical wind standard deviation.
void MarsAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
{
  greal sigmaU = 2.0 + 0.1 * height;
  greal htAboveMolaSurface = height - surfaceHeight;
  // If height is too negative, use user input height.
  if (height < -8.7_km) {
    htAboveMolaSurface = inputParameters.heightAboveSurface;
  }

  // Added contribution to SIGU for near-surface heights.
  if (htAboveMolaSurface >= 0.0 && htAboveMolaSurface <= 4.5) {
    sigmaU += 1.5 * (1.0 - htAboveMolaSurface / 4.5);
  }
  sigmaU = min(25.0, sigmaU);
  ewStdDev = sigmaU * ewWindPerturbationScale;
  ewStdDev = min(50.0, ewStdDev);

  greal sigmaW = ewStdDev / 5.0;
  // Added contribution to SIGW for near-surface heights.
  if (htAboveMolaSurface >= 0.0 && htAboveMolaSurface <= 4.5) {
    sigmaW += 1.5 * (1.0 - htAboveMolaSurface / 4.5)* ewWindPerturbationScale;
  }
  vertStdDev = sigmaW * verticalWindPerturbationScale;
  nsStdDev = ewStdDev;
}

//! \brief Initializes wave perturbation parameters.
//!
//! This routine will initialize the wave perturbation parameters.  If a wave file name is
//! present, then the wave data list will be loaded from the file.
//!
//! \param params The Mars input parameter structure.    
//!
//! \retval #waveData
//! \retval #waveDataList
void MarsAtmosphere::initializeWaveData(const MarsInputParameters& params)
{
  dataLoaded = false;
  // Initialize wave values to the input parameters
  waveData.waveDate = params.waveDate;
  waveData.waveMeanOffset = params.waveMeanOffset;
  waveData.waveAmplitude1 = params.waveAmplitude1;
  waveData.waveAmplitude2 = params.waveAmplitude2;
  waveData.waveAmplitude3 = params.waveAmplitude3;
  waveData.wavePhase1 = params.wavePhase1;
  waveData.wavePhase2 = params.wavePhase2;
  waveData.wavePhase3 = params.wavePhase3;
  waveData.wavePhase1Rate = params.wavePhase1Rate;
  waveData.wavePhase2Rate = params.wavePhase2Rate;
  waveData.wavePhase3Rate = params.wavePhase3Rate;

  // The wave formulas use west positive longitudes.
  // Convert east positive to west positive
  if (inputParameters.isEastLongitudePositiveOnInput) {
    waveData.wavePhase1 = 360.0 - waveData.wavePhase1;
    waveData.wavePhase2 = 360.0 - waveData.wavePhase2;
    waveData.wavePhase3 = 360.0 - waveData.wavePhase3;
    waveData.wavePhase1Rate *= -1.0;
    waveData.wavePhase2Rate *= -1.0;
    waveData.wavePhase3Rate *= -1.0;
  }

  // If a wave file name is present, read in the wave file data.
  if (!params.waveFile.empty()) {
    // Start with a clear list.
    waveDataList.clear();

    // Open the input stream
    ifstream dataFile(params.waveFile);
    if (!dataFile) {
      throw string(FILE_OPEN_ERROR_MESSAGE);
    }

    // Flag for first line of input
    bool first = true;

    // Read to the end of the file.
    string restOfLine;
    while (!dataFile.eof()) {
      string lineBuffer = "";
      getline(dataFile, lineBuffer);

      if (lineBuffer.empty()) {
        continue;
      }

      istringstream lineInput(lineBuffer);
      WaveData entry;
      lineInput >> entry.waveTime
        >> entry.waveMeanOffset
        >> entry.waveDate
        >> entry.waveAmplitude1
        >> entry.wavePhase1
        >> entry.wavePhase1Rate
        >> entry.waveAmplitude2
        >> entry.wavePhase2
        >> entry.wavePhase2Rate
        >> entry.waveAmplitude3
        >> entry.wavePhase3
        >> entry.wavePhase3Rate;

      // Perform error checks
      // First time must be zero.
      if (first && entry.waveTime != 0.0) {
        throw string("Error: First data in wave file must be at time 0.\n"
                     "       Please verify the WaveTime parameter.");
      }
      // Insist on increasing times
      if (!waveDataList.empty() && entry.waveTime <= waveDataList.back().waveTime) {
        throw string("Error: Time must be in increasing order in the wave file.\n"
                     "       Please verify the wave file data.");
      }
      // The mean must be in range.
      if (entry.waveMeanOffset < 0.1 || entry.waveMeanOffset > 12.0) {
        throw string("Error: The wave mean offset from the wave file is out of range (0.1 to 12).\n"
                     "       Please verify the WaveMeanOffset parameter.");
      }
      first = false;

      // Convert east positive to west positive
      if (inputParameters.isEastLongitudePositiveOnInput) {
       entry.wavePhase1 = 360.0 - entry.wavePhase1;
       entry.wavePhase2 = 360.0 - entry.wavePhase2;
       entry.wavePhase3 = 360.0 - entry.wavePhase3;
       entry.wavePhase1Rate *= -1.0;
       entry.wavePhase2Rate *= -1.0;
       entry.wavePhase3Rate *= -1.0;
      }

      // Passed all error checks.  Save the data.
      waveDataList.emplace_back(entry);
    }

    // Close the input stream.
    dataFile.close();
    dataLoaded = true;
    // Set the wave data to the first element.
    waveData = waveDataList[0];
  }
}
  
//! \brief Updates wind values with slope wind adjustments.
//!
//! The method uses the SlopeWindsModel to adjust the east/west and north/south winds and daily mean winds
//! due to thermally induced mesoscale upslope flow.  Values are also derived for vertical winds.
//!
//! \b Inputs
//! \arg meanWindsScale      
//! \arg boundaryLayerWindsScale      
//! \arg #position      
//! \arg #solarTime      
//! \arg #temperature      
//! \arg #ewWind      
//! \arg #nsWind      
//! \arg #ewWindDaily      
//! \arg #nsWindDaily      
//!
//! \retval #ewWind      
//! \retval #nsWind      
//! \retval #verticalWind      
//! \retval #ewWindDaily      
//! \retval #nsWindDaily
void MarsAtmosphere::updateSlopeWinds()
{
  // Set inputs to the slope winds model.
  SlopeWindsModel slopeWinds;
  slopeWinds.setAreoidRadiusCallback(getAreoidRadiusCallback);
  slopeWinds.setTopographicHeightCallback(getTopographicHeightCallback);
  slopeWinds.setCallbackData(callbackDataPointer);
  slopeWinds.setMeanWindsScale(inputParameters.meanWindsScale);
  slopeWinds.setBoundaryLayerWindsScale(inputParameters.boundaryLayerWindsScale);
  slopeWinds.setPosition(position);
  slopeWinds.setSolarTime(solarTime);
  slopeWinds.setTemperature(temperature);
  slopeWinds.setSpeedOfSound(speedOfSound);
  slopeWinds.setWinds(ewWind, nsWind);
  slopeWinds.setDailyAverageWinds(ewWindDaily, nsWindDaily);

  // Update the model
  slopeWinds.update();

  // Get the modified winds.
  slopeWinds.getWinds(ewWind, nsWind, verticalWind);
  slopeWinds.getDailyAverageWinds(ewWindDaily, nsWindDaily);
}

//! \fn  MarsAtmosphere::setAreoidRadiusCallback(TopoCallback callback)
//! \copydoc MOLATopography::setAreoidRadiusCallback()

//! \fn  MarsAtmosphere::setTopographicHeightCallback(TopoCallback callback)
//! \copydoc MOLATopography::setTopographicHeightCallback()

//! \fn  MarsAtmosphere::setCallbackData(void* dataPointer)
//! \copydoc MOLATopography::setCallbackData()


} // namespace
