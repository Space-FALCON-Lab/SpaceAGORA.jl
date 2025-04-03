//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include <algorithm>
#include <cmath>
#include "MarsGCM.h"
#include "Interpolator.h"
#include "SlopeWindsModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MarsGCM::MarsGCM()
  : MarsCommon(this)
{
  atmos.setPlanetSpecificMetrics(marsAtmos);
  stewartPtr = new StewartModel();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
MarsGCM::MarsGCM(const MarsGCM& orig)
  : Atmosphere(orig), MarsCommon(this)
{
  mapYear = orig.mapYear;
  exosphericTemperature = orig.exosphericTemperature;
  wavePerturbation = orig.wavePerturbation;
  thermosphereBaseHeight = orig.thermosphereBaseHeight;
  maxHeight = orig.maxHeight;
  marsAtmos = orig.marsAtmos;
  atmos.setPlanetSpecificMetrics(marsAtmos);
  stewartPtr = new StewartModel(*orig.stewartPtr);
  if (mapYear == 0) {
    // A map year of zero means to use the MGCM model
    gcmModelPtr = new MGCM(*static_cast<const MGCM*>(orig.gcmModelPtr));
  }
  else {
    // Use the TES model if map year is 1 or 2.
    gcmModelPtr = new TesGCM(*static_cast<const TesGCM*>(orig.gcmModelPtr));
  }
}

//! \copydoc Atmosphere::~Atmosphere()
MarsGCM::~MarsGCM()
{
  if (gcmModelPtr != NULL) {
    delete gcmModelPtr;
  }
  if (stewartPtr != NULL) {
    delete stewartPtr;
  }
}

//! \fn  MarsGCM::setWavePerturbation(greal wavePert)
//! \brief Set the wave perturbation factor.
//! \param wavePert The wave perturbation factor.

//! \fn  MarsGCM::getMaxHeight()
//! \brief Returns the maximum height boundary in the upper atmosphere model \units{km}.

//! \fn  MarsGCM::getSurfaceFairingHeight()
//! \brief Returns the lower bound of the surface fairing region \units{km}.

//! \brief Set the TES map year.
//!
//! This method will select the model to be used by choice of the map year.
//! The MGCM model will be used when the map year is zero.  The TES model will
//! be used for map years 1 and 2.  Values outside of this range will be clamped.
//!
//! \param year The selected map year.
void MarsGCM::setMapYear(int year)
{
  // valid years are 0, 1, or 2
  mapYear = clampSize(year, 3);

  // If the model has already been set, then delete it.
  if (gcmModelPtr != NULL) {
    delete gcmModelPtr;
  }

  if (mapYear == 0) {
    // A map year of zero means to use the MGCM model
    gcmModelPtr = new MGCM();
  }
  else {
    // Use the TES model if map year is 1 or 2.
    gcmModelPtr = new TesGCM();
  }
}

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MarsGCM::setInputParameters(const MarsInputParameters& params)
{
  // Set local parameters.
  exosphericTemperatureOffset = params.exosphericTemperatureOffset;

  // Setting the map year also sets gcmModelPtr.
  setMapYear(params.mapYear);

  // Pass parameters to models.
  gcmModelPtr->setInputParameters(params);
  stewartPtr->setInputParameters(params);
}

//! \brief Interface for the primary atmosphere computations.
//!
//! This defines the interface for the primary atmosphere computations. 
void MarsGCM::update()
{
  // Guard against a null model pointer.
  if (gcmModelPtr == NULL) {
    setMapYear(mapYear);
  }
  // Prefer to use dots to arrows.  Set a local reference.
  MarsGCMBase& gcmModel = *gcmModelPtr;

  // Set the current ephemeris in the Mars models.
  gcmModel.setPosition(position);
  gcmModel.setEphemerisState(ephem);

  // Update the dust model.
  gcmModel.updateDustModel();

  // Temporarily replace the height with the surface height and perform a model update.
  greal saveHeight = height;
  height = surfaceHeight;
  gcmModel.setPosition(position);
  gcmModel.update();
  height = saveHeight;

  // Get surface values (pressureAtSurface, groundTemperature)
  const MarsAtmosphereState& marsState = gcmModel.getAtmosphereState().getPlanetSpecificMetrics< MarsAtmosphereState>();
  greal grndTemperature = marsState.groundTemperature;
  greal surfacePressure = gcmModel.getAtmosphereState().pressure;

  // Get the max height boundary in the upper atmosphere model.
  maxHeight = gcmModel.getMaxHeight();

  // Compute the atmosphere state based on the atmosphere layer (lower, upper, exo).
  // Update the atmosphere
  if (height <= 80.0_km) {
    // The MGCM model is used for heights below 80 km.
    updateLowerAtmosphere();
  }
  else if (height <= maxHeight) {
    // The MTGCM model is used between 80 km and maxHeight.
    updateUpperAtmosphere();
  }
  else {
    // The Stewart model is used above maxHeight.
    updateExosphere();
  }

  // Set the saved surface values.
  groundTemperature = grndTemperature;
  pressureAtSurface = surfacePressure;

  // Compute the altitude of the peak F1 ionization.
  greal csza = cos(toRadians(ephem.solarZenithAngle));
  if (csza > 0.0 && thermosphereBaseHeight < 900.0) {
    greal airmass = 1.0 / (csza + 0.15 * pow(93.885 - ephem.solarZenithAngle, -1.253));
    f1PeakHeight = thermosphereBaseHeight + thermosphereBasePressureScaleHeight * log(airmass);
  }
  else {
    f1PeakHeight = 999.9;
  }

  // Compute species concentrations.
  updateGases();
}

//! \brief Updates the AtmosphereState for the homosphere.
//!
//! The AtmosphereState is updated using the appropriate MGCM or TES model.  The constituent
//! gas values from the surface are used throughout the homosphere.  The average molecular weight
//! is held constant throughout the homosphere.  And the thermosphereBasePressureScaleHeight is 
//! set to zero.
//!
//! \b Inputs
//! \arg #position
//!
//! \retval #atmos
void MarsGCM::updateLowerAtmosphere()
{
  // Get the atmosphere state for the current height.
  MarsGCMBase& gcmModel = *gcmModelPtr;
  gcmModel.setPosition(position);
  gcmModel.update();
  atmos = gcmModel.getAtmosphereState();

  // Assign molecular weight (constant in homosphere) 
  averageMolecularWeight = amw;

  // This is zero in the homosphere.
  thermosphereBasePressureScaleHeight = 0.0;
}

//! \brief Updates the AtmosphereState for the thermosphere.
//!
//! The AtmosphereState is updated using the appropriate MGCM or TES model in the thermosphere.  
//! The average molecular weight, thermosphere base values, and the exosphere temperature are computed.
//!
//! \b Inputs
//! \arg #position
//!
//! \retval #atmos
void MarsGCM::updateUpperAtmosphere()
{
  // Get the atmosphere state for the current height
  MarsGCMBase& gcmModel = *gcmModelPtr;
  gcmModel.setPosition(position);
  gcmModel.update();
  atmos = gcmModel.getAtmosphereState();

  // Establish average molecular weight
  averageMolecularWeight = 1000.0 * UNIVERSAL_GAS * density * temperature / pressure;

  // Find metrics at thermosphere Base Height (= altitude of 1.26 nbar level)
  if (thermosphereBaseHeight > 900.0) {
    // Height is too large, set temp to bogus value.
    thermosphereBaseTemperature = 999.9;
  }
  else {
    // Evaluate the atmospheric state at the thermosphere base height
    Position zfpos = position;
    zfpos.height = thermosphereBaseHeight;
    gcmModel.setPosition(zfpos);
    gcmModel.update();
    const AtmosphereState& zfAtmos = gcmModel.getAtmosphereState();

    // Save off T and H
    thermosphereBaseTemperature = zfAtmos.temperature;
    thermosphereBasePressureScaleHeight = zfAtmos.pressureScaleHeight;

    // Adjust the base height due to the wave perturbation.
    thermosphereBaseHeight = thermosphereBaseHeight + thermosphereBasePressureScaleHeight * log(wavePerturbation);
  }

  // Compute the exospheric temperature...
  // First, get values at maxHeight km.
  Position topPos = position;
  topPos.height = maxHeight;
  gcmModel.setPosition(topPos);
  gcmModel.update();
  const AtmosphereState& topAtmos = gcmModel.getAtmosphereState();
  const MarsAtmosphereState& topAtmosMars = topAtmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // Evaluate Stewart model at maxHeight for exosphericTemperatureOffset (deltaTEX) adjustment
  // Prefer to use dots to arrows.  Set a local reference.
  StewartModel& stewart = *stewartPtr;
  stewart.setThermosphereBase(topAtmosMars.thermosphereBaseHeight, thermosphereBaseTemperature);
  stewart.setPosition(topPos);
  stewart.setEphemerisState(ephem);
  stewart.setExosphericTemperatureOffset(exosphericTemperatureOffset);
  stewart.update();

  // Adjust the exospheric Temperature Offset by the difference in the MTGCM and the 
  // Stewart model temperatures at maxHeight.
  greal ddtex = exosphericTemperatureOffset + (topAtmos.temperature - stewart.getAtmosphereState().temperature);

  // Evaluate the Stewart model at the current position using the adjusted offset.
  Position curPos = position;
  curPos.height = max(height, thermosphereBaseHeight);
  stewart.setThermosphereBase(thermosphereBaseHeight, thermosphereBaseTemperature);
  stewart.setPosition(curPos);
  stewart.setExosphericTemperatureOffset(ddtex);
  stewart.update();
  // Get the exospheric temperature.
  exosphericTemperature = stewart.getExosphericTemperature();

  // Get the gas values from the Stewart model.
  FOR_ALL_GASES_PRESENT(
    gas = stewart.getConstituentGas(gasType);
  );
}

//! \brief Updates the AtmosphereState for the exosphere.
//!
//! The thermosphere base values are updated using the appropriate MGCM or TES model.  
//! The AtmosphereState is updated using the Stewart model in the exosphere.  
//! The average molecular weight, thermosphere base values, and the exosphere temperature are computed.
//!
//! \b Inputs
//! \arg #position
//!
//! \retval #atmos
void MarsGCM::updateExosphere()
{
  // Compute the thermosphere base height values...
  // Get values at maxHeight km 
  Position topPos = position;
  topPos.height = maxHeight;
  MarsGCMBase& gcmModel = *gcmModelPtr;
  gcmModel.setPosition(topPos);
  gcmModel.update();
  AtmosphereState topAtmos = gcmModel.getAtmosphereState();
  const MarsAtmosphereState& topAtmosMars = topAtmos.getPlanetSpecificMetrics<MarsAtmosphereState>();
  thermosphereBaseHeight = topAtmosMars.thermosphereBaseHeight;

  // The dust optical depth is from the top of the thermosphere.
  dustOpticalDepth = topAtmosMars.dustOpticalDepth;

  // Evaluate the atmospheric state at the thermosphere base height
  Position zfpos = position;
  zfpos.height = thermosphereBaseHeight;
  gcmModel.setPosition(zfpos);
  gcmModel.update();
  const AtmosphereState& zfAtmos = gcmModel.getAtmosphereState();

  // Get thermosphere base height T and H
  thermosphereBaseTemperature = zfAtmos.temperature;
  const MarsAtmosphereState& zfAtmosMars = zfAtmos.getPlanetSpecificMetrics<MarsAtmosphereState>();
  thermosphereBasePressureScaleHeight = zfAtmos.pressureScaleHeight;

  // The local height offset is from the thermosphere base height
  localHeightOffset = zfAtmosMars.localHeightOffset;
  heightOffset = zfAtmosMars.heightOffset;

  // Evaluate Stewart thermosphere at maxHeight for exosphericTemperatureOffset (deltaTEX) adjustment
  StewartModel& stewart = *stewartPtr;
  stewart.setEphemerisState(ephem);
  stewart.setThermosphereBase(thermosphereBaseHeight, thermosphereBaseTemperature);
  stewart.setExosphericTemperatureOffset(exosphericTemperatureOffset);
  stewart.setPosition(topPos);
  stewart.update();

  // Adjust the exospheric Temperature Offset by the difference in the MTGCM and the 
  // Stewart model temperatures at maxHeight.
  greal ddtex = exosphericTemperatureOffset + (topAtmos.temperature - stewart.getAtmosphereState().temperature);

  // Adjust base height for wave perturbation, using scale height.
  thermosphereBaseHeight = thermosphereBaseHeight + thermosphereBasePressureScaleHeight * log(wavePerturbation);

  // Evaluate the Stewart model at the maxHeight using the adjusted offset.
  stewart.setThermosphereBase(thermosphereBaseHeight, thermosphereBaseTemperature);
  stewart.setExosphericTemperatureOffset(ddtex);
  stewart.setPosition(topPos);
  stewart.update();

  // Get density and pressure values at maxHeight.
  greal densTop = stewart.getAtmosphereState().density;
  greal presTop = stewart.getAtmosphereState().pressure;

  // Finally, use the Stewart model at the current height for the atmospheric state...
  stewart.setThermosphereBase(thermosphereBaseHeight, thermosphereBaseTemperature);
  stewart.setExosphericTemperatureOffset(ddtex);
  stewart.setPosition(position);
  stewart.update();

  // Get the atmospheric values.
  exosphericTemperature = stewart.getExosphericTemperature();
  const AtmosphereState& satmos = stewart.getAtmosphereState();

  // Get the TPD and scale heights.
  temperature = satmos.temperature;
  pressure = satmos.pressure;
  density = satmos.density;
  pressureScaleHeight = satmos.pressureScaleHeight;
  densityScaleHeight = satmos.densityScaleHeight;

  // Get the gas values from the Stewart model.
  FOR_ALL_GASES_PRESENT(
    gas = stewart.getConstituentGas(gasType);
  );

  // Adjust the density and pressure by the differences in the MTGCM and the 
  // Stewart models at maxHeight.
  density *= topAtmos.density / densTop;
  pressure *= topAtmos.pressure / presTop;

  // Establish average molecular weight
  averageMolecularWeight = 1000.0 * UNIVERSAL_GAS * density * temperature / pressure;

  // The winds come from the MTGCM model.
  ewWind = topAtmos.ewWind;
  nsWind = topAtmos.nsWind;

  // Set daily average values to zero above maxHeight
  temperatureDaily = 0.0;
  pressureDaily = 0.0;
  densityDaily = 0.0;
  ewWindDaily = 0.0;
  nsWindDaily = 0.0;

  // Set daily max, min Temp and Density to zero above maxHeight
  temperatureMin = 0.0;
  temperatureMax = 0.0;
  densityMin = 0.0;
  densityMax = 0.0;
}


//!  \brief Computes the final species concentration information.
//!                                                                                               
//!  Computes species concentrations as mole (or volume) fraction, mass fraction, and number       
//!  density \units{\#/m^3}. Average mole fraction and isotope ratio data are taken from Kieffer et al.,
//!  editors, "Mars" (1992) and Tables A - 5 and A - 6 of NASA / TM-2001-210935 (2001).           
/*!  \verbatim
  Notes:                                                                                       
  (1) Below 80 km, Mars MGCM assumes pure CO2 atmosphere (for which molecular weight would be  
      44.01).In this height range, molecular weight computed from the perfect gas law relation 
      M = R0 * rho * T / p would give values close to, but not exactly this value.Deviations   
      would be caused by the fact that the ratio of the averages is not the same as the average
      of the ratio, i.e.  Avg(rho) * Avg(T) / Avg(p) /= Avg(rho * T / p)                       
  (2) Below 80 km, this subroutine computes species concentrations and resultant average       
      molecular weight by assumptions given in note (3).  Therefore average molecular weight   
      given by this subroutine is not exactly the same as average molecular weight computed    
      from the perfect gas law[as given in note(1)].                                           
  (3) Below 80 km, this subroutine computes species concentrations by assuming atmospheric mass
      density from MGCM is correct, but species concentrations are calculated using the        
      following assumptions:                                                                   
      (a) Average dry - atmosphere mole fractions(fbar), are assumed.                          
      (b) Mole fractions are adjusted for seasonal variation of mass of CO2 in the atmosphere, 
          due to freezing and sublimation of CO2 at the poles.                                 
      (c) Only the partial pressure of CO2 is assumed to vary seasonally, not the partial      
          pressures of other constituents, which are assumed not to freeze out of the          
          atmosphere.  However, mole fractions of all species vary seasonally due to this      
          effect.                                                                              
      (d) Seasonal variation of total pressure is taken from subroutine PRSEAS_M10, which was  
          used in the original Stewart thermosphere model, and was employed in Mars - GRAM up  
          through version 3.8.                                                                 
      (e) Water vapor concentration is added to the dry atmosphere by assuming relative        
          humidity = 20 % (computed by function qrhtp, the same as used in marsrad.f           
          calculations).                                                                       
  (4) Between 80 km and the base of the Stewart thermosphere(at altitude zbase), a combination 
      of information is used from TesGCM and the modified Stewart model.TesGCM data used in    
      Mars-GRAM do not include calculated mole fractions or number densities. Mole fractions   
      between 80 km and zbase height are computed from the following assumptions :             
      (a) Mole fractions for N2, Ar, and O2 are assumed to be the same as their value at 80 km 
          (from methods described in note(3).                                                  
      (b) Mole fractions for He, H2, H, and H2O are assumed to be zero.                        
      (c) Mole fractions for N2, Ar, O2, and CO are assumed to vary linearly from their values 
          at 80 km[from method in note(3)] to their values at zbase altitude (from the Stewart 
          model).                                                                              
      (d) Mole fractions for CO2 and O are computed from two constraint conditions             
          (i)  Sum of the mole fractions must be 1.0, and                                      
          (ii) Input molecular weight(AMz, from TesGCM value) is preserved.                    
  (5) From height zbase up until TesGCM data runs out, a combination of information is used from
      TesGCM and the modified Stewart model.  In this height range, the following assumptions  
      are used                                                                                 
      (a) Mole fractions for constituents other than CO2 and O are taken from the modified     
          Stewart model                                                                        
      (b) Molecular weight from TesGCM is assumed                                              
      (c) Mole fractions for CO2 and O are computed from two constraint conditions             
          (i)  Sum of the mole fractions must be 1.0, and                                      
          (ii) Input molecular weight(AMz, from TesGCM value) is preserved.                    
  (6) Above the top altitude for which TesGCM data are available, the same methodology is used 
      as in note(5), except that the input value for molecular weight(AMz) is taken directly   
      from the modified Stewart model.                                                         
   \endverbatim */
void MarsGCM::updateGases()
{
  greal fmoltot = 0.0; 

  // Assume dry atmosphere, unless modified for heights up to 80 km
  water.massFraction = 0.0;
  water.moleFraction = 0.0;
  water.numberDensity = 0.0;

  // Two cases, above and below the thermosphere base height.
  // See comments above for details.
	if (height <= thermosphereBaseHeight) {
    // Store input values of N2, Ar, O2, and CO mole fraction
    // (= values at height zbase if input height is between 80 km and zbase)
    ConstituentGas n2_in = dinitrogen;
    ConstituentGas ar_in = argon;
    ConstituentGas o2_in = dioxygen;
    ConstituentGas co_in = carbonMonoxide;
    
    // Assumed average mole fractions for dry Mars atmosphere below 80 km
    // altitude.  See Kieffer et al., editors, "Mars" (1992) and 
    // Table A - 5 of NASA / TM - 2001 - 210935 (2001) 
    // (CO2, N2, Ar, O2, CO, O, He, H2, H)
    const greal fbar[9] = { 0.9537, 0.0275, 0.0165, 0.0013, 0.0010,
                            0.0000, 0.0000, 0.0000, 0.0000 };
    double pr = getSeasonalPressure(longitudeSun, latitude);

    fmoltot = carbonDioxide.moleFraction = 1.0 - (1.0 - fbar[0]) / pr;
    fmoltot += dinitrogen.moleFraction = fbar[1] / pr;
    fmoltot += argon.moleFraction = fbar[2] / pr;
    fmoltot += dioxygen.moleFraction = fbar[3] / pr;
    fmoltot += carbonMonoxide.moleFraction = fbar[4] / pr;
    // Since fbar = 0 below...
    oxygen.moleFraction = 0.0;
    helium.moleFraction = 0.0;
    dihydrogen.moleFraction = 0.0;
    hydrogen.moleFraction = 0.0;

		if (height <= 80.0) {
      averageMolecularWeight = 0.0;
      // Iterate over the set of all activated gases.
      // This macro defines "gas" and "gasType" for each activated gas.
      FOR_ALL_GASES_PRESENT(
        averageMolecularWeight += gas.moleFraction * gas.averageMolecularWeight;
      );
      // Water vapor
      water.massFraction = getSpecificHumidity(0.2, temperature, pressure / 100.0);
			water.moleFraction = averageMolecularWeight * water.massFraction / water.averageMolecularWeight;
			fmoltot += water.moleFraction;
      water.moleFraction *= fmoltot;

      averageMolecularWeight = 0;
      // Iterate over the set of all activated gases.
      // This macro defines "gas" and "gasType" for each activated gas.
      FOR_ALL_GASES_PRESENT(
        gas.moleFraction /= fmoltot;
        averageMolecularWeight += gas.moleFraction * gas.averageMolecularWeight;
      );
    }
		else {
			double hgtfact = (height - 80.0) / (thermosphereBaseHeight - 80.0);

      dinitrogen.moleFraction = dinitrogen.moleFraction + (n2_in.moleFraction - dinitrogen.moleFraction) * hgtfact;
			argon.moleFraction = argon.moleFraction + (ar_in.moleFraction - argon.moleFraction) * hgtfact;
			dioxygen.moleFraction = dioxygen.moleFraction + (o2_in.moleFraction - dioxygen.moleFraction) * hgtfact;
			carbonMonoxide.moleFraction = carbonMonoxide.moleFraction + (co_in.moleFraction - carbonMonoxide.moleFraction) * hgtfact;

      double fsum = 1.0 
        - dinitrogen.moleFraction - argon.moleFraction - dioxygen.moleFraction - carbonMonoxide.moleFraction
        - dihydrogen.moleFraction - hydrogen.moleFraction - helium.moleFraction;

			double fmsum = averageMolecularWeight
        - dinitrogen.moleFraction * dinitrogen.averageMolecularWeight
        - argon.moleFraction * argon.averageMolecularWeight
				- dioxygen.moleFraction * dioxygen.averageMolecularWeight 
        - carbonMonoxide.moleFraction * carbonMonoxide.averageMolecularWeight 
        - helium.moleFraction * helium.averageMolecularWeight  
				- dihydrogen.moleFraction * dihydrogen.averageMolecularWeight 
        - hydrogen.moleFraction * hydrogen.averageMolecularWeight;

			double dd = oxygen.averageMolecularWeight - carbonDioxide.averageMolecularWeight;

			carbonDioxide.moleFraction = (fsum * oxygen.averageMolecularWeight - fmsum) / dd;
			oxygen.moleFraction = (fmsum - carbonDioxide.averageMolecularWeight * fsum) / dd;

			if (oxygen.moleFraction < 0.0) {
				oxygen.moleFraction = 0.0;
				carbonDioxide.moleFraction = fsum;
        averageMolecularWeight += carbonDioxide.moleFraction * carbonDioxide.averageMolecularWeight - fmsum;
			}
		}
	}
	else {
		double fsum = 1.0 - dinitrogen.moleFraction - argon.moleFraction - dioxygen.moleFraction - carbonMonoxide.moleFraction
			- dihydrogen.moleFraction - hydrogen.moleFraction - helium.moleFraction;
		double fmsum = averageMolecularWeight
      - dinitrogen.moleFraction * dinitrogen.averageMolecularWeight
      - argon.moleFraction * argon.averageMolecularWeight
			- dioxygen.moleFraction * dioxygen.averageMolecularWeight
      - carbonMonoxide.moleFraction * carbonMonoxide.averageMolecularWeight
      - helium.moleFraction * helium.averageMolecularWeight
			- dihydrogen.moleFraction * dihydrogen.averageMolecularWeight
      - hydrogen.moleFraction * hydrogen.averageMolecularWeight;

		double dd = oxygen.averageMolecularWeight - carbonDioxide.averageMolecularWeight;

		carbonDioxide.moleFraction = (fsum * oxygen.averageMolecularWeight - fmsum) / dd;
		oxygen.moleFraction = (fmsum - carbonDioxide.averageMolecularWeight * fsum) / dd;

		if (oxygen.moleFraction < 0.0) {
			oxygen.moleFraction = 0.0;
			carbonDioxide.moleFraction = fsum;
      averageMolecularWeight += carbonDioxide.moleFraction * carbonDioxide.averageMolecularWeight - fmsum;
		}
		if (carbonDioxide.moleFraction < 0.0) {
			carbonDioxide.moleFraction = 0.0;
			oxygen.moleFraction = fsum;
      averageMolecularWeight += oxygen.moleFraction * oxygen.averageMolecularWeight - fmsum;
		}
	}

	totalNumberDensity = 0.0;
	double totMass = 0.0;
  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    gas.numberDensity = gas.moleFraction * density * AVOGADRO * 1.0e3 / averageMolecularWeight;
    totalNumberDensity += gas.numberDensity;
    totMass += gas.numberDensity * gas.averageMolecularWeight;
  );

  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    gas.massFraction = gas.numberDensity * gas.averageMolecularWeight  / totMass;
  );

}

//! \brief Computes relative seasonal variation in pressure on a reference ellipsoid due to latitude and
//! season (Ls) variations.                                                                      
//!                                                                                               
//! The source of this parameterization is currently unknown, but was used in the original       
//! Stewart thermosphere model, and was employed in Mars-GRAM up through version 3.8.  
//!
//! \param ls Longitude of the sun \units{degrees}.
//! \param lat Latitude \units{degrees}.
//!
//! \returns pressure
greal MarsGCM::getSeasonalPressure(greal ls, greal lat)
{ 
  greal lat2 = lat * lat;

  // amplitude coefficients
  greal a1 = 0.0847194 - 0.570405e-5 * lat2;
  greal a2 = 0.0690599 - 0.132689e-5 * lat2;

  // phase coefficients
  greal p1 = 304.041 + 0.0080602 * lat2;
  greal p2 = 61.362 + 0.0016533 * lat2;

  // Return seasonal pressure
  return 1.0 + a1 * cos(toRadians(ls - p1)) + a2 * cos(2.0 * toRadians(ls - p2));
}

//! \brief Computes specific humidity.
//!
//! Computes specific humidity as a function of RH, T and P.
//! Taken from Savijarvi, 1991, Contr. Atmos. Phys., 64, 103.
//!
//! \param rh Relative humidity \units{\%} (0-1).
//! \param t Temperature \units{K}.
//! \param p Pressure \units{mbar}.
//!
//! \returns Specific humidity \units{mg/kg}.
greal MarsGCM::getSpecificHumidity(greal rh, greal t, greal p)
{
  greal es = min(p / 1.59, 6.1135 * exp(22.542 * (t - 273.16) / (t + 0.32)));
	return (rh * 0.407 * es) / (p - 0.59 * es);
} 

//! \brief Returns reference TPD from COSPAR data.
//!
//! Interpolates COSPAR mean NH profile to a given height.  Returns 
//! coresponding values for temperature, pressure, and density.             
//! COSPAR values from Table XI, "The Mars Atmosphere: Observations and Model Profiles for 
//!    Mars Missions", David E. Pitts et al., eds., JSC-24455. Values represent COSPAR NH mean
//!    values.                                                                                
//!                                                                          
//! \param height The height \units{km} of the desired reference values.
//! \param[out] refTemperature A greal.
//! \param[out] refPressure A greal.
//! \param[out] refDensity A greal.
//! \retval refTemperature The reference temperature \units{K}.
//! \retval refPressure The reference pressure \units{mb}. 
//! \retval refDensity The reference density \units{g/m3}.
void MarsGCM::getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity)
{
  // Establish the cospar height index for lower z-bounding level
  // Height may be negative, so use int for the index
  int zx;  
  if (height < 130.0_km) {
    // 1 km COSPAR steps starting at -10 km
    zx = int(floor(height + 10.0)); 
  }
  else {
    // 10 km COSPAR steps starting at 130 km
    zx = 140 + int(floor((height - 130.0) / 10.0));
  }

  // error trap for index out of bounds
  if (zx < 0 || zx > int(COSPAR_SIZE - 1)) {
    refTemperature = 0.0;
    refPressure = 0.0;   
    refDensity = 0.0;    
    return;              
  }
  // Allow extrapolation in the 360 to 370 km range.
  zx = min(zx, int(COSPAR_SIZE - 2));

  // Get the bounding COSPAR data
  const CosparData& low = cosparDataList[zx];
  const CosparData& high = cosparDataList[zx + 1];

  // Set up for height interpolation 
  Interpolator zInterp;
  zInterp.makeFraction(low.height, high.height, height);

  // Linearly interpolate temperature
  refTemperature = zInterp.linear(low.temperature, high.temperature);

  // Get the reference pressure.
  if (abs(high.temperature - low.temperature) > 0.01) { 
    // if non-isothermal perform logarithmic interpolation for pressure
    double expon = log(high.pressure / low.pressure) / log(high.temperature / low.temperature);
    refPressure = 100.0 * low.pressure * pow(refTemperature / low.temperature, expon);
  } 
  else { 
    // else (isothermal case) perform exponential interpolation
    refPressure = 100.0 * zInterp.log(low.pressure, high.pressure);
  }

  // compute gas constant at upper and lower levels
  double r1 = low.pressure / (low.density * low.temperature); 
  double r2 = high.pressure / (high.density * high.temperature);
  // linearly interpolate gas constant with height
  double r = zInterp.linear(r1, r2); 

  // compute density from gas law
  refDensity = 10.0 * refPressure / (r * refTemperature); 
}
  
//! COSPAR values from Table XI, "The Mars Atmosphere: Observations and Model Profiles for 
//!    Mars Missions", David E. Pitts et al., eds., JSC-24455. Values represent COSPAR NH mean
//!    values.                                                                                
const MarsGCM::CosparData MarsGCM::cosparDataList[COSPAR_SIZE] = {
// height, temp, pressure, density
  {-10.0, 214.0, 1.57E+01, 3.85E-05},
  { -9.0, 214.0, 1.44E+01, 3.52E-05},
  { -8.0, 214.0, 1.31E+01, 3.21E-05},
  { -7.0, 214.0, 1.20E+01, 2.94E-05},
  { -6.0, 214.0, 1.09E+01, 2.68E-05},
  { -5.0, 214.0, 1.00E+01, 2.45E-05},
  { -4.0, 214.0, 9.16E+00, 2.24E-05},
  { -3.0, 214.0, 8.36E+00, 2.04E-05},
  { -2.0, 214.0, 7.63E+00, 1.87E-05},
  { -1.0, 214.0, 6.96E+00, 1.70E-05},
  {  0.0, 214.0, 6.36E+00, 1.55E-05},
  {  1.0, 213.9, 5.80E+00, 1.42E-05},
  {  2.0, 213.8, 5.30E+00, 1.30E-05},
  {  3.0, 213.6, 4.84E+00, 1.18E-05},
  {  4.0, 213.4, 4.41E+00, 1.08E-05},
  {  5.0, 212.9, 4.03E+00, 9.90E-06},
  {  6.0, 212.4, 3.68E+00, 9.06E-06},
  {  7.0, 210.8, 3.35E+00, 8.32E-06},
  {  8.0, 209.2, 3.06E+00, 7.65E-06},
  {  9.0, 207.1, 2.79E+00, 7.04E-06},
  { 10.0, 205.0, 2.54E+00, 6.47E-06},
  { 11.0, 203.2, 2.31E+00, 5.94E-06},
  { 12.0, 201.4, 2.09E+00, 5.44E-06},
  { 13.0, 199.6, 1.90E+00, 4.99E-06},
  { 14.0, 197.8, 1.73E+00, 4.56E-06},
  { 15.0, 196.2, 1.56E+00, 4.17E-06},
  { 16.0, 194.6, 1.42E+00, 3.81E-06},
  { 17.0, 193.0, 1.28E+00, 3.48E-06},
  { 18.0, 191.5, 1.16E+00, 3.17E-06},
  { 19.0, 189.9, 1.05E+00, 2.89E-06},
  { 20.0, 188.3, 9.47E-01, 2.63E-06},
  { 21.0, 186.8, 8.54E-01, 2.39E-06},
  { 22.0, 185.2, 7.70E-01, 2.18E-06},
  { 23.0, 183.8, 6.94E-01, 1.97E-06},
  { 24.0, 182.5, 6.25E-01, 1.79E-06},
  { 25.0, 181.2, 5.62E-01, 1.62E-06},
  { 26.0, 180.0, 5.05E-01, 1.47E-06},
  { 27.0, 178.7, 4.54E-01, 1.33E-06},
  { 28.0, 177.5, 4.07E-01, 1.20E-06},
  { 29.0, 176.2, 3.66E-01, 1.09E-06},
  { 30.0, 175.0, 3.28E-01, 9.80E-07},
  { 31.0, 173.7, 2.94E-01, 8.84E-07},
  { 32.0, 172.5, 2.63E-01, 7.97E-07},
  { 33.0, 171.2, 2.35E-01, 7.19E-07},
  { 34.0, 170.0, 2.10E-01, 6.47E-07},
  { 35.0, 168.7, 1.88E-01, 5.82E-07},
  { 36.0, 167.5, 1.68E-01, 5.24E-07},
  { 37.0, 166.1, 1.49E-01, 4.71E-07},
  { 38.0, 164.8, 1.33E-01, 4.23E-07},
  { 39.0, 163.6, 1.19E-01, 3.79E-07},
  { 40.0, 162.4, 1.06E-01, 3.40E-07},
  { 41.0, 161.2, 9.38E-02, 3.04E-07},
  { 42.0, 160.0, 8.33E-02, 2.72E-07},
  { 43.0, 159.0, 7.39E-02, 2.43E-07},
  { 44.0, 158.0, 6.56E-02, 2.17E-07},
  { 45.0, 157.0, 5.81E-02, 1.94E-07},
  { 46.0, 156.0, 5.15E-02, 1.73E-07},
  { 47.0, 155.0, 4.56E-02, 1.54E-07},
  { 48.0, 154.1, 4.03E-02, 1.37E-07},
  { 49.0, 153.1, 3.56E-02, 1.22E-07},
  { 50.0, 152.2, 3.15E-02, 1.08E-07},
  { 51.0, 151.2, 2.78E-02, 9.60E-08},
  { 52.0, 150.3, 2.45E-02, 8.52E-08},
  { 53.0, 149.5, 2.16E-02, 7.55E-08},
  { 54.0, 148.7, 1.90E-02, 6.69E-08},
  { 55.0, 147.9, 1.67E-02, 5.92E-08},
  { 56.0, 147.2, 1.47E-02, 5.23E-08},
  { 57.0, 146.4, 1.30E-02, 4.63E-08},
  { 58.0, 145.7, 1.14E-02, 4.09E-08},
  { 59.0, 144.9, 1.00E-02, 3.61E-08},
  { 60.0, 144.2, 8.78E-03, 3.18E-08},
  { 61.0, 143.6, 7.70E-03, 2.80E-08},
  { 62.0, 143.0, 6.75E-03, 2.47E-08},
  { 63.0, 142.5, 5.92E-03, 2.17E-08},
  { 64.0, 142.0, 5.19E-03, 1.91E-08},
  { 65.0, 141.5, 4.54E-03, 1.68E-08},
  { 66.0, 141.0, 3.98E-03, 1.48E-08},
  { 67.0, 140.5, 3.48E-03, 1.30E-08},
  { 68.0, 140.0, 3.04E-03, 1.14E-08},
  { 69.0, 139.7, 2.66E-03, 9.97E-09},
  { 70.0, 139.5, 2.33E-03, 8.73E-09},
  { 71.0, 139.2, 2.04E-03, 7.65E-09},
  { 72.0, 139.0, 1.78E-03, 6.70E-09},
  { 73.0, 139.0, 1.56E-03, 5.85E-09},
  { 74.0, 139.0, 1.36E-03, 5.12E-09},
  { 75.0, 139.0, 1.19E-03, 4.47E-09},
  { 76.0, 139.0, 1.04E-03, 3.91E-09},
  { 77.0, 139.0, 9.09E-04, 3.42E-09},
  { 78.0, 139.0, 7.95E-04, 2.99E-09},
  { 79.0, 139.0, 6.95E-04, 2.62E-09},
  { 80.0, 139.0, 6.08E-04, 2.29E-09},
  { 81.0, 139.0, 5.32E-04, 2.00E-09},
  { 82.0, 139.0, 4.65E-04, 1.75E-09},
  { 83.0, 139.0, 4.07E-04, 1.53E-09},
  { 84.0, 139.0, 3.56E-04, 1.34E-09},
  { 85.0, 139.0, 3.11E-04, 1.17E-09},
  { 86.0, 139.0, 2.72E-04, 1.03E-09},
  { 87.0, 139.0, 2.38E-04, 8.97E-10},
  { 88.0, 139.0, 2.09E-04, 7.85E-10},
  { 89.0, 139.0, 1.83E-04, 6.87E-10},
  { 90.0, 139.0, 1.60E-04, 6.01E-10},
  { 91.0, 139.0, 1.40E-04, 5.26E-10},
  { 92.0, 139.0, 1.22E-04, 4.61E-10},
  { 93.0, 139.0, 1.07E-04, 4.03E-10},
  { 94.0, 139.0, 9.39E-05, 3.53E-10},
  { 95.0, 139.0, 8.22E-05, 3.09E-10},
  { 96.0, 139.0, 7.20E-05, 2.71E-10},
  { 97.0, 139.0, 6.30E-05, 2.37E-10},
  { 98.0, 139.0, 5.52E-05, 2.08E-10},
  { 99.0, 139.0, 4.83E-05, 1.82E-10},
  {100.0, 139.0, 4.23E-05, 1.59E-10},
  {101.0, 140.0, 3.71E-05, 1.39E-10},
  {102.0, 141.1, 3.26E-05, 1.21E-10},
  {103.0, 142.1, 2.86E-05, 1.05E-10},
  {104.0, 143.2, 2.52E-05, 9.17E-11},
  {105.0, 144.2, 2.21E-05, 8.01E-11},
  {106.0, 145.2, 1.95E-05, 7.01E-11},
  {107.0, 146.3, 1.72E-05, 6.13E-11},
  {108.0, 147.3, 1.52E-05, 5.38E-11},
  {109.0, 148.4, 1.34E-05, 4.72E-11},
  {110.0, 149.4, 1.19E-05, 4.14E-11},
  {111.0, 150.4, 1.05E-05, 3.64E-11},
  {112.0, 151.5, 9.33E-06, 3.21E-11},
  {113.0, 152.5, 8.28E-06, 2.82E-11},
  {114.0, 153.5, 7.35E-06, 2.49E-11},
  {115.0, 154.6, 6.53E-06, 2.20E-11},
  {116.0, 155.6, 5.81E-06, 1.94E-11},
  {117.0, 156.6, 5.18E-06, 1.72E-11},
  {118.0, 157.7, 4.61E-06, 1.52E-11},
  {119.0, 158.7, 4.11E-06, 1.35E-11},
  {120.0, 159.7, 3.67E-06, 1.19E-11},
  {121.0, 160.8, 3.28E-06, 1.06E-11},
  {122.0, 161.8, 2.94E-06, 9.41E-12},
  {123.0, 162.8, 2.63E-06, 8.36E-12},
  {124.0, 163.8, 2.35E-06, 7.44E-12},
  {125.0, 164.9, 2.11E-06, 6.63E-12},
  {126.0, 165.9, 1.89E-06, 5.91E-12},
  {127.0, 166.9, 1.70E-06, 5.27E-12},
  {128.0, 167.9, 1.53E-06, 4.70E-12},
  {129.0, 169.0, 1.37E-06, 4.20E-12},
  {130.0, 170.0, 1.23E-06, 3.76E-12},
  {140.0, 245.1, 5.20E-07, 1.09E-12},
  {150.0, 288.6, 2.70E-07, 4.73E-13},
  {160.0, 314.0, 1.53E-07, 2.43E-13},
  {170.0, 328.8, 9.04E-08, 1.35E-13},
  {180.0, 337.5, 5.52E-08, 7.90E-14},
  {190.0, 342.6, 3.45E-08, 4.74E-14},
  {200.0, 345.6, 2.19E-08, 2.92E-14},
  {210.0, 346.4, 1.43E-08, 1.80E-14},
  {220.0, 347.3, 9.53E-09, 1.14E-14},
  {230.0, 348.1, 6.49E-09, 7.42E-15},
  {240.0, 348.9, 4.51E-09, 4.92E-15},
  {250.0, 349.7, 3.19E-09, 3.33E-15},
  {260.0, 349.8, 2.30E-09, 2.24E-15},
  {270.0, 349.8, 1.71E-09, 1.55E-15},
  {280.0, 349.9, 1.29E-09, 1.10E-15},
  {290.0, 349.9, 9.93E-10, 7.94E-16},
  {300.0, 350.0, 7.77E-10, 5.87E-16},
  {310.0, 350.0, 6.17E-10, 4.39E-16},
  {320.0, 350.0, 4.97E-10, 3.34E-16},
  {330.0, 350.0, 4.06E-10, 2.59E-16},
  {340.0, 350.0, 3.35E-10, 2.03E-16},
  {350.0, 350.0, 2.79E-10, 1.61E-16},
  {360.0, 350.0, 2.35E-10, 1.30E-16}
};

} // namespace
