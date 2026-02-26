//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "unittest_friend.h"
#include "gram.h"
#include "EarthCommon.h"
#include "Atmosphere.h"
#include "NCEP.h"
#include "MERRA2.h"
#include "MSIS.h"
#include "MET.h"
#include "JB2008.h"
#include "MAP.h"
#include "EarthInputParameters.h"
#include "EarthAtmosphereState.h"

namespace GRAM {

struct EarthData;
typedef struct EarthData EarthModelData;

//! \brief The model of a Earth atmosphere.
//!
//! The EarthModel class ties together various models of the regions of the Earth's 
//! atmosphere to create a complete model.  The NCEP class models the lower atmosphere.
//! The MAP class models the middle atmosphere.  And the user may choose from the MET,
//! MSIS, and JB2008 models of the thermosphere.  This class fairs the data in the transition
//! layers.
//!
//! \ingroup EarthGRAM
class EarthModel : public Atmosphere, public EarthCommon
{
public:
  EarthModel();
  EarthModel(const EarthModel& orig);
  virtual ~EarthModel() override = default;

  void setInputParameters(const EarthInputParameters& params);
  void setThermosphereModel(ThermosphereModelType model) { thermosphereModel = model; }
  void setUseNCEP(bool flag) { useMERRA2 = !flag; }
  void setNCEPParameters(int NCEPYear, int NCEPHour) { ncep.setYearAndHour(NCEPYear, NCEPHour); }
  void setMERRA2Parameters(int M2Hour, greal latMin, greal latMax, greal lonMin, greal lonMax)
       { merra2.setHour(M2Hour); merra2.setExtents(latMin, latMax, lonMin, lonMax); }
  void setSolarParameters(const EarthInputParameters& params);
  void setJB2008Parameters(const EarthInputParameters& params);
  void setPerturbationScales(const EarthInputParameters& params);

  void setDayOfYear(greal doy, greal jDay);

  void update() override;

  void getScaleHeights(greal pres, greal dens, greal temp, greal& pressureScaleHeight, greal& densityScaleHeight);
  void getStandardDeviations(AtmosphereState& inState, bool init);
  greal getEpsilon() { return epsilon; }
  greal getPpmToND() { return ppmToND; }

private:
  void updateMAP(AtmosphereState& atmos, greal& dtzX);
  void updateNCEP(AtmosphereState& atmosN, greal &dtzX);
  void updateMERRA2(AtmosphereState& atmosM2, greal& dtdzM2);
  void updateMET(AtmosphereState& atmos, greal& dtzX, greal& dmdzX);
  void updateMSIS(AtmosphereState& atmos, greal& dtzX, greal& dmdzX);
  void updateJB2008(AtmosphereState& atmos, greal& dtzX, greal& dmdzX);
  void updateSpeciesConcentrations();
  greal getMixingRatio(greal tx, greal p);         // was mixrat
  greal getDewPoint(greal t, greal rh, greal p);   // was tdbuck
  void updateSurface(AtmosphereState& inState, bool init);

  // Models
  bool useMERRA2 = true; //!< Set to true to use MERRA-2 model, false for NCEP.
  NCEP ncep;       //!< Lower atmosphere model.
  MERRA2 merra2;   //!< Lower atmosphere model.
  MAP map;         //!< Middle atmosphere model.
  MSIS msis;       //!< Thermosphere model.
  MET met;         //!< Thermosphere model.
  JB2008 jb2008;   //!< Thermosphere model.

  // Input Parameters
  int year = 2010;                                   //!< Current year.
  int month = 0;                                     //!< Current month.
  ThermosphereModelType thermosphereModel = EG_MSIS; //!< Thermosphere model selection.
  greal ewWindPerturbationScale = 1.0;               //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;               //!< Scales North/South wind perturbations.

  EarthPPM ppm;               //!< Parts per million gas concentrations.
  greal ppmToND = 0.0;        //!< Conversion factor from parts per million to number density.
  greal epsilon = 0.0;        //!< Ratio of molecular weights of water to air.
  greal waterSD = 0.0;        //!< Water standard deviation.
  greal ppmWaterLow = 0.0;   //!< Concentraion of water in ppm from the NCEP model.

  greal lowerLowAtmosFairingHeight = 25.0;  //!< 
  greal upperLowAtmosFairingHeight = 30.0;  //!< 
  const greal lowerThermoFairingHeight = 90.0;   //!< Use MAP model below this height. \note was hj1
  const greal upperThermoFairingHeight = 120.0;  //!< Use thermosphere model above this height.  \note was hj2

  greal vaporPressureLow = 0.0;       //!< Vapor pressure from the NCEP model.  \note was vpn
  greal vaporPressureSDLow = 0.0;     //!< Vapor pressure standard deviation from the NCEP model.  \note was svpn
  greal relativeHumidityLow = 0.0;    //!< Relative humidity from the NCEP model.  \note was rhn
  greal relativeHumiditySDLow = 0.0;  //!< Relative humidity standard deviation from the NCEP model.  \note was srhn
  greal temperatureSDLow = 0.0;       //!< Temperature standard deviation from the NCEP model.  \note was stg
  greal pressureSDLow = 0.0;
  greal densitySDLow = 0.0;
  greal ewWindSDLow = 0.0;
  greal nsWindSDLow = 0.0;
  greal windCorrelationLow = 0.0;

  greal dtdz = 0.0;  //!< Temperature gradient.
  greal dmdz = 0.0;  //!< Molecular Weight Gradient

  EarthAtmosphereState earthAtmos;  //!< Earth specific metrics added to the AtmosphereState.

  // References for convenience.
  greal& temperatureStandardDeviation = atmos.temperatureStandardDeviation;   //!< \copydoc AtmosphereState::temperatureStandardDeviation    \note was tstd
  greal& windSpeed = earthAtmos.windSpeed;                                    //!< \copydoc EarthAtmosphereState::windSpeed                  \note was spdavg
  greal& windSpeedStandardDeviation = earthAtmos.windSpeedStandardDeviation;  //!< \copydoc EarthAtmosphereState::windSpeedStandardDeviation \note was spdsd
  greal& dewPoint = earthAtmos.dewPoint;                                      //!< \copydoc EarthAtmosphereState::dewPoint                   \note was tdmean
  greal& dewPointSD = earthAtmos.dewPointSD;                                  //!< \copydoc EarthAtmosphereState::dewPointSD                 \note was stdd
  greal& vaporPressure = earthAtmos.vaporPressure;                            //!< \copydoc EarthAtmosphereState::vaporPressure              \note was eoft
  greal& vaporPressureSD = earthAtmos.vaporPressureSD;                        //!< \copydoc EarthAtmosphereState::vaporPressureSD            \note was seoft
  greal& relativeHumidity = earthAtmos.relativeHumidity;                      //!< \copydoc EarthAtmosphereState::relativeHumidity           \note was rhp
  greal& relativeHumiditySD = earthAtmos.relativeHumiditySD;                  //!< \copydoc EarthAtmosphereState::relativeHumiditySD         \note was srhp
  greal& vaporDensity = earthAtmos.vaporDensity;                              //!< \copydoc EarthAtmosphereState::vaporDensity               \note was rhov
  greal& vaporDensitySD = earthAtmos.vaporDensitySD;                          //!< \copydoc EarthAtmosphereState::vaporDensitySD             \note was srhov

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(EarthModel, initializeData);
  FRIEND_TEST(EarthModel, getStandardDeviations);
  FRIEND_TEST(EarthModel, fair);
  FRIEND_TEST(EarthModel, update);
#endif // GRAM_UNIT_TEST
};

} // namespace

