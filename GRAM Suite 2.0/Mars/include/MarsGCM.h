//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "MarsCommon.h"
#include "Atmosphere.h"
#include "MarsInputParameters.h"
#include "MarsAtmosphereState.h"
#include "MGCM.h"
#include "TesGCM.h"
#include "StewartModel.h"

namespace GRAM {

//! \brief The interface for the MGCM and TesGCM models.
//!
//! This class defines the interface for the Mars General Circulation Model (MGCM), the Thermal
//! Emission Spectrometer General Circulation Model (TesGCM), and the Steward model.  These two models
//! are programmatically similar except for their dust models.  This class manages the choice
//! of the models two GCM models and the blending with the Stewart model used for exospheric states.
//!
//! \ingroup MarsGRAM
class MarsGCM : public Atmosphere, public MarsCommon
{
public:
  MarsGCM();
  MarsGCM(const MarsGCM& orig);
  virtual ~MarsGCM() override;

  void setMapYear(int year);
  void setWavePerturbation(greal wavePert) { wavePerturbation = wavePert; }
  void setInputParameters(const MarsInputParameters& params);

  void update() override;

  void getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity);
  greal getSurfaceFairingHeight() const { return gcmModelPtr->getSurfaceFairingHeight(); }
  greal getMaxHeight() const { return maxHeight; }

private:
  void updateLowerAtmosphere();
  void updateUpperAtmosphere();
  void updateExosphere();
  void updateGases();
  greal getSeasonalPressure(greal ls, greal lat);
  greal getSpecificHumidity(greal rh, greal t, greal p);

  // Input Parameters
  int mapYear = 0;                           //!< The selected map year for TES data. Zero for MGCM.
  greal exosphericTemperatureOffset = 0.0;   //!< User supplied offset for Stewart model \units{K}.
  greal wavePerturbation = 1.0;              //!< A wave perturbation factor.

  // Internal
  MarsGCMBase* gcmModelPtr = NULL;   //!< Pointer to the MGCM or TES model, as appropriate.
  StewartModel* stewartPtr = NULL;   //!< Pointer to the Stewart exosphere model.
  MarsAtmosphereState marsAtmos;     //!< Mars specific metrics added to the AtmosphereState.

  greal thermosphereBasePressureScaleHeight = 0; //!< Scale height \units{km}.
  greal maxHeight = 0.0;  //!< The max height boundary in the upper atmosphere model \units{km}.

  //! COSPAR Northern hemisphere mean reference data, as given in Pitts et al.,  
  //! "The Mars Atmosphere : Observations and Model Profiles for Mars Missions", NASA JSC - 24455, 1990
  struct CosparData {
    greal height;      //!< \units{ km     }.
    greal temperature; //!< \units{ K      }.
    greal pressure;    //!< \units{ mbar   }.
    greal density;     //!< \units{ g/cm^3 }.
  };

  static const size_t COSPAR_SIZE = 164;                //!< Size of the COSPAR data set.
  static const CosparData cosparDataList[COSPAR_SIZE];  //!< COSPAR reference data (-10 to 360 km).

  // Convenience References
  greal& temperatureDaily = marsAtmos.temperatureDaily;                  //!< \copydoc MarsAtmosphereState::temperatureDaily
  greal& pressureDaily = marsAtmos.pressureDaily;                        //!< \copydoc MarsAtmosphereState::pressureDaily
  greal& densityDaily = marsAtmos.densityDaily;                          //!< \copydoc MarsAtmosphereState::densityDaily
  greal& ewWindDaily = marsAtmos.ewWindDaily;                            //!< \copydoc MarsAtmosphereState::ewWindDaily
  greal& nsWindDaily = marsAtmos.nsWindDaily;                            //!< \copydoc MarsAtmosphereState::nsWindDaily
  greal& densityMin = marsAtmos.densityMin;                              //!< \copydoc MarsAtmosphereState::densityMin
  greal& densityMax = marsAtmos.densityMax;                              //!< \copydoc MarsAtmosphereState::densityMax
  greal& temperatureMin = marsAtmos.temperatureMin;                      //!< \copydoc MarsAtmosphereState::temperatureMin
  greal& temperatureMax = marsAtmos.temperatureMax;                      //!< \copydoc MarsAtmosphereState::temperatureMax
  greal& exosphericTemperature = marsAtmos.exosphericTemperature;        //!< \copydoc MarsAtmosphereState::exosphericTemperature
  greal& groundTemperature = marsAtmos.groundTemperature;                //!< \copydoc MarsAtmosphereState::groundTemperature
  greal& thermosphereBaseHeight = marsAtmos.thermosphereBaseHeight;      //!< \copydoc MarsAtmosphereState::thermosphereBaseHeight
  greal& thermosphereBaseTemperature = marsAtmos.thermosphereBaseTemperature; //!< \copydoc MarsAtmosphereState::thermosphereBaseTemperature
  greal& heightOffset = marsAtmos.heightOffset;                          //!< \copydoc MarsAtmosphereState::heightOffset
  greal& localHeightOffset = marsAtmos.localHeightOffset;                //!< \copydoc MarsAtmosphereState::localHeightOffset
  greal& dustOpticalDepth = marsAtmos.dustOpticalDepth;                  //!< \copydoc MarsAtmosphereState::dustOpticalDepth
  greal& f1PeakHeight = marsAtmos.f1PeakHeight;                          //!< \copydoc MarsAtmosphereState::f1PeakHeight

};

} // namespace
