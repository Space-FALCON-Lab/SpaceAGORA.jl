//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"

namespace GRAM {

//! \brief The augmented output state of a Mars atmosphere model.
//!
//! The MarsAtmosphereState provides Mars model specific metrics as an
//! extension to the AtmosphereState class using AtmosphereState::planetSpecificMetrics.
//!
//! \ingroup MarsGRAM Cpp_Mars
class MarsAtmosphereState
{
public:
  MarsAtmosphereState();
  MarsAtmosphereState(const MarsAtmosphereState& orig) = default;
  virtual ~MarsAtmosphereState() = default;

  greal temperatureDaily = 0.0;          //!< Mean daily temperature       \units{ K      }.
  greal pressureDaily = 0.0;             //!< Mean daily pressure          \units{ N/m^2  }.
  greal densityDaily = 0.0;              //!< Mean daily density           \units{ kg/m^3 }.
  greal ewWindDaily = 0.0;               //!< Mean daily east/west winds   \units{ m/s    }. 
  greal nsWindDaily = 0.0;               //!< Mean daily north/south winds \units{ m/s    }.
  greal densityMin = 0.0;                //!< Daily minimun density        \units{ kg/m^3 }.
  greal densityMax = 0.0;                //!< Daily maximun density        \units{ kg/m^3 }.
  greal temperatureMin = 9999.0;         //!< Daily minimun temperature    \units{ K      }.
  greal temperatureMax = -9999.0;        //!< Daily maximun temperature    \units{ K      }.

  greal planetoGraphicHeight = 0.0;      //!< Planetographic height                               \units{ km }.
  greal planetoGraphicLatitude = 0.0;    //!< Planetographic latitude                             \units{ \text{degrees} }.
  greal referenceHeight = 0.0;           //!< Height relative to the reference ellipsoid          \units{ km }.
  greal referenceRadius = 0.0;           //!< Latitude radius relative to the reference ellipsoid \units{ km }.

  greal groundTemperature = 0.0;         //!< Temperature at ground level           \units{ K  }.
  greal thermosphereBaseHeight = 0.0;    //!< Height of 1.26 nbar level             \units{ km }. (was zf)
  greal thermosphereBaseTemperature = 0.0;  //!< Temperature at the 1.26 nbar level \units{ K  }.
  greal exosphericTemperature = 0.0;     //!< Temperature of the exosphere          \units{ K  }.

  greal wavePerturbation = 0.0;          //!< Perturbation factor for traveling or standing waves.
  greal f1PeakHeight = 0.0;              //!< Altitude of the peak F1 ionization       \units{ km }.
  greal albedo = 0.0;                    //!< Surface albedo.
  greal heightOffset = 0.0;              //!< Height offset as selected by offsetModel \units{ km }.
  greal localHeightOffset = 0.0;         //!< Local height offset                      \units{ km }.

  greal dustOpticalDepth = 0.0;          //!< Dust optical depth.
  greal dustColumnArealDensity = 0.0;    //!< Dust column areal density \units{ kg/m^2 }.
  greal dustMixingRatio = 0.0;           //!< Dust mixing ratio         \units{ kg \text{ dust}/kg \text{ air}     }.
  greal dustMassDensity = 0.0;           //!< Dust mass density         \units{ \text{micrograms dust} / m^3       }.
  greal dustNumberDensity = 0.0;         //!< Dust number density       \units{ \text{number dust particles} / m^3 }.
  int   iceIsPresent = 99;               //!< When true (1), CO2 ice is present on the polar surface.
};

} // namespace