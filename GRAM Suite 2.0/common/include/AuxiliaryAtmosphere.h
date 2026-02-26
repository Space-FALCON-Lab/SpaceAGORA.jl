//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <vector>
#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

class AuxiliaryData
{
public:
  AuxiliaryData() {}
  AuxiliaryData(const AuxiliaryData& orig) = default;
  virtual ~AuxiliaryData() = default;

  greal height = 0;      //!< Height Above Reference Ellipsoid \units{km}
  greal latitude = 0;    //!< Geocentric latitude \units{\text{degrees}}.
  greal longitude = 0;   //!< East longitude positive \units{\text{degrees}}.
  greal temperature = 0; //!< \units{K}.
  greal pressure = 0;    //!< \units{Pa}.
  greal density = 0;     //!< \units{kg/m^3}.
  greal nsWind = 0;      //!< East/west winds \units{m/s}.
  greal ewWind = 0;      //!< North/south winds \units{m/s}.
  greal temperatureStandardDeviation = 0.0;     //!< \units{\%}.
  greal pressureStandardDeviation = 0.0;        //!< \units{\%}.
  greal densityStandardDeviation = 0.0;         //!< \units{\%}.
  greal ewStandardDeviation = 0.0;              //!< \units{\%}.
  greal nsStandardDeviation = 0.0;              //!< \units{\%}.

};

//! \brief An AuxiliaryAtmosphere overrides the values of atmospheric model.
//!
//! In the legacy code, an AuxiliaryAtmosphere was referred to as an auxiliary profile.
//! The name change is due to this now being a subclass of Atmosphere.
//! Use the AuxiliaryAdapter class to add AuxiliaryAtmosphere capabilities to 
//! a PerturbedAtmosphere.
//!
//! Data can be loaded into an AuxiliaryAtmosphere by setDataFile() and loadData() or
//! more directly with setData().  The data should be strictly increasing or decreasing
//! in height.  Input longitudes can be east or west positive, but should be declared
//! with setEastLongitudePositive() or setWestLongitudePositive().
//!
//! Before calling update(), the AuxiliaryAtmosphere needs an AtmosphereState from
//! the model via setInputState().  If the position is within the outer bounds set
//! by setBounds(), then the input state will be replaced by the auxiliary values with 
//! a smoothing between the inner and outer bounds.
//! \ingroup CommonGRAM
class AuxiliaryAtmosphere : public Atmosphere
{
public:
  AuxiliaryAtmosphere();
  AuxiliaryAtmosphere(const AuxiliaryAtmosphere& orig);
  virtual ~AuxiliaryAtmosphere() override;

  virtual void loadData();

  void setBounds(greal inner, greal outer);
  void setDataFile(const std::string& fileName);
  void setData(const std::vector<AuxiliaryData>& auxData);
  void setValues(greal dens, greal pres, greal temp, greal ew, greal ns);
  void setEastLongitudePositive(bool flag = true) { eastLongitudePositive = flag; }
  void setWestLongitudePositive() { eastLongitudePositive = false;  }
  void setAuxId(size_t id) { auxId = id; }
  void setSDFlag(bool dataHasStandardDeviations) { hasStandardDeviations = dataHasStandardDeviations; }
 
  void setInputState(const AtmosphereState& state) { inputState = state; }
  virtual void update() override;
  virtual void updateStandardDeviations();

  bool isActive() { return active; }

private:
  greal innerRadius = 0;                //!< (km) The inner bound for fairing.
  greal outerRadius = 0;                //!< (km) The outer bound for fairing.
  bool eastLongitudePositive = true;    //!< Longitude positive convention flag.
  bool active = false;                  //!< True when the AuxiliaryAtmosphere has data.
  size_t auxId = 0;                     //!< Index of the profileWeight in the AtmoshereState.
  greal profileWeight = 0.0;            //!< The fairing factor between 0 and 1.
  bool useValues = false;
  bool hasStandardDeviations = false;

  AtmosphereState inputState;              //!< AtmosphereState from the model.
  std::vector<AuxiliaryData> auxDataTable; //!< A vertical profile of PTD data.
  std::string auxFileName;                 //!< File name of the input data.
  AuxiliaryData auxValues;
};

} // namespace
