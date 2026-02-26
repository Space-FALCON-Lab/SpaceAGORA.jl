//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <map>
#include <array>
#include <vector>
#include <string>
#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief A three-tiered height based atmosphere model base class.
//!
//! The MinMaxModel base class provides a quick way to create an atmospheric
//! model for height based data with minimum, average, and maximum values. 
//! The data for the model should be two-dimensional, depending on height and
//! the three levels.  It should include temperature, pressure, and density
//! values as well as number densities for the appropriate gases.
//!
//! No data members are declared in the MinMaxModel base class so that they may be
//! declared statically in each subclass.  The class definition contains commented lines
//! to be copied into any subclass.  The subclass should also call initializeGases() in
//! its constructor. This class does not implement a winds model.  It is up to each subclass
//! to implement a winds model.
//! \ingroup CommonGRAM
class MinMaxModel : public Atmosphere
{
public:
  MinMaxModel();
  MinMaxModel(const MinMaxModel& orig);
  virtual ~MinMaxModel() override = default;

  virtual void update() override = 0;

  void setMinMaxFactor(greal factor, bool computeFlag);
  greal getMinMaxFactor() { return minMaxFactor; } 

  void getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity);
  greal getPressureAtSurface() const;

protected:
  static const size_t minMaxSize = 3;     //!< \brief Number of levels for the min-max model data.
  
  void initializeGases();
  virtual void updateMinMaxFactor() = 0;
  virtual void updateAtmosphereState();

  // Special inputs
  greal& minMaxFactor = atmos.minMaxFactor;  //!< \brief Min-max scale factor between min (-1), average (0), and max (1)
  greal userMinMaxFactor = 0;                //!< \brief Min-max scale factor between min (-1), average (0), and max (1)
  bool computeMinMaxFactor = true;           //!< \brief Set to false to suppress the calculation of a #minMaxFactor.

  virtual const std::vector<greal>& mdHeight() const = 0;                                     //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdTemperature() const = 0;        //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdPressure() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDensity() const = 0;            //!< \brief Implement this access method in the subclass.

  virtual const std::array<std::vector<greal>, minMaxSize>& mdArgonND() const = 0;            //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdHeliumND() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdHydrogenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDihydrogenND() const = 0;       //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrogenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDinitrogenND() const = 0;       //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdOxygenND() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDioxygenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdMethaneND() const = 0;          //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonMonoxideND() const = 0;   //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonDioxideND() const = 0;    //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdOzoneND() const = 0;            //!< \brief Implement this access method in the subclass.
  virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrousOxideND() const = 0;     //!< \brief Implement this access method in the subclass.

  //! \brief A map of all gas data for easy processing.
  std::map<GasType, std::reference_wrapper<const std::array<std::vector<greal>, minMaxSize>>> allGasData;

  //! \brief MinMax model data indices.
  enum MinMaxIndex {
    MIN_IDX = 0,   //!< The index of the minimum model.
    AVG_IDX = 1,   //!< The index of the average model.
    MAX_IDX = 2    //!< The index of the maximum model.
  };


  //virtual const std::vector<greal>& mdHeight() const override { return mHeight; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdTemperature() const override { return mTemp; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdPressure() const override { return mPres; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdDensity() const override { return mDens; }

  //virtual const std::array<std::vector<greal>, minMaxSize>& mdArgonND() const override { return mArgonND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdHeliumND() const override { return mHeliumND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdHydrogenND() const override { return mHydrogenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdDihydrogenND() const override { return mDihydrogenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrogenND() const override { return mNitrogenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdDinitrogenND() const override { return mDinitrogenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdOxygenND() const override { return mOxygenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdDioxygenND() const override { return mDioxygenND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdMethaneND() const override { return mMethaneND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonMonoxideND() const override { return mCarbonMonoxideND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonDioxideND() const override { return mCarbonDioxideND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdOzoneND() const override { return mOzoneND; }
  //virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrousOxideND() const override { return mNitrousOxideND; }

  //static std::vector<greal> mHeight;                                   //!< \brief Model data. Lookup key is height in km.
  //static std::array<std::vector<greal>, minMaxSize> mTemp;             //!< \brief Model data. Temperature in degrees K.
  //static std::array<std::vector<greal>, minMaxSize> mPres;             //!< \brief Model data. Pressure in \f$N/m^2\f$.
  //static std::array<std::vector<greal>, minMaxSize> mDens;             //!< \brief Model data. Density in \f$kg/m^3\f$.

  //static std::array<std::vector<greal>, minMaxSize> mArgonND;          //!< \brief Model data. Argon number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mHeliumND;         //!< \brief Model data. Helium number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mHydrogenND;       //!< \brief Model data. Hydrogen atom number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mDihydrogenND;     //!< \brief Model data. Hydrogen diatom number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mNitrogenND;       //!< \brief Model data. Nitrogen atom number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mDinitrogenND;     //!< \brief Model data. Nitrogen diatom number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mOxygenND;         //!< \brief Model data. Oxygen number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mDioxygenND;       //!< \brief Model data. Oxygen diatom number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mMethaneND;        //!< \brief Model data. Methane number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mCarbonMonoxideND; //!< \brief Model data. Carbon Monoxide number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mCarbonDioxideND;  //!< \brief Model data. Carbon Dioxide number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mOzoneND;          //!< \brief Model data. Ozone number density in \f$\#/m^3\f$.
  //static std::array<std::vector<greal>, minMaxSize> mNitrousOxideND;   //!< \brief Model data. Nitrous oxide number density in \f$\#/m^3\f$.

};

constexpr bool COMPUTE_MIN_MAX_FACTOR = true;
constexpr bool DO_NOT_COMPUTE_MIN_MAX_FACTOR = false;


} // namespace

