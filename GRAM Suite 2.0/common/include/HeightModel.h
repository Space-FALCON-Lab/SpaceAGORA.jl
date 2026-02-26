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
#include <vector>
#include <string>
#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief A height based atmosphere model base class.
//!
//! The HeightModel base class provides a quick way to create an atmospheric
//! model for height based data. The data for the model should be one-dimensional
//! and height dependent.  It should include temperature, pressure, and density
//! values as well as number densities for the appropriate gases.
//!
//! No data members are declared in the HeightModel base class so that they may be
//! declared statically in each subclass.  The class definition contains commented lines
//! to be copied into any subclass.  The subclass should also call initializeGases() in
//! its constructor. This class does not implement a winds model.  It is up to each subclass
//! to implement a winds model.
//! \ingroup CommonGRAM
class HeightModel : public Atmosphere
{
public:
  HeightModel();
  HeightModel(const HeightModel& orig) = default;
  virtual ~HeightModel() override = default;

  virtual void update() override = 0;

  virtual void getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity);
  virtual greal getPressureAtSurface() const;

protected:
  void initializeGases();
  void updateAtmosphereState();

  virtual const std::vector<greal>& mdHeight() const = 0;             //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdTemperature() const = 0;        //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdPressure() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdDensity() const = 0;            //!< \brief Implement this access method in the subclass.

  virtual const std::vector<greal>& mdArgonND() const = 0;            //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdHeliumND() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdHydrogenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdDihydrogenND() const = 0;       //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdNitrogenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdDinitrogenND() const = 0;       //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdOxygenND() const = 0;           //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdDioxygenND() const = 0;         //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdMethaneND() const = 0;          //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdCarbonMonoxideND() const = 0;   //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdCarbonDioxideND() const = 0;    //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdOzoneND() const = 0;            //!< \brief Implement this access method in the subclass.
  virtual const std::vector<greal>& mdNitrousOxideND() const = 0;     //!< \brief Implement this access method in the subclass.

  std::map<GasType, std::reference_wrapper<const std::vector<greal>>> allGasData;  //!< A map to all gases for easy processing.

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Copy the lines below to subclass this model.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //virtual const std::vector<greal>& mdHeight() const override { return mHeight; }
  //virtual const std::vector<greal>& mdTemperature() const override { return mTemp; }
  //virtual const std::vector<greal>& mdPressure() const override { return mPres; }
  //virtual const std::vector<greal>& mdDensity() const override { return mDens; }

  //virtual const std::vector<greal>& mdArgonND() const override { return mArgonND; }
  //virtual const std::vector<greal>& mdHeliumND() const override { return mHeliumND; }
  //virtual const std::vector<greal>& mdHydrogenND() const override { return mHydrogenND; }
  //virtual const std::vector<greal>& mdDihydrogenND() const override { return mDihydrogenND; }
  //virtual const std::vector<greal>& mdNitrogenND() const override { return mNitrogenND; }
  //virtual const std::vector<greal>& mdDinitrogenND() const override { return mDinitrogenND; }
  //virtual const std::vector<greal>& mdOxygenND() const override { return mOxygenND; }
  //virtual const std::vector<greal>& mdDioxygenND() const override { return mDioxygenND; }
  //virtual const std::vector<greal>& mdMethaneND() const override { return mMethaneND; }
  //virtual const std::vector<greal>& mdCarbonMonoxideND() const override { return mCarbonMonoxideND; }
  //virtual const std::vector<greal>& mdCarbonDioxideND() const override { return mCarbonDioxideND; }
  //virtual const std::vector<greal>& mdOzoneND() const override { return mOzoneND; }
  //virtual const std::vector<greal>& mdNitrousOxideND() const override { return mNitrousOxideND; }

  //static std::vector<greal> mHeight;           //!< \brief Model data. Lookup key is height in km.
  //static std::vector<greal> mTemp;             //!< \brief Model data. Temperature in degrees K.
  //static std::vector<greal> mPres;             //!< \brief Model data. Pressure in \f$N/m^2\f$.
  //static std::vector<greal> mDens;             //!< \brief Model data. Density in \f$kg/m^3\f$.

  //static std::vector<greal> mArgonND;          //!< \brief Model data. Argon number density in \f$\#/m^3\f$.
  //static std::vector<greal> mHeliumND;         //!< \brief Model data. Helium number density in \f$\#/m^3\f$.
  //static std::vector<greal> mHydrogenND;       //!< \brief Model data. Hydrogen atom number density in \f$\#/m^3\f$.
  //static std::vector<greal> mDihydrogenND;     //!< \brief Model data. Hydrogen diatom number density in \f$\#/m^3\f$.
  //static std::vector<greal> mNitrogenND;       //!< \brief Model data. Nitrogen atom number density in \f$\#/m^3\f$.
  //static std::vector<greal> mDinitrogenND;     //!< \brief Model data. Nitrogen diatom number density in \f$\#/m^3\f$.
  //static std::vector<greal> mOxygenND;         //!< \brief Model data. Oxygen number density in \f$\#/m^3\f$.
  //static std::vector<greal> mDioxygenND;       //!< \brief Model data. Oxygen diatom number density in \f$\#/m^3\f$.
  //static std::vector<greal> mMethaneND;        //!< \brief Model data. Methane number density in \f$\#/m^3\f$.
  //static std::vector<greal> mCarbonMonoxideND; //!< \brief Model data. Carbon Monoxide number density in \f$\#/m^3\f$.
  //static std::vector<greal> mCarbonDioxideND;  //!< \brief Model data. Carbon Dioxide number density in \f$\#/m^3\f$.
  //static std::vector<greal> mOzoneND;          //!< \brief Model data. Ozone number density in \f$\#/m^3\f$.
  //static std::vector<greal> mNitrousOxideND;   //!< \brief Model data. Nitrous oxide number density in \f$\#/m^3\f$.

};

} // namespace

