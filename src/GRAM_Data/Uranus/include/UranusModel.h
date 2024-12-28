//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "UranusCommon.h"
#include "HeightModel.h"

namespace GRAM {

struct UraData;
typedef struct UraData UranusModelData;

//! \brief The height based model of a Uranus Atmosphere.
//!
//! The UranusModel utilizes the simple height based data model.
//! The data for the model was provided by Gary Allen (ARC).
//! \ingroup UranusGRAM
class UranusModel : public HeightModel, public UranusCommon
{
public:
  UranusModel();
  UranusModel(const UranusModel& orig);
  virtual ~UranusModel() override = default;

  void update() override;

protected:
  virtual const std::vector<greal>& mdHeight() const override { return mHeight; }
  virtual const std::vector<greal>& mdTemperature() const override { return mTemp; }
  virtual const std::vector<greal>& mdPressure() const override { return mPres; }
  virtual const std::vector<greal>& mdDensity() const override { return mDens; }

  virtual const std::vector<greal>& mdArgonND() const override { return mArgonND; }
  virtual const std::vector<greal>& mdHeliumND() const override { return mHeliumND; }
  virtual const std::vector<greal>& mdHydrogenND() const override { return mHydrogenND; }
  virtual const std::vector<greal>& mdDihydrogenND() const override { return mDihydrogenND; }
  virtual const std::vector<greal>& mdNitrogenND() const override { return mNitrogenND; }
  virtual const std::vector<greal>& mdDinitrogenND() const override { return mDinitrogenND; }
  virtual const std::vector<greal>& mdOxygenND() const override { return mOxygenND; }
  virtual const std::vector<greal>& mdDioxygenND() const override { return mDioxygenND; }
  virtual const std::vector<greal>& mdMethaneND() const override { return mMethaneND; }
  virtual const std::vector<greal>& mdCarbonMonoxideND() const override { return mCarbonMonoxideND; }
  virtual const std::vector<greal>& mdCarbonDioxideND() const override { return mCarbonDioxideND; }
  virtual const std::vector<greal>& mdOzoneND() const override { return mOzoneND; }
  virtual const std::vector<greal>& mdNitrousOxideND() const override { return mNitrousOxideND; }

private:
  void initializeData();
  void updateWinds();

  static bool initialized;                       //!< True if data has been initialized.
  static const size_t heightSize = 213;          //!< The number of model data records.
  static const UranusModelData data[heightSize]; //!< Model data in table format.

  // Table Data based on height levels.  (The m is for model data.)
  static std::vector<greal> mHeight;           //!< Model data. Lookup key is height \units{km}.
  static std::vector<greal> mTemp;             //!< Model data. Temperature \units{K}.
  static std::vector<greal> mPres;             //!< Model data. Pressure \units{Pa}.
  static std::vector<greal> mDens;             //!< Model data. Density \units{kg/m^3}.

  static std::vector<greal> mArgonND;          //!< Model data. Argon number density \units{\#/m^3}.
  static std::vector<greal> mHeliumND;         //!< Model data. Helium number density \units{\#/m^3}.
  static std::vector<greal> mHydrogenND;       //!< Model data. Hydrogen atom number density \units{\#/m^3}.
  static std::vector<greal> mDihydrogenND;     //!< Model data. Hydrogen diatom number density \units{\#/m^3}.
  static std::vector<greal> mNitrogenND;       //!< Model data. Nitrogen atom number density \units{\#/m^3}.
  static std::vector<greal> mDinitrogenND;     //!< Model data. Nitrogen diatom number density \units{\#/m^3}.
  static std::vector<greal> mOxygenND;         //!< Model data. Oxygen number density \units{\#/m^3}.
  static std::vector<greal> mDioxygenND;       //!< Model data. Oxygen diatom number density \units{\#/m^3}.
  static std::vector<greal> mMethaneND;        //!< Model data. Methane number density \units{\#/m^3}.
  static std::vector<greal> mCarbonMonoxideND; //!< Model data. Carbon Monoxide number density \units{\#/m^3}.
  static std::vector<greal> mCarbonDioxideND;  //!< Model data. Carbon Dioxide number density \units{\#/m^3}.
  static std::vector<greal> mOzoneND;          //!< Model data. Ozone number density \units{\#/m^3}.
  static std::vector<greal> mNitrousOxideND;   //!< Model data. Nitrous oxide number density \units{\#/m^3}.

};

} // namespace

