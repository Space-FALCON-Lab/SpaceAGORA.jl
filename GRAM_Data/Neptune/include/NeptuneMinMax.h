//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "unittest_friend.h"
#include "gram.h"
#include "NeptuneCommon.h"
#include "MinMaxModel.h"

namespace GRAM {

struct NepData;
typedef struct NepData NeptuneMinMaxData;

//! \brief The min-max data model of the Neptune atmosphere.
//!
//! The NeptuneMinMax model utilizes the min-max data model.
//! \ingroup NeptuneGRAM
class NeptuneMinMax : public MinMaxModel, public NeptuneCommon
{
public:
  NeptuneMinMax();
  NeptuneMinMax(const NeptuneMinMax& orig);
  virtual ~NeptuneMinMax() override = default;

  void update() override;

  void setDinitrogenMoleFraction(greal n2mf);

protected:
  void updateWinds() override;

  virtual const std::vector<greal>& mdHeight() const override { return mHeight; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdTemperature() const override { return mTemp; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdPressure() const override { return mPres; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDensity() const override { return mDens; }

  virtual const std::array<std::vector<greal>, minMaxSize>& mdArgonND() const override { return mArgonND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdHeliumND() const override { return mHeliumND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdHydrogenND() const override { return mHydrogenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDihydrogenND() const override { return mDihydrogenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrogenND() const override { return mNitrogenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDinitrogenND() const override { return mDinitrogenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdOxygenND() const override { return mOxygenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdDioxygenND() const override { return mDioxygenND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdMethaneND() const override { return mMethaneND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonMonoxideND() const override { return mCarbonMonoxideND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdCarbonDioxideND() const override { return mCarbonDioxideND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdOzoneND() const override { return mOzoneND; }
  virtual const std::array<std::vector<greal>, minMaxSize>& mdNitrousOxideND() const override { return mNitrousOxideND; }

private:
  virtual void updateMoleFractions() override;
  virtual void updateMinMaxFactor();
  void initializeData();
  greal n2mixr(greal ht, greal fmolnitro);

  greal fmolnitro = 0;             //!< \brief User supplied Nitrogen mole fraction (between 0.0 and 0.006).
  
  static bool initialized;                             //!< \brief True if data has been initialized.
  static const int heightSize = 518;                   //!< \brief The number of model data records.
  static const NeptuneMinMaxData minData[heightSize];  //!< \brief Model data in table format.
  static const NeptuneMinMaxData avgData[heightSize];  //!< \brief Model data in table format.
  static const NeptuneMinMaxData maxData[heightSize];  //!< \brief Model data in table format.

  static std::vector<greal> mHeight;                                   //!< \brief Model data. Lookup key is height in km.
  static std::array<std::vector<greal>, minMaxSize> mTemp;             //!< \brief Model data. Temperature in degrees K.
  static std::array<std::vector<greal>, minMaxSize> mPres;             //!< \brief Model data. Pressure in \f$N/m^2\f$.
  static std::array<std::vector<greal>, minMaxSize> mDens;             //!< \brief Model data. Density in \f$kg/m^3\f$.

  static std::array<std::vector<greal>, minMaxSize> mArgonND;          //!< \brief Model data. Argon number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mHeliumND;         //!< \brief Model data. Helium number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mHydrogenND;       //!< \brief Model data. Hydrogen atom number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mDihydrogenND;     //!< \brief Model data. Hydrogen diatom number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mNitrogenND;       //!< \brief Model data. Nitrogen atom number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mDinitrogenND;     //!< \brief Model data. Nitrogen diatom number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mOxygenND;         //!< \brief Model data. Oxygen number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mDioxygenND;       //!< \brief Model data. Oxygen diatom number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mMethaneND;        //!< \brief Model data. Methane number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mCarbonMonoxideND; //!< \brief Model data. Carbon Monoxide number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mCarbonDioxideND;  //!< \brief Model data. Carbon Dioxide number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mOzoneND;          //!< \brief Model data. Ozone number density in \f$\#/m^3\f$.
  static std::array<std::vector<greal>, minMaxSize> mNitrousOxideND;   //!< \brief Model data. Nitrous oxide number density in \f$\#/m^3\f$.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(NeptuneMinMax, updateMoleFractions);
  FRIEND_TEST(NeptuneMinMax, n2mixr);
  FRIEND_TEST(NeptuneMinMax, updateWinds);
#endif // GRAM_UNIT_TEST
};

} // namespace

