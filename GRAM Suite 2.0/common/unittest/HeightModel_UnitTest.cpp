#include "unittest.h"
#include "HeightModel.h"

using namespace GRAM;

// HeightModel is abstract, so make a test class.
//! \cond Hide_this_from_doxygen
class HeightModelTest : public HeightModel
{
public:
  HeightModelTest() { initializeGases();  initializeData(); }
  virtual ~HeightModelTest() override = default;

  // Will use update() to test updateAtmosphereState()
  virtual void update() { updateAtmosphereState(); }
  void initializeData();
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

  std::vector<greal> mHeight;           //!< Model data. Lookup key is height in km.
  std::vector<greal> mTemp;             //!< Model data. Temperature in degrees K.
  std::vector<greal> mPres;             //!< Model data. Pressure in \f$N/m^2\f$.
  std::vector<greal> mDens;             //!< Model data. Density in \f$kg/m^3\f$.

  std::vector<greal> mArgonND;          //!< Model data. Argon number density in \f$#/m^3\f$.
  std::vector<greal> mHeliumND;         //!< Model data. Helium number density in \f$#/m^3\f$.
  std::vector<greal> mHydrogenND;       //!< Model data. Hydrogen atom number density in \f$#/m^3\f$.
  std::vector<greal> mDihydrogenND;     //!< Model data. Hydrogen diatom number density in \f$#/m^3\f$.
  std::vector<greal> mNitrogenND;       //!< Model data. Nitrogen atom number density in \f$#/m^3\f$.
  std::vector<greal> mDinitrogenND;     //!< Model data. Nitrogen diatom number density in \f$#/m^3\f$.
  std::vector<greal> mOxygenND;         //!< Model data. Oxygen number density in \f$#/m^3\f$.
  std::vector<greal> mDioxygenND;       //!< Model data. Oxygen diatom number density in \f$#/m^3\f$.
  std::vector<greal> mMethaneND;        //!< Model data. Methane number density in \f$#/m^3\f$.
  std::vector<greal> mCarbonMonoxideND; //!< Model data. Carbon Monoxide number density in \f$#/m^3\f$.
  std::vector<greal> mCarbonDioxideND;  //!< Model data. Carbon Dioxide number density in \f$#/m^3\f$.
  std::vector<greal> mOzoneND;          //!< \brief Model data. Ozone number density in \f$\#/m^3\f$.
  std::vector<greal> mNitrousOxideND;   //!< \brief Model data. Nitrous oxide number density in \f$\#/m^3\f$.
};

void HeightModelTest::initializeData()
{
  setGasConstant(HELIUM, 4.0);
  setGasConstant(HYDROGEN, 1.0);
  if (mHeight.size() > 0)
    return;
  mHeight = { (greal)-20.0, (greal)-10.0, (greal)0.0, (greal)10.0, (greal)20.0, };
  mPres = { (greal)11.1, (greal)22.2, (greal)33.3, (greal)44.4, (greal)55.5, };
  mDens = { (greal)11.1, (greal)22.2, (greal)33.3, (greal)44.4, (greal)55.5, };
  mTemp = { (greal)11.1, (greal)22.2, (greal)33.3, (greal)44.4, (greal)55.5, };
  mHeliumND = { (greal)11.1, (greal)22.2, (greal)33.3, (greal)44.4, (greal)55.5, };
  mHydrogenND = { (greal)11.1, (greal)22.2, (greal)33.3, (greal)44.4, (greal)55.5, };
}
//! \endcond

TEST(HeightModel, getPressureAtSurface)
{
  // SETUP
  HeightModelTest heightModel;

  // GIVEN (INPUTS)

  // RUN
  greal pas = heightModel.getPressureAtSurface();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(33.3, pas);

  // TEAR-DOWN
}

TEST(HeightModel, getReferenceValues)
{
  // SETUP
  HeightModelTest heightModel;
  greal t, p, d;

  // GIVEN (INPUTS)
  greal height = 15.0;

  // RUN
  heightModel.getReferenceValues(height, t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(49.95, t);
  EXPECT_DOUBLE_EQ(49.64070910049534, p);
  EXPECT_DOUBLE_EQ(49.0278608399954, d);

  // TEAR-DOWN
}

TEST(HeightModel, updateAtmosphereState)
{
  // SETUP
  HeightModelTest heightModel;

  // GIVEN (INPUTS)
  Position pos;
  pos.height = 15.0;
  heightModel.setPosition(pos);

  // RUN
  heightModel.update();
  AtmosphereState atmos = heightModel.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(49.95, atmos.temperature);
  EXPECT_DOUBLE_EQ(49.64070910049534, atmos.pressure);
  EXPECT_DOUBLE_EQ(49.0278608399954, atmos.density);
  EXPECT_DOUBLE_EQ(-44.8142011772455, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(-44.8142011772455, atmos.densityScaleHeight);

  EXPECT_DOUBLE_EQ(49.64070910049534, atmos.helium.numberDensity);
  EXPECT_DOUBLE_EQ(49.64070910049534, atmos.hydrogen.numberDensity);
  // All other gases should be 0, test argon
  EXPECT_DOUBLE_EQ(0.0, atmos.argon.numberDensity);
  EXPECT_DOUBLE_EQ(2.5, atmos.averageMolecularWeight);

  // TEAR-DOWN
}

