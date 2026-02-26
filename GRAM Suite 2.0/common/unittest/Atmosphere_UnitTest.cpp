#include "unittest.h"
#define GRAM_UNIT_TEST
#include "Atmosphere.h"

namespace GRAM {

// Atmosphere is abstract, so make a test class.
//! \cond Hide_this_from_doxygen
class AtmosphereTest : public Atmosphere
{
public:
  AtmosphereTest() { initializeData(); }
  virtual ~AtmosphereTest() override = default;

  virtual void update() override {}
  void initializeData();
};

void AtmosphereTest::initializeData()
{
  if (helium.averageMolecularWeight > 0.0)
    return;
  setGasConstant(HELIUM, 3.4);
  setGasConstant(HYDROGEN, 2.5);
  helium.numberDensity = 123.4;
  hydrogen.numberDensity = 345.6;
  updateTotalNumberDensity();
}
//! \endcond

TEST(Atmosphere, copy)
{
  // SETUP
  AtmosphereTest source;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  source.setPlanetaryConstants(1.1, 2.2, 3.3, 4.4, 5.5, 6.6);
  pos.height = 123.4;
  source.setPosition(pos);
  ephem.longitudeSun = 234.5;
  source.setEphemerisState(ephem);
  source.setGasConstant(WATER, 7.7);

  // RUN
  AtmosphereTest copy(source);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.1, copy.getGravitationalParameter());
  EXPECT_DOUBLE_EQ(2.2, copy.getEquatorialRadius());
  EXPECT_DOUBLE_EQ(3.3, copy.getPolarRadius());
  EXPECT_DOUBLE_EQ(4.4, copy.getJ2());
  EXPECT_DOUBLE_EQ(5.5, copy.getPeriod());
  EXPECT_DOUBLE_EQ(123.4, copy.getPosition().height);
  EXPECT_DOUBLE_EQ(234.5, copy.getEphemerisState().longitudeSun);
  EXPECT_DOUBLE_EQ(7.7, copy.getAtmosphereState().water.averageMolecularWeight);
  // Make sure the map was rebuilt
  EXPECT_DOUBLE_EQ(7.7, copy.getConstituentGas(WATER).averageMolecularWeight);

  // TEAR-DOWN
}
TEST(Atmosphere, getRadius)
{
  // SETUP
  AtmosphereTest atmoTest;

  // GIVEN (INPUTS)
  atmoTest.setPlanetaryConstants(1.0, 5.5, 5.4, 1.0, 1.0, 1.0);

  // RUN
  greal radius = atmoTest.getRadius(32.3_deg);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(5.4708815986142998, radius);

  // TEAR-DOWN
}

TEST(Atmosphere, getGravity)
{
  // SETUP
  AtmosphereTest atmoTest;

  // GIVEN (INPUTS)
  atmoTest.setPlanetaryConstants(99.9, 5.5, 5.4, 2.3, 765.8, 1.4);

  // RUN
  greal gravity = atmoTest.getGravity(32.3_deg, 5.47, 0.43);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.28376397718448243, gravity);

  // TEAR-DOWN
}

TEST(Atmosphere, getTotalMass)
{
  // SETUP
  AtmosphereTest atmoTest;

  // GIVEN (INPUTS)

  // RUN
  greal totalMass = atmoTest.getTotalMass();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1283.56, totalMass);

  // TEAR-DOWN
}

TEST(Atmosphere, updateTotalNumberDensity)
{
  // SETUP
  AtmosphereTest atmoTest;
  atmoTest.totalNumberDensity = 0.0;

  // GIVEN (INPUTS)
  // see AtmosphereTest::initializeData()

  // RUN
  atmoTest.updateTotalNumberDensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(469.0, atmoTest.totalNumberDensity);

  // TEAR-DOWN
}

TEST(Atmosphere, updateMoleFractions)
{
  // SETUP
  AtmosphereTest atmoTest;

  // GIVEN (INPUTS)
  // see AtmosphereTest::initializeData()

  // RUN
  atmoTest.updateMoleFractions();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.2631130063965885, atmoTest.helium.moleFraction);
  EXPECT_DOUBLE_EQ(0.73688699360341147, atmoTest.hydrogen.moleFraction);
  EXPECT_DOUBLE_EQ(0.0, atmoTest.argon.moleFraction);

  // TEAR-DOWN
}

TEST(Atmosphere, updateMassFractions)
{
  // SETUP
  AtmosphereTest atmoTest;

  // GIVEN (INPUTS)
  // see AtmosphereTest::initializeData()

  // RUN
  atmoTest.updateMassFractions();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.326872136869332195, atmoTest.helium.massFraction);
  EXPECT_DOUBLE_EQ(0.67312786313066781, atmoTest.hydrogen.massFraction);
  EXPECT_DOUBLE_EQ(0.0, atmoTest.argon.massFraction);

  // TEAR-DOWN
}

} // namespace