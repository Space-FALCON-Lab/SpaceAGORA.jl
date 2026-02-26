#include "unittest.h"
#define GRAM_UNIT_TEST
#include "gram.h"
#include "Yelle.h"
#include "Position.h"
#include "AtmosphereState.h"

namespace GRAM {

TEST( Yelle, setMinMaxFactor)
{
  // SETUP
  Yelle yelle;

  // GIVEN (INPUTS)
  greal mdf = 0.1234;

  // RUN
  yelle.setMinMaxFactor(mdf, false);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(mdf, yelle.getMinMaxFactor());

  // GIVEN (INPUTS)
  mdf = 2.1234; // Too Large

  // RUN
  yelle.setMinMaxFactor(mdf, false);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, yelle.getMinMaxFactor());

  // GIVEN (INPUTS)
  mdf = -2.1234; // Too Small

  // RUN
  yelle.setMinMaxFactor(mdf, false);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.0, yelle.getMinMaxFactor());

  // TEAR-DOWN
}

TEST( Yelle, update )
{
  // SETUP
  Position p;
  AtmosphereState s;
  Yelle yelle;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 35.0_deg;
  p.height = 100.0;
  p.elapsedTime = 0;
  yelle.setPosition(p);
  yelle.setMinMaxFactor(0.5, false);

  // RUN
  yelle.update();
  s = yelle.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(149.59316201713162, s.temperature);
  EXPECT_DOUBLE_EQ(1130.1318288367959, s.pressure);
  EXPECT_DOUBLE_EQ(34.923656179431475, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(0.025679899333543012, s.density);
  EXPECT_NEAR(30.482946469296145, s.densityScaleHeight, 1e-12);
  EXPECT_DOUBLE_EQ(2.4449358314033463e+022, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(28.170967048855662, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(9.4692072462271287e+021, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(5.0270017762578639e+023, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(5.3661874318604700e+023, s.totalNumberDensity);

  // TEAR-DOWN
}

TEST( Yelle, getReferenceValues )
{
  // SETUP
  greal temperature, pressure, density;
  Yelle yelle;

  // GIVEN (INPUTS)
  greal height = 101.27;

  // RUN
  yelle.getReferenceValues(height, temperature, pressure, density);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ( 142.89185000000001, temperature);
  EXPECT_DOUBLE_EQ( 907.97718960581221, pressure);
  EXPECT_DOUBLE_EQ( 2.1256587442126228E-002, density);

  // GIVEN (INPUTS)
  height = 2500.0;

  // RUN
  yelle.getReferenceValues(height, temperature, pressure, density);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(175.0, temperature);
  EXPECT_DOUBLE_EQ(7.0104799999999944e-10, pressure);
  EXPECT_DOUBLE_EQ(8.5904799999999904e-15, density);

  // TEAR-DOWN
}

TEST( Yelle, getPressureAtSurface )
{
  // SETUP
  Yelle yelle;

  // GIVEN (INPUTS)
  yelle.setMinMaxFactor(0.234, false);

  // RUN
  greal pas = yelle.getPressureAtSurface();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(149769.48642899672, pas);

  // GIVEN (INPUTS)
  yelle.setMinMaxFactor(-0.678, false);

  // RUN
  pas = yelle.getPressureAtSurface();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(138506.99084249517, pas);

  // TEAR-DOWN
}

TEST(Yelle, updateWinds)
{
	// SETUP
	Yelle yelle;

	// GIVEN (INPUTS)
	yelle.height = 100.0_km;
	yelle.latitude = 22.0_deg;

	// RUN
	yelle.updateWinds();

	// EXPECT (OUTPUTS)
	EXPECT_DOUBLE_EQ(58.19215607176653, yelle.ewWind);
	EXPECT_DOUBLE_EQ(0.0000000000, yelle.nsWind);

	// TEAR-DOWN
}

TEST(Yelle, updateMinMaxFactor)
{
	// SETUP
	Yelle yelle;

	// GIVEN (INPUTS)
	yelle.latitude = 22.0_deg;
	yelle.solarTime = 2.26;
	yelle.longitudeSun = 261.92_deg;
	yelle.setMinMaxFactor(0.0, true);
	
  // RUN
  yelle.updateMinMaxFactor();

  // EXPECT (OUTPUTS)
	EXPECT_DOUBLE_EQ(-0.3366108412480984, yelle.minMaxFactor);

  // GIVEN (INPUTS)
  yelle.setMinMaxFactor(0.5, false);

  // RUN
  yelle.updateMinMaxFactor();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.5, yelle.minMaxFactor);

	// GIVEN (INPUTS)
	yelle.latitude = 22.3_deg;
	yelle.solarTime = 2.43;
	yelle.longitudeSun = 265.02_deg;
	yelle.setMinMaxFactor(0.0, true);

  // RUN
  yelle.updateMinMaxFactor();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.3347660453753690, yelle.minMaxFactor);
}

} // namespace
