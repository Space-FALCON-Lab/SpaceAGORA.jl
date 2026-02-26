#include "gtest/gtest.h"
#include "gram.h"
#define GRAM_UNIT_TEST
#include "NeptuneMinMax.h"
#include "Position.h"
#include "AtmosphereState.h"

namespace GRAM {
	TEST(NeptuneMinMax, update)
	{
		// SETUP
		Position p;
		AtmosphereState s;
		NeptuneMinMax minMax;

		// GIVEN (INPUTS)
		p.longitude = 56.0_deg;
		p.latitude = 35.0_deg;
		p.height = 100.0;
		p.elapsedTime = 0;
		minMax.setPosition(p);
		minMax.setMinMaxFactor(0.5, false);

		// RUN
		minMax.update();
		s = minMax.getAtmosphereState();

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(111.74548013539997, s.temperature);
		EXPECT_DOUBLE_EQ(497.16757657755596, s.pressure);
		EXPECT_DOUBLE_EQ(31.585177569271174, s.pressureScaleHeight);
		EXPECT_DOUBLE_EQ(1.3938415061573526E-003, s.density);
		EXPECT_DOUBLE_EQ(21.937205648920138, s.densityScaleHeight);
		EXPECT_DOUBLE_EQ(6.0149989127181186E+022, s.helium.numberDensity);
		EXPECT_DOUBLE_EQ(2.6048923343470167, s.averageMolecularWeight);
		EXPECT_DOUBLE_EQ(4.9174045430491036E+021, s.methane.numberDensity);
		EXPECT_DOUBLE_EQ(2.5642594642508307E+023, s.dihydrogen.numberDensity);
		EXPECT_DOUBLE_EQ(3.2149334009531335E+023, s.totalNumberDensity);

		// TEAR-DOWN
	}

	TEST(NeptuneMinMax, getReferenceValues)
	{
		// SETUP
		greal temperature, pressure, density;
		NeptuneMinMax minMax;

		// GIVEN (INPUTS)
		greal height = 101.27;

		// RUN
		minMax.getReferenceValues(height, temperature, pressure, density);

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(106.51165999999999, temperature);
		EXPECT_DOUBLE_EQ(340.76230280004881, pressure);
		EXPECT_DOUBLE_EQ(1.0022576256282765E-003, density);

		// GIVEN (INPUTS)
		height = 2500.0;

		// RUN
		minMax.getReferenceValues(height, temperature, pressure, density);

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(684.62000000000000, temperature);
		EXPECT_DOUBLE_EQ(3.3952000000000001E-007, pressure);
		EXPECT_DOUBLE_EQ(1.2054000000000001E-013, density);

		// TEAR-DOWN
	}

	TEST(NeptuneMinMax, getPressureAtSurface)
	{
		// SETUP
		NeptuneMinMax minMax;

		// GIVEN (INPUTS)
		minMax.setMinMaxFactor(0.234, false);

		// RUN
		greal pas = minMax.getPressureAtSurface();

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(100010.0, pas);

		// GIVEN (INPUTS)
		minMax.setMinMaxFactor(-0.678, false);

		// RUN
		pas = minMax.getPressureAtSurface();

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(100003.21989084809, pas);

		// TEAR-DOWN
	}

	TEST(NeptuneMinMax, n2mixr)
	{
		// SETUP
		NeptuneMinMax minMax;

		// GIVEN (INPUTS)
		greal res = minMax.n2mixr(100.0, 0.006);

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(5.937392390672726E-03, res);

		// TEAR-DOWN
	}

	TEST(NeptuneMinMax, updateMoleFractions)
	{
		// SETUP
		NeptuneMinMax minMax;
		minMax.height = 100.0;
		//TEST 1: dinitrogen not present
		// GIVEN (INPUTS)
		minMax.fmolnitro = 0.0;
		minMax.setGasConstant(HELIUM, 4.04);
		minMax.setGasConstant(DIHYDROGEN, 2.02);
		minMax.setGasConstant(METHANE, 16.04);
		minMax.setGasConstant(DINITROGEN, 0.0);

		minMax.helium.numberDensity = 6.014998912718088E+22;
		minMax.dihydrogen.numberDensity = 2.564259464250832E+23;
		minMax.methane.numberDensity = 4.917404543049111E+21;
		minMax.updateTotalNumberDensity();

		// RUN
		minMax.updateMoleFractions();
		
		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(0.1870955992722841, minMax.helium.moleFraction);
		EXPECT_DOUBLE_EQ(0.7976088908997634, minMax.dihydrogen.moleFraction);
		EXPECT_DOUBLE_EQ(.01529550982795241, minMax.methane.moleFraction);
		EXPECT_DOUBLE_EQ(0.000000000000000E+00, minMax.dinitrogen.moleFraction);

		
		//TEST 2: dinitrogen is present
		// GIVEN (INPUTS)

		minMax.fmolnitro = 0.006;
		minMax.setGasConstant(HELIUM, 4.04);
		minMax.setGasConstant(DIHYDROGEN, 2.02);
		minMax.setGasConstant(METHANE, 16.04);
		minMax.setGasConstant(DINITROGEN, 28.0);

		minMax.helium.numberDensity = 6.014998912718088E+22;
		minMax.dihydrogen.numberDensity = 2.564259464250832E+23;
		minMax.methane.numberDensity = 4.917404543049111E+21;
		minMax.dinitrogen.numberDensity = 2.564259464250832E+20;
		minMax.updateTotalNumberDensity();

		// RUN
		minMax.updateMoleFractions();
		// EXPECT (OUTPUTS) 
		EXPECT_DOUBLE_EQ(-1.182438414223615, minMax.helium.moleFraction);
		EXPECT_DOUBLE_EQ(2.161217702116678, minMax.dihydrogen.moleFraction);
		EXPECT_DOUBLE_EQ(1.528331971626426E-02, minMax.methane.moleFraction);
		EXPECT_DOUBLE_EQ(5.937392390672726E-03, minMax.dinitrogen.moleFraction);

		// TEAR-DOWN
	}
	
	TEST(NeptuneMinMax, updateMinMaxFactor)
	{
		class NeptuneMinMaxTest : public NeptuneMinMax
		{
		public:
			using NeptuneMinMax::latitude;
			using NeptuneMinMax::solarTime;
			using NeptuneMinMax::longitudeSun;
			using NeptuneMinMax::computeMinMaxFactor;
			using NeptuneMinMax::minMaxFactor;
		};
		// SETUP
		NeptuneMinMaxTest minMax;

		// GIVEN (INPUTS)
		minMax.latitude = 22.0_deg;
		minMax.solarTime = 2.26;
		minMax.longitudeSun = 261.92_deg;
		minMax.setMinMaxFactor(0.0, true);

		// RUN & EXPECT (OUTPUTS)
		minMax.update();
		EXPECT_DOUBLE_EQ(-0.1517857326294690, minMax.minMaxFactor);
		minMax.setMinMaxFactor(0.5, false);
		minMax.update();
		EXPECT_DOUBLE_EQ(0.5, minMax.minMaxFactor);

		// GIVEN (INPUTS)
		minMax.latitude = 22.3_deg;
		minMax.solarTime = 2.43;
		minMax.longitudeSun = 265.02_deg;
		minMax.setMinMaxFactor(0.0, true);
		// RUN & EXPECT (OUTPUTS)
		minMax.update();
		EXPECT_DOUBLE_EQ(-0.1551962348347123, minMax.minMaxFactor);

		// TEAR-DOWN
	}

	TEST(NeptuneMinMax, updateWinds)
	{
		// SETUP
		NeptuneMinMax minMax;
				
		// GIVEN (INPUTS)
		minMax.height = 100.0_km;
		minMax.latitude = 22.0_deg;

		// RUN
		minMax.updateWinds();

		// EXPECT (OUTPUTS)
		EXPECT_DOUBLE_EQ(-244.0808383589788, minMax.ewWind);
		EXPECT_DOUBLE_EQ(0.0000000000, minMax.nsWind);
	}
} //namespace