#include "unittest.h"
#define GRAM_UNIT_TEST
#include "TitanAtmosphere.h"

namespace GRAM {

	TEST(TitanAtmosphere, getPerturbationFactors)
	{
		// SETUP
		TitanAtmosphere a;
		greal plow, phigh;

		// GIVEN (INPUTS)
		a.height = 45.0;
		
		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.8928571428571428, plow);
		EXPECT_DOUBLE_EQ(1.120000000000000, phigh);

		// GIVEN (INPUTS)
		a.height = 100.0;

		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.8866995073891626, plow);
		EXPECT_DOUBLE_EQ(1.127777777777778, phigh);

		// GIVEN (INPUTS)
		a.height = 350.0;

		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.9090909090909090, plow);
		EXPECT_DOUBLE_EQ(1.1000000000000000, phigh);

	}

	TEST(TitanAtmosphere, getScaleParameters)
	{
		// SETUP
		TitanAtmosphere a;
		greal vls, hls;

		// GIVEN (INPUTS)
		a.height = 230.0;

		// RUN & EXPECT (OUTPUTS)
		a.getScaleParameters(vls, hls);
		EXPECT_DOUBLE_EQ(86.00000000000000, vls);
		EXPECT_DOUBLE_EQ(430.0000000000000, hls);

		// GIVEN (INPUTS)
		a.height = 415.0;

		// RUN & EXPECT (OUTPUTS)
		a.getScaleParameters(vls, hls);
		EXPECT_DOUBLE_EQ(110.0000000000000, vls);
		EXPECT_DOUBLE_EQ(550.0000000000000, hls);

		// GIVEN (INPUTS)
		a.height = 750.0;

		// RUN & EXPECT (OUTPUTS)
		a.getScaleParameters(vls, hls);
		EXPECT_DOUBLE_EQ(120.0000000000000, vls);
		EXPECT_DOUBLE_EQ(600.0000000000000, hls);

	}

	TEST(TitanAtmosphere, getWindDeviations)
	{
		// SETUP
		TitanAtmosphere a;
		greal ewStdDev, nsStdDev, vertStdDev;

		// GIVEN (INPUTS)
		a.height = 30.0;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(11.20000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(11.20000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

		// GIVEN (INPUTS)
		a.height = 80.0;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(18.20000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(18.20000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

	}
}