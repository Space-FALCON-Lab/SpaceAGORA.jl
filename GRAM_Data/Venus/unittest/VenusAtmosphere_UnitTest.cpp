#include "gtest/gtest_prod.h"
#include "unittest.h"
#define GRAM_UNIT_TEST
#include "VenusAtmosphere.h"

namespace GRAM {

	TEST(VenusAtmosphere, getPerturbationFactors)
	{
		// SETUP
		VenusAtmosphere a;
		greal plow, phigh;

		// GIVEN (INPUTS)
		a.height = 45.0;
		
		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.9803921568627451, plow);
		EXPECT_DOUBLE_EQ(1.020000000000000, phigh);

		// GIVEN (INPUTS)
		a.height = 100.0;

		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.9216589861751152, plow);
		EXPECT_DOUBLE_EQ(1.085000000000000, phigh);

		// GIVEN (INPUTS)
		a.height = 350.0;

		// RUN & EXPECT (OUTPUTS)
		a.getPerturbationFactors(plow, phigh);
		EXPECT_DOUBLE_EQ(0.8695652173913044, plow);
		EXPECT_DOUBLE_EQ(1.150000000000000, phigh);

	}

	TEST(VenusAtmosphere, getScaleParameters)
	{
		// SETUP
		VenusAtmosphere a;
		greal vls, hls;

		// GIVEN (INPUTS)
		a.height = 50.0;

		// RUN & EXPECT (OUTPUTS)
		a.getScaleParameters(vls, hls);
		EXPECT_DOUBLE_EQ(15.07692307692308, vls);
		EXPECT_DOUBLE_EQ(904.6153846153846, hls);

		// GIVEN (INPUTS)
		a.height = 235.50;

		// RUN & EXPECT (OUTPUTS)
		a.getScaleParameters(vls, hls);
		EXPECT_DOUBLE_EQ(61.30000000000000, vls);
		EXPECT_DOUBLE_EQ(3678.000000000000, hls);
	}

	TEST(VenusAtmosphere, getWindDeviations)
	{
		// SETUP
		VenusAtmosphere a;
		greal ewStdDev, nsStdDev, vertStdDev;

		// GIVEN (INPUTS)
		a.height = 15.0;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(4.500000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(2.250000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

		// GIVEN (INPUTS)
		a.height = 30.0;
		
		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(8.000000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(4.000000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

		// GIVEN (INPUTS)
		a.height = 55.50;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(12.90000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(6.450000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

		// GIVEN (INPUTS)
		a.height = 125.0;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(21.50000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(10.75000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

		// GIVEN (INPUTS)
		a.height = 180.0;

		// RUN & EXPECT (OUTPUTS)
		a.getWindDeviations(ewStdDev, nsStdDev, vertStdDev);
		EXPECT_DOUBLE_EQ(25.00000000000000, ewStdDev);
		EXPECT_DOUBLE_EQ(12.50000000000000, nsStdDev);
		EXPECT_DOUBLE_EQ(0.00000000000000, vertStdDev);

	}
}