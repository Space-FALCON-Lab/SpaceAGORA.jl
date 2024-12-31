#include "gtest/gtest_prod.h"
#include "unittest.h"
#define GRAM_UNIT_TEST
#include "NeptuneAtmosphere.h"

namespace GRAM {

	TEST(NeptuneAtmosphere, getPerturbationFactors)
	{
		// SETUP
		NeptuneAtmosphere a;
		greal plow, phigh;
		
		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[6] = {-1.0, 45.0, 100.0, 350.0, 515.0, 600.0};
					
		// RUN & EXPECT (OUTPUTS)
		greal anslow[6] = { 0.9803921568627451, 0.9148084619782731,
			0.8730158730158731, 0.9090909090909091,
			0.9004952723998199, 0.8810572687224669 };
		greal anshigh[6] = { 1.020000000000000, 1.093125000000000,
			1.145454545454545, 1.100000000000000,
			1.110500000000000, 1.135000000000000 };

		for (int i = 0; i < 6; i++) {
			a.height = h[i];
			a.getPerturbationFactors(plow, phigh);
			EXPECT_DOUBLE_EQ(anslow[i], plow);
			EXPECT_DOUBLE_EQ(anshigh[i], phigh);
		}
		
	}

	TEST(NeptuneAtmosphere, getScaleParameters)
	{
		// SETUP
		NeptuneAtmosphere a;
		greal vls, hls;

		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[8] = { -1.0, 5.0, 80.0, 115.0, 230.0, 415.0, 700.0, 1000.0 };
		
		// RUN & EXPECT (OUTPUTS)
		greal anslow[8] = {90.44705882352940, 82.00000000000000,
			34.00000000000000, 47.80000000000000,
			93.33333333333333, 110.0000000000000,
			141.3953488372093, 200.0000000000000};
		greal anshigh[8] = { 452.2352941176470, 410.0000000000000,
			170.0000000000000, 239.0000000000000,
			466.6666666666666, 550.0000000000000,
			706.9767441860465, 1000.000000000000 };
		
		for (int i = 0; i < 8; i++) {
			a.height = h[i];
			a.getScaleParameters(vls, hls);
			EXPECT_DOUBLE_EQ(anslow[i], vls);
			EXPECT_DOUBLE_EQ(anshigh[i], hls);
		}
	}

	TEST(NeptuneAtmosphere, getWindDeviations)
	{
		// SETUP
		NeptuneAtmosphere a;
		greal sigu, sigv, svd;

		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[8] = { -1.0, 5.0, 80.0, 115.0, 230.0, 415.0, 700.0, 1000.0 };

		// RUN & EXPECT (OUTPUTS)
		greal ansu[8] = { 25.00000000000000, 25.00000000000000, 41.59999942779541,
			52.79999917745590, 89.59999835491180, 112.0000000000000,
			112.0000000000000, 112.0000000000000};
		greal ansv[8] = { 25.00000000000000, 25.00000000000000,
			25.00000000000000, 26.39999958872795,
			44.79999917745590, 56.00000000000000,
			56.00000000000000, 56.00000000000000 };

		for (int i = 0; i < 8; i++) {
			a.height = h[i];
			a.getWindDeviations(sigu, sigv, svd);
			EXPECT_NEAR(ansu[i], sigu, 0.00001);
			EXPECT_NEAR(ansv[i], sigv, 0.00001);
			EXPECT_EQ(svd, 0.0);
		}
	}
}