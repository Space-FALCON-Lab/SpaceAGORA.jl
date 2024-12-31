#include "unittest.h"
#define GRAM_UNIT_TEST
#include "JupiterAtmosphere.h"

namespace GRAM {

	TEST(JupiterAtmosphere, getPerturbationFactors)
	{
		// SETUP
		JupiterAtmosphere a;
		greal plow, phigh;
		
		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[6] = {-1.0, 45.0, 100.0, 350.0, 515.0, 600.0};
					
		// RUN & EXPECT (OUTPUTS)
		greal anslow[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		greal anshigh[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

		for (int i = 0; i < 6; i++) {
			a.height = h[i];
			a.getPerturbationFactors(plow, phigh);
			EXPECT_DOUBLE_EQ(anslow[i], plow);
			EXPECT_DOUBLE_EQ(anshigh[i], phigh);
		}
		
	}

	TEST(JupiterAtmosphere, getScaleParameters)
	{
		// SETUP
		JupiterAtmosphere a;
		greal vls, hls;

		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[8] = { -1.0, 5.0, 80.0, 115.0, 230.0, 415.0, 700.0, 1000.0 };
		
		// RUN & EXPECT (OUTPUTS)
		greal anslow[8] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		greal anshigh[8] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
		
		for (int i = 0; i < 8; i++) {
			a.height = h[i];
			a.getScaleParameters(vls, hls);
			EXPECT_DOUBLE_EQ(anslow[i], vls);
			EXPECT_DOUBLE_EQ(anshigh[i], hls);
		}
	}

	TEST(JupiterAtmosphere, getWindDeviations)
	{
		// SETUP
		JupiterAtmosphere a;
		greal sigu, sigv, svd;

		// GIVEN (INPUTS)
		// Test for the 6 conditional statements in the function
		greal h[8] = { -1.0, 5.0, 80.0, 115.0, 230.0, 415.0, 700.0, 1000.0 };

		// RUN & EXPECT (OUTPUTS)
    // No wind deviations, but we'll leave all this here for when someone
    // make a wind model
		greal ansu[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		greal ansv[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

		for (int i = 0; i < 8; i++) {
			a.height = h[i];
			a.getWindDeviations(sigu, sigv, svd);
			EXPECT_NEAR(ansu[i], sigu, 0.00001);
			EXPECT_NEAR(ansv[i], sigv, 0.00001);
			EXPECT_EQ(svd, 0.0);
		}
	}
}