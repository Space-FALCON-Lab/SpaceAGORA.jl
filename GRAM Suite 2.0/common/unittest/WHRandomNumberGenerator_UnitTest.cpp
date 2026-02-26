#include <vector>
#include "unittest.h"
#include "WHRandomNumberGenerator.h"


using namespace GRAM;


TEST(WHRandomNumberGenerator, setSeed)
{
	// SETUP
  WHRandomNumberGenerator a;

	// GIVEN (INPUTS)
	int s[10] = { 2511, 43, 7811, 167, 29011, 551, 1211, 110, 4511, 10020 };

	// RUN & EXPECT (OUTPUTS)
	int x[10] = { 2511, 43, 7811, 167, 29011, 551, 1211, 110, 4511, 10020 };
	int y[10] = { 7594, 7396, 9984, 28724, 19544, 3851, 26450, 18920, 18217, 26248 };
	int z[10] = { 2348, 7310, 23981, 28390, 19544, 2701, 23932, 18700, 8795, 5312 };
	for (int i = 0; i < 10; i++)
	{
		a.setSeed(s[i]);
		EXPECT_EQ(x[i], a.ix);
		EXPECT_EQ(y[i], a.iy);
		EXPECT_EQ(z[i], a.iz);
	}
}


TEST(WHRandomNumberGenerator, getRandomNumbers)
{
	// SETUP

  WHRandomNumberGenerator a;
	
	// GIVEN (INPUTS)
  a.setSeed(456789);
	
  // RUN 
  std::vector<greal> rvec(10);
  a.getRandomNumbers(rvec);

  // EXPECT (OUTPUTS)
  greal ans[10] = { 0.8011662082163975, 0.63658413664637448,
    0.2584486197932383, 0.23823737693551705, 0.7325462034467809,
    0.57393381427350154, 0.42416062614455496, 0.29142755342489579,
    0.31260833256546139, 0.85738088690688929 };
  for (int i = 0; i < 10; i++)
	{
		EXPECT_DOUBLE_EQ(ans[i], rvec[i]);
	}
}
