#include <vector>
#include "unittest.h"
#include "RandomNumberGenerator.h"


using namespace GRAM;


//TEST(RandomNumberGenerator, setSeed)
//{
//	// SETUP
//  RandomNumberGenerator a;
//
//	// GIVEN (INPUTS)
//	int s[10] = { 2511, 43, 7811, 167, 29011, 551, 1211, 110, 4511, 10020 };
//
//	// RUN & EXPECT (OUTPUTS)
//	int x[10] = { 2511, 43, 7811, 167, 29011, 551, 1211, 110, 4511, 10020 };
//	int y[10] = { 7594, 7396, 9984, 28724, 19544, 3851, 26450, 18920, 18217, 26248 };
//	int z[10] = { 2348, 7310, 23981, 28390, 19544, 2701, 23932, 18700, 8795, 5312 };
//	for (int i = 0; i < 10; i++)
//	{
//		a.setSeed(s[i]);
//		EXPECT_EQ(x[i], a.ix);
//		EXPECT_EQ(y[i], a.iy);
//		EXPECT_EQ(z[i], a.iz);
//	}
//}


TEST(RandomNumberGenerator, getRandomNumbers)
{
  // SETUP

  RandomNumberGenerator a;

  // GIVEN (INPUTS)
  a.setSeed(456789);

  // RUN 
  std::vector<greal> rvec(10);
  a.getRandomNumbers(rvec);

  // EXPECT (OUTPUTS)
  greal ans[10] = { 0.80561107397079468, 0.54926562309265137,
    0.29496359825134277, 0.83761930465698242, 0.42442446947097778,
    0.74129277467727661, 0.41806495189666748, 0.94419044256210327,
    0.66398674249649048, 0.25070881843566895 };
  for (int i = 0; i < 10; i++)
  {
    EXPECT_DOUBLE_EQ(ans[i], rvec[i]);
  }
}

TEST(RandomNumberGenerator, copyConstructor)
{
  // SETUP

  RandomNumberGenerator a;

  // GIVEN (INPUTS)
  a.setSeed(456789);
  std::vector<greal> avec(100);
  std::vector<greal> bvec(100);
  std::vector<greal> cvec(100);
  a.getRandomNumbers(avec);

  // RUN 
  RandomNumberGenerator b(a);
  a.getRandomNumbers(avec);
  b.getRandomNumbers(bvec);

  // EXPECT (OUTPUTS)
  for (int i = 0; i < 100; i++)
  {
    EXPECT_DOUBLE_EQ(avec[i], bvec[i]);
  }

  // RUN 
  RandomNumberGenerator c;
  c = b;
  c.getRandomNumbers(cvec);
  b.getRandomNumbers(bvec);

  // EXPECT (OUTPUTS)
  for (int i = 0; i < 100; i++)
  {
    EXPECT_DOUBLE_EQ(bvec[i], cvec[i]);
  }
}
