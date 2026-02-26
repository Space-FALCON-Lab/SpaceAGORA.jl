#include "unittest.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

TEST(Interpolator, one_d)
{
  // SETUP
  Interpolator interp;

  // GIVEN (INPUTS)
  interp.makeFraction(21.3, 36.9, 27.54); // fraction is 0.4
  greal inputs[2] = { 23.3, 54.8 };

  // RUN
  greal linear_a = interp.linear(inputs[0], inputs[1]);
  greal linear_b = interp.linear(inputs);
  greal linear_c = interp.linear(inputs[1], inputs[0]);
  greal log_a = interp.log(inputs[0], inputs[1]);
  greal log_b = interp.log(inputs);
  greal log_c = interp.log(inputs[1], inputs[0]);
  greal inverse_a = interp.inverse(inputs[0], inputs[1]);
  greal inverse_c = interp.inverse(inputs[1], inputs[0]);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(35.9, linear_a);
  EXPECT_DOUBLE_EQ(35.9, linear_b);
  EXPECT_DOUBLE_EQ(42.2, linear_c);
  EXPECT_DOUBLE_EQ(32.80392241044440, log_a);
  EXPECT_DOUBLE_EQ(32.80392241044440, log_b);
  EXPECT_DOUBLE_EQ(38.92339409976984, log_c);
  EXPECT_DOUBLE_EQ(30.25687203791470, inverse_a);
  EXPECT_DOUBLE_EQ(35.56657381615600, inverse_c);

  // TEAR-DOWN
}

TEST(Interpolator, two_d)
{
  // SETUP
  Interpolator interp;

  // GIVEN (INPUTS)
  vector<greal> fracs = { 0.4, 0.7 };
  interp.setFraction(fracs);
  greal inputs[2][2] =
  { { 23.3, 54.8 },
  { 45.7, 34.9 } };

  // RUN
  greal linear_a = interp.linear(inputs);
  greal log_a = interp.log(inputs);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(42.466, linear_a);
  EXPECT_DOUBLE_EQ(40.51286566027470, log_a);

  // TEAR-DOWN
}

TEST(Interpolator, two_d_inline)
{
  // SETUP

  // GIVEN (INPUTS)
  Interpolator interp(0.4, 0.7);

  // RUN
  greal linear_a = interp.linear(23.3, 54.8, 
                                 45.7, 34.9);
  greal log_a = interp.log(23.3, 54.8, 
                           45.7, 34.9);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(42.466, linear_a);
  EXPECT_DOUBLE_EQ(40.51286566027470, log_a);

  // TEAR-DOWN
}

TEST(Interpolator, three_d)
{
  // SETUP

  // GIVEN (INPUTS)
  Interpolator interp(0.4, 0.7, 0.3);
  greal inputs[2][2][2] =
  { { { 23.3, 54.8 },
      { 45.7, 34.9 } },
    { { 14.7, 31.4 },
      { 45.7, 67.8 } } };

  // RUN
  greal linear_a = interp.linear(inputs);
  greal log_a = interp.log(inputs);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(40.7458, linear_a);
  EXPECT_DOUBLE_EQ(37.99297971645337, log_a);

  // TEAR-DOWN
}

TEST(Interpolator, four_d)
{
  // SETUP

  // GIVEN (INPUTS)
  Interpolator interp(0.4, 0.7, 0.3, 0.9);
  greal inputs[2][2][2][2] =
  { { { { 23.3, 54.8 },
        { 45.7, 34.9 } },
      { { 14.7, 31.4 },
        { 45.7, 67.8 } } },
    { { { 16.5, 18.6 },
        { 54.7, 54.9 } },
      { { 24.7, 35.4 },
        { 49.0, 76.8 } } } };

  // RUN
  greal linear_a = interp.linear(inputs);
  greal log_a = interp.log(inputs);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(41.92258, linear_a);
  EXPECT_DOUBLE_EQ(38.320158977800297, log_a);

  // TEAR-DOWN
}

} // namespace