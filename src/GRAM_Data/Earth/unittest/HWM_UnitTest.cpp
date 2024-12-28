#include "unittest.h"
#include "HWM.h"
 
using namespace GRAM;



TEST(HWM, gws5)
{
	// SETUP
	HWM hwm;

  // GIVEN (INPUTS)
  int iyd = 1111;
  double sec = 12.3;
  double alt = 0;
  double glat = 0;
  double glon = 0;
  double stl = 0;
  double f107a = 0;
  double f107 = 0;
  double ap[7] = { 0 };

  // RUN
  double w[2];
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, w[0]);
  EXPECT_DOUBLE_EQ(-4.8208626164335362, w[1]);

  // GIVEN (INPUTS)
  iyd = 1111;
  sec = 12.3;
  alt = 34.5;
  glat = 56.7;
  glon = 123.4;
  stl = 1.5;
  f107a = 178.9;
  f107 = 189.7;
  ap[0] = 8.0;

  // RUN
  w[0] = w[1] = 0.0;
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.9281228494233549, w[0]);
  EXPECT_DOUBLE_EQ(3.3227688098391885, w[1]);

  // GIVEN (INPUTS)
  iyd = 1111;
  sec = 12.3;
  alt = 81.2;
  glat = 56.7;
  glon = 123.4;
  stl = 1.5;
  f107a = 178.9;
  f107 = 189.7;
  ap[0] = 8.0;

  // RUN
  w[0] = w[1] = 0.0;
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.475255102519201, w[0]);
  EXPECT_DOUBLE_EQ(-9.3681857602121887, w[1]);

  // GIVEN (INPUTS)
  iyd = 1111;
  sec = 12.3;
  alt = 100.0;
  glat = 56.7;
  glon = 123.4;
  stl = 1.5;
  f107a = 178.9;
  f107 = 189.7;
  ap[0] = 8.0;

  // RUN
  w[0] = w[1] = 0.0;
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.48919420023068483, w[0]);
  EXPECT_DOUBLE_EQ(-0.65712249785568133, w[1]);

  // GIVEN (INPUTS)
  iyd = 1111;
  sec = 12.3;
  alt = 123.4;
  glat = 56.7;
  glon = 123.4;
  stl = 1.5;
  f107a = 178.9;
  f107 = 189.7;
  ap[0] = 8.0;

  // RUN
  w[0] = w[1] = 0.0;
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(16.217092531373858, w[0]);
  EXPECT_DOUBLE_EQ(39.707986670905193, w[1]);

  // GIVEN (INPUTS)
  iyd = 1111;
  sec = 12.3;
  alt = 321.6;
  glat = 56.7;
  glon = 123.4;
  stl = 1.5;
  f107a = 178.9;
  f107 = 189.7;
  ap[0] = 8.0;

  // RUN
  w[0] = w[1] = 0.0;
  hwm.gws5(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-188.14004950504486, w[0]);
  EXPECT_DOUBLE_EQ(-50.650549758275744, w[1]);

  // TEAR-DOWN
}

//TEST(HWM, pdtuv)
//{
//	// SETUP
//	HWM hwm;
//
//	// GIVEN (INPUTS)
//
//	// RUN
//
//	// EXPECT (OUTPUTS)
//	EXPECT_DOUBLE_EQ(0.0, x);
//
//	// TEAR-DOWN
//}

