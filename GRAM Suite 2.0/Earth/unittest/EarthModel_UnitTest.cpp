#include "unittest.h"
#include "EarthModel.h"
 
namespace GRAM {


//TEST(EarthModel, fair)
//{
//  // SETUP
//  EarthModel earth;
//  AtmosphereState a, b, c;
//
//  // GIVEN (INPUTS)
//  a.pressure = 0.098;
//  a.density = 1.9e-6;
//  a.temperature = 181.0;
//  a.ewWind = -3.4;
//  a.nsWind = 7.2;
//  a.verticalWind = 0.012;
//  b.pressure = 0.113;
//  b.density = 2.1e-6;
//  b.temperature = 185.0;
//  b.ewWind = -11.2;
//  b.nsWind = 5.7;
//  b.verticalWind = 0.004;
//
//  // RUN
//  earth.fair(1, 3.3, a, 5.4, b, 4.8, c);
//
//  // EXPECT (OUTPUTS)
//  EXPECT_DOUBLE_EQ(0.11001807948163447, c.pressure);
//  EXPECT_DOUBLE_EQ(2.0608038462714444e-06, c.density);
//  EXPECT_DOUBLE_EQ(184.24697960371748, c.temperature);
//  EXPECT_DOUBLE_EQ(-9.7316102272490586, c.ewWind);
//  EXPECT_DOUBLE_EQ(5.9823826486059506, c.nsWind);
//  EXPECT_DOUBLE_EQ(0.0055060407925650676, c.verticalWind);
//
//  // RUN
//  earth.fair(2, 3.3, a, 5.4, b, 4.8, c);
//
//  // EXPECT (OUTPUTS)
//  EXPECT_DOUBLE_EQ(0.1101761735139405, c.pressure);
//  EXPECT_DOUBLE_EQ(2.0623489801858733e-06, c.density);
//  EXPECT_DOUBLE_EQ(184.24697960371748, c.temperature);
//  EXPECT_DOUBLE_EQ(-9.7316102272490586, c.ewWind);
//  EXPECT_DOUBLE_EQ(5.9823826486059506, c.nsWind);
//  EXPECT_DOUBLE_EQ(0.0055060407925650676, c.verticalWind);
//
//  // TEAR-DOWN
//}

//TEST(EarthModel, update)
//{
//  // SETUP
//  EarthModel earth;
//  Position pos;
//  EarthInputParameters params;
//
//  // GIVEN (INPUTS)
//  pos.height = 29.0;
//  pos.latitude = 0.0;
//  earth.setPosition(pos);
//  earth.setInputParameters(params);
//
//  // RUN
//  earth.update();
//
//  // EXPECT (OUTPUTS)
//  EXPECT_DOUBLE_EQ(17.489951040781552, earth.windSpeed);
//  EXPECT_DOUBLE_EQ(9.2288449956827368, earth.windSpeedStandardDeviation);
//
//  // TEAR-DOWN
//}

} // namespace