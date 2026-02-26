//#include "gtest/gtest.h"
//#include "MarsEphemeris.h"
//
//using namespace GRAM;
//
//TEST(UranusEphemeris, setget )
//{
//  // SETUP
//  EphemerisState inState, outState;
//  UranusEphemeris ephem;
//  GramTime tTime;
//
//  // GIVEN (INPUTS)
//  inState.subsolarLongitude = 123.4_deg;
//  inState.subsolarLatitude  = 234.5_deg;
//  inState.longitudeSun      = 345.6_deg;
//  inState.orbitalRadius     = 456.7;
//  inState.oneWayLightTime   = 567.8;
//
//  // RUN
//  tTime.setStartTime(123.4, false, false);
//  ephem.setTime(tTime);
//  ephem.setLongitude(23.5_deg);
//  ephem.setEphemerisState(inState);
//  ephem.update();
//  outState = ephem.getEphemerisState();
//
//  // EXPECT (OUTPUTS)
//  EXPECT_DOUBLE_EQ(inState.subsolarLongitude, outState.subsolarLongitude);
//  EXPECT_DOUBLE_EQ(inState.subsolarLatitude , outState.subsolarLatitude );
//  EXPECT_DOUBLE_EQ(inState.longitudeSun     , outState.longitudeSun     );
//  EXPECT_DOUBLE_EQ(inState.orbitalRadius    , outState.orbitalRadius    );
//  EXPECT_DOUBLE_EQ(inState.oneWayLightTime  , outState.oneWayLightTime  );
//
//  // TEAR-DOWN
//}
//
//TEST( UranusEphemeris, update )
//{
//  // SETUP
//  EphemerisState outState;
//  UranusEphemeris ephem;
//  GramTime tTime;
//
//  // GIVEN (INPUTS)
//  tTime.setStartTime(1234.4 + 2451545.0, false, false);
//  ephem.setTime(tTime);
//  ephem.setLongitude(23.5_deg);
//
//  // RUN
//  ephem.update();
//  outState = ephem.getEphemerisState();
//
//  // EXPECT (OUTPUTS)
//  EXPECT_NEAR( 132.22592483667617_deg, outState.subsolarLongitude, 1.0e-14);
//  EXPECT_DOUBLE_EQ(-28.571949360582355_deg, outState.subsolarLatitude);
//  EXPECT_DOUBLE_EQ( 265.36207910874413_deg , outState.longitudeSun     );
//  EXPECT_DOUBLE_EQ( 30.081210875607056, outState.orbitalRadius);
//  EXPECT_DOUBLE_EQ( 247.80221945754047, outState.oneWayLightTime );
//  EXPECT_NEAR( 1.5994855992491352_deg, outState.equationOfTime, 1.0e-13);
//
//  //  sunlon =  132.22592483667617
//  //  sunlat = -28.571949360582355
//  //sunLsubs =  265.36207910874413
//  //     rpl = 30.081210875607056
//  //    owlt = 247.80221945754047
//  //     EOT = 1.5994855992491352
//
//  // TEAR-DOWN
//}
//
