#include "unittest.h"
#include "JB2008.h"
 
namespace GRAM {


TEST(JB2008, update)
{
  // SETUP
  JB2008 jb;
  EarthInputParameters params;
  Position pos;

  // GIVEN (INPUTS)
  params.minute = 0;
  params.seconds = 0.0;
  params.hour = 0;
  params.month = 1;
  params.day = 1;
  params.year = 2019;
  params.dailyF10 = 230.0;
  params.meanF10 = 230.4;
  params.ap = 16.7;
  params.dailyS10 = 2.3;
  params.meanS10 = 1.4;
  params.dailyXM10 = 3.4;
  params.meanXM10 = 3.1;
  params.dailyY10 = 4.4;
  params.meanY10 = 4.2;
  params.dstdtc = 3.7;

  jb.setInputParameters(params);
  pos.height = 95.6;
  jb.setPosition(pos);
  jb.updatePosition();

  // RUN
  jb.update();
  const AtmosphereState& atmos = jb.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(186.69454413641296, atmos.temperature);
  EXPECT_DOUBLE_EQ(0.066989062694351648, atmos.pressure);
  EXPECT_DOUBLE_EQ(1.2308545057464622e-06, atmos.density);
  EXPECT_DOUBLE_EQ(6.1116923880153751, atmos.ewWind);
  EXPECT_DOUBLE_EQ(14.571726052687104, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-8.2096423860497517e-06, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(28.521254910080312, atmos.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(1.9992586502928781e+19, atmos.dinitrogen.numberDensity);

  // TEAR-DOWN
}

}