#include <string>
#include <fstream>
#include "unittest.h"
#include "RRA.h"
 
using namespace std;

namespace GRAM {

TEST(RRA, parseTopLine_and_readSiteList)
{
  // SETUP
  RRA rra;
  rra.setUseRRA(true);

  for (size_t index = 0; index < rra.siteNameList.size(); ++index) {
    const string& siteName = rra.siteNameList[index];
    rra.siteLatitude = rra.siteGeodeticLatitudeList[index];
    rra.siteLongitude = rra.siteLongitudeList[index];
    bool goodRead = false;
    try {
      // GIVEN (INPUTS)
      string fileName = rra.rraPath + "T3" + siteName + ".txt";
      ifstream rraDataFile;
      rraDataFile.open(fileName);

      // RUN
      rra.parseTopLine(rraDataFile, fileName);
      goodRead = true;
      rraDataFile.close();
    }
    catch (const std::string&) {
      goodRead = false;
    }

    // EXPECT (OUTPUTS)
    EXPECT_EQ(true, goodRead);
  }

  // TEAR-DOWN
}

TEST(RRA, locateNearestSite)
{
  // SETUP
  RRA rra;
  rra.setUseRRA(true);
  rra.setRRAParameters(RRA_2019, 0.0, 1.0);
  Position pos;

  for (size_t index = 0; index < rra.siteNameList.size(); ++index) {
    const string& siteName = rra.siteNameList[index];

    // GIVEN (INPUTS)
    pos.latitude = rra.siteGeocentricLatitudeList[index] + 0.01;
    pos.longitude = rra.siteLongitudeList[index] - 0.01;
    rra.setPosition(pos);

    // RUN
    rra.locateNearestSite();

    // EXPECT (OUTPUTS)
    EXPECT_EQ(siteName, rra.siteName);
  }

  // TEAR-DOWN
}

TEST(RRA, getWeightFactor)
{
  // SETUP
  RRA rra;

  // GIVEN (INPUTS)
  rra.setRRAParameters(RRA_2019, 4.0, 6.0);

  // RUN
  greal weight = rra.getWeightFactor(3.8);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, weight);

  // RUN
  weight = rra.getWeightFactor(6.8);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, weight);

  // RUN
  weight = rra.getWeightFactor(5.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.5, weight);

  // TEAR-DOWN
}

TEST(RRA, findHeightIndex)
{
  // SETUP
  RRA rra;

  // GIVEN (INPUTS)
  std::vector<greal> z = { 0.5, 1.7, 2.2, 5.0 };

  // RUN
  int index = rra.findHeightIndex(0.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0, index);

  // RUN
  index = rra.findHeightIndex(0.8, z);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0, index);

  // RUN
  index = rra.findHeightIndex(2.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1, index);

  // RUN
  index = rra.findHeightIndex(6.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(3, index);

  // TEAR-DOWN
}

TEST(RRA, getHeightWeight)
{
  // SETUP
  RRA rra;

  // GIVEN (INPUTS)
  std::vector<greal> z = { 0.5, 1.7, 3.0, 5.0 };

  // RUN
  greal weight = rra.getHeightWeight(1.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, weight);

  // RUN
  weight = rra.getHeightWeight(0.8, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, weight);

  // RUN
  weight = rra.getHeightWeight(4.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.5, weight);

  // RUN
  weight = rra.getHeightWeight(6.0, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, weight);

  // TEAR-DOWN
}

TEST(RRA, dedt)
{
  // SETUP
  RRA rra;

  // GIVEN (INPUTS)

  // RUN
  double d = rra.dedt(323.4, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(619.44145155896911, d);

  // RUN
  d = rra.dedt(1234.5, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(4110547939.630332, d);

  // RUN
  d = rra.dedt(0.0001, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(5.9831475013085843e-08, d);

  // TEAR-DOWN
}

TEST(RRA, d2edt2)
{
  // SETUP
  RRA rra;

  // GIVEN (INPUTS)

  // RUN
  double d = rra.d2edt2(323.4, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(26.298707614362534, d);

  // RUN
  d = rra.d2edt2(1234.5, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(93782804.262978882, d);

  // RUN
  d = rra.d2edt2(0.0001, 4.56);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(357980540.22235668, d);

  // TEAR-DOWN
}


}