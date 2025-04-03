#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "VenusTopography.h"

namespace GRAM {

TEST(VenusTopography, getTopographicHeight)
{
  // SETUP
  VenusTopography venusTopo;

  // GIVEN (INPUTS)
  greal lat = 89.5_deg;
  greal lon = 240.0_deg;

  // RUN
  greal surfaceHeight = venusTopo.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.291277 + 0.048, surfaceHeight);

  // GIVEN (INPUTS)
  lat = 89.5_deg;
  lon = 250.0_deg;

  // RUN
  surfaceHeight = venusTopo.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.33229 + 0.048, surfaceHeight);

  // GIVEN (INPUTS)
  lat = 89.5_deg;
  lon = 239.0_deg;

  // RUN
  surfaceHeight = venusTopo.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.286478 + 0.048, surfaceHeight);

  // GIVEN (INPUTS)
  lat = -89.5_deg;
  lon = 239.0_deg;

  // RUN
  surfaceHeight = venusTopo.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.722356 + 0.048, surfaceHeight);

  // TEAR-DOWN
}


}