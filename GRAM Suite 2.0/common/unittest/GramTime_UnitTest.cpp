#include "unittest.h"
#include "GramTime.h"

namespace GRAM {

TEST(GramTime, setStartTime)
{
  // SETUP
  GramTime time;

  // GIVEN (INPUTS)

  // RUN
  time.setStartTime(2035, 11, 26, 14, 14, 14.14, TDB, ERT);
  double et1 = time.getSpiceTime();
  time.setStartTime(2035, 11, 26, 14, 14, 14.14, TDT, ERT);
  double et2 = time.getSpiceTime();
  time.setStartTime(2035, 11, 26, 14, 14, 14.14, UTC, ERT);
  double et3 = time.getSpiceTime();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1132971254.14, et1);
  EXPECT_DOUBLE_EQ(1132971254.1389618, et2);
  EXPECT_DOUBLE_EQ(1132971323.3229617, et3);

  // TEAR-DOWN
}

TEST(GramTime, setStartTimeJD)
{
  // SETUP
  GramTime time;

  // GIVEN (INPUTS)

  // RUN
                   // time in seconds/secsperday + J2000inDays
  time.setStartTime(1132971254.14 / 86400.0 + 2451545.0, TDB, ERT);
  double et1 = time.getSpiceTime();
  time.setStartTime(1132971254.14 / 86400.0 + 2451545.0, TDT, ERT);
  double et2 = time.getSpiceTime();
  time.setStartTime(1132971254.14 / 86400.0 + 2451545.0, UTC, ERT);
  double et3 = time.getSpiceTime();

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(1132971254.14, et1, 1.0e-4);
  EXPECT_NEAR(1132971254.1389618, et2, 1.0e-4);
  EXPECT_NEAR(1132971323.3229617, et3, 1.0e-4);

  // TEAR-DOWN
}

TEST(GramTime, getTime)
{
  // SETUP
  GramTime time;

  int year, month, day, hour, minute;
  double seconds = 0;
  year = month = day = hour = minute = 0;

  // GIVEN (INPUTS)

  // RUN
  time.setSpiceTime(1132971254.14);
  time.getTime(TDB, PET, year, month, day, hour, minute, seconds);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(2035, year);
  EXPECT_EQ(11, month);
  EXPECT_EQ(26, day);
  EXPECT_EQ(14, hour);
  EXPECT_EQ(14, minute);
  EXPECT_DOUBLE_EQ(14.14, seconds);

  // RUN
  seconds = 0;
  year = month = day = hour = minute = 0;
  time.setSpiceTime(1132971254.1389618);
  time.getTime(TDT, PET, year, month, day, hour, minute, seconds);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(2035, year);
  EXPECT_EQ(11, month);
  EXPECT_EQ(26, day);
  EXPECT_EQ(14, hour);
  EXPECT_EQ(14, minute);
  EXPECT_DOUBLE_EQ(14.14, seconds);

  // RUN
  seconds = 0;
  year = month = day = hour = minute = 0;
  time.setSpiceTime(1132971323.3229617);
  time.getTime(UTC, PET, year, month, day, hour, minute, seconds);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(2035, year);
  EXPECT_EQ(11, month);
  EXPECT_EQ(26, day);
  EXPECT_EQ(14, hour);
  EXPECT_EQ(14, minute);
  EXPECT_DOUBLE_EQ(14.14, seconds);

  // TEAR-DOWN
}

TEST(GramTime, getTimeJD)
{
  // SETUP
  GramTime time;
  double julianDate = 0;

  // GIVEN (INPUTS)

  // RUN
  time.setSpiceTime(1132971254.14);
  time.getTime(TDB, PET, julianDate);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2464658.0932192, julianDate);

  // RUN
  julianDate = 0;  
  time.setSpiceTime(1132971254.1389618);
  time.getTime(TDT, PET, julianDate);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2464658.0932192, julianDate);

  // RUN
  julianDate = 0;
  time.setSpiceTime(1132971323.3229617);
  time.getTime(UTC, PET, julianDate);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2464658.0932192, julianDate);

  // TEAR-DOWN
}


} // namespace