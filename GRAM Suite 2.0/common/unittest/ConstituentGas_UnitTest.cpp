#include <vector>
#include "unittest.h"
#include "ConstituentGas.h"

using namespace std;

namespace GRAM {

TEST(ConstituentGas, updateSpecificHeatCapacity_argon)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low     center  right    high   interp
  vector<greal> temperature = {   30.0,  150.0,  450.0, 1111.0, 432.1 };
  vector<greal> pressure =    {   0.04,  1.0e4,  3.3e7,  2.2e9, 6.7e5 };
  vector<greal> expected_cp = { 0.5203, 0.5210, 0.5637, 0.5280, 0.52372722};
  vector<greal> expected_cv = { 0.3122, 0.3124, 0.3187, 0.3140, 0.3127253466666667};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(ARGON, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_helium)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high    high   
  vector<greal> temperature = {   30.0, 1111.0,  150.0};
  vector<greal> pressure =    {   0.04,  2.2e9,  1.0e7};
  vector<greal> expected_cp = { 5.1932, 5.1881, 5.2435};
  vector<greal> expected_cv = { 3.1159, 3.1222, 3.1642};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(HELIUM, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_hydrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   
  vector<greal> temperature = {    30.0, 1111.0 };
  vector<greal> pressure =    {    0.04,  2.2e9 };
  vector<greal> expected_cp = { 20.6224, 20.6223 };
  vector<greal> expected_cv = { 12.3734, 12.3733 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(HYDROGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_dihydrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   random   
  vector<greal> temperature = {   30.0, 1111.0, 600.0};
  vector<greal> pressure =    {   0.04,  2.2e9, 100.0};
  vector<greal> expected_cp = { 10.313, 14.919, 14.549};
  vector<greal> expected_cv = { 6.1887, 10.81 , 10.424};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(DIHYDROGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_nitrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   
  vector<greal> temperature = {   30.0, 1111.0};
  vector<greal> pressure =    {   0.04,  2.2e9};
  vector<greal> expected_cp = {  1.484, 1.4849};
  vector<greal> expected_cv = { 0.8904, 0.8913};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(NITROGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_dinitrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   random
  vector<greal> temperature = {   30.0, 1111.0, 400.0 };
  vector<greal> pressure =    {   0.04,  2.2e9, 1.0e4 };
  vector<greal> expected_cp = { 1.0389, 1.1669, 1.0442 };
  vector<greal> expected_cv = {  0.742, 0.8632, 0.7473 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(DINITROGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_oxygen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high    interp   
  vector<greal> temperature = {   30.0, 1111.0,  474.4};
  vector<greal> pressure =    {   0.04,  2.2e9,  2.2e6};
  vector<greal> expected_cp = { 1.2995, 1.2995, 1.332072};
  vector<greal> expected_cv = { 0.7797, 0.7797, 0.812272};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(OXYGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_dioxygen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low    special     special   high   
  vector<greal> temperature = {   30.0,   50.0,       66.0, 1111.0};
  vector<greal> pressure =    {4.0e-12, 1.0e-7,     8.2e-7,  2.2e9};
  vector<greal> expected_cp = { 0.9106, 0.9103,  0.9096824, 1.0918};
  vector<greal> expected_cv = { 0.6508, 0.6505,  0.6498504, 0.8241};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(DIOXYGEN, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_methane)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   random   
  vector<greal> temperature = {   30.0, 1111.0,   40.0};
  vector<greal> pressure =    {   0.04,  2.2e9,  1.0e4};
  vector<greal> expected_cp = { 2.0736333333333334, 4.348, 2.0731};
  vector<greal> expected_cv = { 1.5553333333333335, 3.8297, 1.5548};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(METHANE, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_carbon_monoxide)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low     random   high   
  vector<greal> temperature = {   30.0,   40.0, 1111.0};
  vector<greal> pressure =    {4.0e-12,   0.01,  2.2e9};
  vector<greal> expected_cp = { 1.0389, 1.0389, 1.1736};
  vector<greal> expected_cv = { 0.7421, 0.7421, 0.8768};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(CARBON_MONOXIDE, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_carbon_dioxide)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   random   interp   interp
  vector<greal> temperature = {   50.0, 1111.0,  550.0,    123.0, 432.1 };
  vector<greal> pressure =    { 2.3e-5,  2.2e9,  1.0e5,    0.045, 6.7e5 };
  vector<greal> expected_cp = { 0.6612,  1.243, 1.0472, 0.674477, 0.9785919066666667};
  vector<greal> expected_cv = { 0.4723, 1.0347, 0.8576, 0.485577, 0.7792665666666667};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(CARBON_DIOXIDE, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_water)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low     random   high   special  special     
  vector<greal> temperature = {   30.0,  450.0, 1111.0,  150.0,  166.0};
  vector<greal> pressure =    {   0.04,  1.0e6,  2.2e9, 1.0e-7, 8.2e-7};
  vector<greal> expected_cp = { 1.8461, 4.3924, 2.4369, 1.8497, 1.8494472};
  vector<greal> expected_cv = { 1.3846, 3.4076, 1.8617, 1.3882, 1.3879472};

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(WATER, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_ozone)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   
  vector<greal> temperature = {   30.0, 1111.0 };
  vector<greal> pressure    = {   0.04,  2.2e9 };
  vector<greal> expected_cp = { 0.6929, 1.1384 };
  vector<greal> expected_cv = { 0.5197, 0.9652 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(OZONE, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, updateSpecificHeatCapacity_nitrous_oxide)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  //                            low      high   
  vector<greal> temperature = {   30.0, 1111.0 };
  vector<greal> pressure    = {   0.04,  2.2e9 };
  vector<greal> expected_cp = { 0.6612,  1.234 };
  vector<greal> expected_cv = { 0.4723, 1.0451 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    gas.updateSpecificHeatCapacity(NITROUS_OXIDE, temperature[i], pressure[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_DOUBLE_EQ(expected_cp[i], gas.specificHeatPressure);
    EXPECT_DOUBLE_EQ(expected_cv[i], gas.specificHeatVolume);
  }
  // TEAR DOWN
}

#ifdef OLDCP
TEST(ConstituentGas, getSpecificHeatCapacity_argon)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 500.0, 1000.0, 6000.0 };
  vector<greal> expected = { 20.79, 20.79, 20.79, 20.79 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(ARGON, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_helium)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 500.0, 1000.0, 6000.0 };
  vector<greal> expected = { 20.79, 20.79, 20.79, 20.79 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(HELIUM, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_hydrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 500.0, 1000.0, 6000.0 };
  vector<greal> expected = { 20.79, 20.79, 20.79, 20.79 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(HYDROGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 175.0, 225.0, 298.15, 500.0, 1000.0, 1100.0, 1700.0, 2500.0, 2600.0, 4900.0, 6000.0 };
  vector<greal> expected = { 26.45, 27.88, 28.84, 29.26, 30.20, 30.58, 33.14, 35.84, 36.11, 40.68, 41.97 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(DIHYDROGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_nitrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 500.0, 1000.0, 6000.0 };
  vector<greal> expected = { 20.74, 20.85, 20.79, 25.70 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(NITROGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_dinitrogen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 100.0, 300.0, 500.0, 600.0, 1100.0, 2000.0, 2100.0, 4900.0, 6000.0 };
  vector<greal> expected = { 29.10, 29.12, 29.58, 30.10, 33.24,  35.98,  36.13, 37.88, 38.27 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(DINITROGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_oxygen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 500.0, 1000.0, 6000.0 };
  vector<greal> expected = { 20.81, 20.81, 20.81, 20.81 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(OXYGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_dioxygen)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 100.0, 300.0, 700.0, 800.0, 1100.0, 2000.0, 2100.0, 4900.0, 6000.0 };
  vector<greal> expected = { 29.10, 29.39, 32.98, 33.74, 35.29,  37.75,  37.97,  42.55,  44.39 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(DIOXYGEN, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_methane)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 100.0, 200.0, 298.0, 300.0, 1100.0, 1300.0, 1400.0, 2700.0, 6000.0 };
  vector<greal> expected = { 33.28, 33.51, 35.69, 35.71, 75.54,  81.74,  84.34,  99.91,  106.41 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(METHANE, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_carbon_monoxide)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 175.0, 200.0, 298.0, 300.0, 800.0, 1300.0, 1400.0, 2700.0, 6000.0 };
  vector<greal> expected = { 29.10, 29.10, 29.10, 29.15, 31.88,  34.55,  34.93,  36.98,  38.37 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(CARBON_MONOXIDE, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_carbon_dioxide)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 175.0, 200.0, 275.0, 300.0, 800.0, 1200.0, 1300.0, 2700.0, 6000.0 };
  vector<greal> expected = { 31.20, 32.35, 36.05, 37.22, 51.44,  56.35,  57.14,  61.80,  64.98 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(CARBON_DIOXIDE, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}

TEST(ConstituentGas, getSpecificHeatCapacity_water)
{
  // SETUP
  ConstituentGas gas;

  // GIVEN (INPUTS)
  vector<greal> temperature = { 298.15, 400.0, 500.0, 600.0, 1100.0, 1700.0, 1800.0, 2700.0, 6000.0 };
  vector<greal> expected = { 75.37,  76.74, 83.66, 36.32, 42.52,  48.92,  49.75,  54.71,  60.59 };

  for (size_t i = 0; i < temperature.size(); ++i) {
    // RUN
    greal cp = gas.getSpecificHeatCapacity(WATER, temperature[i]);

    // EXPECTED (OUTPUTS)
    EXPECT_NEAR(expected[i], cp, 5e-3);
  }
  // TEAR DOWN
}
#endif
} //namespace