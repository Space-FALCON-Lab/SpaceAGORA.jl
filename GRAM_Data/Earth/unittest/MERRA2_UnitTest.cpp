#include <math.h>
#include "unittest.h"
#include "MERRA2.h"
 
namespace GRAM {

TEST(MERRA2, getIndex)
{
  // SETUP
  MERRA2* m2 = new MERRA2();
  m2->M2Hour = 1;
  m2->month = 1;
  m2->initialized = false;
  m2->setExtents(0, 10, 0, 10);
  m2->readInfoFile();
  m2->computeArraySizes();

  // GIVEN (INPUTS)
  size_t iPres = 0;
  size_t iLat = 0;
  size_t iLon = 0;
  size_t iLonS = 0;

  // RUN
  size_t index3D = m2->getIndex(iPres, iLat, iLon);
  size_t index2D = m2->getIndex(iLat, iLon);
  size_t index2DS = m2->getSurfaceIndex(iLat, iLonS);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(size_t(0), index3D);
  EXPECT_EQ(size_t(0), index2D);
  EXPECT_EQ(size_t(0), index2DS);

  // GIVEN (INPUTS)
  iPres = m2->presSize - 1;
  iLat = m2->latSize - 1;
  iLon = m2->lonSize - 1;
  iLonS = m2->lonSurfaceSize - 1;

  // RUN
  index3D = m2->getIndex(iPres, iLat, iLon);
  index2D = m2->getIndex(iLat, iLon);
  index2DS = m2->getSurfaceIndex(iLat, iLonS);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(size_t(m2->arraySize3D - 1), index3D);
  EXPECT_EQ(size_t(m2->arraySize2D - 1), index2D);
  EXPECT_EQ(size_t(m2->arraySize2DSurface - 1), index2DS);

  // GIVEN (INPUTS)
  iPres = m2->presSize / 2;
  iLat = 0;
  iLon = 0;

  // RUN
  index3D = m2->getIndex(iPres, iLat, iLon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(size_t(m2->arraySize3D / 2), index3D);

  // TEAR-DOWN
  delete m2;
}

TEST(MERRA2, readM2File)
{
  // SETUP
  MERRA2* m1 = new MERRA2();
  MERRA2* m2 = new MERRA2();

  // GIVEN (INPUTS)
  m1->M2Hour = 1;
  m1->month = 1;
  m1->initialized = false;
  m1->setExtents(-20, 20, 150, 170);

  m2->M2Hour = 1;
  m2->month = 1;
  m2->initialized = false;
  m2->setExtents(-40, 0, 160, 180);

  // RUN
  m1->readM2File();
  m2->readM2File();

  size_t halfLat = (m1->latSize - 1) / 2;
  size_t halfLon = (m1->lonSize - 1) / 2;
  size_t halfSfcLon = (m1->lonSurfaceSize - 1) / 2;
  for (size_t jLat = 0; jLat < halfLat; ++jLat) {
    for (size_t jLon = 0; jLon < halfLon; ++jLon) {
      for (size_t iPres = 0; iPres < m1->presSize; ++iPres) {
        size_t index1 = m1->getIndex(iPres, jLat, jLon + halfLon);
        size_t index2 = m2->getIndex(iPres, jLat + halfLat, jLon);
        if (isnan(m1->temp[index1])) {
          EXPECT_TRUE(isnan(m2->temp[index2]));
        }
        else {
          EXPECT_DOUBLE_EQ((double)m1->temp[index1], (double)m2->temp[index2]);
        }
        EXPECT_DOUBLE_EQ((double)m1->spdsd[index1], (double)m2->spdsd[index2]);
        EXPECT_DOUBLE_EQ((double)m1->dens[index1], (double)m2->dens[index2]);
        EXPECT_DOUBLE_EQ((double)m1->hgt[index1], (double)m2->hgt[index2]);
      }
      size_t index1 = m1->getIndex(jLat, jLon + halfLon);
      size_t index2 = m2->getIndex(jLat + halfLat, jLon);
      EXPECT_DOUBLE_EQ((double)m1->psfc[index1], (double)m2->psfc[index2]);
      EXPECT_DOUBLE_EQ((double)m1->srhsfc[index1], (double)m2->srhsfc[index2]);
    }
    for (size_t jLon = 0; jLon < halfSfcLon; ++jLon) {
      size_t index1 = m1->getSurfaceIndex(jLat, jLon + halfSfcLon);
      size_t index2 = m2->getSurfaceIndex(jLat + halfLat, jLon);
      EXPECT_DOUBLE_EQ((double)m1->usfc[index1], (double)m2->usfc[index2]);
      EXPECT_DOUBLE_EQ((double)m1->sdsfc[index1], (double)m2->sdsfc[index2]);
    }
  }

  // TEAR-DOWN
  delete m1;
  delete m2;

  // SETUP
  m1 = new MERRA2();
  m2 = new MERRA2();

  // GIVEN (INPUTS)
  m1->M2Hour = 1;
  m1->month = 1;
  m1->initialized = false;
  m1->setExtents(-20, 20, -40, 40);

  m2->M2Hour = 1;
  m2->month = 1;
  m2->initialized = false;
  m2->setExtents(-40, 0, 0, 80);

  // RUN
  m1->readM2File();
  m2->readM2File();

  halfLat = (m1->latSize - 1) / 2;
  halfLon = (m1->lonSize - 1) / 2;
  halfSfcLon = (m1->lonSurfaceSize - 1) / 2;
  for (size_t jLat = 0; jLat < halfLat; ++jLat) {
    for (size_t jLon = 0; jLon < halfLon; ++jLon) {
      for (size_t iPres = 0; iPres < m1->presSize; ++iPres) {
        size_t index1 = m1->getIndex(iPres, jLat, jLon + halfLon);
        size_t index2 = m2->getIndex(iPres, jLat + halfLat, jLon);
        if (isnan(m1->temp[index1])) {
          EXPECT_TRUE(isnan(m2->temp[index2]));
        }
        else {
          EXPECT_DOUBLE_EQ((double)m1->temp[index1], (double)m2->temp[index2]);
        }
        if (isnan(m1->spdsd[index1])) {
          EXPECT_TRUE(isnan(m2->spdsd[index2]));
        }
        else {
          EXPECT_DOUBLE_EQ((double)m1->spdsd[index1], (double)m2->spdsd[index2]);
        }
        if (isnan(m1->dens[index1])) {
          EXPECT_TRUE(isnan(m2->dens[index2]));
        }
        else {
          EXPECT_DOUBLE_EQ((double)m1->dens[index1], (double)m2->dens[index2]);
        }
        EXPECT_DOUBLE_EQ((double)m1->hgt[index1], (double)m2->hgt[index2]);
      }
      size_t index1 = m1->getIndex(jLat, jLon + halfLon);
      size_t index2 = m2->getIndex(jLat + halfLat, jLon);
      EXPECT_DOUBLE_EQ((double)m1->psfc[index1], (double)m2->psfc[index2]);
      EXPECT_DOUBLE_EQ((double)m1->srhsfc[index1], (double)m2->srhsfc[index2]);
    }
    for (size_t jLon = 0; jLon < halfSfcLon; ++jLon) {
      size_t index1 = m1->getSurfaceIndex(jLat, jLon + halfSfcLon);
      size_t index2 = m2->getSurfaceIndex(jLat + halfLat, jLon);
      EXPECT_DOUBLE_EQ((double)m1->usfc[index1], (double)m2->usfc[index2]);
      EXPECT_DOUBLE_EQ((double)m1->sdsfc[index1], (double)m2->sdsfc[index2]);
    }
  }

  // TEAR-DOWN
  delete m1;
  delete m2;
}

// Requires large data set.
TEST(MERRA2, DISABLED_readM2File_dataPoints)
{
  // SETUP
  MERRA2* m1 = new MERRA2();
  MERRA2* m2 = new MERRA2();
  MERRA2* m3 = new MERRA2();

  // GIVEN (INPUTS)
  m1->M2Hour = 9;
  m1->month = 1;
  m1->initialized = false;
  m1->setExtents(-90, 90, 0, 360);

  m2->M2Hour = 9;
  m2->month = 1;
  m2->initialized = false;
  m2->setExtents(-10, 10, 170, 190);

  m3->M2Hour = 9;
  m3->month = 5;
  m3->initialized = false;
  m3->setExtents(-10, 10, 170, 190);

  // RUN
  m1->readM2File();
  m2->readM2File();
  m3->readM2File();

  m1->longitude = 180.0;
  m1->updateIndices();
  size_t ilon = m1->lonIndex;

  EXPECT_DOUBLE_EQ(1.321703005176143e+00, m1->dens[m1->getIndex(0, 20, ilon)]);
  EXPECT_DOUBLE_EQ(1.003125446894883e-02, m1->sden[m1->getIndex(0, 20, ilon)]);
  EXPECT_DOUBLE_EQ(1.320758523491196e+00, m1->dens[m1->getIndex(0, 20, ilon + 1)]);
  EXPECT_DOUBLE_EQ(1.028455781812772e-02, m1->sden[m1->getIndex(0, 20, ilon + 1)]);
  EXPECT_DOUBLE_EQ(1.400522537357747e+00, m1->dens[m1->getIndex(0, 359, 570 - ilon)]);
  EXPECT_DOUBLE_EQ(3.486419706413960e-02, m1->sden[m1->getIndex(0, 359, 570 - ilon)]);

  m2->latitude = 0.0;
  m2->longitude = 180.0;
  m2->updateIndices();
  size_t ilat2 = m2->latIndex;
  size_t ilon2 = m2->lonIndex;
  size_t ilons2 = m2->lonSurfaceIndex;

  EXPECT_DOUBLE_EQ(1.163474154353986e+00, m2->dens[m2->getIndex(0, ilat2, ilon2)]);
  EXPECT_DOUBLE_EQ(4.201215118028831e-03, m2->sden[m2->getIndex(0, ilat2, ilon2)]);
  EXPECT_DOUBLE_EQ(1.163629964488590e+00, m2->dens[m2->getIndex(0, ilat2, ilon2 + 1)]);
  EXPECT_DOUBLE_EQ(4.259756392522181e-03, m2->sden[m2->getIndex(0, ilat2, ilon2 + 1)]);
  EXPECT_DOUBLE_EQ(1.162399589938743e+00, m2->dens[m2->getIndex(0, ilat2, ilon2 - 7)]);
  EXPECT_DOUBLE_EQ(3.723538339256247e-03, m2->sden[m2->getIndex(0, ilat2, ilon2 - 7)]);

  EXPECT_DOUBLE_EQ(2.997146167385963e+02,  m2->tsfc[m2->getSurfaceIndex(ilat2, ilons2)]);
  EXPECT_DOUBLE_EQ(9.764457417535012e-01, m2->stsfc[m2->getSurfaceIndex(ilat2, ilons2)]);
  EXPECT_DOUBLE_EQ(2.996886948616274e+02,  m2->tsfc[m2->getSurfaceIndex(ilat2, ilons2 + 1)]);
  EXPECT_DOUBLE_EQ(9.914426857414584e-01, m2->stsfc[m2->getSurfaceIndex(ilat2, ilons2 + 1)]);
  EXPECT_DOUBLE_EQ(2.998818184145035e+02,  m2->tsfc[m2->getSurfaceIndex(ilat2, ilons2 - 7)]);
  EXPECT_DOUBLE_EQ(8.860921790184959e-01, m2->stsfc[m2->getSurfaceIndex(ilat2, ilons2 - 7)]);

  EXPECT_DOUBLE_EQ(8.109996002751666e+01, m2->rhsfc[m2->getIndex(0, ilat2, ilon2)]);
  EXPECT_DOUBLE_EQ(4.199347685975059e+00, m2->srhsfc[m2->getIndex(0, ilat2, ilon2)]);
  EXPECT_DOUBLE_EQ(8.110238654406007e+01, m2->rhsfc[m2->getIndex(0, ilat2, ilon2 + 1)]);
  EXPECT_DOUBLE_EQ(4.207887271034185e+00, m2->srhsfc[m2->getIndex(0, ilat2, ilon2 + 1)]);
  EXPECT_DOUBLE_EQ(8.108741673934475e+01, m2->rhsfc[m2->getIndex(0, ilat2, ilon2 - 7)]);
  EXPECT_DOUBLE_EQ(4.195904188442703e+00, m2->srhsfc[m2->getIndex(0, ilat2, ilon2 - 7)]);

  m3->latitude = 0.0;
  m3->longitude = 180.0;
  m3->updateIndices();
  size_t ilat3 = m3->latIndex;
  size_t ilon3 = m3->lonIndex;

  EXPECT_DOUBLE_EQ(7.126054086383955, m3->spdav[m3->getIndex(19, ilat3, ilon3)]);
  EXPECT_DOUBLE_EQ(3.786852590240998, m3->spdsd[m3->getIndex(19, ilat3, ilon3)]);
  EXPECT_DOUBLE_EQ(7.130901718232131, m3->spdav[m3->getIndex(19, ilat3, ilon3 + 1)]);
  EXPECT_DOUBLE_EQ(3.813703996306534, m3->spdsd[m3->getIndex(19, ilat3, ilon3 + 1)]);
  EXPECT_DOUBLE_EQ(6.968475117365108, m3->spdav[m3->getIndex(19, ilat3, ilon3 - 7)]);
  EXPECT_DOUBLE_EQ(3.802420381955755, m3->spdsd[m3->getIndex(19, ilat3, ilon3 - 7)]);

  // TEAR-DOWN
  delete m1;
}

TEST(MERRA2, getGasConstant)
{
  // SETUP
  MERRA2* m2 = new MERRA2();

  // GIVEN (INPUTS)

  // RUN
  greal con = m2->getGasConstant(0.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05500000000001, con);

  // RUN
  con = m2->getGasConstant(75.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05500000000001, con);

  // RUN
  con = m2->getGasConstant(300.0, 101325.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.09302739172239, con);

  // RUN
  con = m2->getGasConstant(200.0, 10000.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05503609910716, con);

  // TEAR-DOWN
  delete m2;
}

TEST(MERRA2, updateIndices_basic)
{
  // SETUP
  MERRA2* m2 = new MERRA2();
  m2->latFullSize = 181;
  m2->lonFullSize = 360;
  m2->lonSurfaceFullSize = 360;
  m2->minimumLatitude = -90.0_deg;
  m2->maximumLatitude = 90.0_deg;
  m2->minimumLongitude = 0.0_deg;
  m2->maximumLongitude = 360.0_deg;
  m2->computeArraySizes();

  // GIVEN (INPUTS)
  m2->latitude = -90.0_deg;
  m2->longitude = 0.0_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)0, m2->latIndex);
  EXPECT_EQ((size_t)0, m2->lonIndex);
  EXPECT_EQ((size_t)0, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLon);
  EXPECT_DOUBLE_EQ(0.0, m2->latFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonSurfaceFactor);

  // GIVEN (INPUTS)
  m2->latitude = 90.0_deg;
  m2->longitude = 360.0_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)179, m2->latIndex);
  EXPECT_EQ((size_t)360, m2->lonIndex);
  EXPECT_EQ((size_t)360, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLon);
  EXPECT_DOUBLE_EQ(1.0, m2->latFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonSurfaceFactor);

  // GIVEN (INPUTS)
  m2->latitude = 12.5_deg;
  m2->longitude = 187.4_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)102, m2->latIndex);
  EXPECT_EQ((size_t)187, m2->lonIndex);
  EXPECT_EQ((size_t)187, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLon);
  EXPECT_NEAR(0.5, m2->latFactor, 1.0e-12);
  EXPECT_NEAR(0.4, m2->lonFactor, 1.0e-12);
  EXPECT_NEAR(0.4, m2->lonSurfaceFactor, 1.0e-12);

  // GIVEN (INPUTS)
  m2->latitude = -76.2_deg;
  m2->longitude = 359.1_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)13, m2->latIndex);
  EXPECT_EQ((size_t)359, m2->lonIndex);
  EXPECT_EQ((size_t)359, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(1.0, m2->gridLon);
  EXPECT_NEAR(0.8, m2->latFactor, 1.0e-12);
  EXPECT_NEAR(0.1, m2->lonFactor, 1.0e-12);
  EXPECT_NEAR(0.1, m2->lonSurfaceFactor, 1.0e-12);

  // TEAR-DOWN
  delete m2;
}

TEST(MERRA2, updateIndices_complex)
{
  // SETUP
  MERRA2* m2 = new MERRA2();
  m2->latFullSize = 361;
  m2->lonFullSize = 360 * 4;
  m2->lonSurfaceFullSize = 180;
  m2->minimumLatitude = -30.0_deg;
  m2->maximumLatitude = 20.0_deg;
  m2->minimumLongitude = 100.0_deg;
  m2->maximumLongitude = 200.0_deg;
  m2->computeArraySizes();

  // GIVEN (INPUTS)
  m2->latitude = -30.0_deg;
  m2->longitude = 100.0_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)0, m2->latIndex);
  EXPECT_EQ((size_t)0, m2->lonIndex);
  EXPECT_EQ((size_t)0, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(2.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(4.0, m2->gridLon);
  EXPECT_DOUBLE_EQ(0.0, m2->latFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonSurfaceFactor);

  // GIVEN (INPUTS)
  m2->latitude = 20.0_deg;
  m2->longitude = 200.0_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)99, m2->latIndex);
  EXPECT_EQ((size_t)400, m2->lonIndex);
  EXPECT_EQ((size_t)50, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(2.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(4.0, m2->gridLon);
  EXPECT_DOUBLE_EQ(1.0, m2->latFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonFactor);
  EXPECT_DOUBLE_EQ(0.0, m2->lonSurfaceFactor);

  // GIVEN (INPUTS)
  m2->latitude = 12.25_deg;
  m2->longitude = 187.4_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)84, m2->latIndex);
  EXPECT_EQ((size_t)349, m2->lonIndex);
  EXPECT_EQ((size_t)43, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(2.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(4.0, m2->gridLon);
  EXPECT_NEAR(0.5, m2->latFactor, 1.0e-12);
  EXPECT_NEAR(0.6, m2->lonFactor, 1.0e-12);
  EXPECT_NEAR(0.7, m2->lonSurfaceFactor, 1.0e-12);

  // GIVEN (INPUTS)
  m2->latitude = -12.1_deg;
  m2->longitude = 199.9_deg;

  // RUN
  m2->updateIndices();

  // EXPECT (OUTPUTS)
  EXPECT_EQ((size_t)35, m2->latIndex);
  EXPECT_EQ((size_t)399, m2->lonIndex);
  EXPECT_EQ((size_t)49, m2->lonSurfaceIndex);
  EXPECT_DOUBLE_EQ(2.0, m2->gridLat);
  EXPECT_DOUBLE_EQ(4.0, m2->gridLon);
  EXPECT_NEAR(0.8, m2->latFactor, 1.0e-12);
  EXPECT_NEAR(0.6, m2->lonFactor, 1.0e-12);
  EXPECT_NEAR(0.95, m2->lonSurfaceFactor, 1.0e-12);

  // TEAR-DOWN
  delete m2;
}

// Requires large data set.
TEST(MERRA2, DISABLED_getHeights)
{
  // SETUP
  MERRA2* m2 = new MERRA2();
  m2->M2Hour = 1;
  m2->month = 1;
  m2->initialized = false;
  m2->setExtents(-90, 90, 0, 360);

  // GIVEN (INPUTS)

  // RUN
  greal lowest0p3mb, highest0p1mb;
  m2->getHeights(0.0, 0.0, lowest0p3mb, highest0p1mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(56.528831243699599, lowest0p3mb);
  EXPECT_DOUBLE_EQ(64.23528557207662, highest0p1mb);

  // RUN
  m2->getHeights(45.0, 123.0, lowest0p3mb, highest0p1mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(55.707364151965734, lowest0p3mb);
  EXPECT_DOUBLE_EQ(63.041275468750001, highest0p1mb);

  // RUN
  m2->getHeights(90.0, 213.0, lowest0p3mb, highest0p1mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(53.182574640877014, lowest0p3mb);
  EXPECT_DOUBLE_EQ(61.10720078125, highest0p1mb);

  // RUN
  m2->getHeights(-90.0, 321.0, lowest0p3mb, highest0p1mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(59.85945326360887, lowest0p3mb);
  EXPECT_DOUBLE_EQ(67.84506761592742, highest0p1mb);

  // TEAR-DOWN
  delete m2;
}

} // namespace