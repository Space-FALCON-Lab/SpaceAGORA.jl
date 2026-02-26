#include "unittest.h"
#include "MAP.h"
 
namespace GRAM {

// This is for regression purposes.  To be removed.
//TEST(MAP, init1)
//{
//  // SETUP
//  MAP map;
//
//  // GIVEN (INPUTS)
//
//  // RUN
//  map.init1(2);
//
//  // EXPECT (OUTPUTS)
//  EXPECT_DOUBLE_EQ(5358.0, map.pg[0][0]);
//  EXPECT_DOUBLE_EQ(0.002159, map.pg[map.zHeightSize - 1][map.zLatSize - 1]);
//
//  EXPECT_DOUBLE_EQ(0.08046, map.dg[0][0]);
//  EXPECT_DOUBLE_EQ(1.788e-08, map.dg[map.zHeightSize - 1][map.zLatSize - 1]);
//
//  EXPECT_DOUBLE_EQ(231.9, map.tg[0][0]);
//  EXPECT_DOUBLE_EQ(381.4, map.tg[map.zHeightSize - 1][map.zLatSize - 1]);
//
//  EXPECT_DOUBLE_EQ(2.0, map.ug[0][0]);
//  EXPECT_DOUBLE_EQ(-4.5, map.ug[map.zHeightSize - 1][map.zLatSize - 1]);
//
//  EXPECT_DOUBLE_EQ(0.0, map.psp[0][0][0]);
//  EXPECT_DOUBLE_EQ(0.005, map.psp[1][1][1]);
//  EXPECT_DOUBLE_EQ(0.013, map.psp[map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);
//
//  EXPECT_DOUBLE_EQ(0.0, map.dsp[0][0][0]);
//  EXPECT_DOUBLE_EQ(0.007, map.dsp[1][1][1]);
//  EXPECT_DOUBLE_EQ(0.02, map.dsp[map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);
//
//  EXPECT_DOUBLE_EQ(0.0, map.tsp[0][0][0]);
//  EXPECT_DOUBLE_EQ(-0.001, map.tsp[1][1][1]);
//  EXPECT_DOUBLE_EQ(-0.009, map.tsp[map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);
//
//  EXPECT_DOUBLE_EQ(0.0, map.usp[0][0][0]);
//  EXPECT_DOUBLE_EQ(-0.8, map.usp[1][1][1]);
//  EXPECT_DOUBLE_EQ(30.0, map.usp[map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);
//
//  EXPECT_DOUBLE_EQ(0.0, map.vsp[0][0][0]);
//  EXPECT_DOUBLE_EQ(1.6, map.vsp[1][1][1]);
//  EXPECT_DOUBLE_EQ(16.7, map.vsp[map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);
//
//
//  EXPECT_DOUBLE_EQ(8031.4170372082499, map.h2oa[0][0]);
//  EXPECT_DOUBLE_EQ(0.019901108397324299, map.o3[0][0]);
//  EXPECT_DOUBLE_EQ(0.30271310835320675, map.n2o[0][0]);
//  EXPECT_DOUBLE_EQ(0.15, map.co[0][0]);
//  EXPECT_DOUBLE_EQ(1.7, map.ch4[0][0]);
//
//  EXPECT_DOUBLE_EQ(0.2, map.h2oa[map.aLatSize - 1][map.aHeightSize  - 1]);
//  EXPECT_DOUBLE_EQ(0.0005, map.o3[map.aLatSize - 1][map.aHeightSize - 1]);
//  EXPECT_DOUBLE_EQ(0.00018500000000000019, map.n2o[map.aLatSize - 1][map.aHeightSize - 1]);
//  EXPECT_DOUBLE_EQ(49.99999999999995, map.co[map.aLatSize - 1][map.aHeightSize - 1]);
//  EXPECT_DOUBLE_EQ(0.03, map.ch4[map.aLatSize - 1][map.aHeightSize - 1]);
//
//  EXPECT_DOUBLE_EQ(1.4, map.h2omap[0][0]);
//  EXPECT_DOUBLE_EQ(0.6, map.sh2omap[0][0]);
//  EXPECT_DOUBLE_EQ(0.20775632527346732, map.n2omap[0][0]);
//  EXPECT_DOUBLE_EQ(0.14, map.ch4map[0][0]);
//  EXPECT_DOUBLE_EQ(26068213384.351967, map.oxymap[0][0]);
//  EXPECT_DOUBLE_EQ(0.0, map.ph2omap[0]);
//  EXPECT_DOUBLE_EQ(-1.2039728043259361, map.po3map[0]);
//  EXPECT_DOUBLE_EQ(2.3025850929940459, map.pmap31[0]);
//
//  EXPECT_DOUBLE_EQ(6.9421860288189139, map.h2omap[7][18]);
//  EXPECT_DOUBLE_EQ(0.68264829283385986, map.sh2omap[7][18]);
//  EXPECT_DOUBLE_EQ(154.96201849399003, map.n2omap[18][16]);
//  EXPECT_DOUBLE_EQ(1.5245642178118513, map.ch4map[18][16]);
//  EXPECT_DOUBLE_EQ(47260662.748543113, map.oxymap[18][18]);
//  EXPECT_DOUBLE_EQ(9.2103403719761836, map.ph2omap[18]);
//  EXPECT_DOUBLE_EQ(7.6009024595420822, map.po3map[23]);
//  EXPECT_DOUBLE_EQ(7.6009024595420822, map.pmap31[16]);
//
//
//  // TEAR-DOWN
//}


TEST(MAP, initializeData)
{
  // SETUP
  MAP map;

  // GIVEN (INPUTS)
  map.month = 2;

  // RUN
  map.initializeData();

  int m = map.month - 1;
  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(5358.0, map.pg[m][0][0]);
  EXPECT_DOUBLE_EQ(0.002159, map.pg[m][map.zHeightSize - 1][map.zLatSize - 1]);

  EXPECT_DOUBLE_EQ(0.08046, map.dg[m][0][0]);
  EXPECT_DOUBLE_EQ(1.788e-08, map.dg[m][map.zHeightSize - 1][map.zLatSize - 1]);

  EXPECT_DOUBLE_EQ(231.9, map.tg[m][0][0]);
  EXPECT_DOUBLE_EQ(381.4, map.tg[m][map.zHeightSize - 1][map.zLatSize - 1]);

  EXPECT_DOUBLE_EQ(2.0, map.ug[m][0][0]);
  EXPECT_DOUBLE_EQ(-4.5, map.ug[m][map.zHeightSize - 1][map.zLatSize - 1]);


  EXPECT_DOUBLE_EQ(0.0, map.psp[m][0][0][0]);
  EXPECT_DOUBLE_EQ(0.005, map.psp[m][1][1][1]);
  EXPECT_DOUBLE_EQ(0.013, map.psp[m][map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);

  EXPECT_DOUBLE_EQ(0.0, map.dsp[m][0][0][0]);
  EXPECT_DOUBLE_EQ(0.007, map.dsp[m][1][1][1]);
  EXPECT_DOUBLE_EQ(0.02, map.dsp[m][map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);

  EXPECT_DOUBLE_EQ(0.0, map.tsp[m][0][0][0]);
  EXPECT_DOUBLE_EQ(-0.001, map.tsp[m][1][1][1]);
  EXPECT_DOUBLE_EQ(-0.009, map.tsp[m][map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);

  EXPECT_DOUBLE_EQ(0.0, map.usp[m][0][0][0]);
  EXPECT_DOUBLE_EQ(-0.8, map.usp[m][1][1][1]);
  EXPECT_DOUBLE_EQ(30.0, map.usp[m][map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);

  EXPECT_DOUBLE_EQ(0.0, map.vsp[m][0][0][0]);
  EXPECT_DOUBLE_EQ(1.6, map.vsp[m][1][1][1]);
  EXPECT_DOUBLE_EQ(16.7, map.vsp[m][map.sHeightSize - 2][map.sLatSize - 2][map.sLonSize - 2]);

  // GIVEN (INPUTS)
  map.month = 3;

  // RUN
  map.initializeData();

  // EXPECT (OUTPUTS)
  m = map.month - 1;
  EXPECT_DOUBLE_EQ(0.0, map.pr[m][0][0]);
  EXPECT_DOUBLE_EQ(0.0, map.dr[m][0][0]);
  EXPECT_DOUBLE_EQ(0.0, map.tr[m][0][0]);
  EXPECT_DOUBLE_EQ(36.0, map.ur[m][0][0]);
  EXPECT_DOUBLE_EQ(25.0, map.vr[m][0][0]);
  EXPECT_DOUBLE_EQ(0.008836, map.pr[m][28][18]);
  EXPECT_DOUBLE_EQ(0.006400, map.dr[m][28][18]);
  EXPECT_DOUBLE_EQ(0.011236, map.tr[m][28][18]);
  EXPECT_DOUBLE_EQ(7569.000, map.ur[m][28][18]);
  EXPECT_DOUBLE_EQ(7569.000, map.vr[m][28][18]);


  // TEAR-DOWN
}

TEST(MAP, concinit)
{
  // SETUP
  MAP map;
  map.initializeData();

  // GIVEN (INPUTS)

  // RUN
  //map.concinit();

  int m = 2 - 1;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(8031.4170372082499, map.h2oa[m][0][0]);
  EXPECT_DOUBLE_EQ(0.019901108397324299, map.o3[m][0][0]);
  EXPECT_DOUBLE_EQ(0.30271310835320675, map.n2o[m][0][0]);
  EXPECT_DOUBLE_EQ(0.15, map.co[m][0][0]);
  EXPECT_DOUBLE_EQ(1.7, map.ch4[m][0][0]);

  EXPECT_DOUBLE_EQ(0.2, map.h2oa[m][map.aLatSize - 1][map.aHeightSize - 1]);
  EXPECT_DOUBLE_EQ(0.0005, map.o3[m][map.aLatSize - 1][map.aHeightSize - 1]);
  EXPECT_DOUBLE_EQ(0.00018500000000000019, map.n2o[m][map.aLatSize - 1][map.aHeightSize - 1]);
  EXPECT_DOUBLE_EQ(49.99999999999995, map.co[m][map.aLatSize - 1][map.aHeightSize - 1]);
  EXPECT_DOUBLE_EQ(0.03, map.ch4[m][map.aLatSize - 1][map.aHeightSize - 1]);

  // TEAR-DOWN
}

TEST(MAP, mapinit)
{
  // SETUP
  MAP map;
  map.initializeData();

  // GIVEN (INPUTS)

  // RUN
 // map.mapinit();

  int m = 2 - 1;
  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.4, map.h2omap[m][0][0]);
  EXPECT_DOUBLE_EQ(0.6, map.sh2omap[m][0][0]);
  EXPECT_DOUBLE_EQ(0.20775632527346732, map.n2omap[m][0][0]);
  EXPECT_DOUBLE_EQ(0.14, map.ch4map[m][0][0]);
  EXPECT_DOUBLE_EQ(26068213384.351967, map.oxymap[m][0][18]);
  EXPECT_DOUBLE_EQ(0.0, map.ph2omap[0]);
  EXPECT_DOUBLE_EQ(-1.2039728043259361, map.po3map[0]);
  EXPECT_DOUBLE_EQ(2.3025850929940459, map.pmap31[0]);

  EXPECT_DOUBLE_EQ(6.9421860288189139, map.h2omap[m][7][18]);
  EXPECT_DOUBLE_EQ(0.68264829283385986, map.sh2omap[m][7][18]);
  EXPECT_DOUBLE_EQ(154.96201849399003, map.n2omap[m][18][16]);
  EXPECT_DOUBLE_EQ(1.5245642178118513, map.ch4map[m][18][16]);
  EXPECT_DOUBLE_EQ(47260662.748543113, map.oxymap[m][18][0]);
  EXPECT_DOUBLE_EQ(9.2103403719761836, map.ph2omap[18]);
  EXPECT_DOUBLE_EQ(7.6009024595420822, map.po3map[23]);
  EXPECT_DOUBLE_EQ(7.6009024595420822, map.pmap31[16]);

  // TEAR-DOWN
}

//TEST(MAP, gterp)
//{
//	// SETUP
//	MAP map;
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

//TEST(MAP, pdtuv)
//{
//	// SETUP
//	MAP map;
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

TEST(MAP, mapmod)
{
  // SETUP
  MAP map;
  Position pos;

  // GIVEN (INPUTS)
  map.month = 2;
  map.initializeData();
  pos.height = 23.4;
  pos.latitude = 34.5;
  pos.longitude = 123.4;
  map.setPosition(pos);
  map.updatePosition();  // for gravity and radii

  // RUN
  map.update();
  const AtmosphereState& atmos = map.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3184.7106485525219, atmos.pressure);
  EXPECT_DOUBLE_EQ(0.050878287161947151, atmos.density);
  EXPECT_DOUBLE_EQ(217.909680678272, atmos.temperature);
  EXPECT_DOUBLE_EQ(9.4739580000000068, atmos.ewWind);
  EXPECT_DOUBLE_EQ(0.93547399999999759, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.00034952513039676082, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(0.0013290776399999969, map.temperatureGradient);

//  EXPECT_DOUBLE_EQ(-0.00034953006041461051, atmos.verticalWind);

  // Given
  pos.height = 73.4;
  pos.latitude = -34.5;
  pos.longitude = 323.4;
  map.setPosition(pos);
  map.updatePosition();  // for gravity and radii

  // RUN
  map.update();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3.2359171826613928, atmos.pressure);
  EXPECT_DOUBLE_EQ(5.404397835458626e-05, atmos.density);
  EXPECT_DOUBLE_EQ(208.54458633676396, atmos.temperature);
  EXPECT_DOUBLE_EQ(-13.256671999999986, atmos.ewWind);
  EXPECT_DOUBLE_EQ(0.55027600000000088, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.0004016004034657344, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(-0.0020473825605000057, map.temperatureGradient);

//  EXPECT_DOUBLE_EQ(-0.00040160752217206427, atmos.verticalWind);

  // TEAR-DOWN
}

TEST(MAP, heightInterpolation)
{
	// SETUP
	MAP map;

  // GIVEN (INPUTS)
  double p1 = 10.0;
  double d1 = 10.0;
  double t1 = 10.0;
  double z1 = 123.0;
  double p2 = 20.0;
  double d2 = 20.0;
  double t2 = 20.0;
  double z2 = 127.0;
  double z = 124.0;

  // RUN
  double p;
  double d;
  double t;
  map.heightInterpolation(p1, d1, t1, z1, p2, d2, t2, z2, p, d, t, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(12.5, p);
  EXPECT_DOUBLE_EQ(11.428571428571429, d);
  EXPECT_DOUBLE_EQ(12.5, t);

  // GIVEN (INPUTS)
  p1 = 10.0;
  d1 = 10.0;
  t1 = 10.0;
  z1 = 123.0;
  p2 = 20.0;
  d2 = 20.0;
  t2 = 10.0005;
  z2 = 127.0;
  z = 124.0;

  // RUN
  map.heightInterpolation(p1, d1, t1, z1, p2, d2, t2, z2, p, d, t, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(11.89207115002721, p);
  EXPECT_DOUBLE_EQ(11.89207114445308, d);
  EXPECT_DOUBLE_EQ(10.000125, t);

  // GIVEN (INPUTS)
  z1 = 123.0;
  z2 = 123.0005;

  // RUN
  map.heightInterpolation(p1, d1, t1, z1, p2, d2, t2, z2, p, d, t, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(p1, p);
  EXPECT_DOUBLE_EQ(d1, d);
  EXPECT_DOUBLE_EQ(t1, t);

  // GIVEN (INPUTS)
  p1 = 0.0;
  z1 = 123.0;
  z2 = 127.0;

  // RUN
  map.heightInterpolation(p1, d1, t1, z1, p2, d2, t2, z2, p, d, t, z);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, p);
  EXPECT_DOUBLE_EQ(0.0, d);
  EXPECT_DOUBLE_EQ(0.0, t);

  // TEAR-DOWN
}

TEST(MAP, afglconc)
{
	// SETUP
	MAP map;
  Position pos;

  // GIVEN (INPUTS)
  map.month = 2;
  map.initializeData();
  pos.height = 12.3;
  pos.latitude = 32.1;
  map.setPosition(pos);

  // RUN
  EarthPPM ppm;
  double water = map.afglconc(ppm);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(10.426198162869733, water);
  EXPECT_DOUBLE_EQ(0.24624630098568234, ppm.ozone);
  EXPECT_DOUBLE_EQ(0.29873336137184447, ppm.nitrousOxide);
  EXPECT_DOUBLE_EQ(0.07346788768948774, ppm.carbonMonoxide);
  EXPECT_DOUBLE_EQ(1.5661541968508701, ppm.methane);
  EXPECT_DOUBLE_EQ(329.99999999999994, ppm.carbonDioxide);
  EXPECT_DOUBLE_EQ(209469.99999999994, ppm.dioxygen);
  EXPECT_DOUBLE_EQ(780829.99999999988, ppm.dinitrogen);

  // GIVEN (INPUTS)
  pos.height = 114.7;
  pos.latitude = -12.1;
  map.setPosition(pos);

  // RUN
  water = map.afglconc(ppm);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.24223006688176782, water);
  EXPECT_DOUBLE_EQ(0.0057407681074844071, ppm.ozone);
  EXPECT_DOUBLE_EQ(0.00021379576533063887, ppm.nitrousOxide);
  EXPECT_DOUBLE_EQ(40.977510245983609, ppm.carbonMonoxide);
  EXPECT_DOUBLE_EQ(0.061677333799596626, ppm.methane);
  EXPECT_DOUBLE_EQ(40.985049776125493, ppm.carbonDioxide);
  EXPECT_DOUBLE_EQ(95387.410076813016, ppm.dioxygen);
  EXPECT_DOUBLE_EQ(765299.08230784012, ppm.dinitrogen);


  // GIVEN (INPUTS)
  pos.height = 123.0;
  pos.latitude = 32.1;
  map.setPosition(pos);

  // RUN
  water = map.afglconc(ppm);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, water);
  EXPECT_DOUBLE_EQ(0.0, ppm.ozone);
  EXPECT_DOUBLE_EQ(0.0, ppm.nitrousOxide);
  EXPECT_DOUBLE_EQ(0.0, ppm.carbonMonoxide);
  EXPECT_DOUBLE_EQ(0.0, ppm.methane);
  EXPECT_DOUBLE_EQ(0.0, ppm.carbonDioxide);
  EXPECT_DOUBLE_EQ(0.0, ppm.dioxygen);
  EXPECT_DOUBLE_EQ(0.0, ppm.dinitrogen);
                           
  // TEAR-DOWN
}

TEST(MAP, mapconc)
{
  // SETUP
  MAP map;
  Position pos;

  // GIVEN (INPUTS)
  map.month = 8;
  map.initializeData();
  pos.height = 12.3;
  pos.latitude = 23.4;
  double pgh = 34.5;
  map.setPosition(pos);

  // RUN
  double o3m;
  double h2om;
  double sh2om;
  double n2om;
  double ch4m;
  double oxym;
  map.mapconc(pgh, o3m, h2om, sh2om, n2om, ch4m, oxym);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.5726844349094102, o3m);
  EXPECT_DOUBLE_EQ(5.6972549156415191, h2om);
  EXPECT_DOUBLE_EQ(0.73789020648613513, sh2om);
  EXPECT_DOUBLE_EQ(0.0036532456895392021, n2om);
  EXPECT_DOUBLE_EQ(0.22693050757188032, ch4m);
  EXPECT_DOUBLE_EQ(0.0, oxym);

  // GIVEN (INPUTS)
  pos.height = 102.3;
  pos.latitude = -83.4;
  pgh = 9.6;
  map.setPosition(pos);

  // RUN
  map.mapconc(pgh, o3m, h2om, sh2om, n2om, ch4m, oxym);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.50079930421462826, o3m);
  EXPECT_DOUBLE_EQ(4.8791287792294806, h2om);
  EXPECT_DOUBLE_EQ(0.71043763597087795, sh2om);
  EXPECT_DOUBLE_EQ(0.0, n2om);
  EXPECT_DOUBLE_EQ(0.0, ch4m);
//  EXPECT_DOUBLE_EQ(2.3834074974580768e+17, oxym);
  EXPECT_DOUBLE_EQ(2.3834074974580682e+17, oxym);

  // TEAR-DOWN
}

TEST(MAP, larcwat)
{
  // SETUP
  MAP map;
  Position pos;

  // GIVEN (INPUTS)
  map.month = 8;
  map.initializeData();
  pos.height = 12.3;
  pos.latitude = 23.4;
  map.setPosition(pos);

  // RUN
  double wat = map.larcwat();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(34.507558181852623, wat);

  // GIVEN (INPUTS)
  pos.height = 50.3;
  pos.latitude = 23.4;
  map.setPosition(pos);

  // RUN
  wat = map.larcwat();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, wat);

  // GIVEN (INPUTS)
  pos.height = 32.3;
  pos.latitude = 90.0;
  map.setPosition(pos);

  // RUN
  wat = map.larcwat();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(6.0998020558347799, wat);

  // TEAR-DOWN
}

TEST(MAP, getStandardDeviations)
{
  // SETUP
  MAP map;
  Position pos;
  double p, d, t;

  // GIVEN (INPUTS)
  map.initializeData();
  pos.height = 0.0;
  pos.latitude = 0.0;
  map.setPosition(pos);
  map.month = 5;

  // RUN
  map.getStandardDeviations(p, d, t);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, p);
  EXPECT_DOUBLE_EQ(0.0, d);
  EXPECT_DOUBLE_EQ(0.0, t);

  // GIVEN (INPUTS)
  pos.height = 30.0;
  pos.latitude = 30.0;
  map.setPosition(pos);

  // RUN
  map.getStandardDeviations(p, d, t);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.000289, p);
  EXPECT_DOUBLE_EQ(0.000256, d);
  EXPECT_DOUBLE_EQ(0.000144, t);

  // GIVEN (INPUTS)
  pos.height = 134.7;
  pos.latitude = -36.6;
  map.setPosition(pos);

  // RUN
  map.getStandardDeviations(p, d, t);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.02035896, p);
  EXPECT_DOUBLE_EQ(0.005464105, d);
  EXPECT_DOUBLE_EQ(0.027465325, t);

  // TEAR-DOWN
}

} // namespace