#include "unittest.h"
#include "MET.h"
 
namespace GRAM {


TEST(MET, molwt)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)

  // RUN
  double x = met.molwt(90.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(28.878081999999999, x);

  // RUN
  x = met.molwt(101.4);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(28.032413290267716, x);

  // RUN
  x = met.molwt(105.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(27.725942781249998, x);

  // RUN
  x = met.molwt(106.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, x);

  // RUN
  x = met.molwt(89.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(28.907485049985997, x);

  // TEAR-DOWN
}

TEST(MET, temp)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  double alt = 90.0;
  double tx = 231.5;
  double t1 = 23.4;
  double t3 = 0.123;
  double t4 = 0.006;
  double a2 = 35.7;

  // RUN
  double x = met.temp(alt, tx, t1, t3, t4, a2);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3142.625, x);


  // GIVEN (INPUTS)
  alt = 190.0;
  tx = 231.5;
  t1 = 23.4;
  t3 = 0.123;
  t4 = 0.006;
  a2 = 35.7;

  // RUN
  x = met.temp(alt, tx, t1, t3, t4, a2);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(286.85096984199362, x);


  // TEAR-DOWN
}

TEST(MET, gauss)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  int nmin = 2;
  double z2 = 90.0;
  double tx = 231.5;
  double t1 = 23.4;
  double t3 = 0.123;
  double t4 = 0.006;
  double a2 = 35.7;
  double gphi = 43.5;
  double rphi = 27.8;

  // RUN
  double r = 0;
  met.gauss(nmin, z2, tx, t1, t3, t4, a2, gphi, rphi, r);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.9334191346696403, r);

  // GIVEN (INPUTS)
  nmin = 4;
  z2 = 190.0;
  tx = 231.5;
  t1 = 23.4;
  t3 = 0.123;
  t4 = 0.006;
  a2 = 35.7;
  gphi = 43.5;
  rphi = 27.8;

  // RUN
  r = 0;
  met.gauss(nmin, z2, tx, t1, t3, t4, a2, gphi, rphi, r);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.086059520752504698, r);


  // TEAR-DOWN
}

TEST(MET, met_07_jac)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  double z = 95.6;
  double t = 123.4;
  double gphi = 34.5;
  double rphi = 234.5;
  // RUN
  double tz, an, ao2, ao, aa, ahe, ah, em, dens, dl;
  met.met_07_jac(z, t, tz, an, ao2, ao, aa, ahe, ah, em, dens, dl, gphi, rphi);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(182.50131634954715, tz);
  EXPECT_DOUBLE_EQ(12.931659773608937, an);
  EXPECT_DOUBLE_EQ(12.327151973225863, ao2);
  EXPECT_DOUBLE_EQ(11.526728490891804, ao);
  EXPECT_DOUBLE_EQ(11.009439484650922, aa);
  EXPECT_DOUBLE_EQ(8.149206053156675, ahe);
  EXPECT_DOUBLE_EQ(0.0, ah);
  EXPECT_DOUBLE_EQ(28.526377397401134, em);
  EXPECT_DOUBLE_EQ(5.2610722608400882e-10, dens);
  EXPECT_DOUBLE_EQ(-9.2789257331261279, dl);

  // GIVEN (INPUTS)
  z = 195.6;
  t = 123.4;
  gphi = 34.5;
  rphi = 234.5;
  // RUN
  met.met_07_jac(z, t, tz, an, ao2, ao, aa, ahe, ah, em, dens, dl, gphi, rphi);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(124.5568325811972, tz);
  EXPECT_DOUBLE_EQ(-0.88923456335377749, an);
  EXPECT_DOUBLE_EQ(-3.3681672039727135, ao2);
  EXPECT_DOUBLE_EQ(3.6054460154560375, ao);
  EXPECT_DOUBLE_EQ(-8.2133512112298881, aa);
  EXPECT_DOUBLE_EQ(5.1377231880822141, ahe);
  EXPECT_DOUBLE_EQ(0.0, ah);
  EXPECT_DOUBLE_EQ(4.3447523105934787, em);
  EXPECT_DOUBLE_EQ(1.0197814768462893e-18, dens);
  EXPECT_DOUBLE_EQ(-17.991492880754997, dl);

  // GIVEN (INPUTS)
  z = 795.6;
  t = 123.4;
  gphi = 34.5;
  rphi = 234.5;
  // RUN
  met.met_07_jac(z, t, tz, an, ao2, ao, aa, ahe, ah, em, dens, dl, gphi, rphi);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(123.4027164162828, tz);
  EXPECT_DOUBLE_EQ(-31.297024985738101, an);
  EXPECT_DOUBLE_EQ(-38.102577039500481, ao2);
  EXPECT_DOUBLE_EQ(-13.75973748428429, ao);
  EXPECT_DOUBLE_EQ(-51.577550538141161, aa);
  EXPECT_DOUBLE_EQ(0.79493771475380981, ahe);
  EXPECT_DOUBLE_EQ(14.469687337394182, ah);
  EXPECT_DOUBLE_EQ(1.0079700000000635, em);
  EXPECT_DOUBLE_EQ(4.9361010462054316e-10, dens);
  EXPECT_DOUBLE_EQ(-9.3066159585047448, dl);

  // TEAR-DOWN
}

TEST(MET, slv)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  greal alt = 90.0;
  greal lat = 67.8;
  greal day = 45.0;

  // RUN
  double den = 0;
  met.slv(alt, lat, day, den);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, den);

  // GIVEN (INPUTS)
  alt = 125.0;
  lat = 67.8;
  day = 45.0;

  // RUN
  met.slv(alt, lat, day, den);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.051545781341766378, den);

  // GIVEN (INPUTS)
  alt = 125.0;
  lat = -67.8;
  day = 45.0;

  // RUN
  met.slv(alt, lat, day, den);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.051545781341766378, den);

  // GIVEN (INPUTS)
  alt = 185.0;
  lat = 67.8;
  day = 45.0;

  // RUN
  met.slv(alt, lat, day, den);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, den);

  // TEAR-DOWN
}

TEST(MET, slvh)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  greal lat = 27.8;
  greal sda = toRadians(45.0);
  greal den = 1.123;
  greal denhe = 23.4;

  // RUN
  met.slvh(lat, sda, den, denhe);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0970023352653726, den);
  EXPECT_DOUBLE_EQ(23.130770489000394, denhe);

  // TEAR-DOWN
}

TEST(MET, met07)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  MET::InData inData;
  inData.z = 132.0;   // z;
  inData.lat = -16.3;   // phi;
  inData.lng = 145.7;   // thet;
  inData.i1 = 2;   // 2.0;
  inData.f10 = 23.4;   // f10;
  inData.f10b = 12.3;  // f10b;
  inData.gi = 7.9;  // ap;

  // RUN
  MET::OutData outData;
  met.met07(inData, outData);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(527.51988613470928, outData.tex);
  EXPECT_DOUBLE_EQ(382.12583952945079, outData.th);
  EXPECT_DOUBLE_EQ(85850661893319840.0, outData.n2nd);
  EXPECT_DOUBLE_EQ(11627179937341834.0, outData.o2nd);
  EXPECT_DOUBLE_EQ(37417495951954920.0, outData.ond);
  EXPECT_DOUBLE_EQ(270949100742939.06, outData.arnd);
  EXPECT_DOUBLE_EQ(26221207573636.977, outData.hend);
  EXPECT_DOUBLE_EQ(989410.17380044539, outData.hnd);
  EXPECT_DOUBLE_EQ(25.050286421014324, outData.wtmol);
  EXPECT_DOUBLE_EQ(5.6235999533128768e-09, outData.dh);
  EXPECT_DOUBLE_EQ(-8.2499855813467295, outData.unk);
  EXPECT_DOUBLE_EQ(0.00071325087596809864, outData.ph);

  // TEAR-DOWN
}

TEST(MET, tinf)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  int i1 = 1;
  double f10 = 34.5;
  double f10b = 23.4;
  double gi = 9.8;
  double xlat = 2.34;
  double sda = 0.123;
  double sha = 0.234;
  double dy = 43.2;

  // RUN
  double te;
  met.tinf(i1, f10, f10b, gi, xlat, sda, sha, dy, te);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1326.32705066419, te);

  // GIVEN (INPUTS)
  i1 = 2;
  f10 = 34.5;
  f10b = 23.4;
  gi = 9.8;
  xlat = 2.34;
  sda = 0.123;
  sha = 0.234;
  dy = 43.2;

  // RUN
  met.tinf(i1, f10, f10b, gi, xlat, sda, sha, dy, te);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(575.05709786700277, te);

  // TEAR-DOWN
}

TEST(MET, jacch)
{
  //// SETUP
  //Met met;

  //// GIVEN (INPUTS)
  //double xlat = toRadians(27.8);
  //double sda = toRadians(45.0);
  //double den = 1.123;
  //double denhe = 23.4;

  //// RUN
  //met.slvh(xlat, sda, den, denhe);

  //// EXPECT (OUTPUTS)
  //EXPECT_DOUBLE_EQ(1.0970023352653726, den);
  //EXPECT_DOUBLE_EQ(23.130770489000394, denhe);

  //// TEAR-DOWN
}

TEST(MET, update)
{
  // SETUP
  MET met;
  Position pos;

  // GIVEN (INPUTS)
  met.year = 2020;
  met.month = 3;
  met.day = 17;
  met.minute = 14;
  met.seconds = 34.5;
  met.hour = 22;

  met.f10 = 34.5;
  met.f10b = 23.4;
  met.ap = 45.6;
  pos.height = 123.0;
  pos.latitude = -34.5;
  pos.longitude = 258.4;
  pos.elapsedTime = 17.3;
  met.setPosition(pos);
  met.updatePosition();

  // RUN
  met.update();
  const AtmosphereState& atmos = met.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0017004775570408598, atmos.pressure);
  EXPECT_DOUBLE_EQ(1.4848388651386517e-08, atmos.density);
  EXPECT_DOUBLE_EQ(359.83609746286385, atmos.temperature);
  EXPECT_DOUBLE_EQ(60.464620146527103, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-9.0286346871338292, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.004546092681328542, atmos.verticalWind);

  EXPECT_DOUBLE_EQ(1146231053571349.5, atmos.argon.numberDensity);
  EXPECT_DOUBLE_EQ(29537850477800.898, atmos.helium.numberDensity);
  EXPECT_DOUBLE_EQ(921338.54480655491, atmos.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(37007759013073816, atmos.dioxygen.numberDensity);
  EXPECT_DOUBLE_EQ(2.3692513153986326e+17, atmos.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(67172420620019952, atmos.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(26.124460778690445, atmos.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(MET, wind)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  double ph = 0.0;
  double dh = 0.0;
  double th = 0.0;
  double h = 0.0;
  double g = 0.0;
  double phid = 0.0;
  double ri = 0.0;
  double dpx = 0.0;
  double dpy = 0.0;
  double dtx = 0.0;
  double dty = 0.0;
  double dtz = 0.0;

  // RUN
  double ugh;
  double vgh;
  double wgh;
  met.wind(ph, dh, th, h, g, phid, ri, dpx, dpy, dtx, dty, dtz, ugh, vgh, wgh);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, ugh);
  EXPECT_DOUBLE_EQ(0.0, vgh);
  EXPECT_DOUBLE_EQ(0.0, wgh);

  // GIVEN (INPUTS)
  ph = 12.3;
  dh = 23.4;
  th = 34.5;
  h = 6.7;
  g = 9.8;
  phid = 78.9;
  ri = 6371.2;
  dpx = 12.3;
  dpy = 23.4;
  dtx = 34.5;
  dty = 45.6;
  dtz = 56.7;

  // RUN
  met.wind(ph, dh, th, h, g, phid, ri, dpx, dpy, dtx, dty, dtz, ugh, vgh, wgh);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.012601919935762842, ugh);
  EXPECT_DOUBLE_EQ(0.034406910750526168, vgh);
  EXPECT_DOUBLE_EQ(5.1555722556614831e-09, wgh);

  // TEAR-DOWN
}

TEST(MET, fair5)
{
  // SETUP
  MET met;

  // GIVEN (INPUTS)
  double dhel1 = 1.1;
  double dhel2 = 2.2;
  double dlg1 = 123.4;
  double dlg2 = 456.7;
  double h = 440.0;

  // RUN
  double fdhel;
  double fdlg;
  met.fair5(dhel1, dhel2, dlg1, dlg2, h, fdhel, fdlg);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.1, fdhel);
  EXPECT_DOUBLE_EQ(123.4, fdlg);

  // GIVEN (INPUTS)
  dhel1 = 1.1;
  dhel2 = 2.2;
  dlg1 = 123.4;
  dlg2 = 456.7;
  h = 500.0;

  // RUN
  met.fair5(dhel1, dhel2, dlg1, dlg2, h, fdhel, fdlg);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.2, fdhel);
  EXPECT_DOUBLE_EQ(456.7, fdlg);

  // GIVEN (INPUTS)
  dhel1 = 1.1;
  dhel2 = 2.2;
  dlg1 = 123.4;
  dlg2 = 456.7;
  h = 543.2;

  // RUN
  met.fair5(dhel1, dhel2, dlg1, dlg2, h, fdhel, fdlg);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.2994168056382196, fdhel);
  EXPECT_DOUBLE_EQ(183.82329210838049, fdlg);

  // TEAR-DOWN
}

}