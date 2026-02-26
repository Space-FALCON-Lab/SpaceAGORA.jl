//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "EarthAtmosphere.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

bool EarthAtmosphere::initialized = false;
int EarthAtmosphere::landcd[topoLonSize][topoLatSize];
greal EarthAtmosphere::ztopo[topoLonSize][topoLatSize];

greal EarthAtmosphere::xlbar[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::zlbar[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::xsigl[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::zsigl[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::xlmin[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::zlmin[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::xscale[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::zscale[rsHeightSize] = { 0.0 };
greal EarthAtmosphere::wr[rsHeightSize] = { 0.0 };

greal EarthAtmosphere::plp[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::dlp[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::tlp[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::ulp[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::vlp[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::uds[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::vds[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::udl[pcHeightSize][pcLatSize] = { { 0.0 } };
greal EarthAtmosphere::uvt[pcHeightSize][pcLatSize] = { { 0.0 } };

const int EarthAtmosphere::ptData[pcHeightSize][15] = {
  { /* PT   13    0 */ 999,  999,  999,  999,  999,  763,  890,  877,  901,  908,  763,  890,  877,  901,  908 },
  { /* PT   13    5 */ 986,  992,  994,  994,  994,  676,  802,  840,  843,  843,  676,  802,  840,  843,  843 },
  { /* PT   13   10 */ 975,  991,  995,  996,  995,  589,  808,  880,  899,  892,  589,  808,  880,  899,  892 },
  { /* PT   13   15 */ 970,  974,  986,  996,  998,  573,  611,  748,  927,  963,  573,  611,  748,  927,  963 },
  { /* PT   13   20 */ 972,  978,  983,  990,  991,  603,  664,  720,  821,  838,  603,  664,  720,  821,  838 },
  { /* PT   13   25 */ 978,  985,  987,  998,  994,  632,  717,  750,  954,  866,  632,  717,  750,  954,  866 },
  { /* PT   13   30 */ 965,  975,  987,  992,  994,  546,  633,  726,  805,  826,  518,  599,  747,  838,  864 },
  { /* PT   13   35 */ 965,  975,  989,  993,  994,  558,  642,  755,  820,  840,  509,  595,  770,  851,  876 },
  { /* PT   13   40 */ 967,  980,  987,  991,  992,  581,  681,  766,  829,  848,  509,  632,  733,  805,  828 },
  { /* PT   13   45 */ 960,  974,  981,  984,  985,  534,  670,  766,  824,  844,  440,  551,  630,  678,  693 },
  { /* PT   13   50 */ 958,  969,  975,  977,  978,  491,  649,  762,  830,  853,  418,  493,  545,  577,  588 },
  { /* PT   13   55 */ 962,  968,  971,  973,  974,  483,  636,  746,  812,  834,  430,  472,  502,  520,  526 },
  { /* PT   13   60 */ 977,  979,  980,  981,  982,  587,  673,  734,  771,  784,  544,  570,  588,  600,  604 },
  { /* PT   13   65 */ 987,  988,  989,  990,  990,  681,  710,  730,  742,  747,  677,  703,  722,  733,  737 },
  { /* PT   13   70 */ 981,  982,  983,  983,  983,  648,  663,  673,  679,  681,  574,  587,  597,  602,  604 },
  { /* PT   13   75 */ 973,  973,  973,  973,  973,  615,  615,  615,  615,  615,  471,  471,  471,  471,  471 },
  { /* PT   13   80 */ 972,  972,  972,  972,  972,  571,  571,  571,  571,  571,  455,  455,  455,  455,  455 },
  { /* PT   13   85 */ 971,  971,  971,  971,  971,  526,  526,  526,  526,  526,  439,  439,  439,  439,  439 },
  { /* PT   13   90 */ 968,  968,  968,  968,  968,  482,  482,  482,  482,  482,  422,  422,  422,  422,  422 },
  { /* PT   13  100 */ 962,  962,  962,  962,  962,  393,  393,  393,  393,  393,  390,  390,  390,  390,  390 },
  { /* PT   13  120 */ 953,  953,  953,  953,  953,  328,  328,  328,  328,  328,  325,  325,  325,  325,  325 },
  { /* PT   13  140 */ 946,  946,  946,  946,  946,  262,  262,  262,  262,  262,  260,  260,  260,  260,  260 },
  { /* PT   13  160 */ 934,  934,  934,  934,  934,  197,  197,  197,  197,  197,  195,  195,  195,  195,  195 },
  { /* PT   13  180 */ 909,  909,  909,  909,  909,  131,  131,  131,  131,  131,  130,  130,  130,  130,  130 },
  { /* PT   13  200 */ 909,  909,  909,  909,  909,  131,  131,  131,  131,  131,  130,  130,  130,  130,  130 }
};

const int EarthAtmosphere::pwData[pcHeightSize][10] = {
  { /* PW   13    0 */ 961,  911,  866,  873,  828,  320,  606,  921,  693,  733 },
  { /* PW   13    5 */ 941,  904,  869,  847,  840,  321,  599,  892,  721,  759 },
  { /* PW   13   10 */ 922,  897,  873,  857,  852,  323,  591,  864,  748,  785 },
  { /* PW   13   15 */ 903,  890,  876,  867,  864,  324,  583,  835,  775,  810 },
  { /* PW   13   20 */ 884,  883,  880,  877,  876,  327,  575,  806,  803,  836 },
  { /* PW   13   25 */ 864,  878,  887,  892,  894,  344,  575,  739,  839,  872 },
  { /* PW   13   30 */ 854,  872,  885,  892,  895,  324,  552,  751,  832,  859 },
  { /* PW   13   35 */ 821,  858,  885,  901,  906,  304,  545,  786,  890,  918 },
  { /* PW   13   40 */ 794,  854,  896,  921,  930,  328,  550,  706,  944,  974 },
  { /* PW   13   45 */ 798,  853,  899,  927,  937,  349,  537,  618,  920,  943 },
  { /* PW   13   50 */ 768,  831,  876,  903,  912,  339,  489,  650,  856,  889 },
  { /* PW   13   55 */ 734,  761,  781,  793,  797,  366,  450,  649,  771,  811 },
  { /* PW   13   60 */ 721,  686,  662,  647,  642,  469,  550,  674,  749,  773 },
  { /* PW   13   65 */ 743,  720,  704,  694,  690,  601,  713,  794,  842,  858 },
  { /* PW   13   70 */ 765,  753,  745,  740,  739,  694,  750,  790,  814,  822 },
  { /* PW   13   75 */ 787,  787,  787,  787,  787,  787,  787,  787,  787,  787 },
  { /* PW   13   80 */ 760,  760,  760,  760,  760,  760,  760,  760,  760,  760 },
  { /* PW   13   85 */ 734,  734,  734,  734,  734,  734,  734,  734,  734,  734 },
  { /* PW   13   90 */ 707,  707,  707,  707,  707,  707,  707,  707,  707,  707 },
  { /* PW   13  100 */ 654,  654,  654,  654,  654,  654,  654,  654,  654,  654 },
  { /* PW   13  120 */ 615,  615,  615,  615,  615,  615,  615,  615,  615,  615 },
  { /* PW   13  140 */ 575,  575,  575,  575,  575,  575,  575,  575,  575,  575 },
  { /* PW   13  160 */ 536,  536,  536,  536,  536,  536,  536,  536,  536,  536 },
  { /* PW   13  180 */ 496,  496,  496,  496,  496,  496,  496,  496,  496,  496 },
  { /* PW   13  200 */ 457,  457,  457,  457,  457,  457,  457,  457,  457,  457 }
};

const int EarthAtmosphere::csData[pcHeightSize][10] = {
  { /* CS   13    0 */ -55,  -17,  -47,  -65,  -70,  108,  -27, -122, -181, -201 },                                   
  { /* CS   13    5 */ -61,  -22,  -42,  -54,  -58,   82,  -29, -107, -156, -173 },                                   
  { /* CS   13   10 */ -67,  -26,  -38,  -44,  -46,   56,  -32,  -93, -131, -144 },                                   
  { /* CS   13   15 */ -74,  -31,  -33,  -34,  -34,   30,  -34,  -78, -106, -116 },                                   
  { /* CS   13   20 */ -80,  -36,  -28,  -24,  -22,    4,  -36,  -64,  -82,  -88 },                                   
  { /* CS   13   25 */ -86,  -71,  -60,  -54,  -52,    1,  -73, -125, -157, -168 },                                   
  { /* CS   13   30 */-122,  -67,  -28,   -4,    4,  -56,  -75,  -88,  -97, -100 },                                   
  { /* CS   13   35 */-113,  -46,    2,   30,   40,  -75,  -25,   11,   32,   39 },                                   
  { /* CS   13   40 */ -72,   -8,   37,   65,   75, -101,   15,   98,  147,  163 },                                   
  { /* CS   13   45 */ -51,   -9,   22,   40,   46, -152,   -8,   95,  157,  177 },                                   
  { /* CS   13   50 */-107,  -52,  -12,   10,   18, -169,  -38,   55,  111,  129 },                                   
  { /* CS   13   55 */-213,  -89,    0,   52,   70, -183,  -72,    8,   56,   72 },                                   
  { /* CS   13   60 */-102, -108,  -15,   41,   60, -173,  -93,  -37,   -2,    9 },                                   
  { /* CS   13   65 */-137,  -77,   15,   68,   87, -231,  -56,   68,  143,  168 },                                   
  { /* CS   13   70 */-143,  -82,   19,   79,   99, -256,  -58,   82,  168,  196 },                                   
  { /* CS   13   75 */-149,  -86,   24,   89,  112, -283,  -61,   97,  193,  225 },                                   
  { /* CS   13   80 */-155,  -91,   29,   99,  124, -309,  -63,  112,  218,  253 },                                   
  { /* CS   13   85 */-162,  -96,   34,  110,  136, -335,  -65,  126,  243,  281 },                                   
  { /* CS   13   90 */-168, -100,   38,  120,  148, -361,  -67,  141,  268,  310 },                                   
  { /* CS   13  100 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
  { /* CS   13  120 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
  { /* CS   13  140 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
  { /* CS   13  160 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
  { /* CS   13  180 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
  { /* CS   13  200 */-181, -109,   48,  140,  172, -413,  -72,  170,  318,  366 },                                   
};

const int EarthAtmosphere::clData[pcHeightSize][10] = {
  { /* CL   13    0 */ 308, -207, -576, -796, -869,   75,  159,  -61,  376,  313 },
  { /* CL   13    5 */ 289, -202, -553, -763, -832,  -20,  173,   63,   26,    8 },                                   
  { /* CL   13   10 */ 270, -197, -530, -729, -795,  -22,  250,   94,  -21,  -36 },                                   
  { /* CL   13   15 */ 251, -191, -507, -696, -759,  -57,  227,   61,  -94,  -87 },                                   
  { /* CL   13   20 */ 232, -186, -484, -663, -722,  -25,  101,   44,  -57,  -54 },                                   
  { /* CL   13   25 */ 283, -150, -459, -644, -705,  -19,   86,  107,   -8,  -28 },                                   
  { /* CL   13   30 */ 223, -159, -432, -595, -650,  -43,   72,  146,    6,  -25 },                                   
  { /* CL   13   35 */  56, -202, -385, -495, -532,   19,  139,   79,   -6,  -20 },                                   
  { /* CL   13   40 */  96, -183, -382, -502, -542,   42,  115,   79,   -6,  -21 },                                   
  { /* CL   13   45 */ 163, -170, -407, -549, -596,   -4,  127,   62,   -7,  -18 },                                   
  { /* CL   13   50 */ 137, -180, -405, -540, -586,    8,  139,   74,   -7,  -21 },                                   
  { /* CL   13   55 */ 126, -144, -335, -450, -489,   68,  114,   91,   -8,  -24 },                                   
  { /* CL   13   60 */  79, -108, -240, -320, -347,   23,  147,   85,   -8,  -24 },                                   
  { /* CL   13   65 */  60, -138, -277, -361, -390,   32,  130,   81,   -9,  -23 },                                   
  { /* CL   13   70 */  41, -133, -254, -328, -353,  -77,   35,  -21,   -9,   -3 },                                   
  { /* CL   13   75 */  22, -127, -231, -294, -316,    5,   81,   43,  -10,  -16 },                                   
  { /* CL   13   80 */   3, -122, -208, -261, -280,    6,   75,   40,  -10,  -16 },                                   
  { /* CL   13   85 */ -16, -117, -185, -227, -243,    7,   68,   37,  -11,  -16 },                                   
  { /* CL   13   90 */ -36, -111, -162, -194, -206,    7,   61,   34,  -11,  -16 },                                   
  { /* CL   13  100 */ -74, -101, -116, -127, -132,    9,   48,   28,  -12,  -16 },                                   
  { /* CL   13  120 */ -74, -101, -116, -127, -132,   11,   22,   17,  -14,  -15 },                                   
  { /* CL   13  140 */ -74, -101, -116, -127, -132,   14,   -5,    5,  -16,  -14 },                                   
  { /* CL   13  160 */ -74, -101, -116, -127, -132,   17,  -31,   -7,  -18,  -13 },                                   
  { /* CL   13  180 */ -74, -101, -116, -127, -132,   19,  -58,  -19,  -20,  -12 },                                   
  { /* CL   13  200 */ -74, -101, -116, -127, -132,   22,  -84,  -31,  -22,  -12 }                                  
};

EarthAtmosphere::RsData EarthAtmosphere::rsData[rsHeightSize] = {
  { /* RS   0. */   7.8,    2.8,   0.52,  117.5,   0.53,   0.08,   0.32,   2.65,   0.63 },
  { /* RS   5. */  56.9,   20.3,   1.04,  819.9,   2.20,   0.61,   0.52,  10.56,   0.63 },
  { /* RS  10. */  59.5,   20.4,   1.23,  820.8,   2.68,   0.70,   0.69,  12.32,   0.63 },
  { /* RS  15. */ 106.0,   34.3,   3.11, 1397.2,   3.26,   0.75,   1.01,  14.36,   0.51 },
  { /* RS  20. */ 161.9,   48.7,   8.64, 2040.4,   3.81,   0.78,   1.36,  16.01,   0.42 },
  { /* RS  25. */ 177.3,   52.3,  12.00, 2127.5,   4.04,   0.78,   1.56,  16.16,   0.45 },
  { /* RS  30. */ 196.5,   53.8,  28.60, 2240.3,   4.59,   0.80,   2.10,  17.44,   0.53 },
  { /* RS  35. */ 211.1,   58.6,  35.40, 2279.4,   5.09,   0.87,   2.47,  18.33,   0.63 },
  { /* RS  40. */ 226.3,   62.5,  42.60, 2308.6,   5.52,   0.92,   2.82,  18.78,   0.75 },
  { /* RS  45. */ 241.5,   66.5,  50.10, 2318.4,   5.97,   0.96,   3.20,  19.09,   0.89 },
  { /* RS  50. */ 255.6,   69.9,  57.90, 2300.4,   6.44,   1.00,   3.62,  19.32,   1.15 },
  { /* RS  55. */ 266.7,   74.1,  66.00, 2240.5,   6.92,   1.07,   4.03,  19.39,   1.42 },
  { /* RS  60. */ 277.4,   77.0,  74.40, 2164.0,   7.43,   1.13,   4.46,  19.32,   1.61 },
  { /* RS  65. */ 292.7,   84.3,  83.20, 2107.6,   7.95,   1.20,   4.97,  19.08,   1.54 },
  { /* RS  70. */ 308.9,   90.2,  92.30, 2038.9,   8.47,   1.24,   5.49,  18.63,   1.91 },
  { /* RS  75. */ 325.2,   94.5, 102.00, 1951.2,   9.00,   1.25,   6.04,  18.00,   2.51 },
  { /* RS  80. */ 340.5,   97.9, 111.00, 2043.0,   9.52,   1.25,   6.58,  19.04,   3.30 },
  { /* RS  85. */ 356.8,  103.9, 121.00, 2140.8,  10.04,   1.27,   7.15,  20.08,   4.09 },
  { /* RS  90. */ 373.2,  106.2, 132.00, 2239.2,  11.02,   1.30,   8.07,  22.04,   4.75 },
  { /* RS  95. */ 419.2,  120.5, 142.50, 2515.5,  11.99,   1.35,   8.89,  23.99,   5.01 },
  { /* RS 100. */ 465.3,  134.0, 153.00, 2791.8,  12.97,   1.40,   9.70,  25.94,   5.26 },
  { /* RS 105. */ 511.5,  145.3, 164.75, 3068.8,  13.72,   1.47,  10.22,  27.44,   5.61 },
  { /* RS 110. */ 557.7,  156.6, 176.50, 3345.9,  14.47,   1.53,  10.75,  28.95,   5.95 },
  { /* RS 115. */ 603.8,  166.2, 188.25, 3622.9,  15.23,   1.58,  11.27,  30.45,   6.29 },
  { /* RS 120. */ 650.0,  174.8, 200.00, 3900.0,  15.98,   1.63,  11.79,  31.96,   6.63 },
  { /* RS 140. */ 698.2,  181.0, 232.00, 4189.2,  17.96,   1.75,  13.46,  35.92,   8.33 },
  { /* RS 160. */ 747.0,  185.7, 268.92, 4482.0,  20.00,   1.85,  15.23,  40.00,  10.19 },
  { /* RS 180. */ 795.0,  197.6, 286.20, 4770.0,  22.02,   2.04,  16.77,  44.04,  11.99 },
  { /* RS 200. */ 840.0,  209.7, 300.00, 5040.0,  24.03,   2.24,  18.27,  48.06,  12.01 }
};

//! \brief Copy tabular data into height/latitude model data fields.
//!
//! This routine transforms the static tabular data into height and latitude based model data.
//! Initialization of static data need only be performed once.
//!
//! \b Inputs:
//! \arg #ptData
//! \arg #pwData
//! \arg #csData
//! \arg #clData
//! \arg #rsData
//!
//! \returns Static data is populated.
void EarthAtmosphere::initializeData()
{
  if (initialized) {
    return;
  }
  initialized = true;

  for (size_t iHeight = 0; iHeight < pcHeightSize; iHeight++) {
    for (size_t i = 0; i < 5; i++) {
      plp[iHeight][i + 5] = plp[iHeight][4 - i] = greal(ptData[iHeight][i] / 1000.0);
      dlp[iHeight][i + 5] = dlp[iHeight][4 - i] = greal(ptData[iHeight][i + 5] / 1000.0);
      tlp[iHeight][i + 5] = tlp[iHeight][4 - i] = greal(ptData[iHeight][i + 10] / 1000.0);
    }
  }

  for (size_t iHeight = 0; iHeight < pcHeightSize; iHeight++) {
    for (size_t j = 0; j < 5; j++) {
      // for ulp, vlp assume large-scale fraction for u equal that for v
      ulp[iHeight][j + 5] = ulp[iHeight][4 - j] = greal(pwData[iHeight][j + 5]) / 1000.0; // no typo, see comment above
      vlp[iHeight][j + 5] = vlp[iHeight][4 - j] = greal(pwData[iHeight][j + 5]) / 1000.0;

      uds[iHeight][j + 5] = uds[iHeight][4 - j] = greal(csData[iHeight][j]) / 1000.0;
      vds[iHeight][j + 5] = vds[iHeight][4 - j] = greal(csData[iHeight][j + 5]) / 1000.0;

      udl[iHeight][j + 5] = udl[iHeight][4 - j] = greal(clData[iHeight][j]) / 1000.0;
      uvt[iHeight][j + 5] = uvt[iHeight][4 - j] = greal(clData[iHeight][j + 5]) / 1000.0;
    }
  }

  for (size_t i = 0; i < rsHeightSize; ++i) {
    xlbar[i] = rsData[i].xlbar;       // xMeanSmall
    xsigl[i] = rsData[i].xsigl;       // xSDSmall
    xlmin[i] = rsData[i].xlmin;       // xMin
    xscale[i] = rsData[i].xscale;     // xScaleSmall
    zlbar[i] = rsData[i].zlbar;       // zMeanSmall
    zsigl[i] = rsData[i].zsigl;       // zSDSmall
    zlmin[i] = rsData[i].zlmin;       // zMin
    zscale[i] = rsData[i].zscale;     // zScaleSmall
    wr[i] = rsData[i].wr;             // verticalStandardDeviation
  }

  initializeTopographyData();
}

//! \brief Reads topography and land type data and stores in static arrays.
//!
//! This routine read the topogography data file (topo.txt) into longitude/latitude based model data.
//! Initialization of static data need only be performed once.
//!
//! \b Inputs:
//! \arg #atmPath
//!
//! \returns Static data, #landcd and #ztopo, is populated.
void EarthAtmosphere::initializeTopographyData()
{
  // Open topography file
  string fileName = atmPath + "/topo.txt";
  ifstream tops;
  tops.open(fileName);
  if (!tops) {
    throw string(FILE_OPEN_ERROR_MESSAGE + fileName);
  }

  try {
    string line;
    // Read topography and land code data and store in arrays
    for (int latIndex = 179; latIndex >= 0; latIndex--) {
      for (size_t lonIndex = 0; lonIndex < 360; lonIndex++) {
        // Going to optimise this a little.  Each line input is equivalent to:
        // tops >> xlo >> xla >> landcd[i][j] >> ltopom;
        // But we will skip parsing the unused data.

        // Read a line.
        getline(tops, line);

        // Get the land code starting at position 13.
        size_t pos = 13;
        landcd[lonIndex][latIndex] = getNextInt(line, pos);

        // Get the height in meters. And convert to km.
        ztopo[lonIndex][latIndex] = getNextInt(line, pos) / greal(1000.0);
      }
    }
  }
  catch (const ios_base::failure& e) {
    throw(FILE_READ_TEXT_ERROR_MESSAGE + fileName
      + "\n       " + e.what());
  }

  tops.close();
}


} // namespace