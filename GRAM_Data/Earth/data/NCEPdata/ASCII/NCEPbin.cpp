// Program to convert ASCII files to binary values, write out
// data for map plots and/or profile plots, and write binary
// version NCEP data arrays

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

// These need to be global to avoid stack overflow.
double avg[145][73][19][5][10];
double sdv[145][73][19][5][10];
double hatsfc[145][73][5];
double geop[145][73][19];
double sprs[145][73][19];

int getInt(const string& input, size_t width, size_t &pos);
void sortlevel(int lsort[19], double zsort[19]);
void writeBlock3D(ofstream& oFile, double data[145][73][19][5][10], int hr, int var);
void writeBlock3D(ofstream& oFile, double data[145][73][19]);
void writeBlock2D(ofstream& oFile, double data[145][73][19][5][10], int lvl, int hr, int var);

int main(int argc, char** argv)
{
  // Indexes for averages (avg) and standard deviations (sdv):                       
  //    0-143 = Lon: 0E to 357.5E (step 2.5) with 144 = 360E = 0 E     
  //    0-72  = Lat: 90S to 90 N (step 2.5)                                          
  //    0-18  = Pressure level:  0=sea level, 1=surface, 2=1000mb, ... 18=10mb       
  //    0-4   = Code for hour of day: 
  //               User input: 1=0 UTC, 2=6 UTC,  3=12UTC,  4=18 UTC, 5=All-day average                                   
  //    0-9  = Variable Number:   
  //               0=Geopotential height (m) for pressure level, or surface pressure 
  //                   (mb) for surface level - height later converted to km         
  //               1=Temperature (K)                                                 
  //               2=Eastward (U) wind component (m/s)                               
  //               3=Northward (V) wind component (m/s)                              
  //               4=Relative humidity (percent) - computed if above 300mb height    
  //               5=Density (kg/m**3) - computed from gas law                       
  //               6=Dewpoint temperature (K) - computed from tdbuck routine         
  //               7=Vapor pressure (Pa) - computed from wexler routine              
  //               8=Wind speed (m/s) - computed from U and V components             
  //               9=U*V product (for computing U-V cross correlation)               
  //                                                                                 
  //-----------------------------------------------

  double p[19], den[19], temp[19], u[19], v[19], rh[19], sp[19], sden[19], stemp[19], su[19],
    sv[19], srh[19], h[19], zsort[19], tdew[19], stdew[19], vp[19], svp[19], SpA[19], SpS[19], Ruv[19];

  // Pressures (mb) at levels 0-18.  Values for sea level (level=0)
  // and surface (level=1) are determined later.
  double  pmb[19] = { 9999.0, 9999.0, 1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0,
                      250.0, 200.0, 150.0, 100.0, 70.0, 50.0, 30.0, 20.0, 10.0 };

  // Reference densities (kg/m**3). Values for sea level (level=0)
  // and surface (level=1) are determined later.
  double denref[19] = { 9999.0, 9999.0, 1.2121, 1.1381, 1.0626,
                       0.90798, 0.80144, 0.69144, 0.57714, 0.45721,
                       0.39446, 0.32159, 0.24119, 0.16079, 0.11255,
                       0.08018, 0.04739, 0.03122, 0.01529 };

  string por;
  cout << " Enter the POR and month (yyzzmm): ";
  cin >> por;
  string inFileName = "Nf" + por + ".txt";
  cout << " Reading NCEP ASCII file = " << inFileName << endl;

  // Open and read Fixed ASCII file (Nf) of NCEP averages and standard deviations

  ifstream inFile(inFileName);
  string lineInput;
  getline(inFile, lineInput);

  for (int hr = 0; hr < 4; ++hr) {
    for (int lon = 0; lon < 144; ++lon) {
      for (int lat = 0; lat < 73; ++lat) {
        for (int lvl = 1; lvl < 19; ++lvl) {
          getline(inFile, lineInput);
          size_t pos = 0;
          int ihr = getInt(lineInput, 1, pos);
          int ilvl = getInt(lineInput, 2, pos);
          int ilat = getInt(lineInput, 2, pos);
          int ilon = getInt(lineInput, 3, pos);
          int nv = getInt(lineInput, 4, pos);
          int ihavg = getInt(lineInput, 6, pos);
          int ihsig = getInt(lineInput, 5, pos);

          int iuavg = getInt(lineInput, 4, pos);
          int iusig = getInt(lineInput, 3, pos);
          int ivavg = getInt(lineInput, 4, pos);
          int ivsig = getInt(lineInput, 3, pos);
          int itavg = getInt(lineInput, 4, pos);
          int itsig = getInt(lineInput, 3, pos);
          int idavg = getInt(lineInput, 4, pos);
          int idsig = getInt(lineInput, 3, pos);

          int nrho = getInt(lineInput, 1, pos);
          int irhavg = getInt(lineInput, 4, pos);
          int irhsig = getInt(lineInput, 3, pos);
          int itdavg = getInt(lineInput, 4, pos);
          int itdsig = getInt(lineInput, 3, pos);
          int ivpavg = getInt(lineInput, 4, pos);
          int ivpsig = getInt(lineInput, 4, pos);

          int nvp = getInt(lineInput, 1, pos);
          int ISavg = getInt(lineInput, 4, pos);
          int ISsig = getInt(lineInput, 3, pos);
          int IRuv = getInt(lineInput, 4, pos);


          if (ihr != hr + 1) {
            cout << " Bad hour" << ihr << ", " << ilvl << ", " << ilat << ", " << ilon;
            exit(0);
          }
          else if (ilvl != lvl) {
            cout << " Bad level" << ihr << ", " << ilvl << ", " << ilat << ", " << ilon;
            exit(0);
          }
          else if (ilat != lat + 1) {
            cout << " Bad lat" << ihr << ", " << ilvl << ", " << ilat << ", " << ilon;
            exit(0);
          }
          else if (ilon != lon + 1) {
            cout << " Bad lon" << ihr << ", " << ilvl << ", " << ilat << ", " << ilon;
            exit(0);
          }
          else {
            // Convert ASCII numbers to correct units:
            //    H = geopotential height (m)
            //    T = temperature (K)
            //    U, V = eastward, northward wind components (m/s)
            //    RH = relative humidity (%)
            //    D = density (g/m**3)
            //    Td = dewpoint (K)
            //    vp = vapor pressure (N/m**2)
            avg[lon][lat][lvl][hr][0] = ihavg / 10.0;
            sdv[lon][lat][lvl][hr][0] = ihsig / 10.0;
            avg[lon][lat][lvl][hr][1] = itavg / 10.0;
            sdv[lon][lat][lvl][hr][1] = itsig / 10.0;
            avg[lon][lat][lvl][hr][2] = iuavg / 10.0;
            sdv[lon][lat][lvl][hr][2] = iusig / 10.0;
            avg[lon][lat][lvl][hr][3] = ivavg / 10.0;
            sdv[lon][lat][lvl][hr][3] = ivsig / 10.0;
            avg[lon][lat][lvl][hr][4] = irhavg / 10.0;
            sdv[lon][lat][lvl][hr][4] = irhsig / 10.0;
            avg[lon][lat][lvl][hr][5] = idavg * pow(10.0, -nrho);
            sdv[lon][lat][lvl][hr][5] = idsig * pow(10.0, -nrho);
            avg[lon][lat][lvl][hr][6] = itdavg / 10.0;
            sdv[lon][lat][lvl][hr][6] = itdsig / 10.0;
            avg[lon][lat][lvl][hr][7] = ivpavg * pow(10.0, -nvp);
            sdv[lon][lat][lvl][hr][7] = ivpsig * pow(10.0, -nvp);
            avg[lon][lat][lvl][hr][8] = ISavg / 10.0;
            sdv[lon][lat][lvl][hr][8] = ISsig / 10.0;
            avg[lon][lat][lvl][hr][9] = IRuv / 1000.0;
            sdv[lon][lat][lvl][hr][9] = 0.0;
          }
        }
      }
    }
  }
  inFile.close();

  // START OF test and data conversion section
  double rlim = 0.999;
  double gref = 9.80665;

  for (int lon = 0; lon < 144; ++lon) {
    for (int lat = 0; lat < 73; ++lat) {
      for (int hr = 0; hr < 4; ++hr) {
        // Convert units and store profile values in single-dimension arrays
        for (int lvl = 1; lvl < 19; ++lvl) {
          h[lvl]     = avg[lon][lat][lvl][hr][0];
          // sh[lvl] = sdv[lon][lat][lvl][hr][0];
          temp[lvl]  = avg[lon][lat][lvl][hr][1];
          stemp[lvl] = sdv[lon][lat][lvl][hr][1];
          u[lvl]     = avg[lon][lat][lvl][hr][2];
          su[lvl]    = sdv[lon][lat][lvl][hr][2];
          v[lvl]     = avg[lon][lat][lvl][hr][3];
          sv[lvl]    = sdv[lon][lat][lvl][hr][3];
          rh[lvl]    = avg[lon][lat][lvl][hr][4];
          srh[lvl]   = sdv[lon][lat][lvl][hr][4];
          // Convert density to kg/m**3
          den[lvl]   = avg[lon][lat][lvl][hr][5] / 1000.0;
          sden[lvl]  = sdv[lon][lat][lvl][hr][5] / 1000.0;
          tdew[lvl]  = avg[lon][lat][lvl][hr][6];
          stdew[lvl] = sdv[lon][lat][lvl][hr][6];
          vp[lvl]    = avg[lon][lat][lvl][hr][7];
          svp[lvl]   = sdv[lon][lat][lvl][hr][7];
          SpA[lvl]   = avg[lon][lat][lvl][hr][8];
          SpS[lvl]   = sdv[lon][lat][lvl][hr][8];
          Ruv[lvl]   = avg[lon][lat][lvl][hr][9];
        }
        // Find and store level index values for surface (k=1) and sea
        // level (k=0).
        int k1slp, k2slp, k1sfc, k2sfc;
        double sprat;
        for (int lvl = 1; lvl < 19; ++lvl) {
          if (lvl == 1) {
            h[0] = 0.0;
            // Convert pressure to N/m**2
            p[1]  = avg[lon][lat][lvl][hr][0] * 100.0;
            sp[1] = sdv[lon][lat][lvl][hr][0] * 100.0;
            sprat = sp[1] / p[1];
            // Use levels 2 (1000 mb) and 3 (925 mb) for interpolation
            // to sea level, if level 2 height is positive
            if (h[2] > 0.0) {
              k1slp = 2;
              k2slp = 3;
            }
            // Use levels 2 (1000 mb) and 3 (925 mb) for interpolation
            // to surface, if level 1 (surface) pressure is greater than
            // 1000 mb.
            if (p[1] > pmb[2] * 100.0) {
              k1sfc = 2;
              k2sfc = 3;
            }
          }
          else {
            // Convert pressures to N/m**2
            p[lvl] = pmb[lvl] * 100.0;
            // Convert standard deviation in geopotential height to
            // pressure standard deviation
            sp[lvl] = gref * sdv[lon][lat][lvl][hr][0] * den[lvl];
            sprat = sp[lvl] / p[lvl];
            // Find level index values below and above sea level, for
            // interpolation
            if (h[lvl] <= 0.0 && h[lvl + 1] > 0.0) {
              k1slp = lvl;
              k2slp = lvl + 1;
            }
            // Find level index values below and above surface, for
            // interpolation
            if (lvl < 18 && p[1] / 100.0 <= pmb[lvl] && p[1] > pmb[lvl + 1]) {
              k1sfc = lvl;
              k2sfc = lvl + 1;
            }
          }
          // Correct sigmas using correlation constraints, based on the
          // Buell relationships
          double sdrat = sden[lvl] / den[lvl];
          double strat = stemp[lvl] / temp[lvl];
          if (sdrat >= rlim * (sprat + strat)) {
            // Correct density standard deviation
            sden[lvl] = rlim * den[lvl] * (sprat + strat);
          }
          else if (strat >= rlim * (sprat + sdrat)) {
            // Correct temperature standard deviation
            stemp[lvl] = rlim * temp[lvl] * (sprat + sdrat);
          }
          else if (sprat >= rlim * (sdrat + strat)) {
            // Correct height (or surface pressure) standard deviation
            sp[lvl] = rlim * p[lvl] * (sdrat + strat);
          }
        }
        // Interpolate for surface geopotential height
        double hpsfc = (h[k2sfc] - h[k1sfc]) / log(p[k1sfc] / p[k2sfc]);
        double hsfc  = h[k1sfc] + hpsfc * log(p[k1sfc] / p[1]);
        hatsfc[lon][lat][hr] = hsfc / 1000.0;
        // Interpolate (or extrapolate) for sea level conditions
        double hpslp    = (h[k2slp] - h[k1slp]) / log(p[k1slp] / p[k2slp]);
        double pslp     = p[k1slp] * exp(h[k1slp] / hpslp);
        double spslp    = pslp * sp[k1slp] / p[k1slp];
        double factint  = h[k1slp] / (h[k1slp] - h[k2slp]);
        temp[0]         = temp[k1slp] + (temp[k2slp] - temp[k1slp]) * factint;
        stemp[0]        = stemp[k1slp] + (stemp[k2slp] - stemp[k1slp]) * factint;
        u[0]            = u[k1slp] + (u[k2slp] - u[k1slp]) * factint;
        su[0]           = su[k1slp] + (su[k2slp] - su[k1slp]) * factint;
        v[0]            = v[k1slp] + (v[k2slp] - v[k1slp]) * factint;
        sv[0]           = sv[k1slp] + (sv[k2slp] - sv[k1slp]) * factint;
        rh[0]           = rh[k1slp] + (rh[k2slp] - rh[k1slp]) * factint;
        srh[0]          = srh[k1slp] + (srh[k2slp] - srh[k1slp]) * factint;
        double rgas1    = p[k1slp] / (temp[k1slp] * den[k1slp]);
        double rgas2    = p[k2slp] / (temp[k2slp] * den[k2slp]);
        double rgas     = rgas1 + (rgas2 - rgas1) * factint;
        den[0]          = pslp / (rgas * temp[0]);
        sden[0]         = sden[k1slp] + (sden[k2slp] - sden[k1slp]) * factint;
        tdew[0]         = tdew[k1slp] + (tdew[k2slp] - tdew[k1slp]) * factint;
        stdew[0]        = stdew[k1slp] + (stdew[k2slp] - stdew[k1slp]) * factint;
        double hvp;
        if (vp[k2slp] == 0.0 || vp[k2slp] - vp[k1slp] == 0.0) {
          hvp = 1.0e10;
        }
        else {
          hvp = (h[k2slp] - h[k1slp]) / log(vp[k1slp] / vp[k2slp]);
        }
        vp[0]  = vp[k1slp] * exp(h[k1slp] / hvp);
        svp[0] = vp[0] * svp[k1slp] / vp[k1slp];
        SpA[0] = SpA[k1slp] + (SpA[k2slp] - SpA[k1slp]) * factint;
        SpS[0] = SpS[k1slp] + (SpS[k2slp] - SpS[k1slp]) * factint;
        Ruv[0] = Ruv[k1slp] + (Ruv[k2slp] - Ruv[k1slp]) * factint;
        // Store sea level parameters in Avg and Std arrays as level 0
        avg[lon][lat][0][hr][0] = pslp;
        sdv[lon][lat][0][hr][0] = spslp;
        avg[lon][lat][0][hr][1] = temp[0];
        sdv[lon][lat][0][hr][1] = stemp[0];
        avg[lon][lat][0][hr][2] = u[0];
        sdv[lon][lat][0][hr][2] = su[0];
        avg[lon][lat][0][hr][3] = v[0];
        sdv[lon][lat][0][hr][3] = sv[0];
        avg[lon][lat][0][hr][4] = rh[0];
        sdv[lon][lat][0][hr][4] = srh[0];
        avg[lon][lat][0][hr][5] = den[0];
        sdv[lon][lat][0][hr][5] = sden[0];
        avg[lon][lat][0][hr][6] = tdew[0];
        sdv[lon][lat][0][hr][6] = stdew[0];
        avg[lon][lat][0][hr][7] = vp[0];
        sdv[lon][lat][0][hr][7] = svp[0];
        avg[lon][lat][0][hr][8] = SpA[0];
        sdv[lon][lat][0][hr][8] = SpS[0];
        avg[lon][lat][0][hr][9] = Ruv[0];
        sdv[lon][lat][0][hr][9] = 0.0;
        // Store Density as kg/m**3, pressure as N/m**2, height as km
        for (int lvl = 1; lvl < 19; ++lvl) {
          avg[lon][lat][lvl][hr][5] = den[lvl];
          sdv[lon][lat][lvl][hr][5] = sden[lvl];
          if (lvl == 1) {
            avg[lon][lat][lvl][hr][0] = p[1];
            sdv[lon][lat][lvl][hr][0] = sp[1];
          }
          else {
            // Convert geopotential height to km 
            avg[lon][lat][lvl][hr][0] = h[lvl] / 1000.0;
            sdv[lon][lat][lvl][hr][0] = sp[lvl];
          }
        }
      }

      // Compute averages and standard deviations for all day (hr=5)
      double sum1 = 0;
      for (int hr = 0; hr < 4; ++hr) {
        sum1 += hatsfc[lon][lat][hr];  
      }
      hatsfc[lon][lat][4] = sum1 / 4.0;

      for (int lvl = 0; lvl < 19; ++lvl) {
        for (int ivar = 0; ivar < 10; ++ivar) {

          sum1 = 0;
          for (int hr = 0; hr < 4; ++hr) {
            sum1 += avg[lon][lat][lvl][hr][ivar] / 4.0;
          }
          avg[lon][lat][lvl][4][ivar] = sum1;

          sum1 = 0;
          for (int hr = 0; hr < 4; ++hr) {
            sum1 += pow(sdv[lon][lat][lvl][hr][ivar], 2) / 4.0
              + (pow(avg[lon][lat][lvl][hr][ivar], 2)
              - pow(avg[lon][lat][lvl][4][ivar], 2)) / 4.0;
          }

          sdv[lon][lat][lvl][4][ivar] = sqrt(abs(sum1));
        }
      }
    }
  }

  // Copy longitude=0 values to longitude=360
  for (int lat = 0; lat < 73; ++lat) {
    for (int hr = 0; hr < 5; ++hr) {
      hatsfc[144][lat][hr] = hatsfc[0][lat][hr];
    }
  }
  for (int lat = 0; lat < 73; ++lat) {
    for (int lvl = 0; lvl < 19; ++lvl) {
      for (int hr = 0; hr < 5; ++hr) {
        for (int ivar = 0; ivar < 10; ++ivar) {
          avg[144][lat][lvl][hr][ivar] = avg[0][lat][lvl][hr][ivar];
          sdv[144][lat][lvl][hr][ivar] = sdv[0][lat][lvl][hr][ivar];
        }
      }
    }
  }
  // END OF test and data conversion section

  // Begin section to write out lon-lat maps by level
  cout << " Enter 1 for Maps, 2 for Maps & Trajectory files, 0 for none: ";
  int inwr;
  cin >> inwr;
  if (inwr >= 1) {
    while (true) {
      cout << " Enter indexes for lvl(0-18), hr(1-5), or 0s to end: ";
      int lvl, hr;
      cin >> lvl;
      cin.ignore(100, ',');
      cin >> hr;
      if (hr == 0) {
        break;
      }
      hr -= 1;
      if (hr > 5) hr = 5;
      if (lvl > 18) lvl = 18;
      string olvl;
      if (lvl < 10) olvl = "0";
      olvl += to_string(lvl);
      string fileId = por + olvl + to_string(hr+1) + ".txt";
      cout << " Writing M" << fileId << endl;
      ofstream mtFile;
      if (inwr > 1) {
        cout << " Writing MT" << fileId << endl;
        mtFile.open("MT" + fileId);
      }
      ofstream mFile("M" + fileId);
      if (lvl > 1) {
        mFile << "   I  J   Lon   Lat  Havg_m  Psigmb  Tavg  Tsig"
          "  Uavg  Usig  Vavg  Vsig RHavg RHsig Dav76 Dsd76"
          " TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv\n";
      }
      else {
        mFile << "   I  J   Lon   Lat  Pavgmb  Psigmb  Tavg  Tsig"
          "  Uavg  Usig  Vavg  Vsig RHavg RHsig  Davg  Dstd"
          " TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv Hatsfc\n";
      }
      if (lvl > 1) {
        for (int lon = 0; lon < 145; ++lon) {
          for (int lat = 0; lat < 73; ++lat) {
            // Normalize density by reference density if not at
            // surface or sea level
            // If at surface or sea level, write density avg and  
            // sdv in kg/m**3                                    
            char buffer[1000];
            snprintf(buffer, 1000, "%4d%3d%6.1f%6.1f%8.1f%8.2f%6.1f%6.1f%6.1f%6.1f"
              "%6.1f%6.1f%6.1f%6.1f%6.3f%6.3f%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%7.3f\n",
              lon + 1, lat + 1,
              2.5 * lon,
              -90.0 + 2.5 * lat,
              avg[lon][lat][lvl][hr][0] * 1000.0,
              sdv[lon][lat][lvl][hr][0] / 100.0,
              avg[lon][lat][lvl][hr][1],
              sdv[lon][lat][lvl][hr][1],
              avg[lon][lat][lvl][hr][2],
              sdv[lon][lat][lvl][hr][2],
              avg[lon][lat][lvl][hr][3],
              sdv[lon][lat][lvl][hr][3],
              avg[lon][lat][lvl][hr][4],
              sdv[lon][lat][lvl][hr][4],
              avg[lon][lat][lvl][hr][5] / denref[lvl],
              sdv[lon][lat][lvl][hr][5] / denref[lvl],
              avg[lon][lat][lvl][hr][6],
              sdv[lon][lat][lvl][hr][6],
              avg[lon][lat][lvl][hr][7],
              sdv[lon][lat][lvl][hr][7],
              avg[lon][lat][lvl][hr][8],
              sdv[lon][lat][lvl][hr][8],
              avg[lon][lat][lvl][hr][9]);
            mFile << buffer;
            if (inwr > 1) {
              // Approximate geometric height from geopotential height(km)
              double z = avg[lon][lat][lvl][hr][0];
              z = z / (1.0 - z / 6356.0);
              // Write to Map trajectory file 
              buffer[0] = 0;
              snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                0.0, z, -90. + 2.5 * lat, 2.5 * lon);
              mtFile << buffer;
            }
          }
        }
      }
      else {
        for (int lon = 0; lon < 145; ++lon) {
          // If at surface or sea level, write density avg and
          // sdv in kg/m**3
          for (int lat = 0; lat < 73; ++lat) {
            char buffer[1000];
            snprintf(buffer, 1000, "%4d%3d%6.1f%6.1f%8.1f%8.2f%6.1f%6.1f%6.1f%6.1f"
              "%6.1f%6.1f%6.1f%6.1f%6.3f%6.3f%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%7.3f%7.3f\n",
              lon + 1, lat + 1,
              2.5 * lon,
              -90.0 + 2.5 * lat,
              avg[lon][lat][lvl][hr][0] / 100.0,
              sdv[lon][lat][lvl][hr][0] / 100.0,
              avg[lon][lat][lvl][hr][1],
              sdv[lon][lat][lvl][hr][1],
              avg[lon][lat][lvl][hr][2],
              sdv[lon][lat][lvl][hr][2],
              avg[lon][lat][lvl][hr][3],
              sdv[lon][lat][lvl][hr][3],
              avg[lon][lat][lvl][hr][4],
              sdv[lon][lat][lvl][hr][4],
              avg[lon][lat][lvl][hr][5],
              sdv[lon][lat][lvl][hr][5],
              avg[lon][lat][lvl][hr][6],
              sdv[lon][lat][lvl][hr][6],
              avg[lon][lat][lvl][hr][7],
              sdv[lon][lat][lvl][hr][7],
              avg[lon][lat][lvl][hr][8],
              sdv[lon][lat][lvl][hr][8],
              avg[lon][lat][lvl][hr][9],
              hatsfc[lon][lat][hr]);
            mFile << buffer;

            if (inwr > 1) {
              buffer[0] = 0;
              if (lvl == 0) {
                snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                  0.0, 0.0, -90. + 2.5 * lat, 2.5 * lon);
              }
              else {
                snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                  0.0, hatsfc[lon][lat][hr], -90. + 2.5 * lat, 2.5 * lon);
              }
              mtFile << buffer;
            }
          }
        }
      }
      mFile.close();
      if (inwr > 1) {
        mtFile.close();
      }
    }
  }
  // End of section to write maps

  // Begin section to write vertical profiles at specified lon-lat
  cout << " Enter 1 for Profiles, 2 for Profiles & Trajectory files, or 0 for none: ";
  cin >> inwr;
  if (inwr >= 1) {
    while (true) {
      cout << " Enter lon(deg E),lat(deg),hr(1-5), or hr=0 to end: ";
      int hr;
      double xlon, xlat;
      cin >> xlon;
      cin.ignore(100, ',');
      cin >> xlat;
      cin.ignore(100, ',');
      cin >> hr;
      if (hr == 0) break;
      hr -= 1;
      if (hr > 5) hr = 5;
      if (abs(xlat) > 90.0 || abs(xlon) > 360.0) {
        cout << " Bad Lat-Lon" << endl;
        break;
      }
      int lat = (int)round((xlat + 90.0) / 2.5);
      double zlon = xlon;
      if (zlon < 0.0) zlon = zlon + 360.0;
      int lon = (int)round(zlon / 2.5);
      if (hr > 0) {
        string olat = to_string(lat + 1);
        while (olat.length() < 2) {
          olat = "0" + olat;
        }
        string olon = to_string(lon + 1);
        while (olon.length() < 3) {
          olon = "0" + olon;
        }
        string fileId = por + olon + olat + to_string(hr + 1) + ".txt";
        cout << "Writing P" << fileId << endl;
        ofstream pFile("P" + fileId); 
        ofstream ptFile;
        if (inwr > 1) {
          cout << "Writing PT" << fileId << endl;
          ptFile.open("PT" + fileId); 
        }
        char buffer[1000];
        snprintf(buffer, 1000, "  L  Havgkm  Pavgmb  Psigmb  Tavg  Tsig  "
          "Uavg  Usig  Vavg  Vsig RHavg RHsig      Davg      Dsig"
          " TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv I=%3d"
          " J=%2d Lon=%6.1f Lat=%5.1f\n",
          lon + 1, lat + 1, xlon, xlat);
        pFile << buffer;
        zsort[0] = 0.0;
        zsort[1] = hatsfc[lon][lat][hr];
        for (int lvl = 2; lvl < 19; ++lvl) {
          zsort[lvl] = avg[lon][lat][lvl][hr][0];
        }
        // Sort levels to put surface and sea level in right order
        int lsort[19];
        sortlevel(lsort, zsort);
        // Write out profile levels
        for (int k = 0; k < 19; ++k) {
          int lvl = lsort[k];
          // Write sea level values
          if (lvl == 0) {
            // Compute reference density for sea level, and
            // write out sea level values
            double h3 = avg[lon][lat][3][hr][0];
            double h2 = avg[lon][lat][2][hr][0];
            double hden = (h3 - h2) / log(denref[2] / denref[3]);
            double denref0 = denref[2] * exp(h2 / hden);
            snprintf(buffer, 1000, "%3d%8.3f%8.2f%8.2f%6.1f%6.1f%6.1f%6.1f%6.1f"
              "%6.1f%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%7.3f\n",
              0, 0.0,
              avg[lon][lat][0][hr][0] / 100.0,
              sdv[lon][lat][0][hr][0] / 100.,
              avg[lon][lat][0][hr][1],
              sdv[lon][lat][0][hr][1],
              avg[lon][lat][0][hr][2],
              sdv[lon][lat][0][hr][2],
              avg[lon][lat][0][hr][3],
              sdv[lon][lat][0][hr][3],
              avg[lon][lat][0][hr][4],
              sdv[lon][lat][0][hr][4],
              avg[lon][lat][0][hr][5],
              sdv[lon][lat][0][hr][5],
              avg[lon][lat][0][hr][6],
              sdv[lon][lat][0][hr][6],
              avg[lon][lat][0][hr][7],
              sdv[lon][lat][0][hr][7],
              avg[lon][lat][0][hr][8],
              sdv[lon][lat][0][hr][8],
              avg[lon][lat][0][hr][9]);
            pFile << buffer;
            if (inwr > 1) {
              buffer[0] = 0;
              snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                0.0, 0.0, -90. + 2.5 * lat, 2.5 * lon);
              ptFile << buffer;
            }
          }
          else if (lvl == 1) {
            // Compute reference density for surface, and
            // write out surface values
            double h3 = avg[lon][lat][3][hr][0];
            double h2 = avg[lon][lat][2][hr][0];
            double hden = (h3 - h2) / log(denref[2] / denref[3]);
            double denrefs = denref[2] * exp((h2 - hatsfc[lon][lat][hr]) / hden);
            snprintf(buffer, 1000, "%3d%8.3f%8.2f%8.2f%6.1f%6.1f%6.1f%6.1f%6.1f"
              "%6.1f%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%7.3f\n",
              1, hatsfc[lon][lat][hr],
              avg[lon][lat][1][hr][0] / 100.0,
              sdv[lon][lat][1][hr][0] / 100.0,
              avg[lon][lat][1][hr][1],
              sdv[lon][lat][1][hr][1],
              avg[lon][lat][1][hr][2],
              sdv[lon][lat][1][hr][2],
              avg[lon][lat][1][hr][3],
              sdv[lon][lat][1][hr][3],
              avg[lon][lat][1][hr][4],
              sdv[lon][lat][1][hr][4],
              avg[lon][lat][1][hr][5],
              sdv[lon][lat][1][hr][5],
              avg[lon][lat][1][hr][6],
              sdv[lon][lat][1][hr][6],
              avg[lon][lat][1][hr][7],
              sdv[lon][lat][1][hr][7],
              avg[lon][lat][1][hr][8],
              sdv[lon][lat][1][hr][8],
              avg[lon][lat][1][hr][9]);
            pFile << buffer;
            if (inwr > 1) {
              buffer[0] = 0;
              snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                0.0, hatsfc[lon][lat][hr], -90. + 2.5 * lat, 2.5 * lon);
              ptFile << buffer;
            }
          }
          else {
            // Output levels k=2-18
            snprintf(buffer, 1000, "%3d%8.3f%8.2f%8.2f%6.1f%6.1f%6.1f%6.1f%6.1f"
              "%6.1f%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%10.3E%10.3E%6.1f%6.1f%7.3f\n",
              lvl,
              avg[lon][lat][lvl][hr][0],
              pmb[lvl],
              sdv[lon][lat][lvl][hr][0] / 100.0,
              avg[lon][lat][lvl][hr][1],
              sdv[lon][lat][lvl][hr][1],
              avg[lon][lat][lvl][hr][2],
              sdv[lon][lat][lvl][hr][2],
              avg[lon][lat][lvl][hr][3],
              sdv[lon][lat][lvl][hr][3],
              avg[lon][lat][lvl][hr][4],
              sdv[lon][lat][lvl][hr][4],
              avg[lon][lat][lvl][hr][5],
              sdv[lon][lat][lvl][hr][5],
              avg[lon][lat][lvl][hr][6],
              sdv[lon][lat][lvl][hr][6],
              avg[lon][lat][lvl][hr][7],
              sdv[lon][lat][lvl][hr][7],
              avg[lon][lat][lvl][hr][8],
              sdv[lon][lat][lvl][hr][8],
              avg[lon][lat][lvl][hr][9]);
            pFile << buffer;
            if (inwr > 1) {
              // Approximate geometric height from geopotential height (km) 
              double z = avg[lon][lat][lvl][hr][0];
              z = z / (1.0 - z / 6356.0);
              buffer[0] = 0;
              snprintf(buffer, 1000, "%5.1f%7.3f%8.1f%8.1f\n",
                0.0, z, -90. + 2.5 * lat, 2.5 * lon);
              ptFile << buffer;
            }
          }
        }
        pFile.close();
        if (inwr > 1) {
          ptFile.close();
        }
      }
    }
  }
  // End of section to write vertical profiles

  // Begin section to write binary file, for reading by GRAM program
  cout << " Enter 1 to write binary file, 0 otherwise: ";
  cin >> inwr;
  if (inwr == 1) {
    cout << " Writing binary file Nb" << por << ".bin" << endl;
    ofstream nFile("Nb" + por + ".bin", ios::binary);
    // Loop to put hr = 1, 5 data into NCEP - dimensioned arrays
    for (int hr = 0; hr < 5; ++hr) {
      for (int lon = 0; lon < 145; ++lon) {
        for (int lat = 0; lat < 73; ++lat) {
          for (int lvl = 0; lvl < 19; ++lvl) {
            if (lvl == 0) {
              geop[lon][lat][lvl] = 0.0;
              sprs[lon][lat][lvl] = sdv[lon][lat][lvl][hr][0];
            }
            else if (lvl == 1) {
              geop[lon][lat][lvl] = hatsfc[lon][lat][hr];
              sprs[lon][lat][lvl] = sdv[lon][lat][lvl][hr][0];
            }
            else {
              geop[lon][lat][lvl] = avg[lon][lat][lvl][hr][0];
              sprs[lon][lat][lvl] = sdv[lon][lat][lvl][hr][0];
            }
          }
        }
      }
      // Write NCEP - dimensioned arrays, in order for GRAM reading
      writeBlock3D(nFile, avg, hr, 1); 
      writeBlock3D(nFile, avg, hr, 5);
      writeBlock3D(nFile, avg, hr, 6); 
      writeBlock3D(nFile, avg, hr, 2); 
      writeBlock3D(nFile, avg, hr, 3); 
      writeBlock3D(nFile, geop);
      writeBlock2D(nFile, avg, 0, hr, 0);
      writeBlock2D(nFile, avg, 1, hr, 0);
      writeBlock3D(nFile, sdv, hr, 1); 
      writeBlock3D(nFile, sdv, hr, 5); 
      writeBlock3D(nFile, sdv, hr, 6); 
      writeBlock3D(nFile, sdv, hr, 2); 
      writeBlock3D(nFile, sdv, hr, 3);  
      writeBlock3D(nFile, sprs); 
      writeBlock2D(nFile, sdv, 0, hr, 0);
      writeBlock2D(nFile, sdv, 1, hr, 0);
      writeBlock3D(nFile, avg, hr, 4);
      writeBlock3D(nFile, sdv, hr, 4);   
      writeBlock3D(nFile, avg, hr, 7);
      writeBlock3D(nFile, sdv, hr, 7);
      writeBlock3D(nFile, avg, hr, 8);
      writeBlock3D(nFile, sdv, hr, 8);
      writeBlock3D(nFile, avg, hr, 9); 
    }
    nFile.close();
  }
}

//------------------------------------------------------------------------------

int getInt(const string& input, size_t width, size_t &pos)
{
  int value = stoi(input.substr(pos, width));
  pos += width;
  return value;
}

void writeBlock3D(ofstream& oFile, double data[145][73][19][5][10], int hr, int var) {
  int size = 145 * 73 * 19 * sizeof(double);
  oFile.write((char*)(&size), sizeof(int));
  for (int lvl = 0; lvl < 19; ++lvl) {
    for (int lat = 0; lat < 73; ++lat) {
      for (int lon = 0; lon < 145; ++lon) {
        oFile.write((char*)(&data[lon][lat][lvl][hr][var]), sizeof(double));
      }
    }
  }
  oFile.write((char*)(&size), sizeof(int));
}

void writeBlock3D(ofstream& oFile, double data[145][73][19]) {
  int size = 145 * 73 * 19 * sizeof(double);
  oFile.write((char*)(&size), sizeof(int));
  for (int lvl = 0; lvl < 19; ++lvl) {
    for (int lat = 0; lat < 73; ++lat) {
      for (int lon = 0; lon < 145; ++lon) {
        oFile.write((char*)(&data[lon][lat][lvl]), sizeof(double));
      }
    }
  }
  oFile.write((char*)(&size), sizeof(int));
}

void writeBlock2D(ofstream& oFile, double data[145][73][19][5][10], int lvl, int hr, int var) {
  int size = 145 * 73 * sizeof(double);
  oFile.write((char*)(&size), sizeof(int));
  for (int lat = 0; lat < 73; ++lat) {
    for (int lon = 0; lon < 145; ++lon) {
      oFile.write((char*)(&data[lon][lat][lvl][hr][var]), sizeof(double));
    }
  }
  oFile.write((char*)(&size), sizeof(int));
}

void sortlevel(int lsort[19], double zsort[19])
{
  //-----------------------------------------------
  // Sorts level numbers, to put surface and sae level in right
  // order
  // Initialize level numbers
  for (int i = 0; i < 19; ++i) {
    lsort[i] = i;
  }
  //...  Store surface (level 1) and sea level (level 0) altitudes
  double zsort1 = zsort[1];
  double zsort0 = zsort[0];
      //...  Sort to put surface level number into right order
  for (int i = 2; i < 19; ++i) {
    if (zsort1 > zsort[i]) {
      lsort[i - 1] = lsort[i];
      lsort[i] = 1;
    }
    else {
      break;
    }
  }
      //...  Store altitudes with surface altitude in right order
  double qsort[19];
  for (int i = 0; i < 19; ++i) {
    qsort[i] = zsort[lsort[i]];
  }
      //...  Sort to put sea level number into right order
  for (int i = 1; i < 19; ++i) {
    if (zsort0 > qsort[i]) {
      lsort[i - 1] = lsort[i];
      lsort[i] = 0;
    }
    else {
      break;
    }
  }
}
//------------------------------------------------------------------------------
