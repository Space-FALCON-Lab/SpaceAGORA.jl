      MODULE sortlevel_I                                                          !SMOD  1
      INTERFACE                                                                   !SMOD  2
      SUBROUTINE sortlevel (LSORT, ZSORT)                                         !SMOD  3
      integer, DIMENSION(0:18), INTENT(INOUT) :: LSORT                            !SMOD  4
      real, DIMENSION(0:18), INTENT(IN) :: ZSORT                                  !SMOD  5
      END SUBROUTINE                                                              !SMOD  6
      END INTERFACE                                                               !SMOD  7
      END MODULE                                                                  !SMOD  8
!-------------------------------------------------------------------------------- !SMOD  9
      program NCEPbinF                                                            !NBIN  1
!...   Program to convert ASCII files to binary values, write out                 !NBIN  2
!      data for map plots and/or profile plots, and write binary                  !NBIN  3
!      version NCEP data arrays                                                   !NBIN  4
!-----------------------------------------------                                  !NBIN  5
!   I n t e r f a c e   B l o c k s                                               !NBIN  6
!-----------------------------------------------                                  !NBIN  7
      use sortlevel_I                                                             !NBIN  8
      implicit none                                                               !NBIN  9
!-----------------------------------------------                                  !NBIN 10
!   L o c a l   V a r i a b l e s                                                 !NBIN 11
!-----------------------------------------------                                  !NBIN 12
      integer , dimension(0:18) :: lsort                                          !NBIN 13
      integer :: ihr, lon, lat, lvl, iihr, ilvl, ilat, ilon, nv, ihavg, ihsig, &  !NBIN 14
         iuavg, iusig, ivavg, ivsig, itavg, itsig, idavg, idsig, nrho, irhavg, &  !NBIN 15
         irhsig, itdavg, itdsig, ivpavg, ivpsig, nvp, i, j, k, k1slp, k2slp,   &  !NBIN 16
         k1sfc, k2sfc, ivar, inwr, ISavg, ISsig, IRuv                             !NBIN 17
      real(8) , dimension(145,73,0:18,5,10) :: avg, std                              !NBIN 18
      real(8) , dimension(0:18) :: denref                                            !NBIN 19
      real(8) , dimension(145,73,5) :: hatsfc                                        !NBIN 20
      real(8) , dimension(0:18) :: p, den, temp, u, v, rh, sp, sden, stemp, su,  &   !NBIN 21
         sv, srh, h,  pmb, zsort, tdew, stdew, vp, svp, SpA, SpS, Ruv             !NBIN 22
      real(8) , dimension(145,73,0:18) :: geop, sprs                                 !NBIN 23
      real(8) :: rlim, gref, sprat, sdrat, strat, hpsfc, hsfc, hpslp,            &   !NBIN 24
         pslp, spslp, factint, rgas1, rgas2, rgas, hvp, xlon, xlat, zlon, h3, &   !NBIN 25
         h2, hden, denref0, denrefs, z                                            !NBIN 26
      character(LEN = 2) :: imo, iyr, iyr1, dummy, olat, olvl                     !NBIN 27
      character(LEN = 3) :: olon                                                  !NBIN 28
      character          :: ohr                                                   !NBIN 29
!                                                                                 !NBIN 30
! Indexes for averages (avg) and standard deviations (std):                       !NBIN 30a
!    1-145 = Lon: 0E to 357.5E (step 2.5) with 145 = 360E = 0 E (added below)     !NBIN 30b
!    1-73  = Lat: 90S to 90 N (step 2.5)                                          !NBIN 30c
!    0-18  = Pressure level:  0=sea level, 1=surface, 2=1000mb, ... 18=10mb       !NBIN 30d 
!    1-5   = Code for hour of day: 1=0 UTC, 2=6 UTC,  3=12UTC,  4=18 UTC,         !NBIN 30e
!               5=All-day average (added below)                                   !NBIN 30f
!    1-10  = Variable Number:                                                     !NBIN 30g
!               1=Geopotential height (m) for pressure level, or surface pressure !NBIN 30h 
!                   (mb) for surface level - height later converted to km         !NBIN 30i
!               2=Temperature (K)                                                 !NBIN 30j
!               3=Eastward (U) wind component (m/s)                               !NBIN 30k
!               4=Northward (V) wind component (m/s)                              !NBIN 30l
!               5=Relative humidity (percent) - computed if above 300mb height    !NBIN 30m
!               6=Density (kg/m**3) - computed from gas law                       !NBIN 30n
!               7=Dewpoint temperature (K) - computed from tdbuck routine         !NBIN 30o
!               8=Vapor pressure (Pa) - computed from wexler routine              !NBIN 30p
!               9=Wind speed (m/s) - computed from U and V components             !NBIN 30q
!              10=U*V product (for computing U-V cross correlation)               !NBIN 30r  
!                                                                                 !NBIN 30s
!-----------------------------------------------                                  !NBIN 31
!...   Pressures (mb) at levels 0-18.  Values for sea level (level=0)             !NBIN 32
!      and surface (level=1) are determined later.                                !NBIN 33
      data pmb/ 9999.d0, 9999.d0, 1000.d0, 925.d0, 850.d0, 700.d0, 600.d0, 500.d0, 400.d0, 300.d0, &  !NBIN 34
         250.d0, 200.d0, 150.d0, 100.d0, 70.d0, 50.d0, 30.d0, 20.d0, 10.d0/                         !NBIN 35
!...   Reference densities (kg/m**3). Values for sea level (level=0)              !NBIN 36
!      and surface (level=1) are determined later.                                !NBIN 37
      data denref/ 9999.d0, 9999.d0, 1.2121d0, 1.1381d0, 1.0626d0, 0.90798d0, 0.80144d0,   &    !NBIN 38
         0.69144d0, 0.57714d0, 0.45721d0, 0.39446d0, 0.32159d0, 0.24119d0, 0.16079d0,      &    !NBIN 39
         0.11255d0, 0.08018d0, 0.04739d0, 0.03122d0, 0.01529d0/                             !NBIN 40
      write (*, *) ' Enter 1st year of POR (yy)'                                  !NBIN 41
      read (*, 5) iyr1                                                            !NBIN 42
    5 format(a)                                                                   !NBIN 43
      write (*, *) ' Enter last year of POR to process (yy)'                      !NBIN 44
      read (*, 5) iyr                                                             !NBIN 45
      write (*, *) ' Enter month (mm)'                                            !NBIN 46
      read (*, 5) imo                                                             !NBIN 47
      write (*, *) ' Reading NCEP ASCII file = Nf', iyr1, iyr, imo, '.txt'        !NBIN 48
!                                                                                 !NBIN 49
!...   Open and read Fixed ASCII file (Nf) of NCEP averages and standard          !NBIN 50
!      deviations                                                                 !NBIN 51
!                                                                                 !NBIN 52
      open(21, file='Nf'//iyr1//iyr//imo//'.txt', position='asis')                !NBIN 53
      read (21, 10) dummy                                                         !NBIN 54
   10 format(a)                                                                   !NBIN 55
      do ihr = 1, 4                                                               !NBIN 56
         do lon = 1, 144                                                          !NBIN 57
            do lat = 1, 73                                                        !NBIN 58
               do lvl = 1, 18                                                     !NBIN 59
                  read (21, 120) iihr, ilvl, ilat, ilon, nv, ihavg, ihsig,     &  !NBIN 60
                     iuavg, iusig, ivavg, ivsig, itavg, itsig, idavg, idsig,   &  !NBIN 61
                     nrho, irhavg, irhsig, itdavg, itdsig, ivpavg, ivpsig,     &  !NBIN 62
                     nvp, ISavg, ISsig, IRuv                                      !NBIN 62
                  if (iihr /= ihr) then                                           !NBIN 63
                     write (*, *) ' Bad hour', iihr, ilvl, ilat, ilon             !NBIN 64
                     stop                                                         !NBIN 65
                  else if (ilvl /= lvl) then                                      !NBIN 66
                     write (*, *) ' Bad level', iihr, ilvl, ilat, ilon            !NBIN 67
                     stop                                                         !NBIN 68
                  else if (ilat /= lat) then                                      !NBIN 69
                     write (*, *) ' Bad lat', iihr, ilvl, ilat, ilon              !NBIN 70
                     stop                                                         !NBIN 71
                  else if (ilon /= lon) then                                      !NBIN 72
                     write (*, *) ' Bad lon', iihr, ilvl, ilat, ilon              !NBIN 73
                     stop                                                         !NBIN 74
                  else                                                            !NBIN 75
!...                 Convert ASCII numbers to correct units:                      !NBIN 76
!                       H = geopotential height (m)                               !NBIN 77
!                       T = temperature (K)                                       !NBIN 78
!                       U, V = eastward, northward wind components (m/s)          !NBIN 79
!                       RH = relative humidity (%)                                !NBIN 80
!                       D = density (g/m**3)                                      !NBIN 81
!                       Td = dewpoint (K)                                         !NBIN 82
!                       vp = vapor pressure (N/m**2)                              !NBIN 83
                     avg(lon,lat,lvl,ihr,1) = ihavg/10.                           !NBIN 84
                     std(lon,lat,lvl,ihr,1) = ihsig/10.                           !NBIN 85
                     avg(lon,lat,lvl,ihr,2) = itavg/10.                           !NBIN 86
                     std(lon,lat,lvl,ihr,2) = itsig/10.                           !NBIN 87
                     avg(lon,lat,lvl,ihr,3) = iuavg/10.                           !NBIN 88
                     std(lon,lat,lvl,ihr,3) = iusig/10.                           !NBIN 89
                     avg(lon,lat,lvl,ihr,4) = ivavg/10.                           !NBIN 90
                     std(lon,lat,lvl,ihr,4) = ivsig/10.                           !NBIN 91
                     avg(lon,lat,lvl,ihr,5) = irhavg/10.                          !NBIN 92
                     std(lon,lat,lvl,ihr,5) = irhsig/10.                          !NBIN 93
                     avg(lon,lat,lvl,ihr,6) = idavg*10.0**(-nrho)                 !NBIN 94
                     std(lon,lat,lvl,ihr,6) = idsig*10.0**(-nrho)                 !NBIN 95
                     avg(lon,lat,lvl,ihr,7) = itdavg/10.                          !NBIN 96
                     std(lon,lat,lvl,ihr,7) = itdsig/10.                          !NBIN 97
                     avg(lon,lat,lvl,ihr,8) = ivpavg*10.0**(-nvp)                 !NBIN 98
                     std(lon,lat,lvl,ihr,8) = ivpsig*10.0**(-nvp)                 !NBIN 99
                     avg(lon,lat,lvl,ihr,9) = ISavg/10.0                          !NBIN 99a
                     std(lon,lat,lvl,ihr,9) = ISsig/10.0                          !NBIN 99b
                     avg(lon,lat,lvl,ihr,10) = IRuv/1000.0                        !NBIN 99c
                     std(lon,lat,lvl,ihr,10) = 0.0                                !NBIN 99d
                  endif                                                           !NBIN100
  120             format(i1,2i2,i3,i4,i6,i5,4(i4,i3),i1,2(i4,i3),2i4,i1,i4,i3,i4) !NBIN101
               end do                                                             !NBIN102
            end do                                                                !NBIN103
         end do                                                                   !NBIN104
      end do                                                                      !NBIN105
      close(21)                                                                   !NBIN106
!                                                                                 !NBIN107
!...   START OF test and data conversion section                                  !NBIN108
!                                                                                 !NBIN109
      write (*, *) ' Processing data. Please stand by.'                           !NBIN110
      rlim = 0.999d0                                                               !NBIN112
      gref = 9.80665d0                                                              !NBIN113
!                                                                                 !NBIN114
      do lon = 1, 144                                                             !NBIN115
         do lat = 1, 73                                                           !NBIN116
            do ihr = 1, 4                                                         !NBIN117
               i = lon                                                            !NBIN118
               j = lat                                                            !NBIN119
!              Convert units and store profile values in single-dimension arrays  !NBIN120
               h(1:18) = avg(i,j,1:18,ihr,1)                                      !NBIN121
!.....         sh(1:18) = std(i,j,1:18,ihr,1)                                     !NBIN122
               temp(1:18) = avg(i,j,1:18,ihr,2)                                   !NBIN123
               stemp(1:18) = std(i,j,1:18,ihr,2)                                  !NBIN124
               u(1:18) = avg(i,j,1:18,ihr,3)                                      !NBIN125
               su(1:18) = std(i,j,1:18,ihr,3)                                     !NBIN126
               v(1:18) = avg(i,j,1:18,ihr,4)                                      !NBIN127
               sv(1:18) = std(i,j,1:18,ihr,4)                                     !NBIN128
               rh(1:18) = avg(i,j,1:18,ihr,5)                                     !NBIN129
               srh(1:18) = std(i,j,1:18,ihr,5)                                    !NBIN130
!...           Convert density to kg/m**3                                         !NBIN131
               den(1:18) = avg(i,j,1:18,ihr,6)/1000.                              !NBIN132
               sden(1:18) = std(i,j,1:18,ihr,6)/1000.                             !NBIN133
               tdew(1:18) = avg(i,j,1:18,ihr,7)                                   !NBIN134
               stdew(1:18) = std(i,j,1:18,ihr,7)                                  !NBIN135
               vp(1:18) = avg(i,j,1:18,ihr,8)                                     !NBIN136
               svp(1:18) = std(i,j,1:18,ihr,8)                                    !NBIN137
               SpA(1:18) = avg(i,j,1:18,ihr,9)                                    !NBIN137a
               SpS(1:18) = std(i,j,1:18,ihr,9)                                    !NBIN137b
               Ruv(1:18) = avg(i,j,1:18,ihr,10)                                   !NBIN137c
!...           Find and store level index values for surface (k=1) and sea        !NBIN138
!              level (k=0).                                                       !NBIN139
               do lvl = 1, 18                                                     !NBIN140
                  k = lvl                                                         !NBIN141
                  if (k == 1) then                                                !NBIN142
                     h(0) = 0.0                                                   !NBIN143
!...                 Convert pressure to N/m**2                                   !NBIN144
                     p(1) = avg(i,j,k,ihr,1)*100.                                 !NBIN145
                     sp(1) = std(i,j,k,ihr,1)*100.                                !NBIN146
                     sprat = sp(1)/p(1)                                           !NBIN147
!...                 Use levels 2 (1000 mb) and 3 (925 mb) for interpolation      !NBIN148
!                    to sea level, if level 2 height is positive                  !NBIN149
                     if (h(2) > 0.0) then                                         !NBIN150
                        k1slp = 2                                                 !NBIN151
                        k2slp = 3                                                 !NBIN152
                     endif                                                        !NBIN153
!...                 Use levels 2 (1000 mb) and 3 (925 mb) for interpolation      !NBIN154
!                    to surface, if level 1 (surface) pressure is greater than    !NBIN155
!                    1000 mb.                                                     !NBIN156
                     if (p(1) > pmb(2)*100.) then                                 !NBIN157
                        k1sfc = 2                                                 !NBIN158
                        k2sfc = 3                                                 !NBIN159
                     endif                                                        !NBIN160
                  else                                                            !NBIN161
!...                 Convert pressures to N/m**2                                  !NBIN162
                     p(k) = pmb(k)*100.                                           !NBIN163
!...                 Convert standard deviation in geopotential height to         !NBIN164
!                    pressure standard deviation                                  !NBIN165
                     sp(k) = gref*std(i,j,k,ihr,1)*den(k)                         !NBIN166
                     sprat = sp(k)/p(k)                                           !NBIN167
!...                 Find level index values below and above sea level, for       !NBIN168
!                    interpolation                                                !NBIN169
                     if (h(k)<=0.0 .and. h(k+1)>0.0) then                         !NBIN170
                        k1slp = k                                                 !NBIN171
                        k2slp = k + 1                                             !NBIN172
                     endif                                                        !NBIN173
!...                 Find level index values below and above surface, for         !NBIN174
!                    interpolation                                                !NBIN175
                     if (k<18 .and. p(1)/100.0<=pmb(k) .and. p(1)>pmb(k+1)) &     !NBIN176
                        then                                                      !NBIN177
                        k1sfc = k                                                 !NBIN178
                        k2sfc = k + 1                                             !NBIN179
                     endif                                                        !NBIN180
                  endif                                                           !NBIN181
!...              Correct sigmas using correlation constraints, based on the      !NBIN182
!                 Buell relationships                                             !NBIN183
                  sdrat = sden(k)/den(k)                                          !NBIN184
                  strat = stemp(k)/temp(k)                                        !NBIN185
                  if (sdrat >= rlim*(sprat + strat)) then                         !NBIN186
!...                 Correct density standard deviation                           !NBIN187
                     sden(k) = rlim*den(k)*(sprat + strat)                        !NBIN188
                  else if (strat >= rlim*(sprat + sdrat)) then                    !NBIN189
!...                 Correct temperature standard deviation                       !NBIN190
                     stemp(k) = rlim*temp(k)*(sprat + sdrat)                      !NBIN191
                  else if (sprat >= rlim*(sdrat + strat)) then                    !NBIN192
!...                 Correct height (or surface pressure) standard deviation      !NBIN193
                     sp(k) = rlim*p(k)*(sdrat + strat)                            !NBIN194
                  endif                                                           !NBIN195
               end do                                                             !NBIN196
!...           Interpolate for surface geopotential height                        !NBIN197
               hpsfc = (h(k2sfc)-h(k1sfc))/log(p(k1sfc)/p(k2sfc))                !NBIN198
               hsfc = h(k1sfc) + hpsfc*log(p(k1sfc)/p(1))                        !NBIN199
               hatsfc(i,j,ihr) = hsfc/1000.                                       !NBIN200
!...           Interpolate (or extrapolate) for sea level conditions              !NBIN201
               hpslp = (h(k2slp)-h(k1slp))/log(p(k1slp)/p(k2slp))                !NBIN202
               pslp = p(k1slp)*exp(h(k1slp)/hpslp)                                !NBIN203
               spslp = pslp*sp(k1slp)/p(k1slp)                                    !NBIN204
               factint = h(k1slp)/(h(k1slp)-h(k2slp))                             !NBIN205
               temp(0) = temp(k1slp) + (temp(k2slp)-temp(k1slp))*factint          !NBIN206
               stemp(0) = stemp(k1slp) + (stemp(k2slp)-stemp(k1slp))*factint      !NBIN207
               u(0) = u(k1slp) + (u(k2slp)-u(k1slp))*factint                      !NBIN208
               su(0) = su(k1slp) + (su(k2slp)-su(k1slp))*factint                  !NBIN209
               v(0) = v(k1slp) + (v(k2slp)-v(k1slp))*factint                      !NBIN210
               sv(0) = sv(k1slp) + (sv(k2slp)-sv(k1slp))*factint                  !NBIN211
               rh(0) = rh(k1slp) + (rh(k2slp)-rh(k1slp))*factint                  !NBIN212
               srh(0) = srh(k1slp) + (srh(k2slp)-srh(k1slp))*factint              !NBIN213
               rgas1 = p(k1slp)/(temp(k1slp)*den(k1slp))                          !NBIN214
               rgas2 = p(k2slp)/(temp(k2slp)*den(k2slp))                          !NBIN215
               rgas = rgas1 + (rgas2 - rgas1)*factint                             !NBIN216
               den(0) = pslp/(rgas*temp(0))                                       !NBIN217
               sden(0) = sden(k1slp) + (sden(k2slp)-sden(k1slp))*factint          !NBIN218
               tdew(0) = tdew(k1slp) + (tdew(k2slp)-tdew(k1slp))*factint          !NBIN219
               stdew(0) = stdew(k1slp) + (stdew(k2slp)-stdew(k1slp))*factint      !NBIN220
               if (vp(k2slp)==0.0 .or. vp(k2slp)-vp(k1slp)==0.0d0) then           !NBIN221
                  hvp = 1.0E10                                                    !NBIN222
               else                                                               !NBIN223
                  hvp = (h(k2slp)-h(k1slp))/log(vp(k1slp)/vp(k2slp))             !NBIN224
               endif                                                              !NBIN225
               vp(0) = vp(k1slp)*exp(h(k1slp)/hvp)                                !NBIN226
               svp(0) = vp(0)*svp(k1slp)/vp(k1slp)                                !NBIN227
               SpA(0) = SpA(k1slp) + (SpA(k2slp) - SpA(k1slp))*factint            !NBIN227a
               SpS(0) = SpS(k1slp) + (SpS(k2slp) - SpS(k1slp))*factint            !NBIN227b
               Ruv(0) = Ruv(k1slp) + (Ruv(k2slp) - Ruv(k1slp))*factint            !NBIN227c
!...           Store sea level parameters in Avg and Std arrays as level 0        !NBIN228
               avg(i,j,0,ihr,1) = pslp                                            !NBIN229
               std(i,j,0,ihr,1) = spslp                                           !NBIN230
               avg(i,j,0,ihr,2) = temp(0)                                         !NBIN231
               std(i,j,0,ihr,2) = stemp(0)                                        !NBIN232
               avg(i,j,0,ihr,3) = u(0)                                            !NBIN233
               std(i,j,0,ihr,3) = su(0)                                           !NBIN234
               avg(i,j,0,ihr,4) = v(0)                                            !NBIN235
               std(i,j,0,ihr,4) = sv(0)                                           !NBIN236
               avg(i,j,0,ihr,5) = rh(0)                                           !NBIN237
               std(i,j,0,ihr,5) = srh(0)                                          !NBIN238
               avg(i,j,0,ihr,6) = den(0)                                          !NBIN239
               std(i,j,0,ihr,6) = sden(0)                                         !NBIN240
               avg(i,j,0,ihr,7) = tdew(0)                                         !NBIN241
               std(i,j,0,ihr,7) = stdew(0)                                        !NBIN242
               avg(i,j,0,ihr,8) = vp(0)                                           !NBIN243
               std(i,j,0,ihr,8) = svp(0)                                          !NBIN244
               avg(i,j,0,ihr,9) = SpA(0)                                          !NBIN244a
               std(i,j,0,ihr,9) = SpS(0)                                          !NBIN244b
               avg(i,j,0,ihr,10) = Ruv(0)                                         !NBIN244c
               std(i,j,0,ihr,10) = 0.0                                            !NBIN244d
!...           Store Density as kg/m**3, pressure as N/m**2, height as km         !NBIN245
               do k = 1, 18                                                       !NBIN246
                  avg(i,j,k,ihr,6) = den(k)                                       !NBIN247
                  std(i,j,k,ihr,6) = sden(k)                                      !NBIN248
                  if (k == 1) then                                                !NBIN249
                     avg(i,j,k,ihr,1) = p(1)                                      !NBIN250
                     std(i,j,k,ihr,1) = sp(1)                                     !NBIN251
                  else                                                            !NBIN252
!                    Convert geopotential height to km                            !NBIN252a       
                     avg(i,j,k,ihr,1) = h(k)/1000.                                !NBIN253
                     std(i,j,k,ihr,1) = sp(k)                                     !NBIN254
                  endif                                                           !NBIN255
               end do                                                             !NBIN256
            end do                                                                !NBIN257
!                                                                                 !NBIN258
!...        Compute averages and standard deviations for all day (ihr=5)          !NBIN259
            hatsfc(lon,lat,5) = 0.0                                               !NBIN260
            hatsfc(lon,lat,5) = hatsfc(lon,lat,5) + sum(hatsfc(lon,lat,:4)/4.)    !NBIN261
            do lvl = 0, 18                                                        !NBIN262
               do ivar = 1, 10                                                    !NBIN263
                  avg(lon,lat,lvl,5,ivar) = 0.0                                   !NBIN264
                  std(lon,lat,lvl,5,ivar) = 0.0                                   !NBIN265
                  avg(lon,lat,lvl,5,ivar) = avg(lon,lat,lvl,5,ivar) +   &         !NBIN266
                     sum(avg(lon,lat,lvl,:4,ivar)/4.)                             !NBIN267
                  std(lon,lat,lvl,5,ivar) = std(lon,lat,lvl,5,ivar) +   &         !NBIN268
                     sum(std(lon,lat,lvl,:4,ivar)**2/4.  +              &         !NBIN269
                     (avg(lon,lat,lvl,:4,ivar)**2  -                    &         !NBIN270
                     avg(lon,lat,lvl,5,ivar)**2)/4.)                              !NBIN271
                  std(lon,lat,lvl,5,ivar) =                             &         !NBIN272
                     sqrt(abs(std(lon,lat,lvl,5,ivar)))                           !NBIN273
               end do                                                             !NBIN274
            end do                                                                !NBIN275
         end do                                                                   !NBIN276
      end do                                                                      !NBIN277
!                                                                                 !NBIN278
!...   Copy longitude=0 values to longitude=360                                   !NBIN279
      hatsfc(145,:,:) = hatsfc(1,:,:)                                             !NBIN280
      avg(145,:,:18,:,:) = avg(1,:,:18,:,:)                                       !NBIN281
      std(145,:,:18,:,:) = std(1,:,:18,:,:)                                       !NBIN282
!                                                                                 !NBIN283
!...   END OF test and data conversion section                                    !NBIN284
!                                                                                 !NBIN285
!       Begin section to write out lon-lat maps by level                          !NBIN286
!                                                                                 !NBIN287
      write (*, *) ' Enter 1 for Maps, 2 for Maps & Trajectory files, 0 for none' !NBIN288
      read (*, *) inwr                                                            !NBIN289
      if (inwr >= 1) then                                                         !NBIN290
  200    continue                                                                 !NBIN291
         write (*, *) ' Enter indexes for lvl(0-18),ihr(1-5), or ', '0s to end'   !NBIN292
         read (*, *) lvl, ihr                                                     !NBIN293
         If (ihr > 5) ihr = 5                                                     !NBIN293a
         If (lvl > 18) lvl = 18                                                   !NBIN293b         
         if (lvl>=0 .and. ihr>0) then                                             !NBIN294
            write (olvl, '(I2.2)') lvl                                            !NBIN295
            write (ohr, '(I1)') ihr                                               !NBIN296
            write (*, *) ' Writing M', iyr1, iyr, imo, olvl, ohr, '.txt'          !NBIN297
            If (inwr > 1)Then                                                     !NBIN297a
              write (*, *) ' Writing MT', iyr1, iyr, imo, olvl, ohr, '.txt'       !NBIN297b
              open(23, file='MT'//iyr1//iyr//imo//olvl//ohr//'.txt', position= &  !NBIN297c
               'asis')                                                            !NBIN297d
            Endif                                                                 !NBIN297e
            open(22, file='M'//iyr1//iyr//imo//olvl//ohr//'.txt', position=    &  !NBIN298
               'asis')                                                            !NBIN299
            if (lvl > 1) then                                                     !NBIN300
               write (22, 260)                                                    !NBIN301
            else                                                                  !NBIN302
               write (22, 270)                                                    !NBIN303
            endif                                                                 !NBIN304
  260       format('   I  J   Lon   Lat  Havg_m  Psigmb  Tavg  Tsig',        &    !NBIN305
               '  Uavg  Usig  Vavg  Vsig RHavg RHsig Dav76 Dsd76',           &    !NBIN306
               ' TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv')             !NBIN307
  270       format('   I  J   Lon   Lat  Pavgmb  Psigmb  Tavg  Tsig',        &    !NBIN308
               '  Uavg  Usig  Vavg  Vsig RHavg RHsig  Davg  Dstd',           &    !NBIN309
               ' TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv Hatsfc')      !NBIN310
            if (lvl > 1) then                                                     !NBIN311
               do lon = 1, 145                                                    !NBIN312
                  do lat = 1, 73                                                  !NBIN313
!...          Normalize density by reference density if not at                    !NBIN314
!             surface or sea level                                                !NBIN315
!...          If at surface or sea level, write density avg and                   !NBIN315a
!             std in kg/m**3                                                      !NBIN315b
                     write (22, 280) lon, lat, 2.5*(lon - 1.), -90. +          &  !NBIN316
                        2.5*(lat- 1.),1000.*avg(lon,lat,lvl,ihr,1),            &  !NBIN317
                        std(lon,lat,lvl,ihr,1)/100.,(avg(lon,lat,lvl,ihr,i),   &  !NBIN318
                        std(lon,lat,lvl,ihr,i),i=2,5),avg(lon,lat,lvl,ihr,6)/  &  !NBIN319
                        denref(lvl), std(lon,lat,lvl,ihr,6)/denref(lvl),       &  !NBIN320
                        avg(lon,lat,lvl,ihr,7), std(lon,lat,lvl,ihr,7),        &  !NBIN321
                        avg(lon,lat,lvl,ihr,8), std(lon,lat,lvl,ihr,8),        &  !NBIN322
                        avg(lon,lat,lvl,ihr,9), std(lon,lat,lvl,ihr,9),        &  !NBIN322a
                        avg(lon,lat,lvl,ihr,10)                                   !NBIN322b
                     If (inwr > 1)Then                                            !NBIN322c
!                      Approximate geometric height from geopotential height (km) !NBIN323
                       z = avg(lon,lat,lvl,ihr,1)                                 !NBIN324
                       z = z/(1.0 - z/6356.)                                      !NBIN324a
!                      Write to Map trajectory file                               !NBIN324b
                       Write(23,295)0.0,z,-90.+2.5*(lat-1.),2.5*(lon-1.)          !NBIN324c
                     Endif                                                        !NBIN324d
                  end do                                                          !NBIN325
               end do                                                             !NBIN326
            else                                                                  !NBIN327
               do lon = 1, 145                                                    !NBIN328
!...            If at surface or sea level, write density avg and                 !NBIN329
!               std in kg/m**3                                                    !NBIN330
                  do lat = 1, 73                                                  !NBIN331
                     write (22, 290) lon, lat, 2.5*(lon - 1.), -90. +          &  !NBIN334
                        2.5*(lat- 1.),avg(lon,lat,lvl,ihr,1)/100.,             &  !NBIN335
                        std(lon,lat,lvl,ihr,1)/100., (avg(lon,lat,lvl,ihr,i),  &  !NBIN336
                        std(lon,lat,lvl,ihr,i),i=2,5), avg(lon,lat,lvl,ihr,6), &  !NBIN337
                        std(lon,lat,lvl,ihr,6), avg(lon,lat,lvl,ihr,7),        &  !NBIN338
                        std(lon,lat,lvl,ihr,7), avg(lon,lat,lvl,ihr,8),        &  !NBIN339
                        std(lon,lat,lvl,ihr,8), avg(lon,lat,lvl,ihr,9),        &  !NBIN339a
                        std(lon,lat,lvl,ihr,9), avg(lon,lat,lvl,ihr,10),       &  !NBIN339b
                        hatsfc(lon,lat,ihr)                                       !NBIN340
                     If (inwr > 1)Then                                            !NBIN340a
                       If (lvl == 0)Then                                          !NBIN340b
                         Write(23,295)0.0,0.0,-90.+2.5*(lat-1.),2.5*(lon-1.)      !NBIN340c
                       Else                                                       !NBIN340d
                         Write(23,295)0.0,hatsfc(lon,lat,ihr),                 &  !NBIN340e
                           -90.+2.5*(lat-1.),2.5*(lon-1.)                         !NBIN340f
                       Endif                                                      !NBIN340g
                     Endif                                                        !NBIN340
                  end do                                                          !NBIN341
               end do                                                             !NBIN342
            endif                                                                 !NBIN343
  280       format(i4,i3,2f6.1,f8.1,f8.2,8f6.1,2f6.3,2f6.1,1p,2e10.3,0p,       &  !NBIN344
              2f6.1,f7.3)                                                         !NBIN344a
  290       format(i4,i3,2f6.1,f8.1,f8.2,8f6.1,2f6.3,2f6.1,1p,2e10.3,0p,       &  !NBIN345
              2f6.1,2f7.3)                                                        !NBIN345a
  295       Format(F5.1,F7.3,2F8.1)                                               !NBIN345b
            close(22)                                                             !NBIN346
            close(23)                                                             !NBIN346a
            go to 200                                                             !NBIN347
         endif                                                                    !NBIN348
      endif                                                                       !NBIN349
!                                                                                 !NBIN350
!...    End of section to write maps                                              !NBIN351
!                                                                                 !NBIN352
!...    Begin section to write vertical profiles at specified lon-lat             !NBIN353
!                                                                                 !NBIN354
      write (*, *) ' Enter 1 for Profiles, 2 for Profiles & Trajectory ',      &  !NBIN355
        'files, or 0 for none'                                                    !NBIN355a
      read (*, *) inwr                                                            !NBIN356
      if (inwr >= 1) then                                                         !NBIN357
  300    continue                                                                 !NBIN358
         write (*, *) ' Enter lon(deg E),lat(deg),ihr(1-5), or ihr=0 to end'      !NBIN359
  305    read (*, *) xlon, xlat, ihr                                              !NBIN360
         If (ihr > 5)ihr = 5                                                      !NBIN360a
         If (Abs(xlat) > 90.0 .or. Abs(xlon) > 360.0)Then                         !NBIN360b
           Write(*,*)' Bad Lat-Lon'                                               !NBIN360c
           Goto 305                                                               !NBIN360d
         Endif                                                                    !NBIN360e
         lat = 1 + nint((xlat + 90.)/2.5)                                         !NBIN361
         zlon = xlon                                                              !NBIN362
         if (zlon < 0.0) zlon = zlon + 360.                                       !NBIN363
         lon = 1 + nint(zlon/2.5)                                                 !NBIN364
         if (ihr > 0) then                                                        !NBIN365
            write (olat, '(I2.2)') lat                                            !NBIN366
            write (olon, '(I3.3)') lon                                            !NBIN367
            write (ohr, '(I1)') ihr                                               !NBIN368
            write (*, *) 'Writing P', iyr1, iyr, imo, olon, olat, ohr, '.txt'     !NBIN369
            open(22, file='P'//iyr1//iyr//imo//olon//olat//ohr//'.txt',        &  !NBIN370
               position='asis')                                                   !NBIN371
            If (inwr > 1)Then                                                     !NBIN371a
              write (*, *) 'Writing PT', iyr1, iyr, imo, olon, olat, ohr, '.txt'  !NBIN371b
              open(23, file='PT'//iyr1//iyr//imo//olon//olat//ohr//'.txt',     &  !NBIN371c
                 position='asis')                                                 !NBIN371d
            Endif                                                                 !NBIN371e
            write (22, 310) lon, lat, xlon, xlat                                  !NBIN372
  310       format('  L  Havgkm  Pavgmb  Psigmb  Tavg  Tsig  ',                &  !NBIN373
               'Uavg  Usig  Vavg  Vsig RHavg RHsig      Davg      Dsig',       &  !NBIN374
               ' TdAvg TdSig     VpAvg     VpSig SpAvg SpSig    Ruv I=',i3,    &  !NBIN375
               ' J=',i2,' Lon=',f6.1,' Lat=',f5.1)                                !NBIN376
            zsort(0) = 0.0                                                        !NBIN377
            zsort(1) = hatsfc(lon,lat,ihr)                                        !NBIN378
            zsort(2:18) = avg(lon,lat,2:18,ihr,1)                                 !NBIN379
!...        Sort levels to put surface and sea level in right order               !NBIN380
            call sortlevel (lsort, zsort)                                         !NBIN381
!...        Write out profile levels                                              !NBIN382
            do k = 0, 18                                                          !NBIN383
               lvl = lsort(k)                                                     !NBIN384
!...          Write sea level values                                              !NBIN385
               if (lvl == 0) then                                                 !NBIN386
!...              Compute reference density for sea level, and                    !NBIN387
!                 write out sea level values                                      !NBIN388
                  h3 = avg(lon,lat,3,ihr,1)                                       !NBIN389
                  h2 = avg(lon,lat,2,ihr,1)                                       !NBIN390
                  hden = (h3 - h2)/log(denref(2)/denref(3))                      !NBIN391
                  denref0 = denref(2)*exp(h2/hden)                                !NBIN392
                  write (22, 320) 0, 0.0, avg(lon,lat,0,ihr,1)/100.,            & !NBIN393
                     std(lon,lat,0,ihr,1)/100., (avg(lon,lat,0,ihr,i),          & !NBIN394
                     std(lon,lat,0,ihr,i),i=2,5), avg(lon,lat,0,ihr,6),         & !NBIN395
                     std(lon,lat,0,ihr,6), avg(lon,lat,0,ihr,7),                & !NBIN396
                     std(lon,lat,0,ihr,7), avg(lon,lat,0,ihr,8),                & !NBIN397
                     std(lon,lat,0,ihr,8), avg(lon,lat,0,ihr,9),                & !NBIN398
                     std(lon,lat,0,ihr,9), avg(lon,lat,0,ihr,10)                  !NBIN398a
                  If (inwr > 1)Write(23,295)0.0,0.0,-90.+2.5*(lat-1.),2.5*(lon-1.)!NBIN398b
               else if (lvl == 1) then                                            !NBIN399
!...              Compute reference density for surface, and                      !NBIN400
!                 write out surface values                                        !NBIN401
                  h3 = avg(lon,lat,3,ihr,1)                                       !NBIN402
                  h2 = avg(lon,lat,2,ihr,1)                                       !NBIN403
                  hden = (h3 - h2)/log(denref(2)/denref(3))                      !NBIN404
                  denrefs = denref(2)*exp((h2 - hatsfc(lon,lat,ihr))/hden)        !NBIN405
                  write (22, 320) 1, hatsfc(lon,lat,ihr), avg(lon,lat,1,ihr,1)/ & !NBIN406
                     100.,  std(lon,lat,1,ihr,1)/100., (avg(lon,lat,1,ihr,i),   & !NBIN407
                     std(lon,lat,1,ihr,i),i=2,5), avg(lon,lat,1,ihr,6),         & !NBIN408
                     std(lon,lat,1,ihr,6), avg(lon,lat,1,ihr,7),                & !NBIN409
                     std(lon,lat,1,ihr,7), avg(lon,lat,1,ihr,8),                & !NBIN410
                     std(lon,lat,1,ihr,8), avg(lon,lat,1,ihr,9),                & !NBIN411
                     std(lon,lat,1,ihr,9), avg(lon,lat,1,ihr,10)                  !NBIN411a
                     If (inwr > 1)Write(23,295)0.0,hatsfc(lon,lat,ihr),         & !NBIN411b
                       -90.+2.5*(lat-1.),2.5*(lon-1.)                             !NBIN411c
               else                                                               !NBIN412
!...            Output levels k=2-18                                              !NBIN413
                  write (22, 320) lvl, avg(lon,lat,lvl,ihr,1), pmb(lvl),      &   !NBIN414
                     std(lon,lat,lvl,ihr,1)/100., (avg(lon,lat,lvl,ihr,i),    &   !NBIN415
                     std(lon,lat,lvl,ihr,i),i=2,5), avg(lon,lat,lvl,ihr,6),   &   !NBIN416
                     std(lon,lat,lvl,ihr,6),                                  &   !NBIN417
                     avg(lon,lat,lvl,ihr,7),std(lon,lat,lvl,ihr,7),           &   !NBIN418
                     avg(lon,lat,lvl,ihr,8), std(lon,lat,lvl,ihr,8),          &   !NBIN419
                     avg(lon,lat,lvl,ihr,9), std(lon,lat,lvl,ihr,9),          &   !NBIN419a
                     avg(lon,lat,lvl,ihr,10)                                      !NBIN419b
                     If (inwr > 1)Then                                            !NBIN419c
!                      Approximate geometric height from geopotential height (km) !NBIN419d
                       z = avg(lon,lat,lvl,ihr,1)                                 !NBIN419e
                       z = z/(1.0 - z/6356.)                                      !NBIN419f
                       Write(23,295)0.0,z,-90.+2.5*(lat-1.),2.5*(lon-1.)          !NBIN419g
                     Endif                                                        !NBIN419h
               endif                                                              !NBIN420
  320          format(i3,f8.3,2f8.2,8f6.1,1p,2e10.3,0p,2f6.1,1p,2e10.3,0p,    &   !NBIN421
                 2f6.1,f7.3)                                                      !NBIN421a
            end do                                                                !NBIN422
            close(22)                                                             !NBIN423
            close(23)                                                             !NBIN423a
            go to 300                                                             !NBIN424
         endif                                                                    !NBIN425
      endif                                                                       !NBIN426
!                                                                                 !NBIN427
!...    End of section to write vertical profiles                                 !NBIN428
!                                                                                 !NBIN429
!...    Begin section to write binary file, for reading by GRAM program           !NBIN430
!                                                                                 !NBIN431
      write (*, *) ' Enter 1 to write binary file, 0 otherwise'                   !NBIN432
      read (*, *) inwr                                                            !NBIN433
      if (inwr == 1) then                                                         !NBIN434
         write (*, *) ' Writing binary file Nb', iyr1, iyr, imo, '.bin'           !NBIN435
         open(22, file='Nb'//iyr1//iyr//imo//'.bin', form='unformatted',   &      !NBIN436
            position='asis')                                                      !NBIN437
!...      Loop to put hr=1,5 data into NCEP-dimensioned arrays                    !NBIN438
         do ihr = 1, 5                                                            !NBIN439
            do lon = 1, 145                                                       !NBIN440
               do lat = 1, 73                                                     !NBIN441
                  do k = 0, 18                                                    !NBIN442
                     if (k == 0) then                                             !NBIN443
                        geop(lon,lat,0) = 0.0                                     !NBIN444
                        sprs(lon,lat,0) = std(lon,lat,0,ihr,1)                    !NBIN445
                     else if (k == 1) then                                        !NBIN446
                        geop(lon,lat,1) = hatsfc(lon,lat,ihr)                     !NBIN447
                        sprs(lon,lat,1) = std(lon,lat,1,ihr,1)                    !NBIN448
                     else                                                         !NBIN449
                        geop(lon,lat,k) = avg(lon,lat,k,ihr,1)                    !NBIN450
                        sprs(lon,lat,k) = std(lon,lat,k,ihr,1)                    !NBIN451
                     endif                                                        !NBIN452
                  end do                                                          !NBIN453
               end do                                                             !NBIN454
            end do                                                                !NBIN455
!...        Write NCEP-dimensioned arrays, in order for GRAM reading              !NBIN456
            write (22) (((avg(lon,lat,k,ihr,2),lon=1,145),lat=1,73),k=0,18)       !NBIN457
            write (22) (((avg(lon,lat,k,ihr,6),lon=1,145),lat=1,73),k=0,18)       !NBIN458
            write (22) (((avg(lon,lat,k,ihr,7),lon=1,145),lat=1,73),k=0,18)       !NBIN459
            write (22) (((avg(lon,lat,k,ihr,3),lon=1,145),lat=1,73),k=0,18)       !NBIN460
            write (22) (((avg(lon,lat,k,ihr,4),lon=1,145),lat=1,73),k=0,18)       !NBIN461
            write (22) geop                                                       !NBIN462
            write (22) ((avg(lon,lat,0,ihr,1),lon=1,145),lat=1,73)                !NBIN463
            write (22) ((avg(lon,lat,1,ihr,1),lon=1,145),lat=1,73)                !NBIN464
            write (22) (((std(lon,lat,k,ihr,2),lon=1,145),lat=1,73),k=0,18)       !NBIN465
            write (22) (((std(lon,lat,k,ihr,6),lon=1,145),lat=1,73),k=0,18)       !NBIN466
            write (22) (((std(lon,lat,k,ihr,7),lon=1,145),lat=1,73),k=0,18)       !NBIN467
            write (22) (((std(lon,lat,k,ihr,3),lon=1,145),lat=1,73),k=0,18)       !NBIN468
            write (22) (((std(lon,lat,k,ihr,4),lon=1,145),lat=1,73),k=0,18)       !NBIN469
            write (22) sprs                                                       !NBIN470
            write (22) ((std(lon,lat,0,ihr,1),lon=1,145),lat=1,73)                !NBIN471
            write (22) ((std(lon,lat,1,ihr,1),lon=1,145),lat=1,73)                !NBIN472
            write (22) (((avg(lon,lat,k,ihr,5),lon=1,145),lat=1,73),k=0,18)       !NBIN473
            write (22) (((std(lon,lat,k,ihr,5),lon=1,145),lat=1,73),k=0,18)       !NBIN474
            write (22) (((avg(lon,lat,k,ihr,8),lon=1,145),lat=1,73),k=0,18)       !NBIN475
            write (22) (((std(lon,lat,k,ihr,8),lon=1,145),lat=1,73),k=0,18)       !NBIN476
            write (22) (((avg(lon,lat,k,ihr,9),lon=1,145),lat=1,73),k=0,18)       !NBIN476a
            write (22) (((std(lon,lat,k,ihr,9),lon=1,145),lat=1,73),k=0,18)       !NBIN476b
            write (22) (((avg(lon,lat,k,ihr,10),lon=1,145),lat=1,73),k=0,18)      !NBIN476c
         end do                                                                   !NBIN477
         close(22)                                                                !NBIN478
      endif                                                                       !NBIN479
      stop                                                                        !NBIN480
      end program NCEPbinF                                                        !NBIN481
                                                                                  !NBIN482
!------------------------------------------------------------------------------   !NBIN483
      subroutine sortlevel(lsort, zsort)                                          !SRTL  1
      implicit none                                                               !SRTL  2
!-----------------------------------------------                                  !SRTL  3
!   D u m m y   A r g u m e n t s                                                 !SRTL  4
!-----------------------------------------------                                  !SRTL  5
      integer , intent(inout) :: lsort(0:18)                                      !SRTL  6
      real , intent(in) :: zsort(0:18)                                            !SRTL  7
!-----------------------------------------------                                  !SRTL  8
!   L o c a l   V a r i a b l e s                                                 !SRTL  9
!-----------------------------------------------                                  !SRTL 10
      integer :: i                                                                !SRTL 11
      real , dimension(0:18) :: qsort                                             !SRTL 12
      real :: zsort1, zsort0                                                      !SRTL 13
!-----------------------------------------------                                  !SRTL 14
!...  Sorts level numbers, to put surface and sae level in right                  !SRTL 15
!     order                                                                       !SRTL 16
!...  Initialize level numbers                                                    !SRTL 17
      do i = 0, 18                                                                !SRTL 18
         lsort(i) = i                                                             !SRTL 19
      end do                                                                      !SRTL 20
!...  Store surface (level 1) and sea level (level 0) altitudes                   !SRTL 21
      zsort1 = zsort(1)                                                           !SRTL 22
      zsort0 = zsort(0)                                                           !SRTL 23
!...  Sort to put surface level number into right order                           !SRTL 24
      do i = 2, 18                                                                !SRTL 25
         if (zsort1 > zsort(i)) then                                              !SRTL 26
            lsort(i-1) = lsort(i)                                                 !SRTL 27
            lsort(i) = 1                                                          !SRTL 28
         else                                                                     !SRTL 29
            exit                                                                  !SRTL 30
         endif                                                                    !SRTL 31
      end do                                                                      !SRTL 32
!...  Store altitudes with surface altitude in right order                        !SRTL 33
      qsort(:18) = zsort(lsort(:18))                                              !SRTL 34
!...  Sort to put sea level number into right order                               !SRTL 35
      do i = 1, 18                                                                !SRTL 36
         if (zsort0 > qsort(i)) then                                              !SRTL 37
            lsort(i-1) = lsort(i)                                                 !SRTL 38
            lsort(i) = 0                                                          !SRTL 39
         else                                                                     !SRTL 40
            exit                                                                  !SRTL 41
         endif                                                                    !SRTL 42
      end do                                                                      !SRTL 43
      return                                                                      !SRTL 44
      end subroutine sortlevel                                                    !SRTL 45
!------------------------------------------------------------------------------   !SRTL 46
