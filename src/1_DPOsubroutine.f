




C  MPM93 - subroutines adapted by Jeff Haferman (NASA/GSFC 5/97)
C  from Liebe's MPM93 model.  His comments are included below.
C  I've based this adaptation on Frank Evans' MPM92 extraction.
C
C---------------------------------------------------------------------
C       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
C       June 1993
C       Hans J. Liebe     (303-497-3310)    
C       George A. Hufford (       -3457)
C       Michael G. Cotton (       -7346)
C       Institute for Telecommunication Sciences
C       NTIA/ITS.S3 
C       325 BROADWAY
C       Boulder, CO  80303,  USA
C
C       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
C       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV
C
C COMMENTS:
C 
C   The Millimeter-wave Propagation Model (MPM85) was reported in Ref.
C [1]. Molecular absorption by O2, H2O, and N2 is considered, as well as
C dielectric loss for haze and fog/cloud conditions (Rayleigh absorption
C approximation), and dielectric plus scatter losses (aR**b -
C approximation to Mie's theory) under rain conditions. The complex
C atmospheric refractivity N (or path-specific rates of attenuation A and
C delay B) were continued to be upgraded as discussed in [2] - [7].
C 
C   Features of the current version, MPM93, are:
C 
C - Haze model to predict the water droplet density for 
C       U = 80 to 99.95%RH , when a hygroscopic aerosol reference density
C       wa(80%RH) and a climatic code ('A, B, C, or D') are provided [2],[3]   
C 
C - Improved model for the dielectric properties of liquid water to
C       calculate RAYLEIGH absorption and delay by suspended water droplets
C       for haze, fog, and cloud conditions [6],[7]
C 
C - Rain attenuation model for Laws & Parsons drop-sizes by Olsen et al. 
C       [11], and associated dispersive delay, approximated from results 
C       reported by Zuffery [12]
C 
C - New temperature-dependent linewidth data (b3 to b6) for the water
C       vapor lines below 1 THz, and a 5 percent increase in the 
C       strength b1 of the 22-GHz and 183-GHz lines [9]
C 
C - New set of line mixing coefficients (a5, a6) for dry air, and 
C       their improved fit to the extensive 60-GHz lab. data [8],[9]
C 
C - Revised water vapor saturation pressure equation [10] 
C 
C - Approximation for Zeeman (O2) [4] and Doppler (H2O) line-broadening
C       to cover heights up to 100 km.
C 
C - New pseudo-line water vapor continuum formulation [9]   
C 
C - Detailed treatment of the anisotropic, mesospheric Zeeman effect
C   of O2 microwave lines [5]. The ZPM  code [9].
C 
C 
C                                 REFERENCES
C 
C  [1] H. Liebe, "An updated model for millimeter-wave propagation in
C       moist air", Radio Science, vol. 20, no. 5, pp. 1069-1089, 1985.
C 
C  [2] H. Liebe,"A contribution to modeling atmospheric mm-wave properties",
C       FREQUENZ, vol.41, no. 1/2, pp. 31-36, 1987.
C 
C  [3] H. Liebe and D. Layton, "MM-wave Properties of the Atmosphere:
C       Laboratory Studies and Propagation Modeling",
C       NTIA Report 87-224, 80p., Oct. 1987 (NTIS Order No. PB88-164215/AF).
C       
C  [4] H. Liebe,"MPM89 - An atmospheric mm-wave propagation model",
C       Int. J. IR & MM Waves, vol.10, no.6, pp. 631-650, June 1989.
C 
C  [5] G. Hufford and H. Liebe, "MM-Wave Propagation in the Mesosphere",
C       NTIA Report 89-249, 67p., Sept. 1989 (NTIS Order No. PB90-119868/AS).
C 
C  [6] H. Liebe, T. Manabe, and G. Hufford, "Mm-wave attenuation and delay
C       rates due to fog/cloud conditions", IEEE Trans. Ant. Prop.,
C       vol. 37, no. 12, pp. 1617-1623, Dec. 1989.
C 
C  [7] H. Liebe, G. Hufford (ice), and T. Manabe, "A model for the complex
C       refractivity of water (ice) at frequencies below 1 THz",
C       Int. J. IR & MM Waves, vol. 12, no. 7, 659-682, 1991.
C   
C  [8] H. Liebe, P. Rosenkranz, and G. Hufford, "Atmospheric 60-GHz   
C       oxygen spectrum: New laboratory measurements and line parameters", 
C       J. Quant. Spectr. Rad. Transf., vol. 48, no. 5/6, pp. 629-643, 1992.
C 
C  [9] H. Liebe, G. Hufford, and M. Cotton, "Propagation modeling of moist air 
C       and suspended water/ice particles at frequencies below 1000 GHz", 
C       Proc. AGARD Conf. Paper 3/1-10, Palma De Mallorca, Spain, May 1993.
C  
C [10] W. Boegel, "Neue Naeherungsgleichungen fuer den Saettigungsdruck des
C       Wasserdampfes, DFVLR Bericht DLR-FB 77-52, 1977.
C 
C [11] R.L. Olsen, D.V. Rogers, and D.B. Hodge, "The aRb relation in the
C       calculation of rain attenuation",
C       IEEE Trans. Ant. Prop., vol. AP-26, no. 2, pp. 318-329, 1978.
C 
C [12] C.H. Zuffery, "A study of rain effects on EM waves in the
C       1 to 600 GHz range", MS-THesis, Dept. Electrical Eng.,
C       University of Colorado, Boulder,  CO 80309, Feb., 1972.
C-----------------------------------------------------------------------
C 
      SUBROUTINE MPM93(F, Pbkpa, Ekpa, Tc, W, ABSCOF)
C****************************************************************** 
C Computes volume absorption coefficient for an atmospheric
C layer given the meteorological properties. The allowed frequency
C range is from 1 to 1000 GHz.  This routine is hacked from Liebe's
C GAS1 subroutine in his MPM93 model, taking out rain and dispersion
C computations.  Included is dry air attenuation, oxygen and "psuedo"
C water vapor line-continuum absorption, and Rayleigh cloud droplet
C absorption. 
C    Parameters:
C      F       frequency (GHz)
C      Pbkpa   total pressure (kPa)
C      Ekpa    water vapor pressure (kPa)
C      Tc      temperature (C)
C      W       cloud liquid water content (g/m^3)
C      ABSCOF  absorption coefficient (km^-1)
C****************************************************************** 
      IMPLICIT NONE
      REAL*8    F, Pbkpa, Ekpa, Tc, W, ABSCOF
      INTEGER IFIRST, I, ICE

      REAL*8 AT1, AT2, AT3, AT4
      REAL*8 GAMMA, S, DELTA, So, GAMMAo, Sn
      REAL*8 GAMH, GAMD2, DELH
      REAL*8 fD, fS, Eps, Epinf, Eopt
      REAL*8 Ai, Bi, fice
      REAL*8 V, P, Pb, E
      
      COMPLEX ZN, ZNw, ZEp, ZF, ZFo, ZFn

C Common block for oxygen and water vapor lines
      REAL*8    F0O2(44), A(6,44)
      REAL*8    F0H2O(35), B(6,35)
      real*8    A1(44),A2(44),A3(44),A4(44),A5(44),A6(44)
      real*8    B1(35),B2(35),B3(35),B4(35),B5(35),B6(35)
      ! COMMON /MWLINES2/ F0O2,A, F0H2O,B

      DATA IFIRST/0/
      DATA ICE/0/     ! hardcoded by JLH

c
c     The following data was in mwlines.93.data      
      data  F0O2 /  50.474239,   50.987747,   51.503349,   52.021412,
     $              52.542393,   53.066906,   53.595749,   54.130001,   
     $              54.671158,   55.221367,   55.783802,   56.264774,
     $              56.363388,   56.968204,   57.612484,   58.323875,
     $              58.446590,   59.164207,   59.590984,   60.306061,
     $              60.434776,   61.150558,   61.800156,   62.411217,
     $              62.486259,   62.997978,   63.568520,   64.127769,
     $              64.678902,   65.224068,   65.764771,   66.302094,
     $              66.836830,   67.369598,   67.900864,   68.431007,
     $              68.960312,  118.750343,  368.498352,  424.763123,
     $             487.249359,  715.393127,  773.839661,  834.145325 /
      data A1 /      0.094,    0.246,    0.608,    1.414,    3.102,
     $               6.410,   12.470,   22.800,   39.180,   63.160,
     $              95.350,   54.890,  134.400,  176.300,  214.100,
     $             238.600,  145.700,  240.400,  211.200,  212.400,
     $             246.100,  250.400,  229.800,  193.300,  151.700,
     $             150.300,  108.700,   73.350,   46.350,   27.480,
     $              15.300,    8.009,    3.946,    1.832,    0.801,
     $               0.330,    0.128,   94.500,    6.790,   63.800,
     $              23.500,    9.960,   67.100,   18.000  /
      data A2 /      9.694,    8.694,    7.744,    6.844,    6.004,
     $               5.224,    4.484,    3.814,    3.194,    2.624,
     $               2.119,    0.015,    1.660,    1.260,    0.915,
     $               0.626,    0.084,    0.391,    0.212,    0.212,
     $               0.391,    0.626,    0.915,    1.260,    0.083,
     $               1.665,    2.115,    2.620,    3.195,    3.815,
     $               4.485,    5.225,    6.005,    6.845,    7.745,
     $               8.695,    9.695,    0.009,    0.049,    0.044,
     $               0.049,    0.145,    0.130,    0.147 /
      data A3 /      0.890,    0.910,    0.940,    0.970,    0.990,
     $               1.020,    1.050,    1.070,    1.100,    1.130,
     $               1.170,    1.730,    1.200,    1.240,    1.280,
     $               1.330,    1.520,    1.390,    1.430,    1.450,
     $               1.360,    1.310,    1.270,    1.230,    1.540,
     $               1.200,    1.170,    1.130,    1.100,    1.070,
     $               1.050,    1.020,    0.990,    0.970,    0.940,
     $               0.920,    0.900,    1.630,    1.920,    1.930,
     $               1.920,    1.810,    1.820,    1.810 /
      data A4 /      0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.600,    0.600,
     $               0.600,    0.600,    0.600,    0.600 /
      data A5 /      0.240,    0.220,    0.197,    0.166,    0.136,
     $               0.131,    0.230,    0.335,    0.374,    0.258,
     $              -0.166,    0.390,   -0.297,   -0.416,   -0.613,
     $              -0.205,    0.748,   -0.722,    0.765,   -0.705,
     $               0.697,    0.104,    0.570,    0.360,   -0.498,
     $               0.239,    0.108,   -0.311,   -0.421,   -0.375,
     $              -0.267,   -0.168,   -0.169,   -0.200,   -0.228,
     $              -0.240,   -0.250,   -0.036,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000 /
      data A6 /      0.790,    0.780,    0.774,    0.764,    0.751,
     $               0.714,    0.584,    0.431,    0.305,    0.339,
     $               0.705,   -0.113,    0.753,    0.742,    0.697,
     $               0.051,   -0.146,    0.266,   -0.090,    0.081,
     $              -0.324,   -0.067,   -0.761,   -0.777,    0.097,
     $              -0.768,   -0.706,   -0.332,   -0.298,   -0.423,
     $              -0.575,   -0.700,   -0.735,   -0.744,   -0.753,
     $              -0.760,   -0.765,    0.009,    0.000,    0.000,
     $               0.000,    0.000,    0.000,    0.000 /
      data F0H2O / 22.235081,   67.803963,  119.995941,  183.310089,
     $            321.225647,  325.152924,  336.222595,  380.197357,
     $            390.134521,  437.346680,  439.150818,  443.018280,
     $            448.001068,  470.888947,  474.689117,  488.491119,
     $            503.568542,  504.482697,  547.676453,  552.020935,
     $            556.935974,  620.700806,  645.866150,  658.005310,
     $            752.033203,  841.053955,  859.962341,  899.306702,
     $            902.616150,  906.207336,  916.171570,  923.118408,
     $            970.315002,  987.926758, 1780.000000 /
      data B1 /     0.01130,    0.00012,    0.00008,    0.24200,   
     $              0.00483,    0.14990,    0.00011,    1.15200,
     $              0.00046,    0.00650,    0.09218,    0.01976,
     $              1.03200,    0.03297,    0.12620,    0.02520,
     $              0.00390,    0.00130,    0.97010,    1.47700,
     $             48.74000,    0.50120,    0.00713,    0.03022,
     $             23.96000,    0.00140,    0.01472,    0.00605,
     $              0.00426,    0.01876,    0.83400,    0.00869,
     $              0.89720,   13.21000, 2230.00000 /
      data B2 /   2.143,    8.735,    8.356,    0.668,    6.181,    
     $            1.540,    9.829,    1.048,    7.350,    5.050,
     $            3.596,    5.050,    1.405,    3.599,    2.381,
     $            2.853,    6.733,    6.733,    0.114,    0.114,
     $            0.159,    2.200,    8.580,    7.820,    0.396,
     $            8.180,    7.989,    7.917,    8.432,    5.111,
     $            1.442,   10.220,    1.920,    0.258,    0.952 /
      data B3 /   2.811,    2.858,    2.948,    3.050,    2.303,
     $            2.783,    2.693,    2.873,    2.152,    1.845,
     $            2.100,    1.860,    2.632,    2.152,    2.355,
     $            2.602,    1.612,    1.612,    2.600,    2.600,
     $            3.210,    2.438,    1.800,    3.210,    3.060,
     $            1.590,    3.060,    2.985,    2.865,    2.408,
     $            2.670,    2.900,    2.550,    2.985,   17.620 /
      data B4 /   4.80,    4.93,    4.78,    5.30,    4.69,    4.85,
     $            4.74,    5.38,    4.81,    4.23,    4.29,    4.23,
     $            4.84,    4.57,    4.65,    5.04,    3.98,    4.01,
     $            4.50,    4.50,    4.11,    4.68,    4.00,    4.14,
     $            4.09,    5.76,    4.09,    4.53,    5.10,    4.70,
     $            4.78,    5.00,    4.94,    4.55,   30.50 /
      data B5 /   0.69,    0.69,    0.70,    0.64,    0.67,    0.68,
     $            0.69,    0.54,    0.63,    0.60,    0.63,    0.60,
     $            0.66,    0.66,    0.65,    0.69,    0.61,    0.61,
     $            0.70,    0.70,    0.69,    0.71,    0.60,    0.69,
     $            0.68,    0.33,    0.68,    0.68,    0.70,    0.70,
     $            0.70,    0.70,    0.64,    0.68,    2.00 /
      data B6 /   1.00,    0.82,    0.79,    0.85,    0.54,    0.74,
     $            0.61,    0.89,    0.55,    0.48,    0.52,    0.50,
     $            0.67,    0.65,    0.64,    0.72,    0.43,    0.45,
     $            1.00,    1.00,    1.00,    0.68,    0.50,    1.00,
     $            0.84,    0.45,    0.84,    0.90,    0.95,    0.53,
     $            0.78,    0.80,    0.67,    0.90,    5.00 /

C---------------------------------------------------------------------
C
      do i = 1, 44
        A(1,i) = A1(i)
        A(2,i) = A2(i)
        A(3,i) = A3(i)
        A(4,i) = A4(i)
        A(5,i) = A5(i)
        A(6,i) = A6(i)
      end do
      do i = 1, 35
        B(1,i) = B1(i)
        B(2,i) = B2(i)
        B(3,i) = B3(i)
        B(4,i) = B4(i)
        B(5,i) = B5(i)
        B(6,i) = B6(i)
      end do
         
C Only read in line data the first time called
c      IF (IFIRST.EQ.0) THEN
c        IFIRST = 1
c	CALL READLINES2
c      ENDIF

C Relative inverse temperature
      V=300./(Tc+273.15)
C This version inputs E.
C Note MPM93 has pressure in mb, whereas MPM92 uses kPA
      Pb = 10.*Pbkpa
      E  = 10.*Ekpa
      P=Pb-E
      IF(P.LT.0)THEN
        P=0.
        Pb=E
      ENDIF

C For OXYGEN
      ZN=CMPLX(0.,0.)
      DO 10 I=1,44
       GAMMA=0.
       S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
       GAMMA=A(3,I)*(P*V**(0.8-A(4,I))+1.1*E*V)*1.E-3
       GAMMA=(GAMMA**2+(25*0.6E-4)**2)**0.5
       DELTA=(A(5,I)+A(6,I)*V)*(P+E)*(V**0.8)*1.E-3
       ZF=F/F0O2(I)*(CMPLX(1.,-DELTA)/CMPLX(F0O2(I)-F,-GAMMA)-
     +               CMPLX(1.,DELTA)/CMPLX(F0O2(I)+F,GAMMA))
       ZN=ZN+S*ZF
10    CONTINUE

C OXYGEN LINE ABSORPTION  
C Cannot be less than 0.  
      AT1=.182*F*AIMAG(ZN)
      IF(AT1.LT.0.) AT1=0.
C
C DRY AIR CONTINUUM
      ZN=CMPLX(0.,0.)
      So=6.14E-5*P*V**2
      GAMMAo=0.56E-3*(P+E)*V**.8
      ZFo=-F/CMPLX(F,GAMMAo)
      Sn=1.40E-12*p**2*V**3.5
      ZFn=CMPLX(0.,F/(1.93E-5*F**1.5+1.))
      ZN=So*ZFo+Sn*ZFn

C NONRESONAT DRY AIR ABSORPTION 
      AT2=.182*F*AIMAG(ZN)
C 
C WATER VAPOR
      ZN=CMPLX(0.,0.)
      DO 20 I=1,35
       GAMH=0.
       S=B(1,I)*E*V**3.5*EXP(B(2,I)*(1.-V))  
C Doppler approximation.
       GAMH=B(3,I)*(P*V**B(5,I)+B(4,I)*E*V**B(6,I))*1.E-3
       GAMD2=1E-12/V*(1.46*F0H2O(I))**2
       GAMH=0.535*GAMH+(0.217*GAMH**2+GAMD2)**0.5
       DELH=0.
       ZF=F/F0H2O(I)*(CMPLX(1.,-DELH)/CMPLX(F0H2O(I)-F,-GAMH)-
     +             CMPLX(1.,DELH)/CMPLX(F0H2O(I)+F,GAMH))
       ZN=ZN+S*ZF
20           CONTINUE

C WATER VAPOR LINE ABSORPTION 
C SEE LIEBE'S COMMENT REGARDING "PSUEDO-LINE WATER VAPOR CONTINUUM" - JLH
      AT3=.182*F*AIMAG(ZN)

C 
C LIQUID WATER PERMITTIVITY [8]
C Use exponential form for gamma for T<0 extrapolation (a la Frank Evans)
      IF(ICE.EQ.0)THEN
CJLH    fD=20.20-146.4*(V-1)+316*(V-1)**2
        fD=20.1*exp(7.88*(1-V))
        fS=39.8*fD 
        Eps=103.3*(V-1)+77.66
        Epinf=0.0671*Eps
        Eopt=3.52
C Complex Permittivity of water (double-Debye model)
        ZEp=Eps-f*((Eps-Epinf)/CMPLX(f,fD)+(Epinf-Eopt)/CMPLX(f,fS))
C
C ICE PERMITTIVITY [8]
      ELSE
        Ai=(62.*V-11.6)*1.E-4*EXP(-22.1*(V-1.))
        Bi=.542E-6*(-24.17+116.79/V+(V/(V-.9927))**2)
        Eps=3.15
C Complex Permittivity of Ice 
        fice=f
        IF(f.LT..001)fice=.001
        ZEp=CMPLX(3.15,Ai/fice+Bi*fice)
      ENDIF
C SUSPENDED PARTICLE RAYLEIGH APPROXIMATION [6]
      ZNw=1.5*W*((ZEp-1.)/(ZEp+2.)-1.+3./(Eps+2.))
C

C SUSPENDED WATER DROPLET EXTINCTION 
      AT4=.182*F*AIMAG(ZNw)

      ABSCOF=0.23026*(AT1+AT2+AT3+AT4)

      RETURN
      END 
C---------------------------------------------------------------------

          subroutine get_atmosG(temp_lev, press_lev, relhum, 
     +mxlyr, nlyr, avg_pressure, vapor_pressure,freq,ABSCOEF)

c     Calculate average air pressure and vapor pressure in specified
c     layers, given the temperature, pressure, and relative humidity 
c     from cloud-resolving model output.
c     vapor_pressure viene data in mb =hPa
C+-------+---------+---------+---------+---------+---------+---------+-+
C+----------------------------------------------------------------------
      implicit none
      integer mxlyr,nlyr,nz
      real*8       press_lev(0:mxlyr)
      real*8       temp_lev(0:mxlyr)
      real*8       relhum(mxlyr), ABScoef(mxlyr)
      real*8       avg_pressure(mxlyr)
      real*8       vapor_pressure(mxlyr)
     
      real*8   freq,tavg,tc,es,esatvap,a0,a1,a2,a3,a4,a5,a6

      data a0/6.107799961e0/
      data a1/4.436518521e-1/
      data a2/1.428945805e-2/
      data a3/2.650648471e-4/
      data a4/3.031240396e-6/
      data a5/2.034080948e-8/
      data a6/6.136820929e-11/

    
          do nz = 1, nlyr
            tavg = 0.5*(temp_lev(nz-1) + temp_lev(nz))
            avg_pressure(nz) = 
     +          (press_lev(nz) - press_lev(nz-1) )/
     +           log(press_lev(nz)/press_lev(nz-1))
            tc = tavg - 273.15
            ES = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
            if (ES .lt. 0.) ES = 0.
            vapor_pressure(nz) =  relhum(nz)*ES/100.
           call MPM93(FREQ, avg_pressure(nz),0.1*vapor_pressure(nz),Tc, 
     $ 0.0d0, ABSCOEF(nz))
c             write(18,*) 'MPM93',FREQ,avg_pressure(nz),
c     $0.1*vapor_pressure(nz),Tc,ABSCOEF(nz), relhum(nz)
           end do

      return
      end       



          subroutine get_atmosGlev(temp_lev, press_lev, relhum, 
     +mxlyr, nlyr, vapor_pressure,freq,ABSCOEF)

c     Calculate average air pressure and vapor pressure in specified
c     layers, given the temperature, pressure, and relative humidity 
c     from cloud-resolving model output.
c     vapor_pressure viene data in mb =hPa
C+-------+---------+---------+---------+---------+---------+---------+-+
C+----------------------------------------------------------------------
      implicit none
      integer mxlyr,nlyr,nz
      real*8       press_lev(0:mxlyr)
      real*8       temp_lev(0:mxlyr)
      real*8       relhum(0:mxlyr), ABScoef(0:mxlyr)
      real*8       vapor_pressure(0:mxlyr)
     
      real*8   freq,tavg,tc,es,esatvap,a0,a1,a2,a3,a4,a5,a6

      data a0/6.107799961e0/
      data a1/4.436518521e-1/
      data a2/1.428945805e-2/
      data a3/2.650648471e-4/
      data a4/3.031240396e-6/
      data a5/2.034080948e-8/
      data a6/6.136820929e-11/


          do nz = 0, nlyr
            tc = temp_lev(nz) - 273.15
            ES = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
            if (ES .lt. 0.) ES = 0.
            vapor_pressure(nz) =  relhum(nz)*ES/100.
           call MPM93(FREQ, press_lev(nz),0.1*vapor_pressure(nz),Tc, 
     $ 0.0d0, ABSCOEF(nz))
c             write(18,*) 'MPM93',FREQ,avg_pressure(nz),
c     $0.1*vapor_pressure(nz),Tc,ABSCOEF(nz), relhum(nz)
           end do

      return
      end       




