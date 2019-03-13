C********************************************************************
      SUBROUTINE  DIECON(S,T,FREQ,E1,E2)
C***********************************************************************
C compute the dieletric constant (E1,E2) of a SURFACE given 
C INPUT 
c     S=salinity in PPM
C     T=temperature in Celsius
C     freq=frequency in GHz                                 

       implicit none
       real*8  ST,S,T,FREQ,S2,T2,SST,STT,SSTT,ES,TAU,
     $ SIGMA,ZNU,OMEGA,DEN,TWOPI,EINF,SOLD,TOLD,E1,E2

      DATA TWOPI /6.283185307/ , EINF /4.9/ , SOLD /0.0/ , TOLD /-99./
C      save es, tau,sigma
C      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10
      ST = S*T
      S2 = S*S
      T2 = T*T
      SST = S2*T
      STT = T2*S
      SSTT = S2*T2
      ES = 88.-4.339E-01*S+1.71E-03*S2-4.035E-01*T+8.065E-04*T2+6.170
     $ E-03 * ST-8.910E-05*SST-6.934E-05*STT+1.439E-06*SSTT
      TAU = (18.70-7.924E-02*S+6.35E-04*S2-5.489E-01*T+5.758E-03*T2+
     $ 1.889E-03*ST-7.209E-06*SST-5.299E-07*STT-2.101E-07*SSTT)*1.0E-12
      SIGMA = (7.788E-03*S-1.672E-06*S2-8.570E-15*T+2.996E-16*T2+
     $ 4.059E-04 * ST-3.215E-06*SST-1.423E-06*STT+3.229E-08*SSTT)*1.0E11
  10  ZNU = FREQ*1.E09
C      write(18,*)'tau,es,sigma',tau,es,sigma
      OMEGA = TWOPI*ZNU
      DEN = 1. + (OMEGA * TAU) ** 2
      E1 = (ES-EINF)/DEN+EINF
      E2 = (ES-EINF)*OMEGA*TAU/DEN+2.*SIGMA/ZNU
      SOLD = S
      TOLD = T
      RETURN
      END

C********************************************************************
     

C********************************************************************
      SUBROUTINE  EMITTanceSURF(E1,E2,EBAR)
C********************************************************************
C  INPUT 
C        (E1,E2) dielectric constant of the surface
C  OUTPUT
C       EBAR= emissivity of the surface as defined in eq.6 Battaglia&Mantovani
C 
      implicit none
      integer Nang,I
      PARAMETER (NANG = 90)
      CHARACTER POLN*1
      real*8  pi,E1,E2,EBAR,angles,EV,EH,
     $ SUM, EAVG,DMU,AVMU
      REAL*8   ANG(NANG), MU(NANG), ESUM(NANG)
C      write(18,*)'F,poln,Ts,W,UMU',F,poln,Ts,W,UMU
      PI = 2.*ASIN(1.0)
C     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
        ANG(I) = I-1
        MU(I) = COS(ANG(I)*PI/180.)
        ANGLES = ANG(I)
        CALL EMISSIVITY(E1,E2,ANGLES,EV,EH)
        ESUM(I) =  EV + EH
C        write(18,*)'esum',esum(i)
  58  CONTINUE

C *** CALCULATE EBAR
      SUM = 0.0
      DO 59 I = 1,NANG-1
        EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
        DMU = MU(I) - MU(I+1)
        AVMU = 0.5*( MU(I) + MU(I+1) )
        SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
C      write(18,*)'ebar',ebar
c 104  FORMAT(1X,F8.2,F10.2,F9.3,A6,F7.2,2F11.5)
      RETURN
      END 

C********************************************************************
      SUBROUTINE  EMISSivity(E1,E2,ANGLE,EMISV,EMISH)
C********************************************************************
C INPUT 
c     E1,E2 = real and imaginary part of dielectric constant 
C     ANGLE = angle in degrees

C OUTPUT
C      EMISV, EMISH emissivity for V-H polarization at the input angle
C
      implicit none
      real*8  e1,e2,CMHU,CMHU2,FACT1,B, B2, ARG,PSIZ,COSPZ,
     $ F11,PSIY,G,COSPY,F22,EMISH,EMISV,
     $ DTR,theta,angle

      DATA DTR / 0.01745329252 /
      THETA = ANGLE*DTR
      CMHU = COS(THETA)
C      write(18,*)'eps',e1,e2,cmhu
      CMHU2 = CMHU*CMHU
      FACT1 = E1+CMHU2-1.
      B = (FACT1*FACT1+E2*E2)**0.25
      B2 = B*B
      ARG = E2/FACT1
      PSIZ = 0.5*DATAN(ARG)
      COSPZ = COS(PSIZ)*CMHU*B
      F11 = CMHU2+B2+2.*COSPZ
      PSIY = DATAN(E2/E1)
      G = CMHU*SQRT(E1*E1+E2*E2)
      COSPY = COS(PSIY-PSIZ)*G*B
      F22 = G*G+B2+2.*COSPY
C *** FOR SPECULAR SURFACES THE HORIZONTAL AND
C     VERTICAL EMISSIVITY CALCULATION
C
      EMISH = 4.*COSPZ/F11
      EMISV = 4.*COSPY/F22
      RETURN
      END
  




 


C********************************************************************
      SUBROUTINE  Reflect_matrix(E1,E2,ANGLE,
     $ a11,a22,a33,a44,a34,a43)
C********************************************************************
C INPUT 
c     E1,E2 = real and imaginary part of dielectric constant 
C     ANGLE = angle in degrees

C OUTPUT
c  elements of the Fresnel matrix at ANGLE
c   a11=RSH,  a22=RSV
C

      implicit none
      real*8  E1,E2,a11,a22,a33,a44,a34,a43,
     $DTR,theta,angle,cmhu,sincompl
       
       complex  epsilon,Rfrev,Rfreh
        DATA DTR / 0.01745329252 /   

      THETA = ANGLE*DTR
      CMHU = COS(THETA)
C      write(18,*)'eps',e1,e2,cmhu
C
      epsilon=e1*(1.0,0.0)+e2*(0.0,1.0)
      sincompl=(1.0-CMHU*CMHU)
      Rfrev= (epsilon*CMHU-sqrt(epsilon-sincompl))/
     $ (epsilon*CMHU+sqrt(epsilon-sincompl))
      Rfreh=(CMHU*(1.0,0.0)-sqrt(epsilon-sincompl))/
     $ (CMHU*(1.0,0.0)+sqrt(epsilon-sincompl))
      a11=(abs(Rfrev))*(abs(Rfrev))
      a22=(abs(Rfreh))*(abs(Rfreh))
      a33=real(Rfrev*conjg(Rfreh))
      a44=real(Rfrev*conjg(Rfreh))
      a34=-imag(Rfrev*conjg(Rfreh))
      a43=+imag(Rfrev*conjg(Rfreh))  !Tsang convention see eq.7 pg 242 (note
                       ! opposite sign respect to Mischenko, eq. 37 pag 511

      RETURN
      END
