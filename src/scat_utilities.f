

      SUBROUTINE MIE (WAVELENGTH, MINDEX, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .            EXTINCTION, ALBEDO,Back_scatt, NLEGEN, LEGEN,
     $LEGEN2,LEGEN3,LEGEN4, aerodist)
c note that Mindex has the convention with negative imaginary part


C       Computes the Mie scattering properties for a gamma or lognormal
C     distribution of spheres.
      IMPLICIT NONE
      INTEGER     MAXLEG, NLEGEN, NUMRAD
      LOGICAL     PHASEFLAG
      REAL*8        WAVELENGTH, RAD1, RAD2
      REAL*8        AD, BD, ALPHA, GAMMA
      COMPLEX*16     MINDEX
      REAL*8        EXTINCTION, ALBEDO,Back_scatt,LEGEN(200),
     $ LEGEN2(200), LEGEN3(200), LEGEN4(200)
      INTEGER     MAXN
      PARAMETER   (MAXN=5000)
      REAL*8      PI
      PARAMETER   (PI = 3.14159265358979D0)
      INTEGER     NTERMS, NQUAD, NMIE, NLEG
      INTEGER     I, L, M, IR, MIN0
      REAL*8      X, DELRAD, RADIUS, NDENS, TMP
      REAL*8      QEXT, QSCAT,Qback,SCATTER
      REAL*8      DISTRIBUTION
      REAL*8      MU(MAXN), WTS(MAXN)
      REAL*8      P1, PL, PL1, PL2,P2,P3,P4
      REAL*8      SUMQE, SUMQS,  SUMQback,ONE
      REAL*8      SUMP1(MAXN), COEF1(MAXN),SUMP2(MAXN), COEF2(MAXN),
     $            SUMP3(MAXN), COEF3(MAXN),SUMP4(MAXN), COEF4(MAXN)
      COMPLEX*16  A(MAXN), B(MAXN), MSPHERE
      CHARACTER aerodist*1
 
 
C           Find the maximum number of terms required in the Mie series,
      MSPHERE = MINDEX
      X = 2.0D0*PI*RAD2/WAVELENGTH
      NTERMS = 0
      CALL MIECALC (NTERMS, X, MSPHERE, A, B)
      NLEGEN = 2*NTERMS
      NLEGEN = MIN0(MAXLEG, NLEGEN)
      NQUAD  = (NLEGEN + 2*NTERMS + 2)/2
      IF (NQUAD .GT. MAXN)  STOP 'MIE: MAXN exceeded' 
 
C           Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, MU, WTS)
 
      SUMQE = 0.0d0
      SUMQS = 0.0d0
      SUMQback=0.0d0
      DO I = 1, NQUAD
        SUMP1(I) = 0.0d0
        SUMP2(I) = 0.0d0
        SUMP3(I) = 0.0d0
        SUMP4(I) = 0.0d0
      ENDDO
 
C               Integration loop over radius of spheres
      IF (NUMRAD .GT. 0)  DELRAD = (RAD2-RAD1)/NUMRAD
      DO IR = 1, NUMRAD+1
          RADIUS = RAD1 + (IR-1)*DELRAD
          NDENS = DISTRIBUTION (DBLE(AD), DBLE(BD), 
     .                     DBLE(ALPHA), DBLE(GAMMA), RADIUS, aerodist)
          IF ((IR .EQ. 1 .OR. IR .EQ. NUMRAD+1)
     .                      .AND. NUMRAD .GT. 0) THEN
              NDENS = 0.5*NDENS
          ENDIF
          X = 2.0D0*PI*RADIUS/WAVELENGTH
          NMIE = 0
          CALL MIECALC (NMIE, X, MSPHERE, A, B)
          CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT,Qback)
          SUMQE = SUMQE + QEXT*NDENS*RADIUS**2
          SUMQS = SUMQS + QSCAT*NDENS*RADIUS**2
          SUMQback = SUMQback + Qback*NDENS*RADIUS**2
c          write(*,*)'Mie in',Qext,qscat,Ndens,radius,x,nmie
          IF (PHASEFLAG) THEN
            NMIE = MIN0(NMIE, NTERMS)
            DO I = 1, NQUAD
              CALL MIEANGLE (NMIE, A, B, MU(I), P1,P2,P3,P4)
              SUMP1(I) = SUMP1(I) + P1*NDENS 
              SUMP2(I) = SUMP2(I) + P2*NDENS
              SUMP3(I) = SUMP3(I) + P3*NDENS
              SUMP4(I) = SUMP4(I) + P4*NDENS
            ENDDO
          ENDIF
        ENDDO
 
 
C           Multiply the sums by the integration delta and other constants
C             Put quadrature weights in angular array for later
      IF (NUMRAD .EQ. 0) DELRAD = 1.0d0
 
      EXTINCTION = PI*SUMQE*DELRAD
      SCATTER = PI*SUMQS*DELRAD
      Back_scatt= PI*SUMQback*DELRAD
      ALBEDO = SCATTER/EXTINCTION
 
C         If the phase function is not desired then leave now
      IF (.NOT. PHASEFLAG) RETURN
 
      TMP = (WAVELENGTH**2/(PI*SCATTER)) * DELRAD
      DO I = 1, NQUAD
        SUMP1(I) = TMP*SUMP1(I) *WTS(I)
        SUMP2(I) = TMP*SUMP2(I) *WTS(I)
        SUMP3(I) = TMP*SUMP3(I) *WTS(I)
        SUMP4(I) = TMP*SUMP4(I) *WTS(I)
      ENDDO
 
C           Integrate the angular scattering functions times Legendre
C             polynomials to find the Legendre coefficients
      DO M = 1, NLEGEN+1
        COEF1(M) = 0.0d0
        COEF2(M) = 0.0d0
        COEF3(M) = 0.0d0
        COEF4(M) = 0.0d0
      ENDDO
C           Use upward recurrence to find Legendre polynomials
      DO I = 1, NQUAD
          PL1 = 1.0d0
          PL = 1.0d0
          DO L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
              COEF1(M) = COEF1(M) + SUMP1(I)*PL 
              COEF2(M) = COEF2(M) + SUMP2(I)*PL 
              COEF3(M) = COEF3(M) + SUMP3(I)*PL 
              COEF4(M) = COEF4(M) + SUMP4(I)*PL
              PL2 = PL1
              PL1 = PL
          ENDDO
      ENDDO
      NLEG = NLEGEN
      DO L = 0, NLEG
          M = L + 1
          LEGEN(M) = (2*L+1)/2.0 *COEF1(M) 
          LEGEN2(M)= (2*L+1)/2.0 *COEF2(M) 
          LEGEN3(M)= (2*L+1)/2.0 *COEF3(M) 
          LEGEN4(M)= (2*L+1)/2.0 *COEF4(M)
          IF (LEGEN(M) .GT. 1.0E-7)  NLEGEN = L
      ENDDO
       
      RETURN
      END
 
         subroutine density_ice(a_mtox,bcoeff,rad2,density)
c      a_mtox  in kg/m^b
c      bcoeff   adimensional
c      rad2 in mm
c      density in kg/dm^3        
        implicit none
        real*8 a_mtox,bcoeff,rad2,density,pi

        pi=dacos(-1.0d0)
        density=6.0d0/pi*1e-3*a_mtox*(2.0*rad2*1e-3)**(bcoeff-3.0d0)
        density=min(0.917,density)
        return  
        end



cc    computing the scattering properties according to 
cc    soft sphere model, i.e. the electromegnetic properties of the 
c     particle by using the Maxwell Garnett model for refractive index   
       SUBROUTINE MIE_densitysizedep_softsphere(WAVELENGTH, 
     $  M_Ice,m_air,
     $ a_mtox,bcoeff,RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .            EXTINCTION, ALBEDO,Back_scatt, NLEGEN, LEGEN,
     $LEGEN2,LEGEN3,LEGEN4, aerodist)
c note that Mindex has the convention with negative imaginary part


C       Computes the Mie scattering properties for a gamma or lognormal
C     distribution of spheres.
      IMPLICIT NONE
      INTEGER     MAXLEG, NLEGEN, NUMRAD
      LOGICAL     PHASEFLAG
      REAL*8        WAVELENGTH, RAD1, RAD2
      REAL*8        AD, BD, ALPHA, GAMMA
      COMPLEX*16     M_ice,m_air,m_MG
      real*8 a_mtox,bcoeff,dens_graup,fvol_ice
      REAL*8        EXTINCTION, ALBEDO,Back_scatt,LEGEN(200),
     $ LEGEN2(200), LEGEN3(200), LEGEN4(200)
      INTEGER     MAXN
      PARAMETER   (MAXN=5000)
      REAL*8      PI
      PARAMETER   (PI = 3.14159265358979D0)
      INTEGER     NTERMS, NQUAD, NMIE, NLEG
      INTEGER     I, L, M, IR, MIN0
      REAL*8      X, DELRAD, RADIUS, NDENS, TMP
      REAL*8      QEXT, QSCAT,Qback,SCATTER
      REAL*8      DISTRIBUTION
      REAL*8      MU(MAXN), WTS(MAXN)
      REAL*8      P1, PL, PL1, PL2,P2,P3,P4
      REAL*8      SUMQE, SUMQS,  SUMQback,ONE
      REAL*8      SUMP1(MAXN), COEF1(MAXN),SUMP2(MAXN), COEF2(MAXN),
     $            SUMP3(MAXN), COEF3(MAXN),SUMP4(MAXN), COEF4(MAXN)
      COMPLEX*16  A(MAXN), B(MAXN), MSPHERE
      CHARACTER aerodist*1
 
 
C           Find the maximum number of terms required in the Mie series,
      call density_ice(a_mtox,bcoeff,rad2,dens_graup)
      fvol_ice=dens_graup/0.917
      call Max_Garn(m_ice,m_air,1-fvol_ice,m_MG)
      MSPHERE =conjg(m_MG) 
      X = 2.0D0*PI*RAD2/WAVELENGTH
      NTERMS = 0
      CALL MIECALC (NTERMS, X, MSPHERE, A, B)
      NLEGEN = 2*NTERMS
      NLEGEN = MIN0(MAXLEG, NLEGEN)
      NQUAD  = (NLEGEN + 2*NTERMS + 2)/2
      IF (NQUAD .GT. MAXN)  STOP 'MIE: MAXN exceeded' 
 
C           Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, MU, WTS)
 
      SUMQE = 0.0d0
      SUMQS = 0.0d0
      SUMQback=0.0d0
      DO I = 1, NQUAD
        SUMP1(I) = 0.0d0
        SUMP2(I) = 0.0d0
        SUMP3(I) = 0.0d0
        SUMP4(I) = 0.0d0
      ENDDO
 
C               Integration loop over radius of spheres
      IF (NUMRAD .GT. 0)  DELRAD = (RAD2-RAD1)/NUMRAD
      DO IR = 1, NUMRAD+1
          RADIUS = RAD1 + (IR-1)*DELRAD
          NDENS = DISTRIBUTION (DBLE(AD), DBLE(BD), 
     .                     DBLE(ALPHA), DBLE(GAMMA), RADIUS, aerodist)
          IF ((IR .EQ. 1 .OR. IR .EQ. NUMRAD+1)
     .                      .AND. NUMRAD .GT. 0) THEN
              NDENS = 0.5*NDENS
          ENDIF
          X = 2.0D0*PI*RADIUS/WAVELENGTH
          NMIE = 0
          call density_ice(a_mtox,bcoeff,radius,dens_graup)
c         write(18,*)'dens',dens_graup
         fvol_ice=dens_graup/0.917
         call Max_Garn(m_ice,m_air,1-fvol_ice,m_MG)
         MSPHERE =conjg(m_MG) 
          CALL MIECALC (NMIE, X, MSPHERE, A, B)
          CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT,Qback)
          SUMQE = SUMQE + QEXT*NDENS*RADIUS**2
          SUMQS = SUMQS + QSCAT*NDENS*RADIUS**2
          SUMQback = SUMQback + Qback*NDENS*RADIUS**2
c          write(*,*)'Mie in',Qext,qscat,Ndens,radius,x,nmie
          IF (PHASEFLAG) THEN
            NMIE = MIN0(NMIE, NTERMS)
            DO I = 1, NQUAD
              CALL MIEANGLE (NMIE, A, B, MU(I), P1,P2,P3,P4)
              SUMP1(I) = SUMP1(I) + P1*NDENS 
              SUMP2(I) = SUMP2(I) + P2*NDENS
              SUMP3(I) = SUMP3(I) + P3*NDENS
              SUMP4(I) = SUMP4(I) + P4*NDENS
            ENDDO
          ENDIF
        ENDDO
 
 
C           Multiply the sums by the integration delta and other constants
C             Put quadrature weights in angular array for later
      IF (NUMRAD .EQ. 0) DELRAD = 1.0d0
 
      EXTINCTION = PI*SUMQE*DELRAD
      SCATTER = PI*SUMQS*DELRAD
      Back_scatt= PI*SUMQback*DELRAD
      ALBEDO = SCATTER/EXTINCTION
 
C         If the phase function is not desired then leave now
      IF (.NOT. PHASEFLAG) RETURN
 
      TMP = (WAVELENGTH**2/(PI*SCATTER)) * DELRAD
      DO I = 1, NQUAD
        SUMP1(I) = TMP*SUMP1(I) *WTS(I)
        SUMP2(I) = TMP*SUMP2(I) *WTS(I)
        SUMP3(I) = TMP*SUMP3(I) *WTS(I)
        SUMP4(I) = TMP*SUMP4(I) *WTS(I)
      ENDDO
 
C           Integrate the angular scattering functions times Legendre
C             polynomials to find the Legendre coefficients
      DO M = 1, NLEGEN+1
        COEF1(M) = 0.0d0
        COEF2(M) = 0.0d0
        COEF3(M) = 0.0d0
        COEF4(M) = 0.0d0
      ENDDO
C           Use upward recurrence to find Legendre polynomials
      DO I = 1, NQUAD
          PL1 = 1.0d0
          PL = 1.0d0
          DO L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
              COEF1(M) = COEF1(M) + SUMP1(I)*PL 
              COEF2(M) = COEF2(M) + SUMP2(I)*PL 
              COEF3(M) = COEF3(M) + SUMP3(I)*PL 
              COEF4(M) = COEF4(M) + SUMP4(I)*PL
              PL2 = PL1
              PL1 = PL
          ENDDO
      ENDDO
      NLEG = NLEGEN
      DO L = 0, NLEG
          M = L + 1
          LEGEN(M) = (2*L+1)/2.0 *COEF1(M) 
          LEGEN2(M)= (2*L+1)/2.0 *COEF2(M) 
          LEGEN3(M)= (2*L+1)/2.0 *COEF3(M) 
          LEGEN4(M)= (2*L+1)/2.0 *COEF4(M)
          IF (LEGEN(M) .GT. 1.0E-7)  NLEGEN = L
      ENDDO
       
      RETURN
      END

cc    computing the scattering properties according to 
cc    ice sphere model, i.e. the electromegnetic properties of the 
c     particle are computed by assuming that they are the same 
c     as the equivalent mass sphere 
   
       SUBROUTINE MIE_densitysizedep_spheremasseq(WAVELENGTH, 
     $ M_Ice,m_air, a_mtox,bcoeff,RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .            EXTINCTION, ALBEDO,Back_scatt, NLEGEN, LEGEN,
     $LEGEN2,LEGEN3,LEGEN4, aerodist)
c note that Mindex has the convention with negative imaginary part


C       Computes the Mie scattering properties for a gamma or lognormal
C     distribution of spheres.
      IMPLICIT NONE
      INTEGER     MAXLEG, NLEGEN, NUMRAD
      LOGICAL     PHASEFLAG
      REAL*8        WAVELENGTH, RAD1, RAD2
      REAL*8        AD, BD, ALPHA, GAMMA
      COMPLEX*16     M_ice,m_air,m_MG
      real*8 a_mtox,bcoeff,dens_graup,fvol_ice
      REAL*8        EXTINCTION, ALBEDO,Back_scatt,LEGEN(200),
     $ LEGEN2(200), LEGEN3(200), LEGEN4(200)
      INTEGER     MAXN
      PARAMETER   (MAXN=5000)
      REAL*8      PI
      PARAMETER   (PI = 3.14159265358979D0)
      INTEGER     NTERMS, NQUAD, NMIE, NLEG
      INTEGER     I, L, M, IR, MIN0
      REAL*8      X, DELRAD, RADIUS, NDENS, TMP,
     $  radius_ice,rad2_ice
      REAL*8      QEXT, QSCAT,Qback,SCATTER
      REAL*8      DISTRIBUTION
      REAL*8      MU(MAXN), WTS(MAXN)
      REAL*8      P1, PL, PL1, PL2,P2,P3,P4
      REAL*8      SUMQE, SUMQS,  SUMQback,ONE
      REAL*8      SUMP1(MAXN), COEF1(MAXN),SUMP2(MAXN), COEF2(MAXN),
     $            SUMP3(MAXN), COEF3(MAXN),SUMP4(MAXN), COEF4(MAXN)
      COMPLEX*16  A(MAXN), B(MAXN), MSPHERE
      CHARACTER aerodist*1
 
 
C           Find the maximum number of terms required in the Mie series,
      call density_ice(a_mtox,bcoeff,rad2,dens_graup)
      rad2_ice=(dens_graup/0.917)**0.33333333*rad2
      MSPHERE =conjg(m_ice) 
      X = 2.0D0*PI*RAD2_ice/WAVELENGTH
      NTERMS = 0
      CALL MIECALC (NTERMS, X, MSPHERE, A, B)
      NLEGEN = 2*NTERMS
      NLEGEN = MIN0(MAXLEG, NLEGEN)
      NQUAD  = (NLEGEN + 2*NTERMS + 2)/2
      IF (NQUAD .GT. MAXN)  STOP 'MIE: MAXN exceeded' 
 
C           Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, MU, WTS)
 
      SUMQE = 0.0d0
      SUMQS = 0.0d0
      SUMQback=0.0d0
      DO I = 1, NQUAD
        SUMP1(I) = 0.0d0
        SUMP2(I) = 0.0d0
        SUMP3(I) = 0.0d0
        SUMP4(I) = 0.0d0
      ENDDO
 
C               Integration loop over radius of spheres
      IF (NUMRAD .GT. 0)  DELRAD = (RAD2-RAD1)/NUMRAD
      DO IR = 1, NUMRAD+1
          RADIUS = RAD1 + (IR-1)*DELRAD
          NDENS = DISTRIBUTION (DBLE(AD), DBLE(BD), 
     .                     DBLE(ALPHA), DBLE(GAMMA), RADIUS, aerodist)
          IF ((IR .EQ. 1 .OR. IR .EQ. NUMRAD+1)
     .                      .AND. NUMRAD .GT. 0) THEN
              NDENS = 0.5*NDENS
          ENDIF
        
          NMIE = 0
         
          call density_ice(a_mtox,bcoeff,radius,dens_graup)
c         write(18,*)'dens',dens_graup
          radius_ice=(dens_graup/0.917)**0.33333333*radius

          X = 2.0D0*PI*RADIUS_ice/WAVELENGTH

          CALL MIECALC (NMIE, X, MSPHERE, A, B)
          CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT,Qback)
          SUMQE = SUMQE + QEXT*NDENS*RADIUS_ice**2
          SUMQS = SUMQS + QSCAT*NDENS*RADIUS_ice**2
          SUMQback = SUMQback + Qback*NDENS*RADIUS_ice**2
c          write(*,*)'Mie in',Qext,qscat,Ndens,radius,x,nmie
          IF (PHASEFLAG) THEN
            NMIE = MIN0(NMIE, NTERMS)
            DO I = 1, NQUAD
              CALL MIEANGLE (NMIE, A, B, MU(I), P1,P2,P3,P4)
              SUMP1(I) = SUMP1(I) + P1*NDENS 
              SUMP2(I) = SUMP2(I) + P2*NDENS
              SUMP3(I) = SUMP3(I) + P3*NDENS
              SUMP4(I) = SUMP4(I) + P4*NDENS
            ENDDO
          ENDIF
        ENDDO
 
 
C           Multiply the sums by the integration delta and other constants
C             Put quadrature weights in angular array for later
      IF (NUMRAD .EQ. 0) DELRAD = 1.0d0
 
      EXTINCTION = PI*SUMQE*DELRAD
      SCATTER = PI*SUMQS*DELRAD
      Back_scatt= PI*SUMQback*DELRAD
      ALBEDO = SCATTER/EXTINCTION
 
C         If the phase function is not desired then leave now
      IF (.NOT. PHASEFLAG) RETURN
 
      TMP = (WAVELENGTH**2/(PI*SCATTER)) * DELRAD
      DO I = 1, NQUAD
        SUMP1(I) = TMP*SUMP1(I) *WTS(I)
        SUMP2(I) = TMP*SUMP2(I) *WTS(I)
        SUMP3(I) = TMP*SUMP3(I) *WTS(I)
        SUMP4(I) = TMP*SUMP4(I) *WTS(I)
      ENDDO
 
C           Integrate the angular scattering functions times Legendre
C             polynomials to find the Legendre coefficients
      DO M = 1, NLEGEN+1
        COEF1(M) = 0.0d0
        COEF2(M) = 0.0d0
        COEF3(M) = 0.0d0
        COEF4(M) = 0.0d0
      ENDDO
C           Use upward recurrence to find Legendre polynomials
      DO I = 1, NQUAD
          PL1 = 1.0d0
          PL = 1.0d0
          DO L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
              COEF1(M) = COEF1(M) + SUMP1(I)*PL 
              COEF2(M) = COEF2(M) + SUMP2(I)*PL 
              COEF3(M) = COEF3(M) + SUMP3(I)*PL 
              COEF4(M) = COEF4(M) + SUMP4(I)*PL
              PL2 = PL1
              PL1 = PL
          ENDDO
      ENDDO
      NLEG = NLEGEN
      DO L = 0, NLEG
          M = L + 1
          LEGEN(M) = (2*L+1)/2.0 *COEF1(M) 
          LEGEN2(M)= (2*L+1)/2.0 *COEF2(M) 
          LEGEN3(M)= (2*L+1)/2.0 *COEF3(M) 
          LEGEN4(M)= (2*L+1)/2.0 *COEF4(M)
          IF (LEGEN(M) .GT. 1.0E-7)  NLEGEN = L
      ENDDO
       
      RETURN
      END

       SUBROUTINE MIE_densitysizedep_parameterization(WAVELENGTH, 
     $ M_Ice,m_air, a_mtox,bcoeff, type,
     $ RAD1,RAD2,NUMRAD,MAXLEG,AD,BD,ALPHA,GAMMA,PHASEFLAG,
     .  EXTINCTION, ALBEDO,Back_scatt, NLEGEN, LEGEN,
     $LEGEN2,LEGEN3,LEGEN4, aerodist)
c note that Mindex has the convention with negative imaginary part


C       Computes the Mie scattering properties for a gamma or lognormal
C     distribution of spheres.
      IMPLICIT NONE
      INTEGER     MAXLEG, NLEGEN, NUMRAD
      LOGICAL     PHASEFLAG
      REAL*8        WAVELENGTH, RAD1, RAD2
      REAL*8        AD, BD, ALPHA, GAMMA
      COMPLEX*16     M_ice,m_air,m_MG
      real*8 a_mtox,bcoeff,dens_graup,fvol_ice
      REAL*8        EXTINCTION, ALBEDO,Back_scatt,LEGEN(200),
     $ LEGEN2(200), LEGEN3(200), LEGEN4(200)
      INTEGER     MAXN
      PARAMETER   (MAXN=5000)
      REAL*8      PI
      PARAMETER   (PI = 3.14159265358979D0)
      INTEGER     NTERMS, NQUAD, NMIE, NLEG
      INTEGER     I, L, M, IR, MIN0
      REAL*8      X, DELRAD, RADIUS, NDENS, TMP,
     $  radius_ice,rad2_ice
      REAL*8      QEXT, QSCAT,Qback,SCATTER,asym,SCATTERFalse
      REAL*8      DISTRIBUTION
      REAL*8      MU(MAXN), WTS(MAXN)
      REAL*8      P1, PL, PL1, PL2,P2,P3,P4
      REAL*8      SUMQE, SUMQS,  SUMQback,ONE,SUMQSfalse
      REAL*8      SUMP1(MAXN), COEF1(MAXN),SUMP2(MAXN), COEF2(MAXN),
     $            SUMP3(MAXN), COEF3(MAXN),SUMP4(MAXN), COEF4(MAXN)
      COMPLEX*16  A(MAXN), B(MAXN), MSPHERE
      CHARACTER aerodist*1,type*5
 
 
C           Find the maximum number of terms required in the Mie series,
      call density_ice(a_mtox,bcoeff,rad2,dens_graup)
      rad2_ice=(dens_graup/0.917)**0.33333333*rad2
      MSPHERE =conjg(m_ice) 
      X = 2.0D0*PI*RAD2_ice/WAVELENGTH
      NTERMS = 0
      CALL MIECALC (NTERMS, X, MSPHERE, A, B)
      NLEGEN = 2*NTERMS
      NLEGEN = MIN0(MAXLEG, NLEGEN)
      NQUAD  = (NLEGEN + 2*NTERMS + 2)/2
      IF (NQUAD .GT. MAXN)  STOP 'MIE: MAXN exceeded' 
 
C           Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, MU, WTS)
 
      SUMQE = 0.0d0
      SUMQS = 0.0d0 
      SUMQSfalse=0.0d0
      SUMQback=0.0d0
      DO I = 1, NQUAD
        SUMP1(I) = 0.0d0
        SUMP2(I) = 0.0d0
        SUMP3(I) = 0.0d0
        SUMP4(I) = 0.0d0
      ENDDO
 
C               Integration loop over radius of spheres
      IF (NUMRAD .GT. 0)  DELRAD = (RAD2-RAD1)/NUMRAD
      DO IR = 1, NUMRAD+1
          RADIUS = RAD1 + (IR-1)*DELRAD
          NDENS = DISTRIBUTION (DBLE(AD), DBLE(BD), 
     .                     DBLE(ALPHA), DBLE(GAMMA), RADIUS, aerodist)
          IF ((IR .EQ. 1 .OR. IR .EQ. NUMRAD+1)
     .                      .AND. NUMRAD .GT. 0) THEN
              NDENS = 0.5*NDENS
          ENDIF
        
          NMIE = 0
         
          call density_ice(a_mtox,bcoeff,radius,dens_graup)
c         write(18,*)'dens',dens_graup
          radius_ice=(dens_graup/0.917)**0.33333333*radius

          X = 2.0D0*PI*RADIUS_ice/WAVELENGTH
C       these call is used to compute the back and phase function
          CALL MIECALC (NMIE, X, MSPHERE, A, B)
          CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT,Qback)
         SUMQSfalse = SUMQSfalse + QSCAT*NDENS*RADIUS_ice**2
c       scat and extinction are computed according to Parameterizations
          call snow_SS_param(WAVELENGTH, radius_ice,
     $ type,qscat,qext,asym)
          SUMQE = SUMQE + QEXT*NDENS*RADIUS_ice**2
          SUMQS = SUMQS + QSCAT*NDENS*RADIUS_ice**2
          SUMQback = SUMQback + Qback*NDENS*RADIUS_ice**2
c          write(*,*)'Mie in',Qext,qscat,Ndens,radius,x,nmie
          IF (PHASEFLAG) THEN
            NMIE = MIN0(NMIE, NTERMS)
            DO I = 1, NQUAD
              CALL MIEANGLE (NMIE, A, B, MU(I), P1,P2,P3,P4)
              SUMP1(I) = SUMP1(I) + P1*NDENS 
              SUMP2(I) = SUMP2(I) + P2*NDENS
              SUMP3(I) = SUMP3(I) + P3*NDENS
              SUMP4(I) = SUMP4(I) + P4*NDENS
            ENDDO
          ENDIF
        ENDDO
 
 
C           Multiply the sums by the integration delta and other constants
C             Put quadrature weights in angular array for later
      IF (NUMRAD .EQ. 0) DELRAD = 1.0d0
 
      EXTINCTION = PI*SUMQE*DELRAD
      SCATTER = PI*SUMQS*DELRAD
      Back_scatt= PI*SUMQback*DELRAD
      ALBEDO = SCATTER/EXTINCTION
      SCATTERFalse = PI*SUMQSFalse*DELRAD
C         If the phase function is not desired then leave now
      IF (.NOT. PHASEFLAG) RETURN
 
      TMP = (WAVELENGTH**2/(PI*SCATTERFalse)) * DELRAD
      DO I = 1, NQUAD
        SUMP1(I) = TMP*SUMP1(I) *WTS(I)
        SUMP2(I) = TMP*SUMP2(I) *WTS(I)
        SUMP3(I) = TMP*SUMP3(I) *WTS(I)
        SUMP4(I) = TMP*SUMP4(I) *WTS(I)
      ENDDO
 
C           Integrate the angular scattering functions times Legendre
C             polynomials to find the Legendre coefficients
      DO M = 1, NLEGEN+1
        COEF1(M) = 0.0d0
        COEF2(M) = 0.0d0
        COEF3(M) = 0.0d0
        COEF4(M) = 0.0d0
      ENDDO
C           Use upward recurrence to find Legendre polynomials
      DO I = 1, NQUAD
          PL1 = 1.0d0
          PL = 1.0d0
          DO L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
              COEF1(M) = COEF1(M) + SUMP1(I)*PL 
              COEF2(M) = COEF2(M) + SUMP2(I)*PL 
              COEF3(M) = COEF3(M) + SUMP3(I)*PL 
              COEF4(M) = COEF4(M) + SUMP4(I)*PL
              PL2 = PL1
              PL1 = PL
          ENDDO
      ENDDO
      NLEG = NLEGEN
      DO L = 0, NLEG
          M = L + 1
          LEGEN(M) = (2*L+1)/2.0 *COEF1(M) 
          LEGEN2(M)= (2*L+1)/2.0 *COEF2(M) 
          LEGEN3(M)= (2*L+1)/2.0 *COEF3(M) 
          LEGEN4(M)= (2*L+1)/2.0 *COEF4(M)
          IF (LEGEN(M) .GT. 1.0E-7)  NLEGEN = L
      ENDDO
       
      RETURN
      END




 
        real*8  FUNCTION DISTRIBUTION (A, B, ALPHA, GAMMA, R, distflag)
C        DISTRIBUTION returns the particle density for a given radius R
C      for a modified gamma distribution specified by A, B, ALPHA, GAMMA:
C           N(r) = a * r^alpha * exp(-b * r^gamma)     .
C      or a log-normal distribution:
C           N(r) = a/r * exp(- ln(r/b)^2 / (2*alpha^2) )     .
C      depending in distflag.
      IMPLICIT NONE
      REAL*8   A, B, ALPHA, GAMMA, R
      CHARACTER distflag*1
 
      IF (DISTFLAG .EQ. 'G') THEN
C           Modified gamma distibution
        DISTRIBUTION = A* R**ALPHA * DEXP(-B*R**GAMMA)
      ELSEIF (DISTFLAG .EQ. 'L') THEN
C           Log-normal distibution
        DISTRIBUTION = A/R *DEXP(-.5*(DLOG(R/B))**2/ALPHA**2)
      ELSE
        WRITE (*,*) 'Unrecognized distflag in DISTRIBUTION'
      ENDIF
 
      RETURN
      END



     
 

 
 
      SUBROUTINE MIECALC (NTERMS, X, MN, A, B)
C        MIECALC calculates the complex Mie coefficients An and Bn
C      given the dimensionless size parameter X and the complex
C      index of refraction (Mre,Mim).  The number of terms calculated
C      is given by NTERMS unless NTERMS <= 0 or in which case the
C      appropriate number is calculated and returned in NTERMS.
      IMPLICIT NONE
      INTEGER   NTERMS
      REAL*8    X
      COMPLEX*16  MN, A(*), B(*)
      INTEGER     MAXTERMS,  NSTOP, N, NN
      PARAMETER   (MAXTERMS=10000)
      REAL*8      PSIN, PSIM, CHIN, CHIM, TMP
      REAL*8      DCOS, DSIN
      COMPLEX*16  M, Y, D(MAXTERMS+15), XIN, XIM, CTMP
      COMPLEX*16  DCMPLX
 
 
C           If NTERMS is not specified calculate it
      NSTOP = X + 4.0*X**0.3334 + 2
      IF (NTERMS .LE. 0)  NTERMS = NSTOP
      IF (NTERMS .GT. MAXTERMS) THEN
          WRITE (*,*)
     .       'Mie calculation requires more terms than available.'
          STOP
      ENDIF
 
C           Generate the Dn's by down recurrence  D = d(log(PSI(y)))/dy
      M = DCONJG(MN)
      Y = M*X
      NN = NTERMS + 15
      D(NN) = DCMPLX (0.0D0, 0.0D0)
      DO N = NN, 2, -1
          D(N-1) = N/Y - 1.0/ (D(N) + N/Y)
      ENDDO
 
C           Generate the PSIn's and XIn'S by upward recurrence
C             and calculate the An's and Bn's from them.
C             (PSIN = PSI(n), PSIM = PSI(n-1), same for CHI)
      PSIM = DCOS(X)
      PSIN = DSIN(X)
      CHIM = -DSIN(X)
      CHIN = DCOS(X)
      DO N = 1, NTERMS
          TMP = PSIN
          PSIN = (2*N-1)/X *PSIN - PSIM
          PSIM = TMP
          TMP = CHIN
          CHIN = (2*N-1)/X *CHIN - CHIM
          CHIM = TMP
          XIN = DCMPLX (PSIN, -CHIN)
          XIM = DCMPLX (PSIM, -CHIM)
          CTMP = D(N)/M + N/X
          A(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
          CTMP = M*D(N) + N/X
          B(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
      ENDDO
 
      RETURN
      END
 
 
 

      SUBROUTINE MIECROSS (NTERMS, X, A, B, QEXT, QSCAT,QBackSCAT)
C        MIECROSS calculates the extinction, scattering, and
C      backscatter efficiencies given the Mie coefficients An and Bn
C      and the size parameter X.
      IMPLICIT NONE
      INTEGER    NTERMS
      REAL*8     X, QEXT, QSCAT,QBackSCAT
      COMPLEX*16 A(*), B(*), SUM3
      INTEGER    N
      REAL*8     SUM1, SUM2
 
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      SUM3 = 0.0d0
      DO N = 1, NTERMS
          SUM1 = SUM1 + (2*N+1)*( DREAL(A(N)) + DREAL(B(N)) )
          SUM2 = SUM2 + (2*N+1)*( DREAL(A(N)*DCONJG(A(N)))
     .                          + DREAL(B(N)*DCONJG(B(N))) )
         SUM3=SUM3+(2*N+1)*(-1.0)**N*(A(N)-B(N))
      ENDDO
      QEXT = 2.0D0/X**2 * SUM1
      QSCAT = 2.0D0/X**2 * SUM2
      QBackSCAT = 1.0D0/X**2 * SUM3*DCONJG(SUM3)
      RETURN
      END
 
 
 
 
      SUBROUTINE MIEANGLE (NTERMS, A, B, MU, P1,P2,P3,P4)
C        MIEANGLE calculates the intensity scattering matrix elements
C      (P1,P2,P3,P4) for a particular value of MU (cos(theta)) from the
C      Mie coefficients An's and Bn's.  The matrix elements are for the
C      stokes intensity vector (I,Q,U,V) and are calculated from the
C      complex scattering amplitudes S1 and S2.
      IMPLICIT NONE
      INTEGER    NTERMS
      REAL*8     MU, P1,P2,P3,P4
      COMPLEX*16 A(*), B(*)
      INTEGER    N
      REAL*8     TMP, PIN, PIM, TAUN, C
      COMPLEX*16 S1, S2
 
 
      S1 = DCMPLX(0.0,0.0)
      S2 = DCMPLX(0.0,0.0)
C               Sum up the series using the An's and Bn's
      PIN = 1.0d0
      PIM = 0.0d0
      DO N = 1, NTERMS
          TAUN = N*MU*PIN - (N+1)*PIM
C               Calculate the scattering functions at +mu and -mu
C                 using the PIn's and the TAUn's.
          C = (2*N+1) / DFLOAT(N*(N+1))
          S1 = S1 + C*( A(N)*PIN + B(N)*TAUN)
          S2 = S2 + C*( B(N)*PIN + A(N)*TAUN)
C               Calculate the angular function PIn by up recurrence
          TMP = PIN
          PIN = ( (2*N+1)*MU*PIN - (N+1)*PIM ) / N
          PIM = TMP
      ENDDO
C           Calculate the first Stokes parameter scattering matrix element
      P1 = 0.5*( CDABS(S2)**2 + CDABS(S1)**2 )
      P2=  0.5*( CDABS(S2)**2 - CDABS(S1)**2 )
      P3=  dreal(DCONJG(S2)*S1)
      P4=  imag(DCONJG(S1)*S2)
      
      RETURN
      END
 
        SUBROUTINE LEGEndre2PHASEFUNCTION(A1,LMAX,NPNA,NPL,P11,ANG) 

      IMPLICIT NONE
      INTEGER J,NPNA,NPL,LMAX,L1,L1MAX,I1,L,K
      REAL*8   PL(NPL),A1(NPL),ANG(NPNA),P11(NPNA),DN,DA,DB
          REAL*8 D6,TAA,TB,U,F11,DL,DL1
    
     
       
    
C *** 
        DN=1D0/DFLOAT(NPNA-1)
        DA=DACOS(-1D0)*DN
        DB=180D0*DN
        L1MAX=LMAX+1

C *** COMPUTATION OF .F11
        TB=-DB
        TAA=-DA
        D6=DSQRT(6D0)*0.25D0

C *** LOOP ON THE NUMBER OF ANGLES
  
        DO 500 I1=1,NPNA
          TAA=TAA+DA
          TB=TB+DB
          U=COS(TAA)   
C ==> FOR LEGENDRE 
            CALL FLEG2(U,PL,L1MAX)
            F11=0.D0
          
            DO 403 L1=1,L1MAX            
              F11=F11+A1(L1)*PL(L1)
 403        CONTINUE
        
     
              P11(I1)=F11
              ANG(I1)=TB*DACOS(-1.D0)/180.0
  500       CONTINUE
  
      RETURN
      END
C    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&



      SUBROUTINE GAUSQUAD (N, XA, WT)
C        Generates the abscissas (X) and weights (W) for an N point
C      Gauss-Legendre quadrature.  
      IMPLICIT NONE
      INTEGER  N
      REAL*8   XA(*), WT(*)
      INTEGER  K, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-13)

      K = (N+1)/2
      DO 130 J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO 120 L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
120       CONTINUE
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J)     = -X
        XA(N-J+1) = X
        WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
        WT(N-J+1) = WT(J)
130   CONTINUE

      RETURN
      END




      DOUBLE PRECISION FUNCTION gammln(xx)
      DOUBLE PRECISION xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +>k-5V1`..

  








      SUBROUTINE REFWAT(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF)
C
C     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR WATER
C     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM .2 MICRONS TO 10 CM
C     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 0.1 CM
C
C     ERIC A. SMITH
C     DEPT OF ATMOSPHERIC SCIENCE
C     COLORADO STATE UNIVERSITY
C     FORT COLLINS,CO  80523
C     TEL   303-491-8533
C
C     REFERENCES
C
C     0.2 UM - 0.69 UM
C
C     HALE,G., AND M. QUERRY,1972.
C     OPTICAL CONSTANTS OF WATER IN THE 200 NM TO 200 UM WAVELENGTH REGI
C     APPLIED OPTICS,12,3,555-563.
C
C     0.69 UM - 2.0 UM
C
C     PALMER,K.F., AND D. WILLIAMS,1974.
C     OPTICAL PROPERTIES OF WATER IN THE NEAR INFRARED.
C     JOURNAL OF THE OPTICAL SOCIETY OF AMERICA,64,8,1107-1110.
C
C     2.0 UM - 1000.0 UM
C
C     DOWNING,H.D., AND D. WILLIAMS,1975.
C     OPTICAL CONSTANTS OF WATER IN THE INFRARED.
C     JOURNAL OF GEOPHYSICAL REVIEW,80,12,1656-1661.
C
C     1.0 MM - 10.0 CM
C
C     RAY,P.S.,1972.
C     BROADBAND COMPLEX REFRACTIVE INDICES OF ICE AND WATER.
C     APPLIED OPTICS,11,8,1836-1844.
C
C     INPUT PARAMETERS
C
C     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
C           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
C           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
C           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
C     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
C     T = TEMPERATURE ( DEGREES KELVIN )
C
C     OUTPUT PARAMETERS
C
C     RN = REAL PORTION ( SCATTERING )
C     CN = COMPLEX PORTION ( ABSORPTION )
C     ABSIND = ABSORPTIVE INDEX ( CN/RN )
C     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX E,M
      DIMENSION WLTABW(518),RNTABW(518),CNTABW(518)
      DATA NUMWAT/518/
      DATA WLMIN,WLMAX/0.2,100000.0/
      DATA CUTWAT/1000.0/
      DATA (WLTABW(I),I=  1, 66)/
     *    .20000,    .22500,    .25000,    .27500,    .30000,    .32500,
     *    .35001,    .37500,    .40000,    .42501,    .45000,    .47499,
     *    .50000,    .52499,    .54999,    .57501,    .59999,    .62500,
     *    .64998,    .67499,    .68966,    .70175,    .71429,    .72464,
     *    .73529,    .74627,    .75188,    .75758,    .76923,    .78125,
     *    .79365,    .80645,    .81301,    .81967,    .83333,    .84746,
     *    .86207,    .87719,    .89286,    .90909,    .92593,    .93458,
     *    .94340,    .95238,    .96154,    .97276,    .98039,    .99010,
     *   1.00000,   1.01010,   1.02041,   1.03093,   1.04167,   1.05263,
     *   1.06952,   1.08696,   1.09890,   1.11111,   1.12360,   1.13636,
     *   1.14943,   1.16279,   1.17647,   1.19048,   1.20482,   1.21951/
      DATA (WLTABW(I),I= 67,132)/
     *   1.23457,   1.25000,   1.26582,   1.28205,   1.29870,   1.31579,
     *   1.33333,   1.35135,   1.36986,   1.38889,   1.40845,   1.42857,
     *   1.44300,   1.47059,   1.49254,   1.51515,   1.53846,   1.56250,
     *   1.58730,   1.61290,   1.63934,   1.66667,   1.69492,   1.72414,
     *   1.75439,   1.78571,   1.80180,   1.81818,   1.85185,   1.88679,
     *   1.92678,   1.96078,   2.00000,   2.02020,   2.04082,   2.06186,
     *   2.08333,   2.10526,   2.12766,   2.15054,   2.17391,   2.19780,
     *   2.22222,   2.24719,   2.27273,   2.29885,   2.32558,   2.35294,
     *   2.38095,   2.40964,   2.43902,   2.46914,   2.50000,   2.50627,
     *   2.51256,   2.51889,   2.52525,   2.53165,   2.53807,   2.54453,
     *   2.55102,   2.55754,   2.56410,   2.57069,   2.57732,   2.58398/
      DATA (WLTABW(I),I=133,198)/
     *   2.59067,   2.59740,   2.60417,   2.61097,   2.61780,   2.62467,
     *   2.63158,   2.63852,   2.64550,   2.65252,   2.65957,   2.66667,
     *   2.67380,   2.68097,   2.68817,   2.69542,   2.70270,   2.71003,
     *   2.71739,   2.72480,   2.73224,   2.73973,   2.74725,   2.75482,
     *   2.76243,   2.77008,   2.77778,   2.78552,   2.79330,   2.80112,
     *   2.80899,   2.81690,   2.82486,   2.83286,   2.84091,   2.84900,
     *   2.85714,   2.86533,   2.87356,   2.88184,   2.89017,   2.89855,
     *   2.90698,   2.91545,   2.92398,   2.93255,   2.94118,   2.94985,
     *   2.95858,   2.96736,   2.97619,   2.98507,   2.99401,   3.00300,
     *   3.01205,   3.02115,   3.03030,   3.03951,   3.04878,   3.05810,
     *   3.06748,   3.07692,   3.08642,   3.09598,   3.10559,   3.11526/
      DATA (WLTABW(I),I=199,264)/
     *   3.12500,   3.13480,   3.14465,   3.15457,   3.16456,   3.17460,
     *   3.18471,   3.19489,   3.20513,   3.21543,   3.22581,   3.23625,
     *   3.24675,   3.25733,   3.26797,   3.27869,   3.28947,   3.30033,
     *   3.31126,   3.32226,   3.33333,   3.34448,   3.35570,   3.36700,
     *   3.37838,   3.38983,   3.40136,   3.41297,   3.42466,   3.43643,
     *   3.44828,   3.46021,   3.47222,   3.48432,   3.49650,   3.50877,
     *   3.52113,   3.53357,   3.54610,   3.55872,   3.57143,   3.58423,
     *   3.59712,   3.61011,   3.62319,   3.63636,   3.64964,   3.66300,
     *   3.67647,   3.69004,   3.70370,   3.71747,   3.73134,   3.74532,
     *   3.75940,   3.77358,   3.78788,   3.80228,   3.81679,   3.83142,
     *   3.84615,   3.86100,   3.87597,   3.89105,   3.90625,   3.92157/
      DATA (WLTABW(I),I=265,330)/
     *   3.93701,   3.95257,   3.96825,   3.98406,   4.00000,   4.01606,
     *   4.03226,   4.04858,   4.06504,   4.08163,   4.09836,   4.11523,
     *   4.13223,   4.14938,   4.16667,   4.18410,   4.20168,   4.21941,
     *   4.23729,   4.25532,   4.27350,   4.29185,   4.31034,   4.32900,
     *   4.34783,   4.36681,   4.38596,   4.40529,   4.42478,   4.44444,
     *   4.46429,   4.48430,   4.50450,   4.52489,   4.54545,   4.56621,
     *   4.58716,   4.60829,   4.62963,   4.65116,   4.67290,   4.69484,
     *   4.71698,   4.73934,   4.76190,   4.78469,   4.80769,   4.83092,
     *   4.85437,   4.87805,   4.90196,   4.92611,   4.95050,   4.97512,
     *   5.00000,   5.02513,   5.05051,   5.07614,   5.10204,   5.12821,
     *   5.15464,   5.18135,   5.20833,   5.23560,   5.26316,   5.29101/
      DATA (WLTABW(I),I=331,396)/
     *   5.31915,   5.34759,   5.37634,   5.40541,   5.43478,   5.46448,
     *   5.49451,   5.52486,   5.55556,   5.58659,   5.61798,   5.64972,
     *   5.68182,   5.71429,   5.74713,   5.78035,   5.81395,   5.84795,
     *   5.88235,   5.91716,   5.95238,   5.98802,   6.02410,   6.06061,
     *   6.09756,   6.13497,   6.17284,   6.21118,   6.25000,   6.28931,
     *   6.32911,   6.36943,   6.41026,   6.45161,   6.49351,   6.53595,
     *   6.57895,   6.62252,   6.66667,   6.71141,   6.75676,   6.80272,
     *   6.84932,   6.89655,   6.94444,   6.99301,   7.04225,   7.09220,
     *   7.14286,   7.19424,   7.24638,   7.29927,   7.35294,   7.40741,
     *   7.46269,   7.51880,   7.57576,   7.63359,   7.69231,   7.75194,
     *   7.81250,   7.87402,   7.93651,   8.00000,   8.06452,   8.13008/
      DATA (WLTABW(I),I=397,462)/
     *   8.19672,   8.26446,   8.33333,   8.40336,   8.47458,   8.54701,
     *   8.62069,   8.69565,   8.77193,   8.84956,   8.92857,   9.00901,
     *   9.09091,   9.17431,   9.25926,   9.34579,   9.43396,   9.52381,
     *   9.61538,   9.70874,   9.80392,   9.90099,  10.00000,  10.10101,
     *  10.20408,  10.30928,  10.41667,  10.52632,  10.63830,  10.75269,
     *  10.86957,  10.98901,  11.11111,  11.23596,  11.36364,  11.49425,
     *  11.62791,  11.76471,  11.90476,  12.04819,  12.19512,  12.34568,
     *  12.50000,  12.65823,  12.82051,  12.98701,  13.15789,  13.33333,
     *  13.51351,  13.69863,  13.88889,  14.08451,  14.28571,  14.49275,
     *  14.70588,  14.92537,  15.15152,  15.38462,  15.62500,  15.87302,
     *  16.12903,  16.39344,  16.66667,  16.94915,  17.24138,  17.54386/
      DATA (WLTABW(I),I=463,518)/
     *  17.85714,  18.18182,  18.51852,  18.86792,  19.23077,  19.60784,
     *  20.00000,  20.40816,  20.83333,  21.27660,  21.73913,  22.22222,
     *  22.72727,  23.25581,  23.80952,  24.39024,  25.00000,  25.64103,
     *  26.31579,  27.02703,  27.77778,  28.57143,  29.41176,  30.30303,
     *  31.25000,  32.25806,  33.33333,  34.48276,  35.71429,  37.03704,
     *  38.46154,  40.00000,  41.66667,  43.47826,  45.45455,  47.61905,
     *  50.00000,  52.63158,  55.55556,  58.82353,  62.50000,  66.66667,
     *  71.42857,  76.92308,  83.33333,  90.90909, 100.00000, 111.11111,
     * 125.00000, 142.85714, 166.66667, 200.00000, 250.00000, 333.33333,
     * 500.00000,1000.00000/
      DATA (RNTABW(I),I=  1, 66)/
     *1.396,1.373,1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,
     *1.336,1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.332,1.332,
     *1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331,
     *1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.328,
     *1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,
     *1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.326,1.325,1.325,1.325/
      DATA (RNTABW(I),I= 67,132)/
     *1.325,1.325,1.324,1.324,1.324,1.324,1.323,1.323,1.323,1.322,1.322,
     *1.321,1.321,1.321,1.320,1.320,1.319,1.319,1.318,1.318,1.317,1.316,
     *1.315,1.314,1.314,1.313,1.312,1.312,1.311,1.310,1.309,1.307,1.306,
     *1.301,1.301,1.300,1.298,1.298,1.296,1.295,1.294,1.293,1.291,1.289,
     *1.287,1.285,1.282,1.280,1.277,1.274,1.270,1.265,1.261,1.260,1.259,
     *1.257,1.256,1.255,1.254,1.252,1.250,1.249,1.247,1.246,1.243,1.241/
      DATA (RNTABW(I),I=133,198)/
     *1.240,1.238,1.235,1.232,1.230,1.227,1.224,1.221,1.218,1.214,1.210,
     *1.205,1.200,1.195,1.191,1.185,1.179,1.172,1.166,1.157,1.149,1.144,
     *1.139,1.138,1.138,1.139,1.141,1.144,1.149,1.154,1.158,1.161,1.165,
     *1.171,1.177,1.183,1.191,1.199,1.212,1.220,1.233,1.246,1.258,1.271,
     *1.282,1.293,1.305,1.317,1.329,1.342,1.353,1.364,1.376,1.386,1.398,
     *1.407,1.417,1.426,1.434,1.442,1.450,1.457,1.465,1.471,1.476,1.480/
      DATA (RNTABW(I),I=199,264)/
     *1.483,1.486,1.487,1.487,1.487,1.486,1.485,1.482,1.479,1.477,1.474,
     *1.472,1.467,1.464,1.461,1.457,1.454,1.451,1.448,1.444,1.441,1.437,
     *1.434,1.431,1.427,1.425,1.421,1.418,1.415,1.413,1.410,1.407,1.405,
     *1.403,1.400,1.398,1.396,1.394,1.392,1.390,1.388,1.387,1.385,1.383,
     *1.382,1.379,1.378,1.377,1.375,1.374,1.372,1.371,1.370,1.369,1.367,
     *1.366,1.365,1.363,1.361,1.361,1.360,1.358,1.358,1.357,1.355,1.354/
      DATA (RNTABW(I),I=265,330)/
     *1.353,1.352,1.351,1.350,1.349,1.348,1.348,1.347,1.346,1.345,1.344,
     *1.344,1.343,1.342,1.341,1.340,1.340,1.338,1.337,1.337,1.335,1.334,
     *1.334,1.333,1.332,1.332,1.331,1.330,1.330,1.330,1.329,1.329,1.329,
     *1.328,1.328,1.327,1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.325,
     *1.325,1.325,1.325,1.325,1.325,1.324,1.324,1.323,1.322,1.322,1.321,
     *1.320,1.319,1.318,1.318,1.317,1.316,1.314,1.313,1.311,1.310,1.308/
      DATA (RNTABW(I),I=331,396)/
     *1.306,1.304,1.302,1.299,1.297,1.294,1.291,1.288,1.285,1.282,1.278,
     *1.275,1.271,1.267,1.262,1.256,1.251,1.247,1.242,1.241,1.241,1.247,
     *1.265,1.289,1.311,1.332,1.349,1.354,1.356,1.354,1.350,1.345,1.341,
     *1.337,1.333,1.330,1.326,1.324,1.322,1.320,1.319,1.318,1.317,1.316,
     *1.315,1.314,1.313,1.311,1.310,1.309,1.308,1.307,1.306,1.305,1.303,
     *1.302,1.301,1.300,1.298,1.296,1.295,1.294,1.293,1.291,1.288,1.286/
      DATA (RNTABW(I),I=397,462)/
     *1.285,1.283,1.281,1.279,1.276,1.274,1.271,1.269,1.267,1.264,1.261,
     *1.259,1.256,1.253,1.249,1.246,1.242,1.238,1.234,1.230,1.224,1.220,
     *1.214,1.208,1.202,1.194,1.189,1.181,1.174,1.168,1.162,1.156,1.149,
     *1.143,1.139,1.135,1.132,1.132,1.131,1.132,1.130,1.130,1.134,1.138,
     *1.142,1.157,1.171,1.182,1.189,1.201,1.213,1.223,1.236,1.249,1.264,
     *1.277,1.289,1.303,1.313,1.324,1.335,1.348,1.361,1.372,1.385,1.396/
      DATA (RNTABW(I),I=463,518)/
     *1.407,1.419,1.431,1.441,1.451,1.462,1.470,1.480,1.488,1.496,1.504,
     *1.510,1.515,1.521,1.527,1.532,1.537,1.541,1.545,1.549,1.552,1.552,
     *1.552,1.550,1.546,1.543,1.541,1.539,1.537,1.534,1.532,1.529,1.525,
     *1.528,1.542,1.567,1.600,1.640,1.689,1.746,1.801,1.848,1.890,1.929,
     *1.960,1.982,1.997,2.000,2.010,2.020,2.040,2.070,2.110,2.150,2.225,
     *2.481/
      DATA (CNTABW(I),I=  1, 66)/
     *1.1000E-07,4.9000E-08,3.4000E-08,2.4000E-08,1.6000E-08,1.1000E-08,
     *6.5000E-09,3.5000E-09,1.9000E-09,1.3000E-09,1.0000E-09,9.4000E-10,
     *1.0000E-09,1.3000E-09,2.0000E-09,3.6000E-09,1.1000E-08,1.4000E-08,
     *1.6000E-08,2.2000E-08,2.7000E-08,3.8000E-08,5.6000E-08,7.7300E-08,
     *1.3900E-07,1.6300E-07,1.6800E-07,1.6400E-07,1.5400E-07,1.4300E-07,
     *1.3300E-07,1.2500E-07,1.2400E-07,1.3000E-07,2.0400E-07,2.6100E-07,
     *2.9400E-07,3.5300E-07,4.3300E-07,5.4300E-07,8.7700E-07,1.1800E-06,
     *1.6100E-06,2.4400E-06,3.6000E-06,3.9800E-06,3.9200E-06,3.7000E-06,
     *3.3100E-06,2.8200E-06,2.3100E-06,1.9000E-06,1.5700E-06,1.3700E-06,
     *1.2600E-06,1.4400E-06,1.6800E-06,2.0500E-06,2.8900E-06,4.9600E-06,
     *8.8700E-06,1.0900E-05,1.1500E-05,1.1800E-05,1.2000E-05,1.1800E-05/
      DATA (CNTABW(I),I= 67,132)/
     *1.1500E-05,1.1000E-05,1.0800E-05,1.1500E-05,1.3800E-05,1.7500E-05,
     *2.3900E-05,4.1600E-05,5.9400E-05,1.0100E-04,2.4100E-04,3.5200E-04,
     *3.6400E-04,3.3400E-04,2.5800E-04,1.8800E-04,1.4800E-04,1.2000E-04,
     *1.0200E-04,8.7300E-05,7.9200E-05,7.4900E-05,7.6200E-05,8.5500E-05,
     *1.0600E-04,1.3000E-04,1.3600E-04,1.3700E-04,1.5900E-04,8.6300E-04,
     *1.9000E-03,1.7000E-03,1.1000E-03,9.0000E-04,7.3100E-04,6.1700E-04,
     *5.1400E-04,4.5200E-04,4.0000E-04,3.5900E-04,3.4100E-04,3.3800E-04,
     *3.4500E-04,3.7600E-04,4.1600E-04,4.6500E-04,5.4200E-04,6.5200E-04,
     *7.9200E-04,9.6800E-04,1.2300E-03,1.5600E-03,1.9000E-03,1.9500E-03,
     *2.0000E-03,2.0500E-03,2.0700E-03,2.1000E-03,2.1200E-03,2.1500E-03,
     *2.1900E-03,2.2400E-03,2.2700E-03,2.3100E-03,2.3400E-03,2.3900E-03/
      DATA (CNTABW(I),I=133,198)/
     *2.4300E-03,2.4800E-03,2.5700E-03,2.7000E-03,2.9800E-03,3.3000E-03,
     *4.0200E-03,4.3700E-03,4.8200E-03,5.3600E-03,6.2700E-03,7.3200E-03,
     *8.5500E-03,1.0500E-02,1.2700E-02,1.4500E-02,1.6400E-02,1.8600E-02,
     *2.0500E-02,2.8200E-02,3.8000E-02,4.6200E-02,5.4800E-02,6.4900E-02,
     *7.4400E-02,8.3600E-02,9.2700E-02,1.0200E-01,1.1200E-01,1.2100E-01,
     *1.3100E-01,1.4200E-01,1.5400E-01,1.6700E-01,1.8000E-01,1.9400E-01,
     *2.0600E-01,2.1800E-01,2.2900E-01,2.3900E-01,2.4900E-01,2.5800E-01,
     *2.6500E-01,2.7100E-01,2.7600E-01,2.8000E-01,2.8100E-01,2.8200E-01,
     *2.8200E-01,2.7900E-01,2.7600E-01,2.7200E-01,2.6700E-01,2.6200E-01,
     *2.5500E-01,2.5000E-01,2.4300E-01,2.3600E-01,2.2800E-01,2.2000E-01,
     *2.1200E-01,2.0400E-01,1.9500E-01,1.8300E-01,1.7300E-01,1.6300E-01/
      DATA (CNTABW(I),I=199,264)/
     *1.5300E-01,1.4400E-01,1.3400E-01,1.2500E-01,1.1700E-01,1.1000E-01,
     *9.9400E-02,9.2000E-02,8.5500E-02,7.8500E-02,7.1600E-02,6.5300E-02,
     *6.0000E-02,5.5000E-02,5.0400E-02,4.6200E-02,4.2200E-02,3.8500E-02,
     *3.4800E-02,3.1500E-02,2.9700E-02,2.7900E-02,2.6200E-02,2.5000E-02,
     *2.2900E-02,2.1000E-02,1.9300E-02,1.7700E-02,1.6300E-02,1.5100E-02,
     *1.3800E-02,1.2800E-02,1.1800E-02,1.1000E-02,1.0100E-02,9.4100E-03,
     *8.6600E-03,8.0700E-03,7.3700E-03,6.8300E-03,6.2500E-03,5.7900E-03,
     *5.3800E-03,5.0600E-03,4.7300E-03,4.4900E-03,4.2400E-03,4.0500E-03,
     *3.8900E-03,3.7600E-03,3.6300E-03,3.5500E-03,3.4700E-03,3.4000E-03,
     *3.3500E-03,3.3600E-03,3.3500E-03,3.3900E-03,3.4000E-03,3.4800E-03,
     *3.5200E-03,3.6300E-03,3.7000E-03,3.7800E-03,3.8900E-03,3.9900E-03/
      DATA (CNTABW(I),I=265,330)/
     *4.1000E-03,4.2200E-03,4.3300E-03,4.5000E-03,4.6500E-03,4.7900E-03,
     *4.9400E-03,5.1200E-03,5.3100E-03,5.4900E-03,5.6800E-03,5.8600E-03,
     *6.0800E-03,6.3100E-03,6.5300E-03,6.7300E-03,6.9600E-03,7.2200E-03,
     *7.4900E-03,7.7900E-03,8.0600E-03,8.3300E-03,8.6400E-03,8.9600E-03,
     *9.2700E-03,9.6600E-03,1.0000E-02,1.0400E-02,1.0800E-02,1.1200E-02,
     *1.1700E-02,1.2200E-02,1.2600E-02,1.3100E-02,1.3600E-02,1.4000E-02,
     *1.4500E-02,1.4900E-02,1.5200E-02,1.5400E-02,1.5600E-02,1.5700E-02,
     *1.5700E-02,1.5700E-02,1.5500E-02,1.5300E-02,1.5100E-02,1.4800E-02,
     *1.4600E-02,1.4300E-02,1.4000E-02,1.3700E-02,1.3300E-02,1.2900E-02,
     *1.2600E-02,1.2200E-02,1.1800E-02,1.1500E-02,1.1000E-02,1.0800E-02,
     *1.0500E-02,1.0300E-02,1.0100E-02,1.0000E-02,9.9300E-03,9.9000E-03/
      DATA (CNTABW(I),I=331,396)/
     *9.9500E-03,1.0000E-02,1.0200E-02,1.0400E-02,1.0700E-02,1.1000E-02,
     *1.1500E-02,1.2000E-02,1.2800E-02,1.3800E-02,1.5000E-02,1.6600E-02,
     *1.8500E-02,2.0500E-02,2.4200E-02,2.9300E-02,3.3200E-02,4.2900E-02,
     *5.4400E-02,6.8800E-02,8.4000E-02,1.0210E-01,1.1700E-01,1.3000E-01,
     *1.3200E-01,1.2400E-01,1.0600E-01,8.8000E-02,7.4000E-02,6.1800E-02,
     *5.3500E-02,4.8400E-02,4.4700E-02,4.2000E-02,3.9800E-02,3.8300E-02,
     *3.7300E-02,3.7000E-02,3.6600E-02,3.6300E-02,3.6000E-02,3.5700E-02,
     *3.5500E-02,3.5200E-02,3.5000E-02,3.4700E-02,3.4600E-02,3.4300E-02,
     *3.4200E-02,3.4200E-02,3.4200E-02,3.4300E-02,3.4200E-02,3.4200E-02,
     *3.4200E-02,3.4200E-02,3.4200E-02,3.4400E-02,3.4500E-02,3.4600E-02,
     *3.4900E-02,3.5100E-02,3.5100E-02,3.5100E-02,3.5200E-02,3.5600E-02/
      DATA (CNTABW(I),I=397,462)/
     *3.5900E-02,3.6100E-02,3.6200E-02,3.6600E-02,3.7000E-02,3.7400E-02,
     *3.7800E-02,3.8300E-02,3.8700E-02,3.9200E-02,3.9800E-02,4.0500E-02,
     *4.1100E-02,4.1700E-02,4.2400E-02,4.3400E-02,4.4300E-02,4.5300E-02,
     *4.6700E-02,4.8100E-02,4.9700E-02,5.1500E-02,5.3400E-02,5.5700E-02,
     *5.8900E-02,6.2200E-02,6.6100E-02,7.0700E-02,7.6400E-02,8.2800E-02,
     *8.9800E-02,9.7300E-02,1.0700E-01,1.1800E-01,1.3000E-01,1.4400E-01,
     *1.5900E-01,1.7600E-01,1.9200E-01,2.0800E-01,2.2600E-01,2.4300E-01,
     *2.6000E-01,2.7700E-01,2.9200E-01,3.0500E-01,3.1700E-01,3.2800E-01,
     *3.3800E-01,3.4700E-01,3.5600E-01,3.6500E-01,3.7300E-01,3.7900E-01,
     *3.8600E-01,3.9200E-01,3.9700E-01,4.0300E-01,4.0800E-01,4.1200E-01,
     *4.1700E-01,4.2000E-01,4.2300E-01,4.2500E-01,4.2700E-01,4.2800E-01/
      DATA (CNTABW(I),I=463,518)/
     *4.2700E-01,4.2700E-01,4.2600E-01,4.2500E-01,4.2300E-01,4.2100E-01,
     *4.1800E-01,4.1500E-01,4.1100E-01,4.0800E-01,4.0400E-01,4.0100E-01,
     *3.9700E-01,3.9400E-01,3.9000E-01,3.8600E-01,3.8200E-01,3.7700E-01,
     *3.7200E-01,3.6800E-01,3.6300E-01,3.5900E-01,3.5600E-01,3.5200E-01,
     *3.5300E-01,3.5700E-01,3.6100E-01,3.6800E-01,3.7500E-01,3.8500E-01,
     *3.9800E-01,4.1400E-01,4.3600E-01,4.6900E-01,5.0500E-01,5.3900E-01,
     *5.7100E-01,5.9700E-01,6.1800E-01,6.2900E-01,6.2200E-01,6.0800E-01,
     *5.9300E-01,5.7700E-01,5.5700E-01,5.3200E-01,5.0700E-01,4.8700E-01,
     *4.6600E-01,4.5000E-01,4.4400E-01,4.3800E-01,4.6000E-01,5.2700E-01,
     *7.1800E-01,8.4657E-01/
      DATA PI/3.14159265/
C
C     FUNCTION FOR TREATING ABSORPTION BANDS NOT CONSIDERED IN THE
C     DEBYE THEOREY
C
      SUM(WL,WLCEN,BET,DEL,GAM)=BET*
     $ DEXP(-ABS(DLOG10(WL/WLCEN)/DEL)**GAM)
C
C     ZERO PARAMETERS
C
      RN=0.0
      CN=0.0
      ABSIND=0.0
      ABSCOF=0.0
C
C     CONVERT WAVELENGTH TO MICRONS
C
      WL=XLAM
      IF(IUNIT.EQ.1)WL=1000*WL
      IF(IUNIT.EQ.2)WL=10000*WL
      IF(IUNIT.EQ.3)WL=10000*(1.0/WL)
      IF(WL.LT.WLMIN.OR.WL.GT.WLMAX)RETURN
C
C     REGION FROM 0.2 MICRON TO 1000.0 MICRON  -  TABLE LOOKUP
C
      IF(WL.GT.CUTWAT)GO TO 3
      DO 1 I=2,NUMWAT
      IF(WL.GT.WLTABW(I))GO TO 1
      I1=I-1
      I2=I
      GO TO 2
 1    CONTINUE
      I1=NUMWAT-1
      I2=NUMWAT
 2    FAC=(WL-WLTABW(I1))/(WLTABW(I2)-WLTABW(I1))
      RN=RNTABW(I1)+FAC*(RNTABW(I2)-RNTABW(I1))
      CN=CNTABW(I1)+FAC*(CNTABW(I2)-CNTABW(I1))
      GO TO 5
C
C     REGION FROM 0.1 CM TO 10 CM
C
C     EXTENSION OF DEBYE THEOREY BASED ON THE WORK OF
C
C        COLE,K.S.,AND R.H.COLE,1941.JOUR.CHEM.PHYS.,9,P 341.
C
C     DEFINE TEMPERATURE TERMS AND WAVELENGTH IN CM
C
 3    TC=T-273.15
      T1=TC+273.0
      T2=TC-25.0
      XL=WL/10000.0
C
C     DEFINE FREQUENCY INDEPENDENT CONDUCTIVITY(SIGMA) AND
C     SPREAD PARAMETER(ALPHA)
C
C     IN CLASSICAL DEBYE THEOREY THESE TERMS ARE ZERO
C
C     SIGMA GIVEN BY SAXTON,J.A.,1949.WIRELESS ENGINEER,26,P 288.
C     ALPHA GIVEN BY RAY ( EQUATION 7B )
C
      SIGMA=12.5664E8
      ALPHA=-16.8129/T1+0.0609265
C
C     DEFINE STATIC DIELECTRIC CONSTANT(ES) - RAY EQN 4
C            HIGH FREQUENCY DIELECTRIC CONSTANT(E00) - RAY EQN 7A
C            RELAXTION WAVELENGTH IN CM(XLAMS) - RAY EQN 7C
C
C     TEMPERATURE DEPENDENCE OF ES GIVEN BY
C
C        WYMAN,J.,AND E.N.INGALLS,1938.JOUR.AM.CHEM.SOC.,60,P 1182.
C
      ES=78.54*(1.0-4.579E-3*T2+1.19E-5*T2*T2-2.8E-8*T2*T2*T2)
      E00=5.27137+0.0216474*TC-0.00131198*TC*TC
      XLAMS=0.00033836*EXP(2513.98/T1)
C
C     CALCULATE EXPRESSIONS USED FOR DIELECTRIC CONSTANT
C
      TERM=PI*ALPHA/2
      SINT=SIN(TERM)
      COST=COS(TERM)
      XLRAT=XLAMS/XL
      POWTRM=XLRAT**(1-ALPHA)
      DENOM=1.0+2*POWTRM*SINT+XLRAT**(2.0*(1-ALPHA))
C
C     CALCULATION OF DIELECTRIC CONSTANT
C
C     REAL PART - RAY EQN 5
C
      ER=E00+(ES-E00)*(1.0+POWTRM*SINT)/DENOM
C
C     IMAGINARY PART OR LOSS TERM - RAY EQN 6
C
      EI=(SIGMA*XL/18.8496E10)+(ES-E00)*POWTRM*COST/DENOM
C
C     COMPLEX PERMITTIVITY
C
      E=CMPLX(ER,-EI)
C
C     COMPLEX INDEX OF REFRACTION - RAY EQN 1
C
      M=CSQRT(E)
      RN=REAL(M)
      CN=-AIMAG(M)
C
C     CORRECTION TO IMAGINARY INDEX TO ACCOUNT FOR THE
C     REMAINING ABSORPTION BANDS - RAY EQN 8(TABLE 2)
C
      IF(WL.GT.3000.0)GO TO 5
      CN=CN+SUM_function(WL, 17.0,0.39,0.45,1.3)
     *     +SUM_function(WL, 62.0,0.41,0.35,1.7)
     *     +SUM_function(WL,300.0,0.25,0.47,3.0)
C
C     ABSORPTIVE QUANITIES
C
 5    ABSIND=CN/RN
      ABSCOF=4.0*PI*CN/WL
      RETURN
      END

      real*8 FUNCTION SUM_function(WL,WLCEN,BET,DEL,GAM)
      implicit none
      real*8 WL,WLCEN,BET,DEL,GAM !,sum_function  : PSG

      SUM_function=BET*DEXP(-ABS(DLOG10(WL/WLCEN)/DEL)**GAM)
      RETURN
      END



      SUBROUTINE REFICE(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF)
C
C     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
C     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
C     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.
C
C     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
C                               RN  VS.        T
C                           LOG(CN) VS. LOG(XLAM)
C                           LOG(CN) VS.        T
C
C     STEPHEN G. WARREN - 1983
C     DEPT. OF ATMOSPHERIC SCIENCES
C     UNIVERSITY OF WASHINGTON
C     SEATTLE, WA  98195
C
C     BASED ON
C
C        WARREN,S.G.,1984.
C        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
C        APPLIED OPTICS,23,1206-1225
C
C     INPUT PARAMETERS
C
C     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
C           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
C           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
C           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
C     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
C     T = TEMPERATURE ( DEGREES KELVIN )
C
C     OUTPUT PARAMETERS
C
C     RN = REAL PORTION ( SCATTERING )
C     CN = COMPLEX PORTION ( ABSORPTION )
C     ABSIND = ABSORPTIVE INDEX ( CN/RN )
C     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NWL=468   )
      PARAMETER(NWLT=62   )
      DIMENSION WL(NWL),WLT(NWLT)
      DIMENSION TABRE(NWL),TABRET(NWLT,4),TABIM(NWL),TABIMT(NWLT,4)
      DIMENSION TEMREF(4)
C
C     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA
C
      DATA TEMREF/272.16,268.16,253.16,213.16/
C
      DATA WLMIN,WLMAX/0.045,8.6E6/
      DATA CUTICE/167.0/
C
      DATA (WL(I),I=1,114)/
     +0.4430E-01,0.4510E-01,0.4590E-01,0.4680E-01,0.4770E-01,0.4860E-01,
     +0.4960E-01,0.5060E-01,0.5170E-01,0.5280E-01,0.5390E-01,0.5510E-01,
     +0.5640E-01,0.5770E-01,0.5900E-01,0.6050E-01,0.6200E-01,0.6360E-01,
     +0.6530E-01,0.6700E-01,0.6890E-01,0.7080E-01,0.7290E-01,0.7380E-01,
     +0.7510E-01,0.7750E-01,0.8000E-01,0.8270E-01,0.8550E-01,0.8860E-01,
     +0.9180E-01,0.9300E-01,0.9540E-01,0.9920E-01,0.1033E+00,0.1078E+00,
     +0.1100E+00,0.1127E+00,0.1140E+00,0.1181E+00,0.1210E+00,0.1240E+00,
     +0.1272E+00,0.1295E+00,0.1305E+00,0.1319E+00,0.1333E+00,0.1348E+00,
     +0.1362E+00,0.1370E+00,0.1378E+00,0.1387E+00,0.1393E+00,0.1409E+00,
     +0.1425E+00,0.1435E+00,0.1442E+00,0.1450E+00,0.1459E+00,0.1468E+00,
     +0.1476E+00,0.1480E+00,0.1485E+00,0.1494E+00,0.1512E+00,0.1531E+00,
     +0.1540E+00,0.1550E+00,0.1569E+00,0.1580E+00,0.1589E+00,0.1610E+00,
     +0.1625E+00,0.1648E+00,0.1669E+00,0.1692E+00,0.1713E+00,0.1737E+00,
     +0.1757E+00,0.1779E+00,0.1802E+00,0.1809E+00,0.1821E+00,0.1833E+00,
     +0.1843E+00,0.1850E+00,0.1860E+00,0.1870E+00,0.1880E+00,0.1890E+00,
     +0.1900E+00,0.1910E+00,0.1930E+00,0.1950E+00,0.2100E+00,0.2500E+00,
     +0.3000E+00,0.3500E+00,0.4000E+00,0.4100E+00,0.4200E+00,0.4300E+00,
     +0.4400E+00,0.4500E+00,0.4600E+00,0.4700E+00,0.4800E+00,0.4900E+00,
     +0.5000E+00,0.5100E+00,0.5200E+00,0.5300E+00,0.5400E+00,0.5500E+00/
      DATA (WL(I),I=115,228)/
     +0.5600E+00,0.5700E+00,0.5800E+00,0.5900E+00,0.6000E+00,0.6100E+00,
     +0.6200E+00,0.6300E+00,0.6400E+00,0.6500E+00,0.6600E+00,0.6700E+00,
     +0.6800E+00,0.6900E+00,0.7000E+00,0.7100E+00,0.7200E+00,0.7300E+00,
     +0.7400E+00,0.7500E+00,0.7600E+00,0.7700E+00,0.7800E+00,0.7900E+00,
     +0.8000E+00,0.8100E+00,0.8200E+00,0.8300E+00,0.8400E+00,0.8500E+00,
     +0.8600E+00,0.8700E+00,0.8800E+00,0.8900E+00,0.9000E+00,0.9100E+00,
     +0.9200E+00,0.9300E+00,0.9400E+00,0.9500E+00,0.9600E+00,0.9700E+00,
     +0.9800E+00,0.9900E+00,0.1000E+01,0.1010E+01,0.1020E+01,0.1030E+01,
     +0.1040E+01,0.1050E+01,0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,
     +0.1100E+01,0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01,
     +0.1160E+01,0.1170E+01,0.1180E+01,0.1190E+01,0.1200E+01,0.1210E+01,
     +0.1220E+01,0.1230E+01,0.1240E+01,0.1250E+01,0.1260E+01,0.1270E+01,
     +0.1280E+01,0.1290E+01,0.1300E+01,0.1310E+01,0.1320E+01,0.1330E+01,
     +0.1340E+01,0.1350E+01,0.1360E+01,0.1370E+01,0.1380E+01,0.1390E+01,
     +0.1400E+01,0.1410E+01,0.1420E+01,0.1430E+01,0.1440E+01,0.1449E+01,
     +0.1460E+01,0.1471E+01,0.1481E+01,0.1493E+01,0.1504E+01,0.1515E+01,
     +0.1527E+01,0.1538E+01,0.1563E+01,0.1587E+01,0.1613E+01,0.1650E+01,
     +0.1680E+01,0.1700E+01,0.1730E+01,0.1760E+01,0.1800E+01,0.1830E+01,
     +0.1840E+01,0.1850E+01,0.1855E+01,0.1860E+01,0.1870E+01,0.1890E+01/
      DATA (WL(I),I=229,342)/
     +0.1905E+01,0.1923E+01,0.1942E+01,0.1961E+01,0.1980E+01,0.2000E+01,
     +0.2020E+01,0.2041E+01,0.2062E+01,0.2083E+01,0.2105E+01,0.2130E+01,
     +0.2150E+01,0.2170E+01,0.2190E+01,0.2220E+01,0.2240E+01,0.2245E+01,
     +0.2250E+01,0.2260E+01,0.2270E+01,0.2290E+01,0.2310E+01,0.2330E+01,
     +0.2350E+01,0.2370E+01,0.2390E+01,0.2410E+01,0.2430E+01,0.2460E+01,
     +0.2500E+01,0.2520E+01,0.2550E+01,0.2565E+01,0.2580E+01,0.2590E+01,
     +0.2600E+01,0.2620E+01,0.2675E+01,0.2725E+01,0.2778E+01,0.2817E+01,
     +0.2833E+01,0.2849E+01,0.2865E+01,0.2882E+01,0.2899E+01,0.2915E+01,
     +0.2933E+01,0.2950E+01,0.2967E+01,0.2985E+01,0.3003E+01,0.3021E+01,
     +0.3040E+01,0.3058E+01,0.3077E+01,0.3096E+01,0.3115E+01,0.3135E+01,
     +0.3155E+01,0.3175E+01,0.3195E+01,0.3215E+01,0.3236E+01,0.3257E+01,
     +0.3279E+01,0.3300E+01,0.3322E+01,0.3345E+01,0.3367E+01,0.3390E+01,
     +0.3413E+01,0.3436E+01,0.3460E+01,0.3484E+01,0.3509E+01,0.3534E+01,
     +0.3559E+01,0.3624E+01,0.3732E+01,0.3775E+01,0.3847E+01,0.3969E+01,
     +0.4099E+01,0.4239E+01,0.4348E+01,0.4387E+01,0.4444E+01,0.4505E+01,
     +0.4547E+01,0.4560E+01,0.4580E+01,0.4719E+01,0.4904E+01,0.5000E+01,
     +0.5100E+01,0.5200E+01,0.5263E+01,0.5400E+01,0.5556E+01,0.5714E+01,
     +0.5747E+01,0.5780E+01,0.5814E+01,0.5848E+01,0.5882E+01,0.6061E+01,
     +0.6135E+01,0.6250E+01,0.6289E+01,0.6329E+01,0.6369E+01,0.6410E+01/
      DATA (WL(I),I=343,456)/
     +0.6452E+01,0.6494E+01,0.6579E+01,0.6667E+01,0.6757E+01,0.6897E+01,
     +0.7042E+01,0.7143E+01,0.7246E+01,0.7353E+01,0.7463E+01,0.7576E+01,
     +0.7692E+01,0.7812E+01,0.7937E+01,0.8065E+01,0.8197E+01,0.8333E+01,
     +0.8475E+01,0.8696E+01,0.8929E+01,0.9091E+01,0.9259E+01,0.9524E+01,
     +0.9804E+01,0.1000E+02,0.1020E+02,0.1031E+02,0.1042E+02,0.1053E+02,
     +0.1064E+02,0.1075E+02,0.1087E+02,0.1100E+02,0.1111E+02,0.1136E+02,
     +0.1163E+02,0.1190E+02,0.1220E+02,0.1250E+02,0.1282E+02,0.1299E+02,
     +0.1316E+02,0.1333E+02,0.1351E+02,0.1370E+02,0.1389E+02,0.1408E+02,
     +0.1429E+02,0.1471E+02,0.1515E+02,0.1538E+02,0.1563E+02,0.1613E+02,
     +0.1639E+02,0.1667E+02,0.1695E+02,0.1724E+02,0.1818E+02,0.1887E+02,
     +0.1923E+02,0.1961E+02,0.2000E+02,0.2041E+02,0.2083E+02,0.2222E+02,
     +0.2260E+02,0.2305E+02,0.2360E+02,0.2460E+02,0.2500E+02,0.2600E+02,
     +0.2857E+02,0.3100E+02,0.3333E+02,0.3448E+02,0.3564E+02,0.3700E+02,
     +0.3824E+02,0.3960E+02,0.4114E+02,0.4276E+02,0.4358E+02,0.4458E+02,
     +0.4550E+02,0.4615E+02,0.4671E+02,0.4736E+02,0.4800E+02,0.4878E+02,
     +0.5003E+02,0.5128E+02,0.5275E+02,0.5350E+02,0.5424E+02,0.5500E+02,
     +0.5574E+02,0.5640E+02,0.5700E+02,0.5746E+02,0.5840E+02,0.5929E+02,
     +0.6000E+02,0.6100E+02,0.6125E+02,0.6250E+02,0.6378E+02,0.6467E+02,
     +0.6558E+02,0.6655E+02,0.6760E+02,0.6900E+02,0.7053E+02,0.7300E+02/
      DATA (WL(I),I=457,468)/
     +0.7500E+02,0.7629E+02,0.8000E+02,0.8297E+02,0.8500E+02,0.8680E+02,
     +0.9080E+02,0.9517E+02,0.1000E+03,0.1200E+03,0.1500E+03,0.1670E+03/
      DATA  WLT/
     +                                 0.1670E+03,0.1778E+03,0.1884E+03,
     +0.1995E+03,0.2113E+03,0.2239E+03,0.2371E+03,0.2512E+03,0.2661E+03,
     +0.2818E+03,0.2985E+03,0.3162E+03,0.3548E+03,0.3981E+03,0.4467E+03,
     +0.5012E+03,0.5623E+03,0.6310E+03,0.7943E+03,0.1000E+04,0.1259E+04,
     +0.2500E+04,0.5000E+04,0.1000E+05,0.2000E+05,0.3200E+05,0.3500E+05,
     +0.4000E+05,0.4500E+05,0.5000E+05,0.6000E+05,0.7000E+05,0.9000E+05,
     +0.1110E+06,0.1200E+06,0.1300E+06,0.1400E+06,0.1500E+06,0.1600E+06,
     +0.1700E+06,0.1800E+06,0.2000E+06,0.2500E+06,0.2900E+06,0.3200E+06,
     +0.3500E+06,0.3800E+06,0.4000E+06,0.4500E+06,0.5000E+06,0.6000E+06,
     +0.6400E+06,0.6800E+06,0.7200E+06,0.7600E+06,0.8000E+06,0.8400E+06,
     +0.9000E+06,0.1000E+07,0.2000E+07,0.5000E+07,0.8600E+07/
      DATA (TABRE(I),I=1,114)/
     +   0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,
     +   0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,
     +   0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,
     +   0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,
     +   1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,
     +   1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,
     +   1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,
     +   1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,
     +   1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,
     +   1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,
     +   1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,
     +   1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,
     +   1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,
     +   1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,
     +   1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,
     +   1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,
     +   1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,
     +   1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,
     +   1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/
      DATA (TABRE(I),I=115,228)/
     +   1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,
     +   1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,
     +   1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,
     +   1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,
     +   1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,
     +   1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,
     +   1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,
     +   1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,
     +   1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,
     +   1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,
     +   1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,
     +   1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,
     +   1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,
     +   1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,
     +   1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,
     +   1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,
     +   1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,
     +   1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,
     +   1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/
      DATA (TABRE(I),I=229,342)/
     +   1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,
     +   1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,
     +   1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,
     +   1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,
     +   1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,
     +   1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,
     +   1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,
     +   1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,
     +   0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,
     +   1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,
     +   1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,
     +   1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,
     +   1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,
     +   1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,
     +   1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,
     +   1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,
     +   1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,
     +   1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,
     +   1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/
      DATA (TABRE(I),I=343,456)/
     +   1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,
     +   1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,
     +   1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,
     +   1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,
     +   1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,
     +   1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,
     +   1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,
     +   1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,
     +   1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,
     +   1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,
     +   1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,
     +   1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,
     +   1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,
     +   1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,
     +   1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,
     +   1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,
     +   1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,
     +   1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,
     +   1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/
      DATA (TABRE(I),I=457,468)/
     +   1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,
     +   1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/
      DATA (TABRET(I,1),I=1,NWLT)/
     +                                    1.82961,   1.83258,   1.83149,
     +   1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,
     +   1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,
     +   1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,
     +   1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,
     +   1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78800/
      DATA (TABRET(I,2),I=1,NWLT)/
     +                         1.82961,   1.83258,   1.83149,   1.82748,
     +   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,
     +   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,
     +   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,
     +   1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78720/
      DATA(TABRET(I,3),I=1,NWLT)/
     +              1.82961,   1.83258,   1.83149,   1.82748,   1.82224,
     +   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,
     +   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,
     +   1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,
     +   1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,
     +   1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78400,   1.78450/
      DATA (TABRET(I,4),I=1,NWLT)/
     +   1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,
     +   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,
     +   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,
     +   1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77800/
      DATA(TABIM(I),I=1,114)/
     +0.1640E+00,0.1730E+00,0.1830E+00,0.1950E+00,0.2080E+00,0.2230E+00,
     +0.2400E+00,0.2500E+00,0.2590E+00,0.2680E+00,0.2790E+00,0.2970E+00,
     +0.3190E+00,0.3400E+00,0.3660E+00,0.3920E+00,0.4160E+00,0.4400E+00,
     +0.4640E+00,0.4920E+00,0.5170E+00,0.5280E+00,0.5330E+00,0.5340E+00,
     +0.5310E+00,0.5240E+00,0.5100E+00,0.5000E+00,0.4990E+00,0.4680E+00,
     +0.3800E+00,0.3600E+00,0.3390E+00,0.3180E+00,0.2910E+00,0.2510E+00,
     +0.2440E+00,0.2390E+00,0.2390E+00,0.2440E+00,0.2470E+00,0.2240E+00,
     +0.1950E+00,0.1740E+00,0.1720E+00,0.1800E+00,0.1940E+00,0.2130E+00,
     +0.2430E+00,0.2710E+00,0.2890E+00,0.3340E+00,0.3440E+00,0.3820E+00,
     +0.4010E+00,0.4065E+00,0.4050E+00,0.3890E+00,0.3770E+00,0.3450E+00,
     +0.3320E+00,0.3150E+00,0.2980E+00,0.2740E+00,0.2280E+00,0.1980E+00,
     +0.1720E+00,0.1560E+00,0.1100E+00,0.8300E-01,0.5800E-01,0.2200E-01,
     +0.1000E-01,0.3000E-02,0.1000E-02,0.3000E-03,0.1000E-03,0.3000E-04,
     +0.1000E-04,0.3000E-05,0.1000E-05,0.7000E-06,0.4000E-06,0.2000E-06,
     +0.1000E-06,0.6377E-07,0.3750E-07,0.2800E-07,0.2400E-07,0.2200E-07,
     +0.1900E-07,0.1750E-07,0.1640E-07,0.1590E-07,0.1325E-07,0.8623E-08,
     +0.5504E-08,0.3765E-08,0.2710E-08,0.2510E-08,0.2260E-08,0.2080E-08,
     +0.1910E-08,0.1540E-08,0.1530E-08,0.1550E-08,0.1640E-08,0.1780E-08,
     +0.1910E-08,0.2140E-08,0.2260E-08,0.2540E-08,0.2930E-08,0.3110E-08/
      DATA(TABIM(I),I=115,228)/
     +0.3290E-08,0.3520E-08,0.4040E-08,0.4880E-08,0.5730E-08,0.6890E-08,
     +0.8580E-08,0.1040E-07,0.1220E-07,0.1430E-07,0.1660E-07,0.1890E-07,
     +0.2090E-07,0.2400E-07,0.2900E-07,0.3440E-07,0.4030E-07,0.4300E-07,
     +0.4920E-07,0.5870E-07,0.7080E-07,0.8580E-07,0.1020E-06,0.1180E-06,
     +0.1340E-06,0.1400E-06,0.1430E-06,0.1450E-06,0.1510E-06,0.1830E-06,
     +0.2150E-06,0.2650E-06,0.3350E-06,0.3920E-06,0.4200E-06,0.4440E-06,
     +0.4740E-06,0.5110E-06,0.5530E-06,0.6020E-06,0.7550E-06,0.9260E-06,
     +0.1120E-05,0.1330E-05,0.1620E-05,0.2000E-05,0.2250E-05,0.2330E-05,
     +0.2330E-05,0.2170E-05,0.1960E-05,0.1810E-05,0.1740E-05,0.1730E-05,
     +0.1700E-05,0.1760E-05,0.1820E-05,0.2040E-05,0.2250E-05,0.2290E-05,
     +0.3040E-05,0.3840E-05,0.4770E-05,0.5760E-05,0.6710E-05,0.8660E-05,
     +0.1020E-04,0.1130E-04,0.1220E-04,0.1290E-04,0.1320E-04,0.1350E-04,
     +0.1330E-04,0.1320E-04,0.1320E-04,0.1310E-04,0.1320E-04,0.1320E-04,
     +0.1340E-04,0.1390E-04,0.1420E-04,0.1480E-04,0.1580E-04,0.1740E-04,
     +0.1980E-04,0.2500E-04,0.5400E-04,0.1040E-03,0.2030E-03,0.2708E-03,
     +0.3511E-03,0.4299E-03,0.5181E-03,0.5855E-03,0.5899E-03,0.5635E-03,
     +0.5480E-03,0.5266E-03,0.4394E-03,0.3701E-03,0.3372E-03,0.2410E-03,
     +0.1890E-03,0.1660E-03,0.1450E-03,0.1280E-03,0.1030E-03,0.8600E-04,
     +0.8220E-04,0.8030E-04,0.8500E-04,0.9900E-04,0.1500E-03,0.2950E-03/
      DATA(TABIM(I),I=229,342)/
     +0.4687E-03,0.7615E-03,0.1010E-02,0.1313E-02,0.1539E-02,0.1588E-02,
     +0.1540E-02,0.1412E-02,0.1244E-02,0.1068E-02,0.8414E-03,0.5650E-03,
     +0.4320E-03,0.3500E-03,0.2870E-03,0.2210E-03,0.2030E-03,0.2010E-03,
     +0.2030E-03,0.2140E-03,0.2320E-03,0.2890E-03,0.3810E-03,0.4620E-03,
     +0.5480E-03,0.6180E-03,0.6800E-03,0.7300E-03,0.7820E-03,0.8480E-03,
     +0.9250E-03,0.9200E-03,0.8920E-03,0.8700E-03,0.8900E-03,0.9300E-03,
     +0.1010E-02,0.1350E-02,0.3420E-02,0.7920E-02,0.2000E-01,0.3800E-01,
     +0.5200E-01,0.6800E-01,0.9230E-01,0.1270E+00,0.1690E+00,0.2210E+00,
     +0.2760E+00,0.3120E+00,0.3470E+00,0.3880E+00,0.4380E+00,0.4930E+00,
     +0.5540E+00,0.6120E+00,0.6250E+00,0.5930E+00,0.5390E+00,0.4910E+00,
     +0.4380E+00,0.3720E+00,0.3000E+00,0.2380E+00,0.1930E+00,0.1580E+00,
     +0.1210E+00,0.1030E+00,0.8360E-01,0.6680E-01,0.5400E-01,0.4220E-01,
     +0.3420E-01,0.2740E-01,0.2200E-01,0.1860E-01,0.1520E-01,0.1260E-01,
     +0.1060E-01,0.8020E-02,0.6850E-02,0.6600E-02,0.6960E-02,0.9160E-02,
     +0.1110E-01,0.1450E-01,0.2000E-01,0.2300E-01,0.2600E-01,0.2900E-01,
     +0.2930E-01,0.3000E-01,0.2850E-01,0.1730E-01,0.1290E-01,0.1200E-01,
     +0.1250E-01,0.1340E-01,0.1400E-01,0.1750E-01,0.2400E-01,0.3500E-01,
     +0.3800E-01,0.4200E-01,0.4600E-01,0.5200E-01,0.5700E-01,0.6900E-01,
     +0.7000E-01,0.6700E-01,0.6500E-01,0.6400E-01,0.6200E-01,0.5900E-01/
      DATA(TABIM(I),I=343,456)/
     +0.5700E-01,0.5600E-01,0.5500E-01,0.5700E-01,0.5800E-01,0.5700E-01,
     +0.5500E-01,0.5500E-01,0.5400E-01,0.5200E-01,0.5200E-01,0.5200E-01,
     +0.5200E-01,0.5000E-01,0.4700E-01,0.4300E-01,0.3900E-01,0.3700E-01,
     +0.3900E-01,0.4000E-01,0.4200E-01,0.4400E-01,0.4500E-01,0.4600E-01,
     +0.4700E-01,0.5100E-01,0.6500E-01,0.7500E-01,0.8800E-01,0.1080E+00,
     +0.1340E+00,0.1680E+00,0.2040E+00,0.2480E+00,0.2800E+00,0.3410E+00,
     +0.3790E+00,0.4090E+00,0.4220E+00,0.4220E+00,0.4030E+00,0.3890E+00,
     +0.3740E+00,0.3540E+00,0.3350E+00,0.3150E+00,0.2940E+00,0.2710E+00,
     +0.2460E+00,0.1980E+00,0.1640E+00,0.1520E+00,0.1420E+00,0.1280E+00,
     +0.1250E+00,0.1230E+00,0.1160E+00,0.1070E+00,0.7900E-01,0.7200E-01,
     +0.7600E-01,0.7500E-01,0.6700E-01,0.5500E-01,0.4500E-01,0.2900E-01,
     +0.2750E-01,0.2700E-01,0.2730E-01,0.2890E-01,0.3000E-01,0.3400E-01,
     +0.5300E-01,0.7550E-01,0.1060E+00,0.1350E+00,0.1761E+00,0.2229E+00,
     +0.2746E+00,0.3280E+00,0.3906E+00,0.4642E+00,0.5247E+00,0.5731E+00,
     +0.6362E+00,0.6839E+00,0.7091E+00,0.6790E+00,0.6250E+00,0.5654E+00,
     +0.5433E+00,0.5292E+00,0.5070E+00,0.4883E+00,0.4707E+00,0.4203E+00,
     +0.3771E+00,0.3376E+00,0.3056E+00,0.2835E+00,0.3170E+00,0.3517E+00,
     +0.3902E+00,0.4509E+00,0.4671E+00,0.4779E+00,0.4890E+00,0.4899E+00,
     +0.4873E+00,0.4766E+00,0.4508E+00,0.4193E+00,0.3880E+00,0.3433E+00/
      DATA(TABIM(I),I=457,468)/
     +0.3118E+00,0.2935E+00,0.2350E+00,0.1981E+00,0.1865E+00,0.1771E+00,
     +0.1620E+00,0.1490E+00,0.1390E+00,0.1200E+00,0.9620E-01,0.8300E-01/
      DATA(TABIMT(I,1),I=1,NWLT)/
     +                                 0.8300E-01,0.6900E-01,0.5700E-01,
     +0.4560E-01,0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,
     +0.1760E-01,0.1665E-01,0.1620E-01,0.1550E-01,0.1470E-01,0.1390E-01,
     +0.1320E-01,0.1250E-01,0.1180E-01,0.1060E-01,0.9540E-02,0.8560E-02,
     +0.6210E-02,0.4490E-02,0.3240E-02,0.2340E-02,0.1880E-02,0.1740E-02,
     +0.1500E-02,0.1320E-02,0.1160E-02,0.8800E-03,0.6950E-03,0.4640E-03,
     +0.3400E-03,0.3110E-03,0.2940E-03,0.2790E-03,0.2700E-03,0.2640E-03,
     +0.2580E-03,0.2520E-03,0.2490E-03,0.2540E-03,0.2640E-03,0.2740E-03,
     +0.2890E-03,0.3050E-03,0.3150E-03,0.3460E-03,0.3820E-03,0.4620E-03,
     +0.5000E-03,0.5500E-03,0.5950E-03,0.6470E-03,0.6920E-03,0.7420E-03,
     +0.8200E-03,0.9700E-03,0.1950E-02,0.5780E-02,0.9700E-02/
      DATA(TABIMT(I,2),I=1,NWLT)/
     +                      0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,
     +0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,0.1760E-01,
     +0.1665E-01,0.1600E-01,0.1500E-01,0.1400E-01,0.1310E-01,0.1230E-01,
     +0.1150E-01,0.1080E-01,0.9460E-02,0.8290E-02,0.7270E-02,0.4910E-02,
     +0.3300E-02,0.2220E-02,0.1490E-02,0.1140E-02,0.1060E-02,0.9480E-03,
     +0.8500E-03,0.7660E-03,0.6300E-03,0.5200E-03,0.3840E-03,0.2960E-03,
     +0.2700E-03,0.2520E-03,0.2440E-03,0.2360E-03,0.2300E-03,0.2280E-03,
     +0.2250E-03,0.2200E-03,0.2160E-03,0.2170E-03,0.2200E-03,0.2250E-03,
     +0.2320E-03,0.2390E-03,0.2600E-03,0.2860E-03,0.3560E-03,0.3830E-03,
     +0.4150E-03,0.4450E-03,0.4760E-03,0.5080E-03,0.5400E-03,0.5860E-03,
     +0.6780E-03,0.1280E-02,0.3550E-02,0.5600E-02/
      DATA(TABIMT(I,3),I=1,NWLT)/
     +           0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,0.3790E-01,
     +0.3140E-01,0.2620E-01,0.2190E-01,0.1880E-01,0.1660E-01,0.1540E-01,
     +0.1470E-01,0.1350E-01,0.1250E-01,0.1150E-01,0.1060E-01,0.9770E-02,
     +0.9010E-02,0.7660E-02,0.6520E-02,0.5540E-02,0.3420E-02,0.2100E-02,
     +0.1290E-02,0.7930E-03,0.5700E-03,0.5350E-03,0.4820E-03,0.4380E-03,
     +0.4080E-03,0.3500E-03,0.3200E-03,0.2550E-03,0.2120E-03,0.2000E-03,
     +0.1860E-03,0.1750E-03,0.1660E-03,0.1560E-03,0.1490E-03,0.1440E-03,
     +0.1350E-03,0.1210E-03,0.1160E-03,0.1160E-03,0.1170E-03,0.1200E-03,
     +0.1230E-03,0.1320E-03,0.1440E-03,0.1680E-03,0.1800E-03,0.1900E-03,
     +0.2090E-03,0.2160E-03,0.2290E-03,0.2400E-03,0.2600E-03,0.2920E-03,
     +0.6100E-03,0.1020E-02,0.1810E-02/
      DATA(TABIMT(I,4),I=1,NWLT)/
     +0.8300E-01,0.6900E-01,0.5700E-01,0.4450E-01,0.3550E-01,0.2910E-01,
     +0.2440E-01,0.1970E-01,0.1670E-01,0.1400E-01,0.1235E-01,0.1080E-01,
     +0.8900E-02,0.7340E-02,0.6400E-02,0.5600E-02,0.5000E-02,0.4520E-02,
     +0.3680E-02,0.2990E-02,0.2490E-02,0.1550E-02,0.9610E-03,0.5950E-03,
     +0.3690E-03,0.2670E-03,0.2510E-03,0.2290E-03,0.2110E-03,0.1960E-03,
     +0.1730E-03,0.1550E-03,0.1310E-03,0.1130E-03,0.1060E-03,0.9900E-04,
     +0.9300E-04,0.8730E-04,0.8300E-04,0.7870E-04,0.7500E-04,0.6830E-04,
     +0.5600E-04,0.4960E-04,0.4550E-04,0.4210E-04,0.3910E-04,0.3760E-04,
     +0.3400E-04,0.3100E-04,0.2640E-04,0.2510E-04,0.2430E-04,0.2390E-04,
     +0.2370E-04,0.2380E-04,0.2400E-04,0.2460E-04,0.2660E-04,0.4450E-04,
     +0.8700E-04,0.1320E-03/
C
      DATA PI/3.14159265/
C
C     ZERO PARAMETERS
C
      RN=0.0
      CN=0.0
      ABSIND=0.0
      ABSCOF=0.0
C
C     CONVERT WAVELENGTH TO MICRONS
C
      ALAM=XLAM
      IF(IUNIT.EQ.1)ALAM=1000*ALAM
      IF(IUNIT.EQ.2)ALAM=10000*ALAM
      IF(IUNIT.EQ.3)ALAM=10000*(1.0/ALAM)
      IF(ALAM.LT.WLMIN.OR.ALAM.GT.WLMAX)RETURN
      IF(ALAM.GT.CUTICE)GO TO 10
C
C     REGION FROM 0.045 MICRONS TO 167.0 MICRONS - NO TEMPERATURE DEPEND
C
      DO 1 I=2,NWL
      IF(ALAM.LT.WL(I)) GO TO 2
 1    CONTINUE
 2    X1=DLOG(WL(I-1))
      X2=DLOG(WL(I))
      Y1=TABRE(I-1)
      Y2=TABRE(I)
      X=DLOG(ALAM)
      Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      RN=Y
      Y1=DLOG(ABS(TABIM(I-1)))
      Y2=DLOG(ABS(TABIM(I)))
      Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      CN=EXP(Y)
      GO TO 20
C
C     REGION FROM 167.0 MICRONS TO 8.6 METERS - TEMPERATURE DEPENDENCE
C
 10   TK=T
      IF(TK.GT.TEMREF(1))TK=TEMREF(1)
      IF(TK.LT.TEMREF(4))TK=TEMREF(4)
      DO 11 I=2,4
      IF(TK.GE.TEMREF(I)) GO TO 12
 11   CONTINUE
 12   LT1=I
      LT2=I-1
      DO 13 I=2,NWLT
      IF(ALAM.LE.WLT(I)) GO TO 14
 13   CONTINUE
 14   X1=DLOG(WLT(I-1))
      X2=DLOG(WLT(I))
      Y1=TABRET(I-1,LT1)
      Y2=TABRET(I,LT1)
      X=DLOG(ALAM)
      YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y1=TABRET(I-1,LT2)
      Y2=TABRET(I,LT2)
      YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      T1=TEMREF(LT1)
      T2=TEMREF(LT2)
      Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
      RN=Y
      Y1=DLOG(ABS(TABIMT(I-1,LT1)))
      Y2=DLOG(ABS(TABIMT(I,LT1)))
      YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y1=DLOG(ABS(TABIMT(I-1,LT2)))
      Y2=DLOG(ABS(TABIMT(I,LT2)))
      YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
      CN=EXP(Y)
C
C     ABSORPTIVE QUANITIES
C
 20   ABSIND=CN/RN
      ABSCOF=4.0*PI*CN/ALAM
      RETURN
      END


      subroutine Max_Garn(m_matrix,m_inclusion, f_vol_incl,m_MG)
      implicit none
C     Purpose:
C     Compute the Maxweel and Garnett refarctive index
C
      
      complex*16  m_matrix, m_inclusion,m_MG,beta,Im
      real*8 f_vol_incl,f_v 
       
         f_v=f_vol_incl
c         beta=2*m_matrix**2/(m_inclusion**2-m_matrix**2)*
c     $  (m_inclusion**2/(m_inclusion**2-m_matrix**2)*
c     $   log(m_inclusion**2/m_matrix**2)-1)
         beta=3.0*m_matrix**2/(m_inclusion**2+2.0*m_matrix**2)
         m_MG=sqrt((m_matrix**2.*(1-f_v)+f_v*beta*m_inclusion**2.)
     $/(1-f_v+f_v*beta))


      return
      end



        subroutine refr_indx_liq(freq, tav, denliq, refre, refim)
c   Tav in K
c  freq in GHz
C     THIS PROGRAM CALCULATES THE COMPLEX INDEX OF REFRACTION
C     FOR WATER DROPLETS, FOR ANY FREQUENCY AND TEPERATURE.  No
c     idea where this code comes from.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DCN,SQ
c
c     Compute some useful constants
      PI = 2.0*asin( 1.0 )
      TEMP = TAV - 273.
      WAVE = 30./FREQ
C
C                ** CALCULATE THE REAL AND IM. PART OF REFR. INDEX
C
      EINF = 5.27137 + 0.021647*TEMP - 0.00131198*TEMP*TEMP
      ALFA = -16.8129/(TEMP + 273) + 0.0609265
      WAVS = 0.00033836*EXP(2513.98/(TEMP + 273))
      SIGM = 12.5664E+08
      EGAM = 78.54*(1.0 - 4.579E-03*(TEMP - 25.0) + 1.19E-05*
     &       (TEMP - 25.0)**2 - 2.8E-08*(TEMP - 25.0)**3)
C     CALCULATE SOME INTERMEDIATE PRODUCTS
      A1 = EGAM - EINF
      A2 = (WAVS/WAVE)**(1-ALFA)
      A3 = SIN(ALFA*PI/2)
      A4 = COS(ALFA*PI/2)
      A5 = 1 + 2*A2*A3 + A2**2
C     CALCULATE DCNR,DCNI
      DCNR = EINF + A1*(1+A2*A3)/A5
      DCNI = A1*A2*A4/A5 + SIGM*WAVE/18.8496E+10
      DCN = CMPLX(DCNR,DCNI)
      SQ = CDSQRT(DCN)
      REFRE = DREAL(SQ)
      REFIM = DIMAG(SQ)
    
C
C                  ** IF THE DENSITY OF LIQUID IS NOT EQUAL TO 1.0
C                  ** THEN REFRACTIVE INDEX IS AVERAGED WITH AIR
C                  ** REFREAL=1.0, REFRIMAG=0.0 TO OBTAIN PROPER
C                  ** DENSITY.
C
      REFRE = (DENLIQ/1.0d0)*REFRE + (1.0d0-DENLIQ)
      REFIM = (DENLIQ/1.0d0)*REFIM 
c       write(18,*)'retet',sq,refre,refim
C
      return
      end      
c

    
      subroutine refr_indx_ice(freq, tav, refre, refim)
c    TAV in K
c    freq in GHz
c
C     THIS PROGRAM CALCULATES THE COMPLEX INDEX OF REFRACTION
C     FOR ice DROPLETS, FOR ANY FREQUENCY AND TEPERATURE

      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX DCN,SQ
C
C     COMPUTE SOME USEFUL CONSTANTS
      !TEMP = TAV - 273.
      !WAVE = 30./FREQ
      !PI = 2.0*ASIN( 1.0 )
C
C                    ** CALCULATE THE REAL AND IM. PART OF REFR. INDEX
C
      REFRE = 1.780
C     DETERMINE THE IMAGINARY INDEX OF REFRACTION BASED ON FREQ. AND TEMP.
      IF (TAV .GT. 264.0) ICET = 0
      IF ((TAV .LE. 264.0) .AND. (TAV .GT. 254.0)) ICET = 10
      IF ((TAV .LE. 254.0) .AND. (TAV .GT. 242.0)) ICET = 20
      IF ((TAV .LE. 242.0) .AND. (TAV .GT. 222.0)) ICET = 40
      IF (TAV .LE. 222.0) ICET = 60

      IF (ICET .EQ.  0) GOTO 100
      IF (ICET .EQ. 10) GOTO 200
      IF (ICET .EQ. 20) GOTO 300
      IF (ICET .EQ. 40) GOTO 400
      IF (ICET .EQ. 60) GOTO 500

100   IF (FREQ .LT. 30.0)                           REFIM = 0.00265
      IF ((FREQ .GE. 30.0) .AND. (FREQ .LT. 60.0))  REFIM = 0.0036
      IF ((FREQ .GE. 60.0) .AND. (FREQ .LT. 120.0)) REFIM = 0.0054
      IF (FREQ .GT. 120.0)                          REFIM = 0.0074
      GOTO 999

200   IF (FREQ .LT. 30.0)                           REFIM = 0.0016
      IF ((FREQ .GE. 30.0) .AND. (FREQ .LT. 60.0))  REFIM = 0.0023
      IF ((FREQ .GE. 60.0) .AND. (FREQ .LT. 120.0)) REFIM = 0.0038
      IF (FREQ .GT. 120.0)                          REFIM = 0.0058
      GOTO 999

300   IF (FREQ .LT. 30.0)                           REFIM = 0.00095
      IF ((FREQ .GE. 30.0) .AND. (FREQ .LT. 60.0))  REFIM = 0.0015
      IF ((FREQ .GE. 60.0) .AND. (FREQ .LT. 120.0)) REFIM = 0.0027
      IF (FREQ .GT. 120.0)                          REFIM = 0.0046
      GOTO 999

400   IF (FREQ .LT. 30.0)                           REFIM = 0.0006
      IF ((FREQ .GE. 30.0) .AND. (FREQ .LT. 60.0))  REFIM = 0.0010
      IF ((FREQ .GE. 60.0) .AND. (FREQ .LT. 120.0)) REFIM = 0.0019
      IF (FREQ .GT. 120.0)                          REFIM = 0.0030
      GOTO 999

500   IF (FREQ .LT. 30.0)                           REFIM = 0.0004
      IF ((FREQ .GE. 30.0) .AND. (FREQ .LT. 60.0))  REFIM = 0.0007
      IF ((FREQ .GE. 60.0) .AND. (FREQ .LT. 120.0)) REFIM = 0.0013
      IF (FREQ .GT. 120.0)                          REFIM = 0.0020
      GOTO 999

 999  CONTINUE
      RETURN
      END





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




       SUBROUTINE fleg2(x,pl,nl)
      implicit none
      INTEGER nl
      REAL*8 x,pl(nl)
      INTEGER j
      REAL*8 d,f1,f2,twox
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.d0
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.d0
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +>k-5V1`..
     






       subroutine snow_SS_param(wave,rad,type,qsca,qext,asym)
c     Compute the extinction, scattering efficinecy, asymmetry parameter
c     for a given radius of equimass ice

c     Input:
c     wave 	       wavelength [mm]
c     rad               radius of equivalent mass ice particles [mm]
c     type		type of particles 
c  possible option 'snowA', 'snowB','rosettes','KimC1','KimC2',
c        'KimC3','KimC4'  
c   NOTe that au upper boundary of x_max=4.0d0 is applied for Liu
c and an upper boundary of  x_max=5.0d0 is applied for Kim,
c   i.e. if x is larger the efficiencies and asym parameter are evaluated
c   for this maximum value.
c     
c     
c     Output:
c     qsca		scattering efficiency  [adimensional number]
c     qabs		absorption efficiency
c     asym		asymmetry factor 
c   
c

      implicit none



      integer i,jj,n_max
      
      real *8 freq,rad,x_max
   
      real *8 wave,pi, num,Deltar,x_snow,b(3),a(3),c(0:3),
     $aK(8),f(8),bK(6),sum1,sum2
      real *8 qext,qabs,qsca, asym
  
      CHARACTER type*5
 

          pi=dacos(-1.d0)   
        
        x_max=4.0d0  !upper threshold (see fig12 Liu)
        if(type(1:3).eq.'Kim')  x_max=5.0d0
   
        x_snow=2.0*pi*rad/wave !size parameter adim
      
        if(x_snow.gt.x_max) x_snow=x_max  !putting a threshold to x_max

         b(1)=0.7446*1e-2    !parameter for absorption
         b(2)=0.010607d0
         b(3)=-0.14505*1e-2
        
                             !parameters for scattering cross section
          if(type.eq.'snowA') then           
             if(x_snow.le.1.4) then
             a(1)=-0.036379
             a(2)=0.11716 
             a(3)=0.18637
             else
             a(1)=-0.1622
             a(2)=0.56253 
             a(3)=-0.066369
             endif

           elseif(type.eq.'snowB') then
             if(x_snow.le.0.5) then
             a(1)=-0.036379d0
             a(2)=0.11716d0 
             a(3)=0.18637d0
             else
             a(1)=-0.0096948d0
             a(2)=0.15898d0 
             a(3)=0.01078d0
             endif
          ELSEIF(type.eq.'roset') then
             if(x_snow.le.2.2) then
             a(1)=-0.036379d0
             a(2)=0.11716d0 
             a(3)=0.18637d0
             else
             a(1)=-0.60643d0
             a(2)=1.0934d0 
             a(3)=-0.14630d0
             endif
          ENDIF
            !parameters for g parameter
          if(x_snow.le.0.825d0) then  
              c(0)=0.d0 
              c(1)=-0.077361 
              c(2)=0.59902 
              c(3)=-0.18825*1e-2 
        else
         IF(type.eq.'roset') then
              c(0)=0.30617 
              c(1)=0.019795 
              c(2)=0.029307 
              c(3)=-0.29968*1e-3 
              ELSE
                 if(x_snow.le.1.0d0) then  
              c(0)=0.d0 
              c(1)=-0.077361 
              c(2)=0.59902 
              c(3)=-0.18825*1e-2
               else 
              c(0)=0.42725 
              c(1)=0.062429 
              c(2)=0.028416 
              c(3)=-0.42245*1e-2
               endif 
             ENDIF
         endif
        



c        write(*,*)'a,b,c',a,b,c   
        qabs=0.0d0
        qsca=0.0d0
        asym=c(0)
        do jj=1,3
          qsca=qsca+a(jj)*x_snow**jj  !formula(9) Liu
          qabs=qabs+b(jj)*x_snow**jj
          asym=asym+c(jj)*x_snow**jj
        enddo

          if(type(1:3).ne.'Kim') goto 300
          if(type.eq.'KimC1') then
          aK(1)=-0.3353d0
          aK(2)= 3.3177d0
          aK(3)= -1.7217d0
          aK(4)= -1.7254d0
          aK(5)= -0.1953d0
          aK(6)= 0.7358d0
          aK(7)= 0.4084d0
          aK(8)= 0.0544d0
          f(1)=-0.6304d0
          f(2)= 1.5281d0
          f(3)= -0.2125d0
          f(4)= -0.9502d0
          f(5)= -1.7090d0
          f(6)= 0.1557d0
          f(7)= 1.4016d0
          f(8)= 0.5477d0
            elseif(type.eq.'KimC2') then
          aK(1)=-0.3533d0
          aK(2)=  3.3295d0
          aK(3)=  -1.6769d0
          aK(4)= -1.9710d0 
          aK(5)= -0.5256d0
          aK(6)=  1.1379d0 
          aK(7)=1.1043d0
          aK(8)= 0.2963d0 
          f(1)=-0.5673d0
          f(2)=  1.5418d0 
          f(3)= -1.0410d0 
          f(4)= -1.0442d0  
          f(5)= -0.0600d0
          f(6)= 0.8422d0
          f(7)=0.6686d0
          f(8)= 0.1597d0
                elseif(type.eq.'KimC3') then
          aK(1)=-0.3597d0    
          aK(2)=3.3643d0
          aK(3)=-1.5013d0
          aK(4)=-2.0822d0 
          aK(5)= -1.2714d0
          aK(6)= 0.9382d0 
          aK(7)= 1.6981d0
          aK(8)=0.6088d0 
          f(1)=-0.5832d0
          f(2)= 1.6818d0 
          f(3)=-1.0855d0 
          f(4)= -1.4262d0  
          f(5)=-0.2155d0
          f(6)=1.0944d0
          f(7)=0.8690d0
          f(8)=0.1937d0
           elseif(type.eq.'KimC4') then
          aK(1)=-0.3432d0    
          aK(2)= 3.4542d0
          aK(3)=-1.4338d0
          aK(4)= -2.6021d0 
          aK(5)= -2.2706d0
          aK(6)=1.1111d0 
          aK(7)=2.8529d0
          aK(8)=1.1258d0 
           If(x_snow.le.1.0d0) then 
          f(1)=-0.6122d0
          f(2)= 2.3329d0 
          f(3)=3.6036d0 
          f(4)=13.9784d0  
          f(5)=26.3336d0
          f(6)=26.3125d0
          f(7)=13.4166d0
          f(8)=2.7443d0
           Else
          f(1)=-0.4654d0     
          f(2)= -3.9724d0 
          f(3)=81.0301d0 
          f(4)=-504.904d0  
          f(5)=1569.3d0
          f(6)=-2620.1d0
          f(7)=2230.9d0
          f(8)= -757.586d0
          Endif
            endif         
            IF (freq.lt.100.0d0) then
           bK(1)=1.508*1e-4 
           bK(2)=0.0021d0
           bK(3)= 0.0081d0
           bK(4)= -0.0051d0
           bK(5)= 0.002d0
           bK(6)= -2.59*1e-4
           ELSEIF (freq.lt.200.0d0) then
           bK(1)=1.122*1e-4
           bK(2)= 0.0061d0
           bK(3)= 0.0086d0 
           bK(4)=-0.0022d0
           bK(5)= 5.35*1e-4 
           bK(6)= -4.82*1e-5
           ENDIF
         
           sum1=0.0d0  
           sum2=0.0d0

           qabs=0.0d0 
          do jj=1,6       
          qabs=qabs+bK(jj)*x_snow**(jj-1)
        enddo
         do jj=1,8
          sum1=sum1+aK(jj)*(log10(x_snow))**(jj-1) 
          sum2=sum2+f(jj)*(log10(x_snow))**(jj-1) 
        enddo
          
           qsca=10.d0**(sum1)
           asym=10.d0**(sum2)
c           write(*,*)'x',x_snow,qsca,qabs,asym
 300        continue

        qabs=max(0.0d0,qabs)
        qsca=max(0.0d0,qsca)
c        qsca=min(6.0,qsca)     !to cut high values at large x
c        asym=min(0.9,asym)
        qext=qabs+qsca


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



