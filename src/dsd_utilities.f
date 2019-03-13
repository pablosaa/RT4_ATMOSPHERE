C
C All subroutines for the Program: data_base_forward_matrix
C_____________________________________________________________________
C
C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
C    N - NUMBER OF POINTS                                             *
C    Z - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************
 
      SUBROUTINE GAUSS_L(N,IND1,IND2,Z,W)
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
C 
C**********************************************************************
C    CALCULATION OF POINTS                                            *
C    Gaussia theory: define  X(i), ALPHA_G(i), BETA_G(i)              *
C    N=NGAUSS - NUMBER OF POINTS                                      *
C**********************************************************************
C
      SUBROUTINE GAUSS_point(N,Y,Dmax_min,Dmin,ALPHAmax_min,ALPHAmin,
     &                       BETAmax_min,BETAmin,
     &                       X,ALPHA,BETA)
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Y(N),X(N),ALPHA(N),BETA(N)
      DO 1 i=1,N
            Y1i=Y(i)
            X(i)=Y1i*Dmax_min+Dmin
            ALPHA(i)=Y1i*ALPHAmax_min+ALPHAmin
            BETA(i)=Y1i*BETAmax_min+BETAmin
    1 CONTINUE
      RETURN
      END
C 
C**********************************************************************
C    CALCULATION OF POINTS ONLY FOR DSD                               *
C    Gaussia theory: define  X(i)                                     *
C    N=NGAUSS - NUMBER OF POINTS                                      *
C**********************************************************************
C
      SUBROUTINE GAUSS_point_for_dsd(N,Y,Dmax_min,Dmin,X)
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Y(N),X(N)
      DO 1 i=1,N
            Y1i=Y(i)
            X(i)=Y1i*Dmax_min+Dmin
    1 CONTINUE
      RETURN
      END     
C 

C 
C**********************************************************************
C    CALCULATION location for interpolation                           *
C    Copr. 1986-92 Numerical Recipes Software +>k-5V1`                *
C                                                                     *
C**********************************************************************
C      
      SUBROUTINE locate(xx,n,x,j)
c    Given an array xx(1:n) and given a value x, returns a value j such 
c    that x is between xx(j) and xx(j+1) xx(1:n) must be monotonic, either
c  increasing or decreasing. j=0 or j=n is returned to indicate 
c    that x is out of range

      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
       jl=0     !initialise lower and
      ju=n+1    ! upper boundaries
10    if(ju-jl.gt.1)then   !if we are nmot yet done
        jm=(ju+jl)/2       !compute a midpoint
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm    !and replace either the lower
        else
          ju=jm     !or the upper limit
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
C
C**********************************************************************
C This program computes the ntegration over a Gamma DSD expressed in 
c terms of Diameters using the Gaussian integration
C   INPUT parameters:
c   
c                                                                  *
cC**********************************************************************
C
      SUBROUTINE DSD_gamma(NGAUSS,Y,W,Dmax,Dmin,NPNG2,
     &                     D,ngridy,mu,N_w,D_0,N_elev,
     &                     av_ab_ReS11,av_ab_ImS11,
     +                     av_ab_ReS22,av_ab_ImS22,
     &                     av_abD_ReS11,av_abD_ImS11,
     &                     av_abD_ReS22,av_abD_ImS22,R,Z) 
      IMPLICIT none
      integer Ngauss,ngridy,j,NPNG2,jj,N_elev,kk
c     REAL*8 (A-H,P-Z)
       REAL*8 N_D,N_w,mu,gammln,X(NGAUSS),Y(NGAUSS),W(NGAUSS),
     &       D(NPNG2),  Gamma,R,Z, 
     $ av_abD_ImS11(N_elev), av_abD_ImS22(N_elev),
     $  av_abD_ReS11(N_elev),av_abD_ReS22(N_elev),
     +    av_ab_ReS11(NPNG2,N_elev),av_ab_ImS11(NPNG2,N_elev),
     +    av_ab_ReS22(NPNG2,N_elev),av_ab_ImS22(NPNG2,N_elev),pi,
     $ av_ab_ReS11j(N_elev),av_ab_ImS11j(N_elev),
     $ av_ab_ReS22j(N_elev),D_0,
     $ av_ab_ImS22j(N_elev),Xj,Wj,f_mu,Dmax_min,Dmin,Dmax,v_D
      

       pi=dacos(-1.0d0)
C
C N_D is gamma distribution
C    N_D is in m^-3mm^-1
C N_w is the normalized number concentration parameter
C    N_w is in m^-3mm^-1
C v_D is the fall speed of rain drops, m s^-1
C R is the rain rate in mm h^-1
C Z is the radar reflectivity factor in dBZ    
C
C calculate Gamma function
      Gamma=DEXP(gammln(mu+4D0))
c        write(19,*)'gamma',gamma,mu
C
      f_mu=(6D0/(3.67D0**4))*((3.67D0+mu)**(mu+4D0))/Gamma
      do kk=1,N_elev
      av_abD_ImS11(kk)=0.D0
      av_abD_ImS22(kk)=0.D0
      av_abD_ReS11(kk)=0.D0
      av_abD_ReS22(kk)=0.D0
      enddo
      R=0D0
      Z=0D0
      Dmax_min=min(Dmax-Dmin,7.0*D_0)  
C Gaussia theory: define  X(i)
      call GAUSS_point_for_dsd(NGAUSS,Y,Dmax_min,Dmin,X)
C
C calculation rain rate R, radar reflectivity factor Zi
        
      DO 2 j=1,NGAUSS
            Xj=X(j)
C Xj is diameter
            Wj=W(j)
            N_D=N_w*f_mu*(Xj/D_0)**mu*dexp(-(3.67D0+mu)*Xj/D_0)  !1/(m^3 mm)
            v_D=9.65D0-10.3D0*DEXP(-0.6D0*Xj)
c            write(19,*)'vel',Dmax_min,v_D,f_mu,Xj,N_D,Wj
            R=R+(pi/6)*N_D*Wj*Xj**3*v_D
c            write(19,*)'R',R
            Z=Z+N_D*Wj*Xj**6
C            
C location for interpolation
C
      call locate(D,ngridy,Xj,jj)
c      write(14,*)'xj',xj,jj 
C
C interpolation
C
      do kk=1,N_elev
      av_ab_ReS11j(kk)=av_ab_ReS11(jj,kk)+((Xj-D(jj))/(D(jj+1)-D(jj)))
     &              *(av_ab_ReS11(jj+1,kk)-av_ab_ReS11(jj,kk))
      av_ab_ImS11j(kk)=av_ab_ImS11(jj,kk)+((Xj-D(jj))/(D(jj+1)-D(jj)))
     &              *(av_ab_ImS11(jj+1,kk)-av_ab_ImS11(jj,kk))
      av_ab_ReS22j(kk)=av_ab_ReS22(jj,kk)+((Xj-D(jj))/(D(jj+1)-D(jj)))
     &              *(av_ab_ReS22(jj+1,kk)-av_ab_ReS22(jj,kk))
      av_ab_ImS22j(kk)=av_ab_ImS22(jj,kk)+((Xj-D(jj))/(D(jj+1)-D(jj)))
     &              *(av_ab_ImS22(jj+1,kk)-av_ab_ImS22(jj,kk))         
C
C average of Im(S11), Im(S22), Re(S11),Re(S22) in beta, alpha, D                    
C
       av_abD_ImS11(kk)=av_abD_ImS11(kk)+N_D*Wj*av_ab_ImS11j(kk)
       av_abD_ImS22(kk)=av_abD_ImS22(kk)+N_D*Wj*av_ab_ImS22j(kk)
       av_abD_ReS11(kk)=av_abD_ReS11(kk)+N_D*Wj*av_ab_ReS11j(kk)
       av_abD_ReS22(kk)=av_abD_ReS22(kk)+N_D*Wj*av_ab_ReS22j(kk)
       enddo
    2       CONTINUE
C
            R=Dmax_min*R *3.6D-03  !mm/h
            Z=10D0*log10(Dmax_min*Z)   !dbZ 
              do kk=1,N_elev
            av_abD_ImS11(kk)=Dmax_min*av_abD_ImS11(kk)
            av_abD_ImS22(kk)=Dmax_min*av_abD_ImS22(kk)
            av_abD_ReS11(kk)=Dmax_min*av_abD_ReS11(kk)
            av_abD_ReS22(kk)=Dmax_min*av_abD_ReS22(kk)
              enddo
      RETURN
      END  
C 

C**********************************************************************
C This program computes the ntegration over a Gamma DSD expressed in 
c terms of Diameters using the Gaussian integration
C   INPUT parameters:
c   NGAUSS number of Gaussian points
c   Y gaussian points between 0 and 1
c   W weights for Gaussian integration
c   Dmax and Dmin min and maximum diamtere where to integrate
c   
cC**********************************************************************
C
      SUBROUTINE DSD_gamma_field(NGAUSS,Y,W,Dmax,Dmin,N_LUT,
     & D_LUT,field_LUT,mu,N_w,D_0,av_DSD_field) 
      IMPLICIT none
      integer Ngauss,j,N_LUT,jj
c     REAL*8 (A-H,P-Z)
       REAL*8 N_D,N_w,mu,gammln,X(NGAUSS),Y(NGAUSS),W(NGAUSS),
     & D_LUT(N_LUT),Gamma,
     $ av_DSD_field,field_LUT(N_LUT),pi,
     $ av_ab_ReS11j,
     $D_0,Xj,Wj,f_mu,Dmax_min,Dmin,Dmax
      

       pi=dacos(-1.0d0)
C
C N_D is gamma distribution
C    N_D is in m^-3mm^-1
C N_w is the normalized number concentration parameter
C    N_w is in m^-3mm^-1
C
C calculate Gamma function
      Gamma=DEXP(gammln(mu+4D0))
c        write(19,*)'gamma',gamma,mu
C
      f_mu=(6D0/(3.67D0**4))*((3.67D0+mu)**(mu+4D0))/Gamma
      av_DSD_field=0.D0

      Dmax_min=min(Dmax-Dmin,7.0*D_0-Dmin)  
C Gaussia theory: define  X(i)
c      CALL GAUSS_L(NGAUSS,0,0,Y,W)
      call GAUSS_point_for_dsd(NGAUSS,Y,Dmax_min,Dmin,X)
C
C calculation rain rate R, radar reflectivity factor Zi
        
      DO 2 j=1,NGAUSS
            Xj=X(j)
C Xj is diameter
            Wj=W(j)
            N_D=N_w*f_mu*(Xj/D_0)**mu*dexp(-(3.67D0+mu)*Xj/D_0)  !1/(m^3 mm)
C            
C location for interpolation
C
      call locate(D_LUT,N_LUT,Xj,jj)
c      write(14,*)'xj',xj,jj 
C
C interpolation
C
      av_ab_ReS11j=field_LUT(jj)+((Xj-D_LUT(jj))/
     $ (D_LUT(jj+1)-D_LUT(jj)))*(field_LUT(jj+1)-field_LUT(jj))       
         av_DSD_field=av_DSD_field+N_D*Wj*av_ab_ReS11j    
    2       CONTINUE
            av_DSD_field=Dmax_min*av_DSD_field
      RETURN
      END  

      SUBROUTINE DSD_gamma_vector(NGAUSS,Y,W,Dmax,Dmin,N_vec,
     $ N_LUT,N_temp,j_temp,D_LUT,field_LUT,
     $ mu,N_w,D_0,av_DSD_field) 
      IMPLICIT none
      integer Ngauss,j,N_LUT,jj,kk,N_vec,N_temp,j_temp
c     REAL*8 (A-H,P-Z)
       REAL*8 N_D,N_w,mu,gammln,X(NGAUSS),Y(NGAUSS),W(NGAUSS),
     & D_LUT(N_LUT),Gamma,
     $ av_DSD_field(N_vec),field_LUT(N_LUT,N_vec,N_temp),pi,
     $ av_ab_ReS11j,
     $D_0,Xj,Wj,f_mu,Dmax_min,Dmin,Dmax
      

       pi=dacos(-1.0d0)
C
C N_D is gamma distribution
C    N_D is in m^-3mm^-1
C N_w is the normalized number concentration parameter
C    N_w is in m^-3mm^-1
C
C calculate Gamma function
      Gamma=DEXP(gammln(mu+4D0))
c        write(26,*)'gamma',gamma,mu,D_0,Dmax,Dmin,N_w
C
      f_mu=(6D0/(3.67D0**4))*((3.67D0+mu)**(mu+4D0))/Gamma
      do kk=1,N_vec
      av_DSD_field(kk)=0.D0
      enddo
      Dmax_min=min(Dmax-Dmin,7.0*D_0-Dmin)  
C Gaussia theory: define  X(i)
c      CALL GAUSS_L(NGAUSS,0,0,Y,W)
      call GAUSS_point_for_dsd(NGAUSS,Y,Dmax_min,Dmin,X)
C
C calculation rain rate R, radar reflectivity factor Zi
c       write(19,*)NGAUSS,X
      DO 2 j=1,NGAUSS
            Xj=X(j)
C Xj is diameter
            Wj=W(j)
            N_D=N_w*f_mu*(Xj/D_0)**mu*dexp(-(3.67D0+mu)*Xj/D_0)  !1/(m^3 mm)
C            
C location for interpolation
C
      call locate(D_LUT,N_LUT,Xj,jj)
c      write(14,*)'xj',xj,jj 
C
C interpolation
C
       do kk=1,N_vec
      av_ab_ReS11j=field_LUT(jj,kk,j_temp)+((Xj-D_LUT(jj))/
     $ (D_LUT(jj+1)-D_LUT(jj)))*
     $ (field_LUT(jj+1,kk,j_temp)-field_LUT(jj,kk,j_temp))       
         av_DSD_field(kk)=av_DSD_field(kk)+N_D*Wj*av_ab_ReS11j
       enddo    
    2       CONTINUE
            do kk=1,N_vec
            av_DSD_field(kk)=Dmax_min*av_DSD_field(kk)
            enddo
      RETURN
      END  


       SUBROUTINE DSD_gamma_matrix4(NGAUSS,Y,W,Dmax,Dmin,N_1,
     $ N_2,N_3,N_4,N_LUT,N_temp,j_temp,
     $ D_LUT,field_LUT,mu,N_w,D_0,av_DSD_field) 
      IMPLICIT none
      integer Ngauss,j,N_LUT,jj,k1,k2,k3,k4,N_1,N_2,N_3,N_4,
     $N_temp,j_temp
c     REAL*8 (A-H,P-Z)
       REAL*8 N_D,N_w,mu,gammln,X(NGAUSS),Y(NGAUSS),W(NGAUSS),
     & D_LUT(N_LUT),Gamma,
     $ av_DSD_field(N_1,N_2,N_3,N_4),
     $ field_LUT(N_LUT,N_1,N_2,N_3,N_4,N_temp),pi,
     $ av_ab_ReS11j,
     $D_0,Xj,Wj,f_mu,Dmax_min,Dmin,Dmax
      

       pi=dacos(-1.0d0)
C
C N_D is gamma distribution
C    N_D is in m^-3mm^-1
C N_w is the normalized number concentration parameter
C    N_w is in m^-3mm^-1
C
C calculate Gamma function
      Gamma=DEXP(gammln(mu+4D0))
c        write(19,*)'gamma',gamma,mu,D_0,Dmax,Dmin
C
      f_mu=(6D0/(3.67D0**4))*((3.67D0+mu)**(mu+4D0))/Gamma
       do k1=1,N_1  
       do k2=1,N_2  
       do k3=1,N_3  
       do k4=1,N_4 

      av_DSD_field(k1,k2,k3,k4)=0.D0
      enddo  
      enddo  
      enddo  
      enddo
      Dmax_min=min(Dmax-Dmin,7.0*D_0-Dmin)  
C Gaussia theory: define  X(i)
c      CALL GAUSS_L(NGAUSS,0,0,Y,W)
      call GAUSS_point_for_dsd(NGAUSS,Y,Dmax_min,Dmin,X)
C
C calculation rain rate R, radar reflectivity factor Zi
c       write(19,*)NGAUSS,X
      DO 2 j=1,NGAUSS
            Xj=X(j)
C Xj is diameter
            Wj=W(j)
            N_D=N_w*f_mu*(Xj/D_0)**mu*dexp(-(3.67D0+mu)*Xj/D_0)  !1/(m^3 mm)
C            
C location for interpolation
C
      call locate(D_LUT,N_LUT,Xj,jj)
c      write(14,*)'xj',xj,jj 
C
C interpolation
C
       do k1=1,N_1  
       do k2=1,N_2  
       do k3=1,N_3  
       do k4=1,N_4       
      av_ab_ReS11j=field_LUT(jj,k1,k2,k3,k4,j_temp)
     $ +((Xj-D_LUT(jj))/
     $ (D_LUT(jj+1)-D_LUT(jj)))*
     $ (field_LUT(jj+1,k1,k2,k3,k4,j_temp)-
     $  field_LUT(jj,k1,k2,k3,k4,j_temp))       
         av_DSD_field(k1,k2,k3,k4)=
     $   av_DSD_field(k1,k2,k3,k4)+N_D*Wj*av_ab_ReS11j
       enddo 
      enddo  
      enddo  
      enddo   
    2       CONTINUE
       do k1=1,N_1  
       do k2=1,N_2  
       do k3=1,N_3  
       do k4=1,N_4       
      av_DSD_field(k1,k2,k3,k4)=Dmax_min*av_DSD_field(k1,k2,k3,k4)
      enddo
      enddo  
      enddo  
      enddo   

      RETURN
      END  


C 


      SUBROUTINE vect_lin_interp(x,y,N,x1,y1,N1) 
      IMPLICIT none
      integer N,N1,ii,jj
      REAL*8 x(N),y(N),x1(N1),y1(N1)

         do ii=1,N1
         call locate(x,N,x1(ii),jj)        
         y1(ii)=y(jj)+((x1(ii)-x(jj))/(x(jj+1)-x(jj)))
     &              *(y(jj+1)-y(jj))

         enddo
         return
         end


