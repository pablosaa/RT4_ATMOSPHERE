program interp1
  implicit none

  real(kind=8), dimension(10) :: elevations, elvmu
  integer, parameter :: nelv = 10, NUMMU = 32, NOUTLEVELS=2, NSTOKES=2, NANGLES=42
  real(kind=8), dimension(32) :: MU_VALUES
  !real(kind=8), dimension(NANGLES) :: ZZ_THTA
  real(kind=8), dimension(NUMMU,NSTOKES,2,NOUTLEVELS) :: OUTVAR
  !real(kind=8), dimension(NANGLES,NSTOKES,2,NOUTLEVELS) :: TB_THTA
  real(kind=8), allocatable, dimension(:,:,:,:) :: TB_THTA
  real(kind=8), allocatable, dimension(:) :: ZZ_THTA
  integer :: i, j

  MU_VALUES = (/(i,i=1,32)/)
  MU_VALUES = MU_VALUES/real(33)
  elevations = 90-(/90.0,30.0, 19.2, 14.4, 11.4, 8.4, 6.6, 5.4, 4.8, 4.2/)
  elvmu = cos(elevations*3.14159/180.)

  OUTVAR = -999.
  OUTVAR(:,1,1,1) = (/(10*i,i=32,1,-1)/)
  OUTVAR(:,1,1,2) = (/(15*i,i=32,1,-1)/)
  OUTVAR(:,1,2,1) = (/(2*i,i=32,1,-1)/)
  OUTVAR(:,1,2,2) = (/(3*i,i=32,1,-1)/)

  allocate(ZZ_THTA(NANGLES))
  allocate(TB_THTA(NANGLES,NSTOKES,2,NOUTLEVELS))
  !print*,MU_VALUES
  !print*,elevations
  print*,OUTVAR(:,1,1,1)
  call interp1_1D(MU_VALUES,OUTVAR,NUMMU,elvmu,nelv,ZZ_THTA,NANGLES,TB_THTA, NSTOKES,2,NOUTLEVELS)

  !print*,acos(MU_VALUES)*180./3.14159
  !print*,acos(ZZ_THTA)*180./3.14159
  write(*,'(I6F8.3F10.2F10.2)') (i,ZZ_THTA(i),TB_THTA(i,1,1,:), i=1,NANGLES)
  if(allocated(ZZ_THTA)) deallocate(ZZ_THTA)
  if(allocated(TB_THTA)) deallocate(TB_THTA)
  stop
end program interp1

subroutine interp1_1D(Xin, Yin, Nin, X0, N0, Xout, Nout, Yout, Nstk, Nflx, Nlev)
  
  ! Interpolated abssice values must increase monotonically, 
  ! the ordenate values can be decreasing or increasing

  implicit none
  integer, intent(in) :: Nin, N0, Nstk, Nflx, Nlev, Nout
  real(kind=8), intent(in) :: Xin(Nin), X0(N0), Yin(Nin,Nstk,Nflx,Nlev)
  real(kind=8), intent(out):: Xout(Nout)
  real(kind=8), intent(out):: Yout(Nout,Nstk,Nflx,Nlev)
  integer :: h, idx, idn
  real(kind=8) :: a, b, Fa(Nstk, Nflx, Nlev), Fb(Nstk, Nflx, Nlev), F0(Nstk, Nflx, Nlev)
  real(kind=8) :: eps = 1.0E-4
  
  ! union of Xin and X0 sets
  Yout = -69.8
  Yout(1:Nin,:,:,:) = Yin
  
  Xout = -69.0
  Xout(1:Nin) = Xin
  do h = 1, N0
     if(X0(h).GT.maxval(Xout)) then
        idx = maxloc(Xout,1)+1
        a = Xout(idx-2)
        b = Xout(idx-1)
        Fa = Yout(idx-2,:,:,:)
        Fb = Yout(idx-1,:,:,:)
     else if(X0(h).LT.minval(Xout)) then
        idx = 1
        a = Xout(idx)
        b = Xout(idx+1)
        Fa = Yout(idx,:,:,:)
        Fb = Yout(idx+1,:,:,:)
     else
        idx = minloc(Xout,DIM=1,MASK=X0(h)<Xout)
        a = Xout(idx-1)
        b = Xout(idx)
        Fa = Yout(idx-1,:,:,:)
        Fb = Yout(idx,:,:,:)
     end if
     
     if(abs(X0(h)-b).LT.eps) then
        F0 = Fb
     else if(abs(X0(h)-a).LT.eps) then
        F0 = Fa
     else
        F0 = Fa + (X0(h)-a)*(Fb-Fa)/(b-a)
     end if
     idn = Nin+h-1
     Xout(idx+1:) = Xout(idx:idn)
     Xout(idx) = X0(h)

     Yout(idx+1:, :, :, :) = Yout(idx:idn,:,:,:)
     Yout(idx, :, :, :)   = F0
     
  end do

end subroutine interp1_1D
