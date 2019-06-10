program liebe_test
	implicit none
	real(kind=8),dimension(4) :: T, P, RH
	real(kind=8),dimension(3) :: KEXTATMO, AVGPRE, VAPPRE

	T = (/272.99, 273.1233, 273.2567,273.39/)
	P = (/10296, 10283, 10269, 10256/)
	RH = (/-57.2, -24.5,8.133,40.8/)

	call get_atmosg(T,P,RH,4,3,AVGPRE,VAPPRE,22.24,KEXTATMO)
	
	print*, T, P, RH
	print*, AVGPRE
	print*, VAPPRE
	print('(E9.4)'), KEXTATMO

	stop
end program
