program sample
  use, intrinsic :: iso_c_binding

  implicit none

  interface
     function mitt_unixt(datum) result(val) bind(c, name='unixtime')
       use, intrinsic :: iso_c_binding
       integer(kind=c_int) :: datum(6)
       real(kind=c_double) :: val
     end function mitt_unixt
  end interface

  integer(c_int) :: date(6)
  real(c_double) :: unixtt

  date = (/2004,4,9,12,0,0/)
  unixtt =  mitt_unixt(date)

  print*, date
  print*, unixtt

  stop
end program sample

