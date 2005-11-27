module tic_toc_class
implicit none

  integer :: time_last_tic   ! time: the last time tic was called
  integer :: time_counts_per_sec
  
  ! this is a modified version of the tic-toc model in Matlab and
  ! Figure 4.10 in "Object-oriented programming via Fortran 90/95"
  ! by Ed Akin
contains

  subroutine tic
    call system_clock(time_last_tic, time_counts_per_sec)
  end subroutine tic

  function toc() result(sec)
    real :: sec
    
    integer :: finish
    
    call system_clock(finish)
    sec = 0.0
    if (finish >= time_last_tic) sec = float(finish - time_last_tic) / &
                                       float(time_counts_per_sec)
  end function toc

end module tic_toc_class