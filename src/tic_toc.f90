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
  
  
  function get_current_date_and_time() result(date_time)
    character(len=50) :: date_time
    
    integer :: today(3), now(3)
    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
    write (date_time, '(a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i2.2)')  &
         "Date ", today(1), "/", today(2), "/", today(3)+2000, "; time ", &
         now(1), ":", now(2), ":", now(3)
  end function get_current_date_and_time

end module tic_toc_class