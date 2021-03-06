! FOM function which compare intermediate scattering function S(q,t) data 
! with calculated S(q,t) equivalent where
!
!   S(q,t) = integral ( exp(i*q*r)*G(r,t)*dr )
!
! G(r,t) is the space-time correlation function (van Hove correlation function). It
! assumes here units of AA^-3, and it follows that S(q,t) is dimensionless. 
! Time is measured in units of 10^-13s=0.1ps and q is AA^-1.
!
!   FOM = sum_i ( S_i^{data} - scale_factor * S_i^{cal} )^2
!
! where S_i is short for the i'th S(q,t) array point
!
    
module s_q_time_fom_class
use various_constants_class
use structure_class
use func_param_class
use s_q_time_class

implicit none

  public :: s_qt_fom_val

  type s_qt_fom_container
    real(db), dimension(:,:), allocatable :: obs  ! holds the data
    
    real(db), dimension(:), allocatable :: q 
    real(db) :: t_bin ! lenght of a time bin
    integer :: n_t_bin ! # of time bins
    
    real(db) :: scale_factor = 1.0 ! apply scale factor to calculated data
    
    character(len=120) :: title = " "
    
    ! stuff which is in common for all function containers
    
    type (func_params) :: params
  end type s_qt_fom_container


contains

  function s_qt_fom_val(s_qt_cal, s_qt_data) result(val)
    type (s_q_time), intent(in) :: s_qt_cal
  	type (s_qt_fom_container), intent(in) :: s_qt_data
    real (db) :: val
    
    integer :: i_q, i_t, n_q
      
    ! First check if dimension of cal and obs fom containers fit.
    
    if ( maxval(s_qt_data%q-s_qt_cal%q) > 0.0001 .or. &
         s_qt_data%n_t_bin /= get_s_q_time_n_t_bin(s_qt_cal) .or. &
         (s_qt_data%t_bin-s_qt_cal%t_bin) > 0.0001 ) then
      write(*,*) " "
      write(*,*) "ERROR in s_qt_fom_val"
      write(*,*) "mis-match between cal and obs s(q,t) data"
      stop
    end if

    val = 0.0
    n_q = size(s_qt_data%q)

    do i_t = 1, s_qt_data%n_t_bin     
      do i_q = 1, n_q 
         val = val + (s_qt_data%obs(i_q, i_t) - s_qt_data%scale_factor*(s_qt_cal%self(i_q, i_t) + s_qt_cal%diff(i_q,i_t)))**2  
      end do 
    end do
  
  end function s_qt_fom_val
  
end module s_q_time_fom_class
