module s_q_omega_fom_class
use various_constants_class
use structure_class
use func_param_class
use s_q_omega_class

implicit none

  public :: s_qo_fom_val


  type s_qo_fom_container
    real(db), dimension(:,:), allocatable :: obs  ! holds the data
    real(db), dimension(:,:), allocatable :: esd  ! holds esds for the data
    
    real(db), dimension(:), allocatable :: q 
    real(db), dimension(:), allocatable :: omega

    
    real(db) :: scale_factor = 1.0
    !real(db) :: weight = 1.0
    
    character(len=120) :: title = " "
    
    ! stuff which is in common for all function containers
    type (func_params) :: params
  end type s_qo_fom_container


contains

  function s_qo_fom_val(s_qo_cal, s_qo_data) result(val)
    type (s_q_omega), intent(in) :: s_qo_cal
  	type (s_qo_fom_container), intent(in) :: s_qo_data
    real (db) :: val
    
    integer :: i_q, i_o, n_q, n_omega
    
    
    ! First check if dimension of cal and obs fom containers fit.
    ! For now makes the assumption that they are the same!
    !if ( maxval(s_qo_data%q-s_qo_cal%q) > 0.0001 .or. &
    !     maxval(s_qo_data%omega-s_qo_cal%omega) > 0.0001 ) then
    !  write(*,*) " "
    !  write(*,*) "ERROR in s_qo_fom_val"
    !  write(*,*) "mis-match between cal and obs s(q,omega) data"
    !  stop
    !end if

    val = 0.0
    n_q = size(s_qo_data%q)
    n_omega = size(s_qo_data%omega)

    do i_o = 1, n_omega     
      do i_q = 1, n_q 
         if ( s_qo_data%obs(i_q, i_o) /= no_datapoint_available ) then
           val = val + &
           (s_qo_data%obs(i_q, i_o) - s_qo_data%scale_factor*(s_qo_cal%self(i_q, i_o) + s_qo_cal%diff(i_q,i_o)))**2 &
           / (s_qo_data%esd(i_q, i_o)**2)
         end if
      end do 
    end do
  
  end function s_qo_fom_val
  
    
end module s_q_omega_fom_class
