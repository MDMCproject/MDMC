! FOM function which compare dynamical structure factor S(q,omega) data 
! with calculated S(q,omega) equivalent where
!
!   S(q,omega) = 1/(2*pi) * integral ( exp(i*[q*r-omega*t])*G(r,t)*dr*dt )
!
! G(r,t) is the space-time correlation function (van Hove correlation function). It
! assumes here units of AA^-3. Time is measured in units of 10^-13s=0.1ps and 
! therefore the unit of S(q,omega) is assumed to 0.1ps. A unit of time is
! also a unit of inverse frequency. Further the unit of q is AA^-1 and unit
! of omega 1/(0.1ps).
!
    
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
    
    ! optional set to true to ignore data esds when calculating this fom
    logical :: ignore_errors = .false.    
    
    character(len=120) :: title = " "
    
    ! stuff which is in common for all function containers
    type (func_params) :: params
  end type s_qo_fom_container


contains

    ! Calculated and observed data arrays must have the same dimensions

    function s_qo_fom_val(s_qo_cal, s_qo_data) result(val)
    type (s_q_omega), intent(in) :: s_qo_cal
  	type (s_qo_fom_container), intent(in) :: s_qo_data
    real (db) :: val
    
    integer :: i_q, i_o, n_q, n_omega

    val = 0.0
    n_q = size(s_qo_data%q)
    n_omega = size(s_qo_data%omega)

    
    ! Check if the array dimensions are the same
    
    if ( n_q /= size(s_qo_cal%q) .or. &
         n_omega /= size(s_qo_cal%omega) ) then
      write(*,*) " "
      write(*,*) "ERROR in s_qo_fom_val"
      write(*,*) "mis-match between calculated and observed s(q,omega) data arrays"
      stop
    end if    
    
         
    ! Note it is currently allowed to have obs data with 'no data available', a bit 
    ! like non-a-number data     
         
    if (s_qo_data%ignore_errors == .true.) then
      do i_o = 1, n_omega     
        do i_q = 1, n_q 
          ! only sum over obs data points
          if ( abs(s_qo_data%obs(i_q, i_o) - no_datapoint_available) > 0.0001 ) then
            val = val + &
                 (s_qo_data%obs(i_q, i_o) - s_qo_data%scale_factor*(s_qo_cal%self(i_q, i_o) + s_qo_cal%diff(i_q,i_o)))**2        
          end if
        end do 
      end do
    else
      do i_o = 1, n_omega     
        do i_q = 1, n_q 
          ! only sum over obs data points
          if ( abs(s_qo_data%obs(i_q, i_o) - no_datapoint_available) > 0.0001 ) then
            val = val + &
                  (s_qo_data%obs(i_q, i_o) - s_qo_data%scale_factor*(s_qo_cal%self(i_q, i_o) + s_qo_cal%diff(i_q,i_o)))**2 &
                  / (s_qo_data%esd(i_q, i_o)**2)          
          end if
        end do 
      end do      
    end if
  
  end function s_qo_fom_val
  
    
end module s_q_omega_fom_class
