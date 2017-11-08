! FOM function which compare "normalised" g_d(r,t) data with calculated g_d(r,t) equivalent where
!
!   (g^norm)_d(r,t) =           V*hist_d(r,t)
!                     --------------------------------
!                     N(N-1)*volume_of_spherical_shell(r)
!
! See also Eq. (29) page 30 in my handwritten notes
!
! g_d(r,t) is the distinct part of the space-time pair correlation function
! and is related very simply to the space-time correlation function G_d(r,t)
! (also called van Hove correlation function) as G_d(r,t) = rho * g_d(r,t) 
! where rho is N/V i.e. the density. 
!
! (g^norm)_d(r,t) = N * g_d(r,t) / (N-1) which has the visual convenient 
! property that -> 1 when t->infinity and r->infinity
!
! In practice not sure if you would ever have real data in this format,
! but this FOM may be useful nevertheless in some cases
!
!   FOM = sum_i ( g_i^{data} - g_i^{cal} )^2
!
! where g_i is short for the i'th (g^norm)_d(r,t) array point
!

module g_d_rt_fom_class
use various_constants_class
use structure_class
use func_param_class
use time_corr_hist_container_class

implicit none

  public :: g_d_rt_fom_val


  type g_d_rt_fom_container
    real(db), dimension(:,:), allocatable :: obs  ! holds the data
    
    real(db) :: t_bin ! lenght of a time bin
    real(db) :: r_bin ! lenght of a bin in the r dimension
    integer :: n_t_bin ! # of time bins
    integer :: n_r_bin ! # of r bins
    
    ! which part of time and r axes of the g_d_data to 
    ! include for calculation fom 
    integer, dimension(2) :: t_cut, r_cut 
    
    ! both copied from initial structure generated by 
    ! handler if <g-d-rt-fom> present in XML input file
    ! Used to normalise time_corr_container values to g_d
    ! - not sure whether this is done in the best possible way.
    ! Perhaps include these in time_corr_container or
    ! something instead if at all.
    integer :: n_atom 
    real(db) :: density
    
    !real(db) :: scale_factor = 1.0
    
    character(len=120) :: title = " "
    
    ! stuff which is in common for all function containers
    type (func_params) :: params
  end type g_d_rt_fom_container


contains

  function g_d_rt_fom_val(time_corr, g_d_data) result(val)
    type (time_corr_hist_container), intent(in) :: time_corr
  	type (g_d_rt_fom_container), intent(in) :: g_d_data
    real (db) :: val
    
    integer :: i_t ! iterator in time direction
    integer :: i_r ! iterator in r direction
    
    real(db) :: g_d_cal
    real(db), dimension(:), allocatable :: prefac
    
    ! First check if dimension of cal and obs fom containers fit.
    ! For now makes the assumption that they are the same!
    if ( g_d_data%n_t_bin /= get_time_corr_hist_n_time_bin(time_corr) .or. &
         g_d_data%n_r_bin /= get_time_corr_hist_n_r_bin(time_corr) .or. &
         g_d_data%r_bin /= get_time_corr_hist_r_bin(time_corr) ) then
      write(*,*) " "
      write(*,*) "ERROR in g_d_rt_fom_val"
      write(*,*) "mis-match between cal and obs g_d_data"
      stop
    end if
    
    
    allocate(prefac(size(time_corr%volume_prefac)))
    

    val = 0.0
    
    prefac = time_corr%volume_prefac / (time_corr%n_accum*g_d_data%density*(g_d_data%n_atom-1))

    do i_t = 1, g_d_data%n_t_bin   
      do i_r = 1, g_d_data%n_r_bin
         g_d_cal = time_corr%g_d_hists_sum(i_t)%val(i_r)*prefac(i_r)
         val = val + (g_d_data%obs(i_r, i_t) - g_d_cal)**2
      end do 
    end do
    
    deallocate(prefac) 
  
  end function g_d_rt_fom_val
  
    
end module g_d_rt_fom_class
