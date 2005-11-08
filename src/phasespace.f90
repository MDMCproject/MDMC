module phasespace_class
use structure_class
use near_neighb_class
use common_potential_block_class

  implicit none


  ! This type is only include to speed up the calculation. Most importantly by
  ! replacing the double summation over atoms by a sum over nearest neighboors
  ! and secondly by storing distances that may be used at a later stage.
  ! Notice that this distances are here stored as REAL rather then REAL (DB)
  ! to save memory, since this type will be the single most memory hungry 
  ! component of this code and knowing these distances to 6 digits accuracy
  ! is e.g. definitely good enough for calculating a histogram, pdf etc from

  type phasespace
    type (structure) :: str
    real(db), dimension(:,:), allocatable :: p  ! momenta
    real(db), dimension(:,:), allocatable :: deriv  ! derivs(=- force = -dp/dt)    
  
    ! the position coordinates are here stored in the same format as p and
    ! deriv for coding convenience and probably a slight speed gain
    real(db), dimension(:,:), allocatable :: r  ! atom positions
    
    ! this one is again added for speed gain and convenience
    real(db), dimension(:,:), allocatable :: inv_mass 
    
    ! this one to e.g. later speed up the calculation of the kinetic energy
    ! p_div_mass stands for momentum divided by mass
    real(db), dimension(:,:), allocatable :: p_div_mass
    
    type (near_neighb_list) :: neighb_list
  end type phasespace


contains

  subroutine trajectory_in_phasespace(ps, n_md, delta_t)
    type (phasespace), intent(inout) :: ps
    integer, intent(in) :: n_md
    real(db), intent(in) :: delta_t   
    
    real(db) :: pot_energy, delta_t_half
    integer :: i
    
    delta_t_half = delta_t / 2.0
    
    ! leap-frog algorithm version used in Rapaport, where it is assumed
    ! that the derivatives are up-to-date
    
    do i = 1, n_md-1
      ! update momenta to time t=t+h/2 assuming we have up-to-date derivs 
      ! for time t=t
      ps%p = ps%p - delta_t_half * ps%deriv
      
      ! update positions to time t=t+h and store p_div_mass which can be used 
      ! for e.g. calculating the kinetic energy a bit faster!
      ps%p_div_mass = ps%p * ps%inv_mass
      ps%r = ps%r + delta_t * ps%p_div_mass
      
      ! calculate derivatives for time t=t+h
      
      ! now get momenta in sync with positions and ready for the next loop
      ps%p = ps%p - delta_t_half * ps%deriv      
    end do
    
    ! do exactly the same as in the loop above but this time calculate
    ! the derivatives + core-loop-md-extra
    
    ps%p = ps%p - delta_t_half * ps%deriv    
    ps%p_div_mass = ps%p * ps%inv_mass
    ps%r = ps%r + delta_t * ps%p_div_mass    
    
    ps%p = ps%p - delta_t_half * ps%deriv
    
    pot_energy = gpe_val(ps%str, common_gpe)
                         
    write (*, '(a,f10.3)') "Epot = ", pot_energy
  
  end subroutine trajectory_in_phasespace

  function make_phasespace(str, r_cut, delta_r) result (ps)
    type (structure), intent(in) :: str
    real (db), intent(in) :: r_cut, delta_r
    type (phasespace) :: ps

    integer :: n_tot
    
    n_tot = size(str%atoms)
    
    ps%neighb_list = make_near_neighb_list(str, r_cut, delta_r)
    
    ps%str = copy_structure(str)
    
    allocate(ps%p(n_tot,3))
    allocate(ps%deriv(n_tot,3))
    allocate(ps%r(n_tot,3))
    allocate(ps%inv_mass(n_tot,3))
    allocate(ps%p_div_mass(n_tot,3))
    
    ! add coordinates to 
    
  end function make_phasespace


end module phasespace_class
