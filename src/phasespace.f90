module phasespace_class
use structure_class
use near_neighb_class
use common_potential_block_class

  implicit none

  public :: trajectory_in_phasespace
  public :: make_phasespace
  

  type phasespace
    type (structure) :: str
    real(db), dimension(:,:), allocatable :: p  ! momenta
    real(db), dimension(:,:), allocatable :: deriv  ! derivs(=- force = -dp/dt)    
  
    ! the position coordinates are here stored in the same format as p and
    ! deriv for coding convenience and probably a slight speed gain
    ! NOTICE I SHOULD PROBABLY CHANGE THE STRUCTURE TYPE TO STORE THE ATOMIC
    ! COORS IN THIS FORMAT! (as a consequence of storing the coordinate I
    ! also have to copy forward in trajectory_in_phasespace()
    real(db), dimension(:,:), allocatable :: r  ! atom positions
    
    ! this one is again added for speed gain and convenience
    real(db), dimension(:,:), allocatable :: inv_mass
    
    ! this one to e.g. later speed up the calculation of the kinetic energy
    ! p_div_mass stands for momentum divided by mass
    real(db), dimension(:,:), allocatable :: p_div_mass
    
    ! to calculate max velocity and kinetic energy
    real(db), dimension(:), allocatable :: v2
    real(db), dimension(:), allocatable :: mass    
    
    type (near_neighb_list) :: neighb_list
  end type phasespace


contains

  subroutine trajectory_in_phasespace(ps, n_md, delta_t)
    type (phasespace), intent(inout) :: ps
    integer, intent(in) :: n_md
    real(db), intent(in) :: delta_t   
    
    real(db) :: pot_energy, delta_t_half
    integer :: i, j, n_tot
    
    ! for testing stuff
    real(db) :: sum_r_start, sum_r_end
    real(db) :: sum_v_start, sum_v_intermediate, sum_v_end
    real(db) :: sum_a_start, sum_a_end
    
    ! for testing
    !real(db), dimension(512) :: dummy

    write(*,'(a,3f12.6)') "first atom ", ps%r(1,:)
    write(*,'(a,f12.6)') "sum_r_start = ", sum(ps%r*ps%r)
    write(*,'(a,f12.6)') "sum_v_start = ", sum(ps%p*ps%p)
    write(*,'(a,f12.6)') "sum_a_start = ", sum(ps%deriv*ps%deriv)    
    
    !dummy = sum(ps%r,2)
    !write(*,'(2f12.6)') dummy(1), sum(ps%r(1,:))
    
    n_tot = size(ps%str%atoms)
    
    delta_t_half = delta_t / 2.0
    
    ! leap-frog algorithm version used in Rapaport, where it is assumed
    ! that the derivatives are up-to-date
    
    do i = 1, n_md
      ! update momenta to time t=t+h/2 assuming we have up-to-date derivs 
      ! for time t=t
      ps%p = ps%p - delta_t_half * ps%deriv
      
      write(*,'(a,f12.6)') "sum_v_intermediate = ", sum(ps%p*ps%p)
    
      ! update positions to time t=t+h and store p_div_mass which can be used 
      ! for e.g. calculating the kinetic energy a bit faster!
      ps%p_div_mass = ps%p * ps%inv_mass
      ps%r = ps%r + delta_t * ps%p_div_mass
      
      write(*,'(a,3f12.6)') "first atom after move ", ps%r(1,:)
      write(*,'(a,f12.6)') "sum_r_end = ", sum(ps%r*ps%r)
      
      ! calculate derivatives for time t=t+h
      do j = 1, n_tot
        ps%str%atoms(j)%r = ps%r(j,:)
      end do
      
      if (ps%neighb_list%ignore_list == .true.) then
        call gpe_deriv(ps%str, ps%deriv, common_gpe)
      else
        ! stored the velocities at times t=t+h/2. Strictly speaking
        ! these should not be used for calculating the total energy at 
        ! t=t+h but could but reused for an adjusting-temperature-function
        ps%v2 = sum(ps%r,2)
        
        ! update nearest neighbour list flags
        call update_nn_list_flags(sqrt(maxval(ps%v2))*delta_t, ps%neighb_list)
        
        if (ps%neighb_list%needs_updating == .true.) then
          call build_near_neighb(ps%str, ps%neighb_list)
        else
          call update_stored_nn_values(ps%str, ps%neighb_list)
        end if
        
        call gpe_deriv_nn(ps%str, ps%deriv, common_gpe, ps%neighb_list)
        
      end if
      
      ! now get momenta in sync with positions and ready for the next loop
      ps%p = ps%p - delta_t_half * ps%deriv     
      
      write(*,'(a,f12.6)') "sum_v_end = ", sum(ps%p*ps%p)
      write(*,'(a,f12.6)') "sum_a_end = ", sum(ps%deriv*ps%deriv)
    end do

    
    ! just some checking
    !pot_energy = gpe_val(ps%str, common_gpe)                      
    !write (*, '(a,f12.6)') "Epot1 = ", pot_energy
    
    !call build_near_neighb_with_cell(ps%str, ps%neighb_list)
    !pot_energy = gpe_val_nn(ps%str, common_gpe, ps%neighb_list)
    !write (*, '(a,f12.6)') "Epot2 = ", pot_energy
    
    !call build_near_neighb_without_cell(ps%str, ps%neighb_list)
    !call update_stored_nn_values(ps%str, ps%neighb_list)
    !pot_energy = gpe_val_nn(ps%str, common_gpe, ps%neighb_list)
    !write (*, '(a,f12.6)') "Epot3 = ", pot_energy
        
  end subroutine trajectory_in_phasespace

  
  function make_phasespace(str, r_cut, delta_r, temperature) result (ps)
    type (structure), intent(in) :: str
    real (db), intent(in) :: r_cut, delta_r, temperature
    type (phasespace) :: ps

    integer :: n_tot, i
    real(db), dimension(ndim) :: momentum_sum  ! mainly for debugging
    
    n_tot = size(str%atoms)
    
    ps%neighb_list = make_near_neighb_list(str, r_cut, delta_r)
    
    ps%str = copy_structure(str)
    
    !ps%delta_r = delta_r
    !ps%r_cut = r_cut
    
    allocate(ps%p(n_tot,3))
    allocate(ps%deriv(n_tot,3))
    allocate(ps%r(n_tot,3))
    allocate(ps%inv_mass(n_tot,3))
    allocate(ps%p_div_mass(n_tot,3))
    allocate(ps%mass(n_tot))
    allocate(ps%v2(n_tot))
    
    ! copy coordinates to r - notice should probably change type of structure
    do i = 1, n_tot
      ps%r(i,:) = str%atoms(i)%r
      ps%inv_mass(i,:) = str%atoms(i)%inv_mass
      ps%mass(i) = str%atoms(i)%mass
    end do
    
    ! initial derivatives set to zero. According to Rapaport this is good
    ! enough, although perhaps I simply calculate the true derivatives here!?
    ps%deriv = 0.0
    
    ! determine initial momenta
    
    if (modulo(n_tot,2) == 1) then
     ps%p(n_tot,:) = 0.0
    end if
       
    do i = 1, floor(n_tot / 2.0)
      ps%p(2*i-1,:) = 1.0
      ps%p(2*i,:) = -1.0
    end do 
    
    ! what is the total momentum

    momentum_sum = 0.0
    do i = 1, n_tot
      momentum_sum = momentum_sum + ps%p(i, :)
    end do
    
    write(*,'(a,3f12.6)') "total momentum in make_phasespace ", momentum_sum
    
    !ps%p = 1.0
    
  end function make_phasespace


end module phasespace_class
