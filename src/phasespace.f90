module phasespace_class
use structure_class
use near_neighb_class
use function_class

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
    !real(db), dimension(:,:), allocatable :: r  ! atom positions
    
    ! this one is again added for speed gain and convenience
    real(db), dimension(:,:), allocatable :: inv_mass
    
    ! this one is at present only used to speedup the calculation of ps%v2 in
    ! this module once
     real(db), dimension(:,:), allocatable :: p_div_mass
    
    ! to calculate max velocity and kinetic energy
    real(db), dimension(:), allocatable :: v2
    real(db), dimension(:), allocatable :: mass    
    
    type (near_neighb_list) :: neighb_list
  end type phasespace


contains

  subroutine trajectory_in_phasespace(ps, pe_list, n_md, delta_t, pressure_comp, pot_energy)
    type (phasespace), intent(inout) :: ps
    type (func_list), intent(in) :: pe_list
    integer, intent(in) :: n_md
    real(db), intent(in) :: delta_t   
 		real (db), optional, intent(out) :: pressure_comp, pot_energy 
 		    
    real(db) :: pot_energy_not_pass, delta_t_half
    integer :: i, j, n_tot
    logical :: extra_args
    
    
    ! test for optional arguments and set extra_args
    
    if (present(pressure_comp) .and. present(pot_energy)) then
      extra_args = .true.
    else if (present(pressure_comp)==.false. .and. present(pot_energy)==.false.) then
      extra_args = .false.
    else
      write(*,*) "ERROR in trajectory_in_phasespace"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if

    
    n_tot = size(ps%str%atoms)
    
    delta_t_half = delta_t / 2.0
    
    ! leap-frog algorithm version used in Rapaport, where it is assumed
    ! that the derivatives are up-to-date
    
    do i = 1, n_md
      ! update momenta to time t=t+h/2 assuming we have up-to-date derivs 
      ! for time t=t

      ps%p = ps%p - delta_t_half * ps%deriv
      
    
      ! update positions to time t=t+h. And store p_div_mass array here. Note it 
      ! correspond to time t=t+h/2 and is at present only used to calculate ps%v2 below 

      ps%p_div_mass = ps%p * ps%inv_mass
      
      ps%str%r = ps%str%r + delta_t * ps%p * ps%inv_mass
 
      
      ! calculate derivatives for time t=t+h
      ! but wrap around first
      
      call apply_boundary_condition(ps%str)   
      
      if (ps%neighb_list%ignore_list == .true.) then
        if (extra_args) then
          call func_deriv(ps%str, ps%deriv, pe_list, pressure_comp, pot_energy)
        else
          call func_deriv(ps%str, ps%deriv, pe_list)
        end if
      else
        ! stored the velocities at times t=t+h/2. Strictly speaking
        ! these should not be used for calculating the total energy at 
        ! t=t+h but could but reused for an adjusting-temperature-function

        ps%v2 = sum(ps%p_div_mass*ps%p_div_mass,2)
        
        
        ! update nearest neighbour list flags

        call update_nn_list_flags(sqrt(maxval(ps%v2))*delta_t, ps%neighb_list)
        
        if (ps%neighb_list%needs_updating == .true.) then
          call build_near_neighb(ps%str, ps%neighb_list)
        else
          call update_stored_nn_values(ps%str, ps%neighb_list)
        end if
        
        if (extra_args) then
          call func_deriv_nn(ps%str, ps%deriv, pe_list, ps%neighb_list, pressure_comp, pot_energy)
        else
          call func_deriv_nn(ps%str, ps%deriv, pe_list, ps%neighb_list)
        end if
      end if
      
      ! now get momenta in sync with positions and ready for the next loop
      ps%p = ps%p - delta_t_half * ps%deriv     
      
    end do
        
  end subroutine trajectory_in_phasespace
  
    
  function make_phasespace(str, r_cut, delta_r, temperature) result (ps)
    type (structure), intent(in) :: str
    real (db), intent(in) :: r_cut, delta_r, temperature
    type (phasespace) :: ps

    integer :: n_tot, i, n_seeds
    !real(db), dimension(ndim) :: momentum_sum  ! mainly for debugging    
    real(db) :: momentum_scale
    real(db), dimension(ndim) :: dummy_vec
    
    n_tot = size(str%atoms)
    
    ps%neighb_list = make_near_neighb_list(str, r_cut, delta_r)
    
    ps%str = copy_structure(str)
    
    
    allocate(ps%p(n_tot,3))
    allocate(ps%deriv(n_tot,3))
    allocate(ps%inv_mass(n_tot,3))
    allocate(ps%p_div_mass(n_tot,3))
    allocate(ps%mass(n_tot))
    allocate(ps%v2(n_tot))
    
    ! copy coordinates to r - notice should probably change type of structure
    do i = 1, n_tot
      ps%inv_mass(i,:) = str%atoms(i)%inv_mass
      ps%mass(i) = str%atoms(i)%mass
    end do
    
    ! initial derivatives set to zero. According to Rapaport this is good
    ! enough, although perhaps I simply calculate the true derivatives here!?
    ps%deriv = 0.0
    
    ! determine initial momenta
    ! momentum_scale = temperature*ndim*(n_tot-1)/sum(ps%inv_mass(:,1))
    momentum_scale = sqrt( temperature*ndim*(n_tot-1)/sum(ps%inv_mass(:,1)) ) / sqrt(3.0)
    !write(*,'(a,f12.6)') "momentum_scale ", momentum_scale * sqrt(3.0)
    
    if (modulo(n_tot,2) == 1) then
     ps%p(n_tot,:) = 0.0
    end if
       
    do i = 1, floor(n_tot / 2.0)
      ps%p(2*i-1,:) = momentum_scale
      ps%p(2*i,:) = - momentum_scale
    end do 
    
    ! write(*,'(a,f12.6)') "a check ", sum(ps%p*ps%p) / ((n_tot-1)*3.0)
  
    call random_seed(size=n_seeds)
    call random_seed(put=(/(i, i = 1, n_seeds)/))

    !write(*,'(a,f)') "length to aim for = ", momentum_scale
    !stop
    do i = 1, n_tot
      call random_number(ps%p(i,:))
      !write(*,'(a,f)') "length before rescaling = ", sqrt(sum(ps%p(i,:)*ps%p(i,:)))
      ps%p(i,:) = ps%p(i,:) * momentum_scale / sqrt(sum(ps%p(i,:)*ps%p(i,:)))
      !write(*,'(a,f)') "length after rescaling = ", sqrt(sum(ps%p(i,:)*ps%p(i,:)))
    end do 
  
    !stop
    
    dummy_vec = sum(ps%p,1) / n_tot
    
    do i = 1, n_tot
      ps%p(i,:) = ps%p(i,:) - dummy_vec
    end do
    
    ! what is the total momentum

    !write(*,'(a,3f12.6)') "total momentum in make_phasespace ", sum(ps%p,1)

    !stop
    
  end function make_phasespace


end module phasespace_class
