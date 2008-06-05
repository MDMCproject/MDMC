module lennard_jones_class
use various_constants_class
use structure_class
use func_param_class
use histogram_class

implicit none

  public :: lj_val
  public :: lj_deriv

  private :: lj_val_nn, lj_val_without_nn

  type lj_pe_container
    ! stuff which is in common for all function containers
    type (func_params) :: params
    type (histogram) :: hist
  end type lj_pe_container


contains

  function lj_val(str, container) result (val)
    type (structure), intent(in) :: str
    type (lj_pe_container), intent(in) :: container
    real (db) :: val    
	
    if (str%nn_list%ignore_list == .true.) then
      val = lj_val_without_nn(str, container)
    else
      ! will remove unnessary 3rd arg later
      val = lj_val_nn(str, container, str%nn_list)
    end if
		
  end function lj_val
  

  subroutine lj_deriv(str, deriv, container, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real (db), dimension(:,:), intent(out) :: deriv
    type (lj_pe_container), intent(in) :: container
    real (db), optional, intent(out) :: pressure_comp, pot_energy
		  
    integer :: i1, i2, n_tot, i
    real (db) :: r_cut, rr_cut 
    real (db) :: sigma, epsilon 
    real (db) :: epsilon_times4    ! 4 * epsilon
    real (db) :: sigma2    ! = sigma*sigma
    real (db) :: epsi48sigma2   ! 48*epsilon/(sigma^2)
    real (db) :: rr, rri, rri3
    real (db), dimension(ndim) :: diff_vec
    real (db) :: prefac  ! used for calculating derivatives
    logical :: extra_args
    
    integer :: j ! used only in nearest neighbour summation loop
    
    
    ! test for optional arguments and set extra_args
    
    if (present(pressure_comp) .and. present(pot_energy)) then
      extra_args = .true.
      pressure_comp = 0.0
      pot_energy = 0.0
    else if (present(pressure_comp)==.false. .and. present(pot_energy)==.false.) then
      extra_args = .false.
    else
      write(*,*) "ERROR in lj_deriv"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if    
    
    
    sigma = get_func_param_val(container%params, "sigma")  !container%params%p(1)%val
    epsilon = get_func_param_val(container%params, "epsilon") !container%params%p(2)%val
    epsilon_times4 = 4*epsilon
    sigma2 = sigma * sigma
    epsi48sigma2 = 48.0 * epsilon / sigma2
        
    
    if (does_func_param_exist(container%params, "r-cut")) then
      r_cut = get_func_param_val(container%params, "r-cut")
    else
      ! otherwise put r_cut to a value higher than any possible distance
      r_cut = sqrt(sum(str%box_edges*str%box_edges))
    end if
    
    
    ! when using neighest neighbour method all distances less than
    ! nn_list%r_cut are 'guarentied' to be in the list but not 
    ! higher distances
    
    if (str%nn_list%ignore_list == .false.) then
      if (str%nn_list%r_cut < r_cut) then
        r_cut = str%nn_list%r_cut
      end if    
    end if
    
    
    rr_cut = r_cut * r_cut
    
    n_tot = size(str%atoms)
    
    ! initiate 
    deriv = 0.0
    
    
    ! Do either single summation loop (nn-method) or double summation
    ! Notice quite a lot of code is repeated twice below. Also, I hope
    ! apply_periodic_box_condition is called 'inline'!?
    
    if (str%nn_list%ignore_list == .true.) then
    
      do i1 = 1, n_tot
        do i2 = i1+1, n_tot
          diff_vec = str%r(i1,:) - str%r(i2,:)
          
          call apply_boundary_condition_to_vector(diff_vec, str%box_edges)
          
          rr = sum(diff_vec*diff_vec)
          
          if (rr < rr_cut) then
            rri = sigma2 / rr
            rri3 = rri * rri * rri
            prefac = epsi48sigma2 * rri3 * (rri3 - 0.5) * rri
            
            ! be aware that derivatives are calculated which are 
            ! equal to - force
            deriv(i1,:) = deriv(i1,:) - prefac * diff_vec
            deriv(i2,:) = deriv(i2,:) + prefac * diff_vec
            
            if (extra_args) then
              ! calculate potential energy
              pot_energy = pot_energy + epsilon_times4 * rri3 * (rri3 - 1.0)
            
              ! calculate sum_i=1^n-1 sum_j>i f_ij |r_ij|^2
              pressure_comp = pressure_comp + prefac * rr     
            end if
                    
          end if
  			
        end do
      end do
    
    else
    
      do j = 1, str%nn_list%n_pairs

  	    i1 = str%nn_list%pairs(2*j-1)
  	    i2 = str%nn_list%pairs(2*j)
        
        
        ! Notice I could speed up the calculation by moving the calculation of
        ! diff_vec into type near_neighb_list with the penalty that this would
        ! take up more memory = nn_list%n_pairs * 3 * 8 bytes
        
        diff_vec = str%r(i1,:) - str%r(i2,:)
          
        call apply_boundary_condition_to_vector(diff_vec, str%box_edges)
        
        rr = str%nn_list%dists(j)
          
        if (rr < rr_cut) then
          rri = sigma2 / rr
          rri3 = rri * rri * rri
          prefac = epsi48sigma2 * rri3 * (rri3 - 0.5) * rri        
          
          ! be aware that derivatives are calculated which are 
          ! equal to - force
          deriv(i1,:) = deriv(i1,:) - prefac * diff_vec
          deriv(i2,:) = deriv(i2,:) + prefac * diff_vec
        
          if (extra_args) then
            ! calculate potential energy
            pot_energy = pot_energy + epsilon_times4 * rri3 * (rri3 - 1.0)
          
            ! calculate sum_i=1^n-1 sum_j>i f_ij |r_ij|^2
            pressure_comp = pressure_comp + prefac * rr     
          end if    
        end if
    
  	  end do    
    
    end if
    
  end subroutine lj_deriv




!!!!!!!!!!!!!!!!!!!!!!! private functions/subroutines !!!!!!!!!!!!!!!!!!!!!!! 

  function lj_val_without_nn(str, container) result (val)
    type (structure), intent(in) :: str
		type (lj_pe_container), intent(in) :: container
    real (db) :: val
  
    integer :: i1, i2, n_tot, i
    real (db) :: r_cut, rr_cut 
    real (db) :: sigma, epsilon 
    real (db) :: epsilon_times4    ! 4 * epsilon
    real (db) :: sigma2    ! = sigma*sigma
    real (db) :: rr, rri, rri3
    real (db), dimension(ndim) :: diff_vec
    
    
    sigma = get_func_param_val(container%params, "sigma")  !container%params%p(1)%val
    epsilon = get_func_param_val(container%params, "epsilon") !container%params%p(2)%val
    epsilon_times4 = 4*epsilon
    sigma2 = sigma * sigma
    
    
    if (does_func_param_exist(container%params, "r-cut")) then
      r_cut = get_func_param_val(container%params, "r-cut")
    else
      ! otherwise put r_cut to a value higher than any possible distance
      r_cut = sqrt(sum(str%box_edges*str%box_edges))
    end if
    
    rr_cut = r_cut * r_cut
    
		n_tot = size(str%atoms)
    
    val = 0.0
    
    do i1 = 1, n_tot
      do i2 = i1+1, n_tot
        diff_vec = str%r(i1,:) - str%r(i2,:)
        
        call apply_boundary_condition_to_vector(diff_vec, str%box_edges)
        
        rr = sum(diff_vec*diff_vec)
        
        if (rr < rr_cut) then
          rri = sigma2 / rr
          rri3 = rri * rri * rri
          val = val + epsilon_times4 * rri3 * (rri3 - 1.0)
          
        end if
			
      end do
    end do
 
  end function lj_val_without_nn


  function lj_val_nn(str, container, nn_list) result (val)
    type (structure), intent(in) :: str
  	type (lj_pe_container), intent(in) :: container
    type (near_neighb_list), intent(in) :: nn_list
    real (db) :: val
  
    integer :: n_tot, j
    real (db) :: r_cut, rr_cut 
    real (db) :: sigma, epsilon 
    real (db) :: epsilon_times4    ! 4 * epsilon
    real (db) :: sigma2    ! = sigma*sigma
    real (db) :: rr, rri, rri3
    
    
    sigma = get_func_param_val(container%params, "sigma")  
    epsilon = get_func_param_val(container%params, "epsilon") 
    epsilon_times4 = 4*epsilon
    sigma2 = sigma * sigma
    
    if (does_func_param_exist(container%params, "r-cut")) then
      r_cut = get_func_param_val(container%params, "r-cut")
    else
      ! otherwise put r_cut to a value higher than any possible distance
      r_cut = sqrt(sum(str%box_edges*str%box_edges))
    end if
    
    
    ! when using neighest neighbour method all distances less than
    ! nn_list%r_cut are 'guarentied' to be in the list but not 
    ! higher distances
        
    if (nn_list%r_cut < r_cut) then
      r_cut = nn_list%r_cut
    end if
    
    
    rr_cut = r_cut * r_cut
    
  	n_tot = size(str%atoms)
    
    val = 0.0
    
    do j = 1, nn_list%n_pairs

  	  !i1 = nn_list%pairs(2*j-1)
  	  !i2 = nn_list%pairs(2*j)
      
      rr = nn_list%dists(j)
      
        
      if (rr < rr_cut) then
        rri = sigma2 / rr
        rri3 = rri * rri * rri
        val = val + epsilon_times4 * rri3 * (rri3 - 1.0)
        
      end if
  
  	end do

  end function lj_val_nn
  
end module lennard_jones_class    