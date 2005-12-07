module rdf_fom_class
use various_constants_class
use structure_class
use func_param_class
use near_neighb_class
use rdf_class

implicit none

  public :: rdf_fom_val_nn, rdf_fom_val
 ! public :: lj_deriv, lj_deriv_nn


  type rdf_fom_container
    type (rdf) :: rdf_data
    type (rdf) :: rdf_cal
    
    real(db) :: scale_factor = 1.0
    real(db) :: weight = 1.0
    
    character(len=120) :: title = " "
    
    type (func_params) :: params
  end type rdf_fom_container


contains

  function rdf_fom_val_nn(str, c, nn_list) result (val)
    type (structure), intent(in) :: str
  	type (rdf_fom_container), intent(in) :: c
    type (near_neighb_list), intent(in) :: nn_list
    real (db) :: val
  
    real(db) :: n_bin
  
    if (nn_list%histogram_needs_updating == .true.) then
      ! for now assume that histogram does not need updating
      
      print *, "ERROR in rdf_fom_val_nn"
      stop
    end if    

    !c%rdf_cal%g_of_r = c%rdf_cal%prefac * nn_list%hist%h
    
    n_bin = size(nn_list%hist%h)
    
    val = c%weight / n_bin * &
          sum((c%rdf_data%g_of_r - c%scale_factor*c%rdf_cal%prefac*nn_list%hist%h)**2)


  end function rdf_fom_val_nn


  function rdf_fom_val(str, container) result (val)
    type (structure), intent(in) :: str
		type (rdf_fom_container), intent(in) :: container
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
    
    ! r_cut = 2.0_db**(1.0_db / 6.0_db) * sigma
    
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
        
        do i = 1, ndim
          if (diff_vec(i) >= 0.5_db * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) - str%box_edges(i)
          end if
          if (diff_vec(i) < -0.5_db * str%box_edges(i)) then
            diff_vec(i) = diff_vec(i) + str%box_edges(i)
          end if       
        end do
        
        rr = sum(diff_vec*diff_vec)
        
        if (rr < rr_cut) then
          rri = sigma2 / rr
          rri3 = rri * rri * rri
          val = val + epsilon_times4 * rri3 * (rri3 - 1.0)
          
        end if
			
			end do
		end do
 
  end function rdf_fom_val
  
    
end module rdf_fom_class
