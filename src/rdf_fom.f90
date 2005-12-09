module rdf_fom_class
use various_constants_class
use structure_class
use func_param_class
use near_neighb_class
use rdf_class

implicit none

  public :: rdf_fom_val_nn


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


!  function rdf_fom_val(str, container) result (val)
!    type (structure), intent(in) :: str
!		type (rdf_fom_container), intent(in) :: container
!    real (db) :: val
  
 
!  end function rdf_fom_val
  
    
end module rdf_fom_class
