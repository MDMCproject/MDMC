module func_params_wrapper_class
use function_class

  implicit none

  public :: move_random_func_params
  public :: print_func_param_vals
  

contains

  subroutine move_random_func_params(list)
    type (func_list), intent(inout) :: list

    !if ( associated (list%pt_lj_pe) ) then
    !  val = lj_val_nn(str, list%pt_lj_pe, nn_list)
    !else if ( associated (list%pt_rdf_fom) ) then
    !  val = rdf_fom_val_nn(str, list%pt_rdf_fom, nn_list)  
    !else
    !  print *, "ERROR in func_val_nn"
    !  stop
    !end if    
  
  end subroutine move_random_func_params
  
  
  subroutine print_func_param_vals(list)
    type (func_list), intent(inout) :: list
    
  
  end subroutine print_func_param_vals

end module func_params_wrapper_class
