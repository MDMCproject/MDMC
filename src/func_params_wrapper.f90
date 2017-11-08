module func_params_wrapper_class
use flib_wxml, only : xmlf_t
use function_class
    
  implicit none

  public :: move_random_func_params
  
  public :: add_xml_attribute_func_params
  
  public :: backup_func_params
  public :: restore_func_params
  
  type (func_params), private :: param_backup

  type (func_params), private :: best_param_backup

contains

  subroutine backup_best_func_params(list)
    type (func_list), intent(in) :: list
    
    best_param_backup = list%pt_lj_pe%params
  end subroutine backup_best_func_params
  

  subroutine restore_best_func_params(list)
    type (func_list), intent(inout) :: list
    
    list%pt_lj_pe%params = best_param_backup  
  end subroutine restore_best_func_params
       

  subroutine backup_func_params(list)
    type (func_list), intent(in) :: list
    
    param_backup = list%pt_lj_pe%params
  end subroutine backup_func_params

  
  subroutine restore_func_params(list)
    type (func_list), intent(inout) :: list
    
    list%pt_lj_pe%params = param_backup
  end subroutine restore_func_params


  subroutine move_random_func_params(list)
    type (func_list), intent(inout) :: list

    if ( associated (list%pt_lj_pe) ) call random_move_func_params(list%pt_lj_pe%params)
    if ( associated (list%pt_rdf_fom) ) call random_move_func_params(list%pt_rdf_fom%params)
  
  end subroutine move_random_func_params
  
  
  ! adds an xml attribute for each param defined in list
  !
  subroutine add_xml_attribute_func_params(xf, list)
    use flib_wxml 
    type (xmlf_t), intent(inout) :: xf      
    type (func_list), intent(in) :: list

    if ( associated (list%pt_lj_pe) ) call xml_add_attribute_func_params(xf, list%pt_lj_pe%params)
    if ( associated (list%pt_rdf_fom) ) call xml_add_attribute_func_params(xf, list%pt_rdf_fom%params)

  end subroutine add_xml_attribute_func_params
    
 
  subroutine print_all_func_params(file_pointer, list)
    type (func_list), intent(in) :: list
    integer, intent(in) :: file_pointer   
        
    if ( associated (list%pt_lj_pe) ) call print_func_params(file_pointer, list%pt_lj_pe%params)
    if ( associated (list%pt_rdf_fom) ) call print_func_params(file_pointer, list%pt_rdf_fom%params)

  end subroutine print_all_func_params

end module func_params_wrapper_class
