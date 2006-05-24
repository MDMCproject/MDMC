module func_params_wrapper_class
use flib_wxml, only : xmlf_t
use function_class
    
  implicit none

  public :: move_random_func_params
  
  public :: add_xml_attribute_func_params
  
  ! THESE TWO BELOW NEED TO BE IMPROVED
  public :: backup_func_params
  public :: restore_func_params
  
  
!  public :: save_func_param_vals
!  public :: open_file_func_param_vals
  
!  type (xmlf_t), private :: xf
!  logical, private :: is_file_open = .false. 
  
  type (func_params), private :: param_backup

contains


  ! NEED TO IMPROVE ON THIS LATER 
  subroutine backup_func_params(list)
    type (func_list), intent(in) :: list
    
    ! VERY PRIMITIVE HERE
    param_backup = list%pt_lj_pe%params
  
  end subroutine backup_func_params

  ! SUPER PRIMITIVE HERE
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
  
  subroutine add_xml_attribute_func_params(xf, list)
    use flib_wxml 
    type (xmlf_t), intent(inout) :: xf      
    type (func_list), intent(in) :: list


    if ( associated (list%pt_lj_pe) ) call xml_add_attribute_func_params(xf, list%pt_lj_pe%params)
    if ( associated (list%pt_rdf_fom) ) call xml_add_attribute_func_params(xf, list%pt_rdf_fom%params)

  
  end subroutine add_xml_attribute_func_params
    
 

    
  
!  subroutine open_file_func_param_vals(filename)
!    use flib_wxml
!    character(len=*), intent(in) :: filename  
!  
!    if (is_file_open) then
!      print *, "ERROR: a func param file is already open"
!      print *, "Cannot open another such file."
!      stop
!    else
!      call xml_OpenFile(filename, xf, indent=.true.)
!      call xml_NewElement(xf, "func-params")
!      is_file_open = .true.
!    end if
!    
!  end subroutine open_file_func_param_vals
!  
!  
!  subroutine close_file_func_param_vals()
!    use flib_wxml 
!  
!    if (is_file_open == .false.) then
!      print *, "ERROR: cannot close a func param file since"
!      print *, "such a file is currently not opened"
!      stop
!    else
!      call xml_EndElement(xf, "func-params")
!      call xml_Close(xf)
!      is_file_open = .false.
!    end if
!    
!  end subroutine close_file_func_param_vals
!  
!  
!  subroutine save_func_param_vals(list)
!    use flib_wxml
!    type (func_list), intent(inout) :: list
!    
!    if ( associated (list%pt_lj_pe) ) then
!      call write_func_param_to_file(list%pt_lj_pe%params, xf)
!    end if
!     
!    !if ( associated (list%pt_rdf_fom) ) then 
!    
!    !end if   
!
!  end subroutine save_func_param_vals

end module func_params_wrapper_class
