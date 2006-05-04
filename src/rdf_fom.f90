module rdf_fom_class
use various_constants_class
use structure_class
use func_param_class
use rdf_class

implicit none

  public :: rdf_fom_val


  type rdf_fom_container
    type (rdf) :: rdf_data
    type (rdf) :: rdf_cal
    
    real(db) :: scale_factor = 1.0
    real(db) :: weight = 1.0
    
    character(len=120) :: title = " "
    
    type (func_params) :: params
  end type rdf_fom_container


contains

  function rdf_fom_val(str, c) result (val)
    type (structure), intent(in) :: str
  	type (rdf_fom_container), intent(in) :: c
    real (db) :: val
  
    real(db) :: n_bin  

    
    n_bin = size(c%rdf_cal%hist%h)
    
    val = c%weight / n_bin * &
          sum((c%rdf_data%g_of_r - c%scale_factor*c%rdf_cal%g_of_r)**2)


  end function rdf_fom_val
  
    
end module rdf_fom_class
