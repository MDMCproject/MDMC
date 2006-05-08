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
    
    ! stuff which is in common for all function containers
    type (func_params) :: params
    type (histogram) :: hist ! used for calculating rdf_cal
  end type rdf_fom_container


contains

  function rdf_fom_val(str, c) result (val)
    type (structure), intent(in) :: str
  	type (rdf_fom_container), intent(inout) :: c
    real (db) :: val
  
    real(db) :: n_bin  
   
    n_bin = size(c%rdf_cal%val)
    
    ! in case hist%n_accum > 0 then this implies that this histogram
    ! has been calculated by some other means (that is in particular
    ! hist%sum)
     
    if (c%hist%n_accum == 0) then
      call cal_histogram(c%hist, str)
    end if
    
    call cal_rdf(c%rdf_cal, c%hist) 
    
    val = c%weight / n_bin * &
      sum((c%rdf_data%val(1:n_bin) - c%scale_factor*c%rdf_cal%val)**2)

  end function rdf_fom_val
  
    
end module rdf_fom_class
