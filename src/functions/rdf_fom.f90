module rdf_fom_class
use various_constants_class
use structure_class
use func_param_class
use rdf_class

implicit none

  public :: rdf_fom_val


  type rdf_fom_container
    type (rdf) :: rdf_data
    
    ! Notice the size of rdf_cal%val <= the size of rdf_data%val 
    type (rdf) :: rdf_cal
    
    real(db) :: scale_factor = 1.0
    real(db) :: weight = 1.0
    
    character(len=120) :: title = " "
    
    ! stuff which is in common for all function containers
    type (func_params) :: params
    !type (histogram) :: hist ! used for calculating rdf_cal
  end type rdf_fom_container

  interface rdf_fom_val
    module procedure rdf_fom_val_structure
    module procedure rdf_fom_val_histogram
  end interface 

contains

  function rdf_fom_val_structure(str, c) result (val)
    type (structure), intent(in) :: str
  	type (rdf_fom_container), intent(inout) :: c
    real (db) :: val
  
    real(db) :: n_bin  
    type (histogram) :: l_hist
   
    n_bin = size(c%rdf_cal%val)
    
    l_hist = make_histogram(n_bin, c%rdf_cal%bin_length)
    call cal_rdf(c%rdf_cal, l_hist) 
    call deallocate_histogram(l_hist)
    
    val = c%weight / n_bin * &
      sum( (c%rdf_data%val(1:n_bin) - c%scale_factor*c%rdf_cal%val)**2 )

  end function rdf_fom_val_structure
  
  function rdf_fom_val_histogram(hist, c) result (val)
    type (histogram), intent(in) :: hist
  	type (rdf_fom_container), intent(inout) :: c
    real (db) :: val
  
    real(db) :: n_bin  
   
    n_bin = size(c%rdf_cal%val)
    
    
    call cal_rdf(c%rdf_cal, hist) 
    
    val = c%weight / n_bin * &
      sum( (c%rdf_data%val(1:n_bin) - c%scale_factor*c%rdf_cal%val)**2 )

  end function rdf_fom_val_histogram
  
    
end module rdf_fom_class
