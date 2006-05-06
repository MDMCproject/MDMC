module function_class
use various_constants_class
use structure_class
use func_param_class
use lennard_jones_class
use rdf_fom_class

implicit none

  public :: func_val
  public :: func_deriv
  
  ! for accummulating and clearing the common histogram attribute
  public :: func_accum_histogram
  public :: func_clear_histogram

  private :: add_lj_pe_container
  private :: add_rdf_fom_container

  type func_list
    !private
    type (lj_pe_container), pointer :: pt_lj_pe => null()
    type (rdf_fom_container), pointer :: pt_rdf_fom => null()
  end type func_list

  interface add_function
    module procedure add_lj_pe_container
    module procedure add_rdf_fom_container
  end interface
  
contains

  subroutine add_lj_pe_container(list, pe_target)
    type (func_list), intent(inout) :: list
    type (lj_pe_container), target, intent(in) :: pe_target
    
    list%pt_lj_pe => pe_target
  end subroutine add_lj_pe_container
  
  
  subroutine add_rdf_fom_container(list, fom_target)
    type (func_list), intent(inout) :: list
    type (rdf_fom_container), target, intent(in) :: fom_target
    
    list%pt_rdf_fom => fom_target
  end subroutine add_rdf_fom_container  

  
  function func_val(str, list) result (val)
    type (structure), intent(in) :: str
    type (func_list), intent(inout) :: list
    real (db) :: val

    if ( associated(list%pt_lj_pe) ) then
      val = lj_val(str, list%pt_lj_pe)
    else if ( associated(list%pt_rdf_fom) ) then
      val = rdf_fom_val(str, list%pt_rdf_fom)
    else
      print *, "Error in func_val"
      stop
    end if
  end function func_val
  
  
  subroutine func_accum_histogram(str, list)
    type (structure), intent(in) :: str
    type (func_list), intent(inout) :: list

    if ( associated(list%pt_lj_pe) ) then
      call accum_histogram(list%pt_lj_pe%hist, str)
    else if ( associated(list%pt_rdf_fom) ) then
      call accum_histogram(list%pt_rdf_fom%hist, str)
    else
      print *, "Error in func_val"
      stop
    end if
  end subroutine func_accum_histogram
  
  
  subroutine func_clear_histogram(str, list)
    type (structure), intent(in) :: str
    type (func_list), intent(inout) :: list

    if ( associated(list%pt_lj_pe) ) then
      call clear_histogram(list%pt_lj_pe%hist)
    else if ( associated(list%pt_rdf_fom) ) then
      call clear_histogram(list%pt_rdf_fom%hist)
    else
      print *, "Error in func_val"
      stop
    end if
  end subroutine func_clear_histogram
  
  
  subroutine func_deriv(str, deriv, list, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (func_list), intent(in) :: list
 		real (db), optional, intent(out) :: pressure_comp, pot_energy
 		
 		! test for optional arguments
    
    if (present(pressure_comp) .and. present(pot_energy)) then
      call lj_deriv(str, deriv, list%pt_lj_pe, pressure_comp, pot_energy)
    else if (present(pressure_comp)==.false. .and. present(pot_energy)==.false.) then
      call lj_deriv(str, deriv, list%pt_lj_pe)
    else
      write(*,*) "ERROR in func_deriv"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if 
 		
  end subroutine func_deriv

      
end module function_class
