module function_class
use various_constants_class
use structure_class
use func_param_class
use lennard_jones_class
use rdf_fom_class

implicit none

  public :: func_val, func_val_nn
  public :: func_deriv, func_deriv_nn


  private ::  add_lj_pe_container

  type func_list
    !private
    type (lj_pe_container), pointer :: pt_lj_pe => null()
    type (rdf_fom_container), pointer :: pt_rdf_fom => null()
  end type func_list

  interface add_potential
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
    type (func_list), intent(in) :: list
    real (db) :: val

    val = lj_val(str, list%pt_lj_pe)
  end function func_val
  
  
  function func_val_nn(str, list, nn_list) result (val)
    type (structure), intent(in) :: str
    type (func_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list    
    real (db) :: val
 
 
    if ( associated (list%pt_lj_pe) ) then
      val = lj_val_nn(str, list%pt_lj_pe, nn_list)
    else if ( associated (list%pt_rdf_fom) ) then
      val = rdf_fom_val_nn(str, list%pt_rdf_fom, nn_list)  
    else
      print *, "ERROR in func_val_nn"
      stop
    end if
    
  end function func_val_nn  
  
  
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


  subroutine func_deriv_nn(str, deriv, list, nn_list, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (func_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list
 		real (db), optional, intent(out) :: pressure_comp, pot_energy
 		    
 		! test for optional arguments
    
    if (present(pressure_comp) .and. present(pot_energy)) then
      call lj_deriv_nn(str, deriv, list%pt_lj_pe, nn_list, pressure_comp, pot_energy)
    else if (present(pressure_comp)==.false. .and. present(pot_energy)==.false.) then
      call lj_deriv_nn(str, deriv, list%pt_lj_pe, nn_list)
    else
      write(*,*) "ERROR in func_deriv_nn"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if 		    
 		     
  end subroutine func_deriv_nn   
      
end module function_class
