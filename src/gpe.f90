module gpe_class
use various_constants_class
use structure_class
use variable_class
use lennard_jones_class

implicit none

  public :: gpe_val, gpe_val_nn
  public :: gpe_deriv, gpe_deriv_nn


  private ::  add_lj_pe_container

  type func_list
    private
    type (lj_pe_container), pointer :: pt_lj_pe => null()
  end type func_list

  interface add_potential
    module procedure add_lj_pe_container
  end interface
  
contains

  subroutine add_lj_pe_container(list, pe_target)
    type (func_list), intent(inout) :: list
    type (lj_pe_container), target, intent(in) :: pe_target
    
    list%pt_lj_pe => pe_target
  end subroutine add_lj_pe_container

  
  function gpe_val(str, list) result (val)
    type (structure), intent(in) :: str
    type (func_list), intent(in) :: list
    real (db) :: val
 
    val = lj_val(str, list%pt_lj_pe)
  end function gpe_val
  
  
  function gpe_val_nn(str, list, nn_list) result (val)
    type (structure), intent(in) :: str
    type (func_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list    
    real (db) :: val
 
    val = lj_val_nn(str, list%pt_lj_pe, nn_list)
  end function gpe_val_nn  
  
  
  subroutine gpe_deriv(str, deriv, list, pressure_comp, pot_energy)
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
      write(*,*) "ERROR in gpe_deriv"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if 
 		
  end subroutine gpe_deriv


  subroutine gpe_deriv_nn(str, deriv, list, nn_list, pressure_comp, pot_energy)
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
      write(*,*) "ERROR in gpe_deriv_nn"
      write(*,*) "Either both pressure_comp and pot_energy must be present"
      write(*,*) "none of them."
      stop    
    end if 		    
 		     
  end subroutine gpe_deriv_nn   
      
end module gpe_class
