module gpe_class
use various_constants_class
use structure_class
use variable_class
use lennard_jones_class

implicit none

  public :: gpe_val, gpe_val_nn
  public :: gpe_deriv, gpe_deriv_nn
  public :: gpe_md_extra, gpe_md_extra_nn


  private ::  add_lj_pe_container

  type pe_list
    private
    type (lj_pe_container), pointer :: pt_lj_pe => null()
  end type pe_list

  interface add_potential
    module procedure add_lj_pe_container
  end interface
  
contains

  subroutine add_lj_pe_container(list, pe_target)
    type (pe_list), intent(inout) :: list
    type (lj_pe_container), target, intent(in) :: pe_target
    
    list%pt_lj_pe => pe_target
  end subroutine add_lj_pe_container

  
  function gpe_val(str, list) result (val)
    type (structure), intent(in) :: str
    type (pe_list), intent(in) :: list
    real (db) :: val
 
    val = lj_val(str, list%pt_lj_pe)
  end function gpe_val
  
  
  function gpe_val_nn(str, list, nn_list) result (val)
    type (structure), intent(in) :: str
    type (pe_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list    
    real (db) :: val
 
    val = lj_val_nn(str, list%pt_lj_pe, nn_list)
  end function gpe_val_nn  
  
  
  subroutine gpe_deriv(str, deriv, list)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (pe_list), intent(in) :: list

    call lj_deriv(str, deriv, list%pt_lj_pe)
  end subroutine gpe_deriv
  
  
  subroutine gpe_md_extra(str, deriv, list, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (pe_list), intent(in) :: list
 		real (db), intent(out) :: pressure_comp, pot_energy    

    call lj_md_extra(str, deriv, list%pt_lj_pe, pressure_comp, pot_energy)
  end subroutine gpe_md_extra  


  subroutine gpe_deriv_nn(str, deriv, list, nn_list)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (pe_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list
    
    call lj_deriv_nn(str, deriv, list%pt_lj_pe, nn_list)
  end subroutine gpe_deriv_nn   
    
  subroutine gpe_md_extra_nn(str, deriv, list, nn_list, pressure_comp, pot_energy)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    type (pe_list), intent(in) :: list
    type (near_neighb_list), intent(in) :: nn_list    
 		real (db), intent(out) :: pressure_comp, pot_energy    

    call lj_md_extra_nn(str, deriv, list%pt_lj_pe, nn_list, pressure_comp, pot_energy)
  end subroutine gpe_md_extra_nn  
      
end module gpe_class
