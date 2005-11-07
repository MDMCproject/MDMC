module gpe_class
use various_constants_class
use structure_class
use variable_class
use lennard_jones_class

implicit none

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
  
  
  subroutine gpe_deriv(str,deriv)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    
    integer :: i
    

  end subroutine gpe_deriv
    
    
end module gpe_class
