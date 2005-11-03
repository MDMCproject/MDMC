module gpe_class
use various_constants_class
use structure_class
use variable_class

implicit none

  type potential
    character(len=99)  :: name = "fisse"
    type (variable), dimension(:), allocatable :: vars
  end type potential

	integer, parameter, private    :: MAX_ITEMS = 64
	type (potential), dimension(MAX_ITEMS) :: common_pe_list

contains

  function gpe_val(str) result (val)
    type (structure), intent(in) :: str
    real (db) :: val
    
    integer :: i
    
    
    
    
  end function gpe_val
  
  
  subroutine gpe_deriv(str,deriv)
    type (structure), intent(in) :: str
    real(db), dimension(:,:), intent(out) :: deriv
    
    integer :: i
    

  end subroutine gpe_deriv
    
    
end module gpe_class
