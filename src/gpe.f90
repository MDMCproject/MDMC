module gpe_class
use various_constants_class

implicit none

  ! originally I imagined calling this type parameter rather then variable
  ! but I think parameter might be a reserved keyword in fortran
  type variable
    character(len=99)  :: name
    real(db) :: val  
    real(db) :: limit_min
    real(db) :: limit_max
  end type variable


contains

  function gpe_val(str) result (val)
    type (structure), intent(in) :: str
    real (db) :: val
    
    integer :: i
    
    
  end function gpe_val
  
  
  subroutine gpe_deriv(str,deriv)
    type (structure), intent(in) :: str
    real(db), dimension(3,:), intent(out) :: deriv
    
    integer :: i
    

  end subroutine gpe_deriv
    
    
end module gpe_class
