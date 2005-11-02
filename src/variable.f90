module variable_class
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

  function get_variable_val(vars, name) result (val)
    type (variable), dimension(:), intent(in) :: vars
    character(len=99), intent(in) :: name
    real (db) :: val
    
    integer :: i
    
    do i = 1, size(vars)
      if (vars(i).name == name) then
        val = vars(i).val
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown
    
  end function get_variable_val
  
  
  subroutine update_variable(vars, name, val)
    type (variable), dimension(:), intent(inout) :: vars
    character(len=99), intent(in) :: name
    real (db), intent(in) :: val
    
     integer :: i
    
    do i = 1, size(vars)
      if (vars(i).name == name) then
        vars(i).val = val
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown    

  end subroutine update_variable
    
    
end module variable_class
