module func_param_class
use various_constants_class

implicit none

  public :: get_func_param_val, update_func_param

  type func_param
    character(len=99) :: name
    real(db) :: val  
    real(db) :: limit_min
    real(db) :: limit_max
  end type func_param


contains

  function get_func_param_val(vars, name) result (val)
    type (func_param), dimension(:), intent(in) :: vars
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
    
  end function get_func_param_val
  
  
  subroutine update_func_param(vars, name, val)
    type (func_param), dimension(:), intent(inout) :: vars
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

  end subroutine update_func_param
    
    
end module func_param_class
