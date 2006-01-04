module func_param_class
use various_constants_class

implicit none

  public :: get_func_param_val, update_func_param
  public :: does_func_param_exist, add_func_param

  type func_param
    private
    character(len=99) :: name
    real(db) :: val  
    logical :: fixed = .true.
    real(db) :: min_limit
    real(db) :: max_limit
    real(db) :: random_max_move
  end type func_param


  integer, parameter, private    :: MAX_ITEMS = 64
  type func_params
    integer :: number_of_items = 0
    type (func_param), dimension(MAX_ITEMS) :: p
  end type func_params


contains

  function does_func_param_exist(params, name) result (val)
    type (func_params), intent(in) :: params
    character(len=*), intent(in) :: name
    logical :: val
    
    integer :: i
    
    val = .false.
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        val = .true.
        return
      end if
    end do
  
  end function does_func_param_exist


  function get_func_param_val(params, name) result (val)
    type (func_params), intent(in) :: params
    character(len=*), intent(in) :: name
    real (db) :: val
    
    integer :: i
   
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        val = params%p(i)%val
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in get_func_param_val"
    write (*,'(3a)') "name = ", name, " not found"
    stop
    
  end function get_func_param_val


  subroutine add_func_param(params, name, val)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    real (db), intent(in) :: val
    
    integer :: n
    
    n = params%number_of_items
    
    if (n == MAX_ITEMS) then
      write(*,*) "Function parameter capacity exceeded"
      stop
    end if   

    n = n + 1
    params%p(n)%name = name
    params%p(n)%val = val
    params%number_of_items = n

  end subroutine add_func_param
  
  
  subroutine update_func_param(params, name, val)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    real (db), intent(in) :: val
    
    integer :: i
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        params%p(i)%val = val
        return
      end if
    end do
    
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in get_func_param_val"
    write (*,'(3a)') "name = ", name, " not found"
    stop  

  end subroutine update_func_param
    
    
end module func_param_class
