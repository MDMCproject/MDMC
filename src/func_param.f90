module func_param_class
use various_constants_class

implicit none

  public :: get_func_param_val
  public :: set_func_param_min, set_func_param_max
  public :: set_func_param_fixed, set_func_param_max_move
  public :: set_func_param_val
  public :: does_func_param_exist, add_func_param
  public :: random_move_func_params
  
  public :: xml_add_attribute_func_params  ! used by add_xml_attribute_func_params

  type func_param
    private
    character(len=99) :: name
    real(db) :: val  
    logical :: fixed = .true.
    real(db) :: min
    real(db) :: max
    real(db) :: max_move   ! max length this parameter is allowed to move
  end type func_param

  integer, parameter, private    :: MAX_ITEMS = 64
  type func_params
    integer :: number_of_items = 0
    type (func_param), dimension(MAX_ITEMS) :: p
  end type func_params

contains

  subroutine print_func_params(file_pointer, params)
    type (func_params), intent(in) :: params
    integer, intent(in) :: file_pointer 
    
    integer :: i
    
    write(file_pointer, *) " "
    
    do i = 1, params%number_of_items
      write(file_pointer,'(a,a,f12.6)') trim(params%p(i)%name), " = ", params%p(i)%val 
    end do
        
  end subroutine print_func_params


  ! used by add_xml_attribute_func_params()
  !
  subroutine xml_add_attribute_func_params(xf, params)
    use flib_wxml 
    type (xmlf_t), intent(inout) :: xf      
    type (func_params), intent(in) :: params

    integer :: i
       
    do i = 1, params%number_of_items
      call xml_AddAttribute(xf, trim(params%p(i)%name), str(params%p(i)%val, format="(f10.5)"))
    end do

  end subroutine xml_add_attribute_func_params


  subroutine write_func_param_to_file(params, xf)
    use flib_wxml
    type (func_params), intent(inout) :: params
    type (xmlf_t), intent(inout) :: xf   
    
    integer :: i
    
    do i = 1, params%number_of_items     
      call xml_NewElement(xf, trim(params%p(i)%name)) 
      call xml_AddAttribute(xf, "val", str(params%p(i)%val, format="(es12.5)"))
      call xml_EndElement(xf, trim(params%p(i)%name))
    end do  
  end subroutine write_func_param_to_file


  subroutine random_move_func_params(params)
    type (func_params), intent(inout) :: params
    
    integer :: i
   
    real (db) :: ran_num
    
    do i = 1, params%number_of_items
      if (.not. params%p(i)%fixed) then
        ! get a random number between 0 and 1
          
        call random_number(ran_num)  
  	    
        params%p(i)%val = params%p(i)%val + params%p(i)%max_move*(ran_num-0.5)
        
        if (params%p(i)%val < params%p(i)%min) params%p(i)%val = params%p(i)%min
        if (params%p(i)%val > params%p(i)%max) params%p(i)%val = params%p(i)%max
      end if
    end do
    
  end subroutine random_move_func_params


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


  subroutine set_func_param_min(params, name, min)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    real (db), intent(in) :: min
    
    integer :: i
   
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        params%p(i)%min = min
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in set_func_param_min"
    write (*,'(3a)') "name = ", name, " not found"
    stop
    
  end subroutine set_func_param_min
  
  
  subroutine set_func_param_max(params, name, max)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    real (db), intent(in) :: max
    
    integer :: i
   
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        params%p(i)%max = max
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in set_func_param_max"
    write (*,'(3a)') "name = ", name, " not found"
    stop
    
  end subroutine set_func_param_max
  
  
  subroutine set_func_param_fixed(params, name, fixed)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    logical, intent(in) :: fixed
    
    integer :: i
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        params%p(i)%fixed = fixed
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in set_func_param_fixed"
    write (*,'(3a)') "name = ", name, " not found"
    stop
    
  end subroutine set_func_param_fixed


  subroutine set_func_param_max_move(params, name, max_move)
    type (func_params), intent(inout) :: params
    character(len=*), intent(in) :: name
    real (db), intent(in) :: max_move
    
    integer :: i
    
    do i = 1, params%number_of_items
      if (trim(params%p(i)%name) == trim(name)) then
        params%p(i)%max_move = max_move
        return
      end if
    end do
    
    ! if hasn't returned by this stage then name was not found
    ! and a error should be thrown but for now brute force
    
    write (*,*) "ERROR in set_func_param_max_move"
    write (*,'(3a)') "name = ", name, " not found"
    stop
    
  end subroutine set_func_param_max_move
  

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
  
  
  subroutine set_func_param_val(params, name, val)
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
    
    write (*,*) "ERROR in set_func_param_val"
    write (*,'(3a)') "name = ", name, " not found"
    stop  

  end subroutine set_func_param_val
    
end module func_param_class
