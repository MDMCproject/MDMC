module handler_class
use flib_sax
use common_block_class, only : my_configuration => common_configuration
use moldyn_control_class
use structure_reader_class

  implicit none
  
  public :: begin_element, end_element, pcdata_chunk
  public :: start_document, end_document
  

  !type (configuration)  :: my_configuration


  ! Used to 'simulate' dynamical polymorphism, see module
  ! configuration_class. I may not need to use targets here
  ! but could instead allocate memery directly from the
  ! constraint pointer in my_configuration rather than using
  ! target, but not convinced that this would work so for
  ! now do it using targets (I would be easy to change later)
  type (fnc_constraint), target, private :: my_fnc_constraint


  logical, private  :: in_constraints

contains

  !START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status
    character(len=40) :: control_object_name
    character(len=120) :: filename

    select case(name)
      case("structure")
        call get_value(attributes,"filename",filename,status)
        call make_structure(filename)
        ! get attributes
        ! my_configuration%my_structure = make_structure()

      case("constraints")
        in_constraints = .true.

      case("control-object")
	call get_value(attributes,"name",control_object_name,status)
	
	select case(control_object_name)
	  case("moldyn_control")
            call run_moldyn_control(my_configuration)	  

	end select
	
    end select


    ! if in constraints

    if (in_constraints) then
      select case(name)
        case("fnc_constraint")
          !perhaps call get_value(attributes, ....)
          !and something like my_fnc_constraint = make_fnc_constraint(....)
          call add_constraint(my_configuration, my_fnc_constraint)
      end select
    end if

  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    character(len=*), intent(in)   :: name

    select case(name)
      case("constraints")
        in_constraints = .false.

    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module handler_class
