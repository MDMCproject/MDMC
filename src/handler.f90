module handler_class
use flib_sax
use common_block_class 
use md_control_class
use structure_reader_class

  implicit none
  
  public :: begin_element, end_element, pcdata_chunk
  public :: start_document, end_document
  

  logical, private  :: in_constraints

contains

  ! START_DOCUMENT
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


      case("constraints")
        in_constraints = .true.

      case("control-object")
	call get_value(attributes,"name",control_object_name,status)
	
	select case(control_object_name)
	  case("moldyn_control")
            call run_md_control(common_configuration)	  

	end select
	
    end select


    ! if in constraints

    if (in_constraints) then
      select case(name)
        case("fnc_constraint")
          call get_value(attributes,"filename",filename,status)
          call make_fnc_constraint(filename)
          call add_constraint(common_configuration, my_fnc_constraint)
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
