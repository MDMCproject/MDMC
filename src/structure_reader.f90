module structure_reader_class
use flib_sax
use common_block_class ! , only : common_configuration
!use moldyn_control_class
!use structure_reader_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  

contains

  subroutine make_structure(filename)
    character(len=*), intent(in) :: filename

    write(*,*) "fisse"
  end subroutine make_structure

  !START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status
    character(len=40) :: control_object_name, read_number_atoms
    character(len=120) :: filename

    select case(name)
      case("molecule")
        call get_value(attributes,"number", read_number_atoms,status)
        !allocate(common_configuration%cf_structure%atoms(number_atoms))        
        !call make_structure(filename)
        ! get attributes
        ! my_configuration%my_structure = make_structure()
	
    end select


  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    character(len=*), intent(in)   :: name



  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
  end subroutine end_document

end module structure_reader_class
