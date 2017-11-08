! Code meant to read fnc constraints
! Current not implemented

module fnc_constraint_reader_class
use flib_sax
use common_block_class, only : common_config
use various_constants_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms

contains

  subroutine make_fnc_constraint(filename)
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat


    call open_xmlfile(filename,fxml,iostat)
    if (iostat /= 0) stop "Cannot open file."

    call xml_parse(fxml,           &
                 begin_element_handler=begin_element, &
                 end_element_handler=end_element,     &
                 pcdata_chunk_handler=pcdata_chunk, &
                 start_document=start_document, &
                 end_document=end_document)

  end subroutine make_fnc_constraint

  ! START_DOCUMENT
  !
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  !
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status
    character(len=40) :: control_object_name, read_number_atoms, read_dp
    integer :: number_atoms;
    character(len=120) :: filename
    real(db) :: number_dp(1)

    select case(name)
  
    end select

  end subroutine begin_element


  ! PCDATA_CHUNK
  !
  subroutine pcdata_chunk(chunk)
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  !
  subroutine end_element(name)
    character(len=*), intent(in)   :: name

  end subroutine end_element


  ! END_DOCUMENT
  !
  subroutine end_document()
  end subroutine end_document

end module fnc_constraint_reader_class
