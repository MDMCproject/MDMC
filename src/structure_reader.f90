module structure_reader_class
use flib_sax
use common_block_class, only : common_configuration


  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms

contains

  subroutine make_structure(filename)
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

  end subroutine make_structure

  !START_DOCUMENT
  subroutine start_document()
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, ndata=0
    character(len=40) :: control_object_name, read_number_atoms, read_dp
    integer :: number_atoms(1);
    character(len=120) :: filename
    double precision :: number_dp(1)

    select case(name)
      case("molecule")
        call get_value(attributes,"number", read_number_atoms,status)
        call build_data_array(read_number_atoms, number_atoms, ndata)
        write(*,*) number_atoms
        allocate(common_configuration%cf_structure%atoms(number_atoms(1)))
        count_number_atoms = 1

      case("atom")
        call get_value(attributes,"x3", read_dp,status)
        ndata = 0
        call build_data_array(read_dp, number_dp, ndata)
        common_configuration%cf_structure%atoms(count_number_atoms)%position(1) = number_dp(1)
        call get_value(attributes,"y3", read_dp,status)
        ndata = 0
        call build_data_array(read_dp, number_dp, ndata)
        common_configuration%cf_structure%atoms(count_number_atoms)%position(2) = number_dp(1)
        call get_value(attributes,"z3", read_dp,status)
        ndata = 0
        call build_data_array(read_dp, number_dp, ndata)
        common_configuration%cf_structure%atoms(count_number_atoms)%position(3) = number_dp(1)
        write (*,*) count_number_atoms, &
          common_configuration%cf_structure%atoms(count_number_atoms)%position(1), &
          common_configuration%cf_structure%atoms(count_number_atoms)%position(2), &
          common_configuration%cf_structure%atoms(count_number_atoms)%position(3)
        count_number_atoms = count_number_atoms + 1        

	
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
