module fom_readers_class
use common_block_class
use various_constants_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms

contains

  subroutine make_rdf_fom(filename)
    use flib_sax  
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

    
    call close_xmlfile(fxml)                             

  end subroutine make_rdf_fom

  
  
  !START_DOCUMENT
  subroutine start_document()
    use flib_sax
  end subroutine start_document


  ! BEGIN_ELEMENT
  subroutine begin_element(name,attributes)
    use flib_sax  
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes
    
    integer :: status, ndata=0
    real(db) :: number_db(1)
    integer :: number_int(1)
    character(len=40) :: read_db, read_int    
    

    select case(name)
      case("rdf")
        call get_value(attributes,"title", target_rdf_fom%title, status)
        
        call get_value(attributes,"r-max", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata) 
        
        call get_value(attributes,"n-bin", read_int,status)
        ndata = 0
        call build_data_array(read_int, number_int, ndata)

        
        target_rdf_fom%rdf_data = make_rdf(product(common_config%str%box_edges), &
                                  size(common_config%str%atoms), number_db(1), &
                                  number_int(1))
        target_rdf_fom%rdf_cal = make_rdf(product(common_config%str%box_edges), &
                                  size(common_config%str%atoms), number_db(1), &
                                  number_int(1))        
                        
        count_number_atoms = 1

                        
      case("g-of-r")

        call get_value(attributes,"g", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        target_rdf_fom%rdf_data%g_of_r(count_number_atoms) = number_db(1)
   
        count_number_atoms = count_number_atoms + 1        


      case("scale-factor")
        call get_value(attributes,"val", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        target_rdf_fom%scale_factor = number_db(1)
        
              
      case("sigma")
        call get_value(attributes,"val", read_db,status)
        ndata = 0
        call build_data_array(read_db, number_db, ndata)
        target_rdf_fom%weight = 1/(number_db(1)*number_db(1))
        
        	
    end select


  end subroutine begin_element


  ! PCDATA_CHUNK
  subroutine pcdata_chunk(chunk)
    use flib_sax
    character(len=*), intent(in) :: chunk

  end subroutine pcdata_chunk


  ! END_ELEMENT
  subroutine end_element(name)
    use flib_sax  
    character(len=*), intent(in)   :: name

    select case(name)
      case("atomArray")
        if (count_number_atoms /= size(common_config%str%atoms)+1) then
          write (*,*) "ERROR: number of atoms read in not equal to element attribute"
          write (*,*) count_number_atoms
          write (*,*) size(common_config%str%atoms)+1
          stop
        end if

    end select

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
    use flib_sax  
  end subroutine end_document

end module fom_readers_class
