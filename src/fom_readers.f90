module fom_readers_class
use common_block_class
use various_constants_class
use converters_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document
  
  integer, private :: count_number_atoms    ! keeps track of how many g-of-r
                                            ! elements have been encounted as the
                                            ! file is read

  integer, private :: number_of_g_of_r      ! used to determine at beginning 
                                            ! # of g-of-r elements in file

contains

  subroutine make_rdf_fom(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    
    ! before reading the rdf data file from a-to-z determine
    ! separately how many g-or-r elements it contains
    
    call open_xmlfile(trim(filename),fxml,iostat)
    
    number_of_g_of_r = 0
    do 
      call get_node(fxml, path="//rdf/g-of-r",status=iostat)
      if (iostat < 0) exit
      number_of_g_of_r = number_of_g_of_r + 1
    end do
  
    call close_xmlfile(fxml)  


    ! now begin reading rdf data file from a-to-z

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
    
    integer :: status
    real(db) :: number_db, number_db2
    character(len=40) :: read_db, read_int    
    

    select case(name)
      case("rdf")
        call get_value(attributes,"title", target_rdf_fom%title, status)
        
        !call get_value(attributes,"r-max", read_db,status)
        !number_db = string_to_db(read_db) 
        
        call get_value(attributes,"bin-length", read_db,status)
        number_db2 = string_to_db(read_db)

        
        target_rdf_fom%rdf_data = make_rdf(product(common_config%str%box_edges), &
                                  size(common_config%str%atoms), number_of_g_of_r, &
                                  number_db2)
                        
        count_number_atoms = 1

                        
      case("g-of-r")

        call get_value(attributes,"g", read_db,status)
        target_rdf_fom%rdf_data%val(count_number_atoms) = string_to_db(read_db)
   
        count_number_atoms = count_number_atoms + 1        


      case("scale-factor")
        call get_value(attributes,"val", read_db,status)
        target_rdf_fom%scale_factor = string_to_db(read_db)
        
              
      case("sigma")
        call get_value(attributes,"val", read_db,status)
        number_db = string_to_db(read_db)
        target_rdf_fom%weight = 1/(number_db*number_db)
        
        	
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

  end subroutine end_element


  !END_DOCUMENT
  subroutine end_document()
    use flib_sax  
  end subroutine end_document

end module fom_readers_class
