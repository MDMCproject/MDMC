module g_d_readers_class
use common_block_class
use various_constants_class
use converters_class
use control_containers_class
use time_correlation_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document                               


contains


  ! read data from filename into public g_d_data array in time_correlation  

  subroutine make_g_d_data_array(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    character(len=40) :: read_db
    type (dictionary_t) :: structure_attributes
    
    real(db) :: r_first, r_last, r_max
    real(db) :: t_last
    real(db) :: bin_length, delta_t
    
    integer :: n_bin, n_time_evals
    integer :: i, j
    
    integer :: number_of_g_d  ! counts # of G-d elements
    
    ! before reading the rdf data file from a-to-z determine
    ! separately how many g-d elements it contains
    
    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! get bin_length
    
    call get_node(fxml, path="//G_d-space-time-pair-correlation-function", attributes=structure_attributes, status=iostat)    
   
    call get_value(structure_attributes,"bin-length",read_db,status=iostat)
    bin_length = string_to_db(read_db)    
    
    ! count first g_d element
    
    number_of_g_d = 1
    
    call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", attributes=structure_attributes, status=iostat)    
    
    call get_value(structure_attributes,"r",read_db,status=iostat)
    r_first = string_to_db(read_db)
    
    
    ! count remainder g_d element
    
    do 
      call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", &
        attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
      
      call get_value(structure_attributes,"r",read_db,status=iostat)
      r_last = string_to_db(read_db)
      
      call get_value(structure_attributes,"t",read_db,status=iostat)
      t_last = string_to_db(read_db)            
      
      number_of_g_d = number_of_g_d + 1
    end do


    r_max = r_last+r_first
    n_bin = nint( (r_last+r_first) / bin_length )
    n_time_evals = number_of_g_d / n_bin
    delta_t = t_last / (n_time_evals-1)
    
    
    ! VERY BAD WHAT I AM DOING HERE
    
    setup_mdmc_control_params%r_max = r_max
    setup_mdmc_control_params%bin_length = bin_length
    setup_mdmc_control_params%n_time_evals = n_time_evals
    setup_mdmc_control_params%g_d_data_time_step = delta_t


    allocate(g_d_data(floor(r_max/bin_length), n_time_evals))
  
  
    call close_xmlfile(fxml)  


    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    do j = 1, n_time_evals
      do i = 1, n_bin
      
        call get_node(fxml, path="//G_d-space-time-pair-correlation-function/G-d", &
          attributes=structure_attributes, status=iostat) 
        call get_value(structure_attributes,"G",read_db,status=iostat)
        g_d_data(i, j) = string_to_db(read_db)             
        
      end do
    end do



!    ! now begin reading rdf data file from a-to-z
!
!    call open_xmlfile(filename,fxml,iostat)
!    if (iostat /= 0) stop "Cannot open file."
!
!    call xml_parse(fxml,           &
!                 begin_element_handler=begin_element, &
!                 end_element_handler=end_element,     &
!                 pcdata_chunk_handler=pcdata_chunk, &
!                 start_document=start_document, &
!                 end_document=end_document)

    
    call close_xmlfile(fxml)                             

  end subroutine make_g_d_data_array

  
  
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
!        call get_value(attributes,"title", target_rdf_fom%title, status)
!        
!        !call get_value(attributes,"r-max", read_db,status)
!        !number_db = string_to_db(read_db) 
!        
!        call get_value(attributes,"bin-length", read_db,status)
!        number_db2 = string_to_db(read_db)
!
!        
!        target_rdf_fom%rdf_data = make_rdf(product(common_config%str%box_edges), &
!                                  size(common_config%str%atoms), number_of_g_d, &
!                                  number_db2)
!                        
!        count_number_atoms = 1

                        
      case("G-d")

        
!
!        call get_value(attributes,"g", read_db,status)
!        target_rdf_fom%rdf_data%val(count_number_atoms) = string_to_db(read_db)
!   
!        count_number_atoms = count_number_atoms + 1        


!      case("scale-factor")
!        call get_value(attributes,"val", read_db,status)
!        target_rdf_fom%scale_factor = string_to_db(read_db)
!        
!              
!      case("sigma")
!        call get_value(attributes,"val", read_db,status)
!        number_db = string_to_db(read_db)
!        target_rdf_fom%weight = 1/(number_db*number_db)
        
        	
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

end module g_d_readers_class
