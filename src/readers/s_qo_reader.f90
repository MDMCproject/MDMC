module s_qo_reader_class
use common_block_class
use various_constants_class
use converters_class
use control_containers_class
use time_corr_algorithm_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document                               


contains

  ! read S(q,omega) data from file and store in target_s_qo_fom
  !
  subroutine make_s_qo_fom_container(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    character(len=40) :: read_db, read_int
    type (dictionary_t) :: structure_attributes
    
    
    integer :: n_q_points, n_omega
    integer :: i, j
    
    integer :: number_of_SQomega  
    
    
    call open_xmlfile(trim(filename),fxml,iostat)
    if (iostat /= 0) stop "Cannot open S(q,omega) data file."
    
    
    ! n_omega and n_q_points
    
    call get_node(fxml, path="//s-q-omega", attributes=structure_attributes, status=iostat)       
    call get_value(structure_attributes,"n-omega-points",read_int,status=iostat)
    if (iostat /= 0) stop "No n-omega-points attribute in S(q,omega) data file."
    n_omega = string_to_int(read_int)
    
    call get_value(structure_attributes,"n-q-points",read_int,status=iostat)
    if (iostat /= 0) stop "No n-q-points attribute in S(q,omega) data file."
    n_q_points = string_to_int(read_int)
    
    
    ! count SQomega elements and populate q array
    
    allocate(target_s_qo_fom%q(n_q_points))
    
    number_of_SQomega = 0
    do 
      call get_node(fxml, path="//s-q-omega/SQomega", attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
                        
      number_of_SQomega = number_of_SQomega + 1
      
      if (number_of_SQomega <= n_q_points) then
        call get_value(structure_attributes,"q",read_db,status=iostat)
        target_s_qo_fom%q(number_of_SQomega) = string_to_db(read_db)      
      end if      
    end do

    if ( number_of_SQomega /= n_q_points*n_omega ) then
      stop "Number of data point in S(q,omega) inconsistent"
    end if
    

    allocate(target_s_qo_fom%omega(n_omega))
    allocate(target_s_qo_fom%obs(n_q_points, n_omega))
  
  
    call close_xmlfile(fxml)  
    

    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! finally populate obs data and omega array note here this equals S-self+S-diff !!
    
    do j = 1, n_omega
      do i = 1, n_q_points
      
        call get_node(fxml, path="//s-q-omega/SQomega", attributes=structure_attributes, status=iostat) 
        
        call get_value(structure_attributes,"S-self",read_db,status=iostat)
        target_s_qo_fom%obs(i, j) = string_to_db(read_db)       
          
        call get_value(structure_attributes,"S-diff",read_db,status=iostat)
        target_s_qo_fom%obs(i, j) = target_s_qo_fom%obs(i, j) + string_to_db(read_db)     
        
        if (i == 1) then
          call get_value(structure_attributes,"omega",read_db,status=iostat)
          target_s_qo_fom%omega(j) = string_to_db(read_db)         
        end if      
        
      end do
    end do

    
    call close_xmlfile(fxml)                             

  end subroutine make_s_qo_fom_container

  
  
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
      !case("G-d")
	
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

end module s_qo_reader_class
