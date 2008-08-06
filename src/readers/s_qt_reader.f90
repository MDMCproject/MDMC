! Please read Umbrello for how this code perhaps should be 
! changed in the future

module s_qt_reader_class
use common_block_class
use various_constants_class
use converters_class
use control_containers_class
use time_corr_algorithm_class

  implicit none
  
  private :: begin_element, end_element, pcdata_chunk
  private :: start_document, end_document                               


contains

  ! read S(q,t) data from file and store in target_s_qt_fom
  !
  subroutine make_s_qt_fom_container(filename)
    use flib_sax  
    use flib_xpath    
    character(len=*), intent(in) :: filename

    type(xml_t) :: fxml
    integer     :: iostat
    
    character(len=40) :: read_db, read_int
    type (dictionary_t) :: structure_attributes
    
    real(db) :: r_first, r_last, r_max
    real(db) :: t_last
    real(db) :: t_bin
    
    integer :: n_q_points, n_t_bin
    integer :: i, j
    
    integer :: number_of_SQt  ! counts # of G-d elements
    
    ! before reading the rdf data file from a-to-z determine
    ! separately how many g-d elements it contains
    
    call open_xmlfile(trim(filename),fxml,iostat)
    if (iostat /= 0) stop "Cannot open S(q,t) data file."
    
    
    ! get r_bin, density and number of atoms
    
    call get_node(fxml, path="//s-q-time", attributes=structure_attributes, status=iostat)       
    call get_value(structure_attributes,"time-bin",read_db,status=iostat)
    if (iostat /= 0) stop "No time-bin attribute in S(q,t) data file."
    t_bin = string_to_db(read_db)   
    
    !print *, "r_bin = ", r_bin
    
    call get_value(structure_attributes,"n-time-bin",read_int,status=iostat)
    if (iostat /= 0) stop "No n-time-bin attribute in S(q,t) data file."
    n_t_bin = string_to_int(read_int)
    
    call get_value(structure_attributes,"n-q-points",read_int,status=iostat)
    if (iostat /= 0) stop "No n-q-points attribute in S(q,t) data file."
    n_q_points = string_to_int(read_int)
        
    
    ! count SQt elements and populate q array
    
    allocate(target_s_qt_fom%q(n_q_points))
    number_of_SQt = 0
    do 
      call get_node(fxml, path="//s-q-time/SQt", attributes=structure_attributes, status=iostat)
      if (iostat < 0) exit
            
      call get_value(structure_attributes,"t",read_db,status=iostat)
      t_last = string_to_db(read_db)            
      
      number_of_SQt = number_of_SQt + 1
      
      if (number_of_SQt <= n_q_points) then
        call get_value(structure_attributes,"q",read_db,status=iostat)
        target_s_qt_fom%q(number_of_SQt) = string_to_db(read_db)      
      end if      
    end do

    if ( number_of_SQt /= n_q_points*n_t_bin ) then
      stop "Number of data point in S(q,t) inconsistent"
    end if
    
    
    target_s_qt_fom%n_t_bin = n_t_bin
    target_s_qt_fom%t_bin = t_bin


    allocate(target_s_qt_fom%obs(n_q_points, n_t_bin))
  
  
    call close_xmlfile(fxml)  
    

    call open_xmlfile(trim(filename),fxml,iostat)
    
    
    ! finally populate obs data note here this equals S-self+S-diff !!
    
    do j = 1, n_t_bin
      do i = 1, n_q_points
      
        call get_node(fxml, path="//s-q-time/SQt", attributes=structure_attributes, status=iostat) 
        
        call get_value(structure_attributes,"S-self",read_db,status=iostat)
        target_s_qt_fom%obs(i, j) = string_to_db(read_db)       
          
        call get_value(structure_attributes,"S-diff",read_db,status=iostat)
        target_s_qt_fom%obs(i, j) = target_s_qt_fom%obs(i, j) + string_to_db(read_db)          
        
      end do
    end do

    
    call close_xmlfile(fxml)                             

  end subroutine make_s_qt_fom_container

  
  
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

end module s_qt_reader_class
